#!/usr/bin/env python3
"""
ISMIP7 control experiment: constant 2015 climate.

Runs a forward simulation with constant SMB (2000-2029 climatology)
and no ocean melt changes. This establishes the model drift baseline
that projections are compared against.

Time-stepping: split diagnostic (dual-form momentum) + prognostic
(DG0 upwind implicit Euler for thickness).

Usage:
    mpiexec -n 12 python scripts/run_control.py
"""

import numpy as np
import os, sys, glob
from time import perf_counter

import firedrake as fd
from firedrake import (
    Constant,
    Function,
    max_value,
    sqrt,
    derivative,
    dx,
    dS,
    ds,
    split,
    assemble,
    Mesh,
    FunctionSpace,
    VectorFunctionSpace,
    TensorFunctionSpace,
    FiniteElement,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    exp,
)
from firedrake.petsc import PETSc

import rasterio, icepack
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
    glen_flow_law as n_glen_val,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
RESULTS_DIR = os.path.join(_ROOT, "results")

sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.mpi_stats import global_mean
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.runconfig import lc as _lc, lc_coarse as _lc_coarse
from mesh_naming import get_buffer_m, mesh_filename

lc = _lc()
lc_coarse = _lc_coarse()
buffer_m = get_buffer_m()

# Simulation parameters
T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2100"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))  # years
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))  # steps between saves


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    mesh_fn = os.environ.get("ISMIP7_MESH", mesh_filename(lc_coarse, lc, buffer_m))
    PETSc.Sys.Print(f"Loading mesh: {mesh_fn}")
    mesh = Mesh(mesh_fn)
    PETSc.Sys.Print(f"  {mesh.num_vertices()} vertices, {mesh.num_cells()} cells")

    # Sidecar resolved (per-mesh preferred, parametric fallback) and
    # HARD-CHECKED against this mesh: an id absent from the mesh makes
    # ds(id) integrate to zero, i.e. silently wrong physics with no crash.
    use_calving_terminus = os.environ.get("ISMIP7_NO_CALVING_TERMINUS") is None
    bnd_ids, calving_ids, bndids_fn = load_boundary_ids(
        mesh, MESH_DIR, mesh_hint=mesh_fn,
        print_coverage=use_calving_terminus,
    )

    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    # ── Load data ──
    PETSc.Sys.Print("Loading data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:bed"), Q)
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "10.0"))
    H = Function(Q, name="thickness").interpolate(
        max_value(
            icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q),
            Constant(h_clamp),
        )
    )
    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q, name="surface").interpolate(
        max_value(b + H, (Constant(1.0) - rho_ratio) * H)
    )

    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V,
        fillvalue=0.0,
    )

    # ── Rheology from inversion ──
    inv_fn = os.path.join(MESH_DIR, f"inversion_icepack2_{lc}.h5")
    PETSc.Sys.Print(f"Loading MAP: {inv_fn}")
    with fd.CheckpointFile(inv_fn, "r") as chk:
        chk_mesh = chk.load_mesh()
        theta = chk.load_function(chk_mesh, name="log_friction")
        phi = chk.load_function(chk_mesh, name="log_fluidity")
    theta_f = Function(Q, name="theta")
    theta_f.dat.data[:] = theta.dat.data_ro
    phi_f = Function(Q, name="phi")
    phi_f.dat.data[:] = phi.dat.data_ro

    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    n = Constant(n_glen_val)
    tau_c = Constant(0.1)
    # global_mean: .dat.data_ro.mean() is the rank-local owned slice.
    # NOTE this driver has no use_residual branch, so unlike the other two
    # u_c feeds friction UNCONDITIONALLY. See icepack2_tools/mpi_stats.
    u_c = Constant(global_mean(Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    )))

    phi_eff = Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rho_W * g * max_value(Constant(0.0), -b) / (rho_I * g * max_value(H, Constant(1.0))),
            Constant(0.01),
        )
    )
    K_base = u_c / (phi_eff * tau_c) ** n
    A_map = A0 * exp(phi_f)
    K_map = K_base * exp(-n * theta_f)

    sparams = {
        "snes_type": "newtonls",
        "snes_max_it": 200,
        "snes_linesearch_type": "nleqerr",
        "snes_divergence_tolerance": -1,
        "snes_stol": 0.0,
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 400,
        "mat_mumps_icntl_24": 1,
        "mat_mumps_cntl_3": 1e-12,
    }
    fc_params = {"quadrature_degree": 4}

    # ── Accumulation (constant control climate) ──
    # For control run: zero accumulation anomaly (SMB = reference state)
    # TODO: load ISMIP7 2000-2029 climatology when available
    accum = Function(Q, name="accumulation").assign(0.0)
    PETSc.Sys.Print("  Accumulation: zero (constant control climate)")

    # ── Initial diagnostic solve ──
    PETSc.Sys.Print("Initial diagnostic solve...")
    z = Function(Z)
    z.sub(0).interpolate(Constant(0.1) * u_obs)

    h = H.copy(deepcopy=True)  # evolving thickness
    h.rename("thickness")

    u_s, M_s, tau_s = split(z)
    fields = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": h,
        "surface": s,
    }
    rheology = {
        "flow_law_exponent": n,
        "flow_law_coefficient": A_map,
        "sliding_exponent": n,
        "sliding_coefficient": K_map,
    }

    L = (
        model.minimization.viscous_power(**fields, **rheology)
        + model.minimization.friction_power(**fields, **rheology)
        + model.minimization.momentum_balance(**fields)
    )
    if use_calving_terminus:
        L += model.minimization.calving_terminus(**fields, outflow_ids=calving_ids)

    prob = NonlinearVariationalProblem(
        derivative(L, z), z, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)

    for exponent in np.linspace(1.0, n_glen_val, 5):
        n.assign(exponent)
        slvr.solve()
    PETSc.Sys.Print("  Done")

    # ── Time-stepping loop ──
    nsteps = int((T_END - T_START) / DT)
    dt_c = Constant(DT)
    PETSc.Sys.Print(
        f"\nTime-stepping: {T_START}→{T_END}, dt={DT}yr, {nsteps} steps"
    )

    # DG0 thickness for upwind transport
    Q_dg = FunctionSpace(mesh, "DG", 0)
    h_dg = Function(Q_dg, name="h_dg")
    h_dg_old = Function(Q_dg)
    phi_dg = fd.TestFunction(Q_dg)
    h_dg_trial = fd.TrialFunction(Q_dg)

    n_facet = fd.FacetNormal(mesh)

    # Volume above flotation for tracking
    s_float = Function(Q).interpolate(
        b + (rho_W / rho_I) * max_value(-b, Constant(0.0))
    )

    results = []

    for k in range(1, nsteps + 1):
        t_start = perf_counter()
        t_yr = T_START + k * DT

        # 1. Diagnostic: solve for velocity given current thickness
        try:
            slvr.solve()
        except fd.ConvergenceError:
            PETSc.Sys.Print(f"  Step {k}: re-doing continuation...")
            for exponent in np.linspace(1.0, n_glen_val, 5):
                n.assign(exponent)
                slvr.solve()

        u_vel = z.subfunctions[0]

        # 2. Prognostic: DG0 upwind implicit Euler for thickness
        h_dg.project(h)
        h_dg_old.assign(h_dg)

        un = fd.dot(u_vel, n_facet)
        un_plus = (un + abs(un)) / 2
        un_minus = (un - abs(un)) / 2

        F_prog = (
            (h_dg_trial - h_dg_old) / dt_c * phi_dg * dx
            - h_dg_trial * fd.div(u_vel * phi_dg) * dx
            + (un_plus("+") * h_dg_trial("+") - un_plus("-") * h_dg_trial("-"))
            * fd.jump(phi_dg)
            * dS
            + un_plus * h_dg_trial * phi_dg * ds
            - accum * phi_dg * dx
        )
        fd.solve(
            fd.lhs(F_prog) == fd.rhs(F_prog),
            h_dg,
            solver_parameters={
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
            },
        )

        # Clamp and project back to CG1
        h_dg.interpolate(max_value(h_dg, Constant(h_clamp)))
        h.project(h_dg)
        h.interpolate(max_value(h, Constant(h_clamp)))

        # Update surface
        s.interpolate(max_value(b + h, (Constant(1.0) - rho_ratio) * h))

        # Update phi_eff for next diagnostic solve
        phi_eff.interpolate(
            max_value(
                Constant(1.0)
                - rho_W * g * max_value(Constant(0.0), -b)
                / (rho_I * g * max_value(h, Constant(1.0))),
                Constant(0.01),
            )
        )

        t_elapsed = perf_counter() - t_start

        # Compute diagnostics
        haf = Function(Q).interpolate(max_value(s - s_float, Constant(0.0)))
        vaf = float(assemble(haf * dx)) * rho_I / 1e12 / 362.5  # mm SLE
        total_mass = float(assemble(h * dx)) * rho_I / 1e12  # Gt

        results.append((t_yr, vaf, total_mass))

        if k % OUTPUT_INTERVAL == 0 or k == 1:
            PETSc.Sys.Print(
                f"  t={t_yr:.1f}  VAF={vaf:.4f} mm SLE  "
                f"mass={total_mass:.1f} Gt  [{t_elapsed:.1f}s]"
            )

        # Periodic checkpoint
        if k % (OUTPUT_INTERVAL * 10) == 0:
            chk_fn = os.path.join(RESULTS_DIR, f"control_{lc}_t{t_yr:.0f}.h5")
            with fd.CheckpointFile(chk_fn, "w") as chk:
                chk.save_mesh(mesh)
                chk.save_function(h, name="thickness")
                chk.save_function(z.subfunctions[0], name="velocity")
            PETSc.Sys.Print(f"    [checkpoint: {chk_fn}]")

    # ── Save final state and time series ──
    PETSc.Sys.Print("\nSimulation complete.")

    chk_fn = os.path.join(RESULTS_DIR, f"control_{lc}_final.h5")
    with fd.CheckpointFile(chk_fn, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(h, name="thickness")
        chk.save_function(s, name="surface")
        chk.save_function(z.subfunctions[0], name="velocity")
    PETSc.Sys.Print(f"Saved: {chk_fn}")

    # Time series CSV
    csv_fn = os.path.join(RESULTS_DIR, f"control_{lc}_timeseries.csv")
    with open(csv_fn, "w") as f:
        f.write("year,vaf_mm_sle,mass_gt\n")
        for t, vaf, mass in results:
            f.write(f"{t:.1f},{vaf:.6f},{mass:.2f}\n")
    PETSc.Sys.Print(f"Saved: {csv_fn}")


if __name__ == "__main__":
    main()
