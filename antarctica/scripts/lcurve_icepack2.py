#!/usr/bin/env python3
"""
L-curve analysis for icepack2 dual-control inversion.

Sweeps regularization strength gamma and records misfit vs regularization
for both friction (theta) and fluidity (phi) to identify optimal trade-off.

Usage:
    # Serial:
    python scripts/lcurve_icepack2.py
    # Parallel:
    mpiexec -n 16 python scripts/lcurve_icepack2.py
"""

import numpy as np
import firedrake as fd
from firedrake import (
    Constant,
    Function,
    max_value,
    inner,
    grad,
    derivative,
    dx,
    split,
    assemble,
    Mesh,
    FunctionSpace,
    VectorFunctionSpace,
    TensorFunctionSpace,
    FiniteElement,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
)
from tlm_adjoint.firedrake import (
    reset_manager,
    start_manager,
    stop_manager,
    clear_caches,
    compute_gradient,
    Functional,
    EquationSolver,
)
from firedrake.petsc import PETSc
from scipy.optimize import minimize as scipy_minimize

import rasterio, icepack, glob, os
import matplotlib.pyplot as plt
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as grav,
    glen_flow_law,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")

import sys
sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.runconfig import lc as _lc, lc_coarse as _lc_coarse
from mesh_naming import get_buffer_m, mesh_filename

lc = _lc()
lc_coarse = _lc_coarse()
buffer_m = get_buffer_m()

# Gamma values to sweep
GAMMAS = [1e-2, 3e-2, 1e-1, 3e-1, 1, 3, 1e1, 3e1, 1e2]
L_REG = 7.5e3
MAX_ITER = int(os.environ.get("ISMIP7_MAXITER", "500"))


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

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

    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:bed"), Q)
    H = Function(Q).interpolate(
        max_value(
            icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q),
            Constant(10.0),
        )
    )
    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q).interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))
    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V,
        fillvalue=0.0,
    )
    phi_eff = Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rho_W * grav * max_value(Constant(0.0), -b) / (rho_I * grav * H),
            Constant(0.01),
        )
    )
    A0 = Function(Q).interpolate(Constant(icepack.rate_factor(Constant(260.0))))

    n_glen = Constant(glen_flow_law)
    tau_c = Constant(0.1)
    u_c = Constant(100.0)
    # Coulomb-regularized sliding: K = u_c / (phi_eff * tau_c)^n
    # In the dual form: friction_power uses |τ|^{1+n} / K
    # phi_eff transitions smoothly from Coulomb (grounded) to free-slip (floating)
    K_base = u_c / (phi_eff * tau_c) ** n_glen

    sparams = {
        "snes_type": "newtonls",
        "snes_max_it": 200,
        "snes_linesearch_type": "nleqerr",
        "snes_divergence_tolerance": -1,
        "snes_stol": 0.0,
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 200,
        "mat_mumps_icntl_24": 1,
        "mat_mumps_cntl_3": 1e-12,
    }
    fc_params = {"quadrature_degree": 4}

    # Warm start
    PETSc.Sys.Print("Warm start...")
    stop_manager()
    z = Function(Z)
    z.sub(0).interpolate(Constant(0.1) * u_obs)

    u_s, M_s, tau_s = split(z)
    flds = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": H,
        "surface": s,
    }
    rh = {
        "flow_law_exponent": n_glen,
        "flow_law_coefficient": A0,
        "sliding_exponent": n_glen,
        "sliding_coefficient": K_base,
    }
    L_init = (
        model.minimization.viscous_power(**flds, **rh)
        + model.minimization.friction_power(**flds, **rh)
        + model.minimization.momentum_balance(**flds)
    )
    if use_calving_terminus:
        L_init += model.minimization.calving_terminus(**flds, outflow_ids=calving_ids)
    PETSc.Sys.Print(
        f"  calving_terminus: {'yes' if use_calving_terminus else 'no (buffered)'}"
    )
    prob = NonlinearVariationalProblem(
        derivative(L_init, z), z, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)
    for exp in np.linspace(1.0, glen_flow_law, 5):
        n_glen.assign(exp)
        slvr.solve()
    PETSc.Sys.Print("  Done")

    # Forward function
    area_val = assemble(Constant(1.0) * dx(mesh))

    first_call = [True]  # mutable flag for closure

    def forward(theta, phi):
        clear_caches()
        K = K_base * fd.exp(-n_glen * theta)
        A = A0 * fd.exp(phi)
        u_s, M_s, tau_s = split(z)
        flds = {
            "velocity": u_s,
            "membrane_stress": M_s,
            "basal_stress": tau_s,
            "thickness": H,
            "surface": s,
        }
        rh = {
            "flow_law_exponent": n_glen,
            "flow_law_coefficient": A,
            "sliding_exponent": n_glen,
            "sliding_coefficient": K,
        }
        L = (
            model.minimization.viscous_power(**flds, **rh)
            + model.minimization.friction_power(**flds, **rh)
            + model.minimization.momentum_balance(**flds)
        )
        if use_calving_terminus:
            L += model.minimization.calving_terminus(**flds, outflow_ids=calving_ids)
        F_fwd = derivative(L, z)
        if first_call[0]:
            # First call per gamma: use continuation for robustness
            for exp in np.linspace(1.0, glen_flow_law, 5):
                n_glen.assign(exp)
                EquationSolver(
                    F_fwd == 0,
                    z,
                    solver_parameters=sparams,
                    form_compiler_parameters=fc_params,
                ).solve()
            first_call[0] = False
        else:
            # Subsequent calls: single solve (warm-started from previous)
            EquationSolver(
                F_fwd == 0,
                z,
                solver_parameters=sparams,
                form_compiler_parameters=fc_params,
            ).solve()
        u_sol, _, _ = split(z)
        J = Functional(name="J")
        J.assign(
            0.5
            / area_val
            * ((u_sol[0] - u_obs[0]) ** 2 + (u_sol[1] - u_obs[1]) ** 2)
            * dx
        )
        return J

    # MPI helpers (same as inversion_icepack2.py)
    def func_to_global(f):
        with f.dat.vec_ro as v:
            scatter, x_seq = PETSc.Scatter.toAll(v)
            scatter.scatter(v, x_seq, mode=PETSc.Scatter.Mode.FORWARD)
            result = x_seq.array.copy()
            scatter.destroy()
            x_seq.destroy()
            return result

    def global_to_func(arr, f):
        with f.dat.vec_wo as v:
            x_seq = PETSc.Vec().createSeq(len(arr), comm=PETSc.COMM_SELF)
            x_seq.array[:] = arr
            scatter, _ = PETSc.Scatter.toAll(v)
            scatter.scatter(x_seq, v, mode=PETSc.Scatter.Mode.REVERSE)
            scatter.destroy()
            x_seq.destroy()

    # ── L-curve sweep ──
    points = []

    for gamma_val in GAMMAS:
        PETSc.Sys.Print(f"\n{'='*50}")
        PETSc.Sys.Print(f"gamma = {gamma_val:.1e}")
        PETSc.Sys.Print(f"{'='*50}")

        theta = Function(Q, name="theta")
        phi = Function(Q, name="phi")
        first_call[0] = True  # reset: first forward call uses continuation
        gamma_theta = Constant(gamma_val)
        gamma_phi = Constant(gamma_val)
        L_reg_c = Constant(L_REG)

        global_ndof = len(func_to_global(theta))

        def obj_grad(x_vec):
            global_to_func(x_vec[:global_ndof], theta)
            global_to_func(x_vec[global_ndof:], phi)

            reset_manager()
            start_manager()
            J = forward(theta, phi)
            stop_manager()
            J_val = float(J)

            dJ_dtheta, dJ_dphi = compute_gradient(J, [theta, phi])

            reg_theta = float(
                assemble(
                    0.5
                    / area_val
                    * gamma_theta
                    * L_reg_c**2
                    * inner(grad(theta), grad(theta))
                    * dx
                )
            )
            reg_phi = float(
                assemble(
                    0.5
                    / area_val
                    * gamma_phi
                    * L_reg_c**2
                    * inner(grad(phi), grad(phi))
                    * dx
                )
            )

            dR_theta = func_to_global(
                assemble(
                    1.0
                    / area_val
                    * gamma_theta
                    * L_reg_c**2
                    * inner(grad(theta), grad(fd.TestFunction(Q)))
                    * dx
                )
            )
            dR_phi = func_to_global(
                assemble(
                    1.0
                    / area_val
                    * gamma_phi
                    * L_reg_c**2
                    * inner(grad(phi), grad(fd.TestFunction(Q)))
                    * dx
                )
            )

            g_theta = func_to_global(dJ_dtheta) + dR_theta
            g_phi = func_to_global(dJ_dphi) + dR_phi

            total = J_val + reg_theta + reg_phi
            total_grad = np.concatenate([g_theta, g_phi])
            return total, total_grad

        x0 = np.concatenate([func_to_global(theta), func_to_global(phi)])
        result = scipy_minimize(
            obj_grad,
            x0,
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": MAX_ITER, "ftol": 1e-12, "gtol": 1e-8, "disp": False},
        )

        global_to_func(result.x[:global_ndof], theta)
        global_to_func(result.x[global_ndof:], phi)

        # Evaluate final misfit and regularization separately
        reset_manager()
        start_manager()
        J = forward(theta, phi)
        stop_manager()
        mis = float(J)

        reg_theta = float(
            assemble(
                0.5
                / area_val
                * gamma_theta
                * L_reg_c**2
                * inner(grad(theta), grad(theta))
                * dx
            )
        )
        reg_phi = float(
            assemble(
                0.5
                / area_val
                * gamma_phi
                * L_reg_c**2
                * inner(grad(phi), grad(phi))
                * dx
            )
        )

        points.append((gamma_val, mis, reg_theta, reg_phi))
        PETSc.Sys.Print(
            f"  misfit={mis:.6e}, reg_θ={reg_theta:.4e}, reg_φ={reg_phi:.4e}"
        )

        # Save checkpoint for this gamma
        chk_fn = os.path.join(MESH_DIR, f"lcurve_icepack2_{lc}.h5")
        mode = "a" if gamma_val != GAMMAS[0] else "w"
        with fd.CheckpointFile(chk_fn, mode) as chk:
            chk.save_mesh(mesh)
            chk.save_function(theta, name=f"gamma{gamma_val:.2e}_log_friction")
            chk.save_function(phi, name=f"gamma{gamma_val:.2e}_log_fluidity")

    # ── Plot ──
    PETSc.Sys.Print("\nL-curve results:")
    PETSc.Sys.Print(
        f"{'gamma':>10s}  {'misfit':>12s}  {'reg_theta':>12s}  {'reg_phi':>12s}"
    )
    for g, m, rt, rp in points:
        PETSc.Sys.Print(f"{g:10.2e}  {m:12.6e}  {rt:12.6e}  {rp:12.6e}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    xs = [p[1] for p in points]
    gammas = [p[0] for p in points]

    for ax, reg_idx, title in [
        (axes[0], 2, "Friction (log C)"),
        (axes[1], 3, "Fluidity (log A)"),
    ]:
        ys = [p[reg_idx] for p in points]
        ax.loglog(xs, ys, "o-", markersize=8, linewidth=1.5)
        for g, x, y in zip(gammas, xs, ys):
            ax.annotate(
                f"{g:.0e}",
                (x, y),
                textcoords="offset points",
                xytext=(6, 6),
                fontsize=8,
            )
        ax.set_xlabel("Velocity misfit")
        ax.set_ylabel("Regularization")
        ax.set_title(f"L-curve: {title}")
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out_fn = os.path.join(FIG_DIR, f"lcurve_icepack2_{lc}.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    PETSc.Sys.Print(f"Saved: {out_fn}")

    csv_fn = os.path.join(FIG_DIR, f"lcurve_icepack2_{lc}.csv")
    with open(csv_fn, "w") as f:
        f.write("gamma,misfit,reg_theta,reg_phi\n")
        for g, m, rt, rp in points:
            f.write(f"{g},{m},{rt},{rp}\n")
    PETSc.Sys.Print(f"Saved: {csv_fn}")


if __name__ == "__main__":
    main()
