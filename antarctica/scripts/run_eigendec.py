#!/usr/bin/env python3
r"""
Gauss-Newton Hessian eigendecomposition for dual-control inversion.

Computes leading eigenmodes of A^{-1} H_GN where:
  H_GN = Gauss-Newton Hessian (PSD via zero-residual trick)
  A    = delta*M + gamma*K  (Laplacian prior, fenics_ice convention)

Follows Recinos et al. (2023) / fenics_ice UQ framework:
  - Prior operator A = delta*M + gamma*K per control
  - Prior covariance Gamma = A^{-1} M A^{-1} (Isaac et al. 2015)
  - DOF-vector Euclidean inner product for eigenvector projections

Modes are stored as MixedFunctions on [Q, Q] (theta, phi components).

Usage:
    python scripts/run_eigendec.py
"""

import numpy as np
import firedrake as fd
from firedrake import (
    Constant,
    Function,
    max_value,
    sqrt,
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
    MixedFunctionSpace,
)
from tlm_adjoint.firedrake import (
    reset_manager,
    start_manager,
    stop_manager,
    clear_caches,
    Functional,
    EquationSolver,
    CachedHessian,
)
from firedrake.petsc import PETSc
from scipy.sparse.linalg import LinearOperator, eigsh

import rasterio, icepack, glob, os, json
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
    glen_flow_law,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
RESULTS_DIR = os.path.join(_ROOT, "results")

import sys
sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.mpi_stats import global_max
from icepack2_tools.runconfig import lc as _lc, lc_coarse as _lc_coarse
from mesh_naming import get_buffer_m, mesh_filename

lc = _lc()
lc_coarse = _lc_coarse()
buffer_m = get_buffer_m()
K_LEADING = 40

# Prior hyperparameters (must match inversion regularization)
# Inversion uses: R = 0.5/A * gamma * ell^2 * |grad(theta)|^2 dx
# In fenics_ice convention A = delta*M + gamma*K:
#   delta = 1/A (mass weight from area normalization)
#   gamma_eff = GAMMA * ELL^2 / A (stiffness weight)
# These are computed at runtime from the mesh area.
GAMMA_THETA = 1.0
GAMMA_PHI = 1.0
ELL = 7.5e3


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # ── Load mesh + data ──
    mesh_fn = os.environ.get("ISMIP7_MESH", mesh_filename(lc_coarse, lc, buffer_m))
    PETSc.Sys.Print(f"Loading mesh: {mesh_fn}")
    mesh = Mesh(mesh_fn)

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
            Constant(1.0) - rho_W * g * max_value(Constant(0.0), -b) / (rho_I * g * H),
            Constant(0.01),
        )
    )
    A0 = Function(Q).interpolate(Constant(icepack.rate_factor(Constant(260.0))))

    n_glen = Constant(glen_flow_law)
    tau_c = Constant(0.1)
    u_c = Constant(100.0)
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
    }
    fc_params = {"quadrature_degree": 4}

    # ── Load MAP ──
    chk_fn = os.path.join(MESH_DIR, f"inversion_icepack2_{lc}.h5")
    PETSc.Sys.Print(f"Loading MAP: {chk_fn}")
    with fd.CheckpointFile(chk_fn, "r") as chk:
        chk_mesh = chk.load_mesh()
        theta_chk = chk.load_function(chk_mesh, name="log_friction")
        phi_chk = chk.load_function(chk_mesh, name="log_fluidity")
    # Project onto our mesh's function space
    theta_map = Function(Q, name="theta_map")
    theta_map.dat.data[:] = theta_chk.dat.data_ro
    phi_map = Function(Q, name="phi_map")
    phi_map.dat.data[:] = phi_chk.dat.data_ro

    # ── Warm start at MAP ──
    PETSc.Sys.Print("Warm start at MAP...")
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
        "flow_law_coefficient": A0 * fd.exp(phi_map),
        "sliding_exponent": n_glen,
        "sliding_coefficient": K_base * fd.exp(-n_glen * theta_map),
    }
    L_map = (
        model.minimization.viscous_power(**flds, **rh)
        + model.minimization.friction_power(**flds, **rh)
        + model.minimization.momentum_balance(**flds)
    )
    if use_calving_terminus:
        L_map += model.minimization.calving_terminus(**flds, outflow_ids=calving_ids)
    prob = NonlinearVariationalProblem(
        derivative(L_map, z), z, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)
    for exp in np.linspace(1.0, glen_flow_law, 5):
        n_glen.assign(exp)
        slvr.solve()
    PETSc.Sys.Print("  Done")

    # u_MAP comes from the warm start (already converged via continuation)
    u_MAP = z.subfunctions[0].copy(deepcopy=True)
    u_MAP_mag = Function(Q).interpolate(sqrt(inner(u_MAP, u_MAP)))
    PETSc.Sys.Print(f"  u_MAP: max={global_max(u_MAP_mag):.0f} m/yr")

    area_val = assemble(Constant(1.0) * dx(mesh))
    invA = Constant(1.0 / area_val)

    # ── GN Hessian via zero-residual trick ──
    # The EquationSolver needs z to already be near the solution.
    # We ensure this by doing a NonlinearVariationalSolver continuation
    # BEFORE the annotated solve. The EquationSolver then converges in ~0 steps.
    def forward_gn(theta, phi):
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
        F = derivative(L, z)
        # Continuation INSIDE annotation — each step recorded on tape
        for exp in np.linspace(1.0, glen_flow_law, 5):
            n_glen.assign(exp)
            EquationSolver(
                F == 0, z, solver_parameters=sparams, form_compiler_parameters=fc_params
            ).solve()
        u_sol, _, _ = split(z)
        J = Functional(name="J")
        J.assign(
            0.5 * invA * ((u_sol[0] - u_MAP[0]) ** 2 + (u_sol[1] - u_MAP[1]) ** 2) * dx
        )
        return J

    # Ensure z is at the MAP solution before annotated solve
    # (warm start already did this, but re-confirm with the UFL expression form)
    PETSc.Sys.Print("Pre-solving with continuation for EquationSolver...")
    stop_manager()
    K_pre = K_base * fd.exp(-n_glen * theta_map)
    A_pre = A0 * fd.exp(phi_map)
    u_s, M_s, tau_s = split(z)
    flds_pre = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": H,
        "surface": s,
    }
    rh_pre = {
        "flow_law_exponent": n_glen,
        "flow_law_coefficient": A_pre,
        "sliding_exponent": n_glen,
        "sliding_coefficient": K_pre,
    }
    L_pre = (
        model.minimization.viscous_power(**flds_pre, **rh_pre)
        + model.minimization.friction_power(**flds_pre, **rh_pre)
        + model.minimization.momentum_balance(**flds_pre)
    )
    if use_calving_terminus:
        L_pre += model.minimization.calving_terminus(
            **flds_pre, outflow_ids=calving_ids
        )
    prob_pre = NonlinearVariationalProblem(
        derivative(L_pre, z), z, form_compiler_parameters=fc_params
    )
    slvr_pre = NonlinearVariationalSolver(prob_pre, solver_parameters=sparams)
    for exp in np.linspace(1.0, glen_flow_law, 5):
        n_glen.assign(exp)
        slvr_pre.solve()
    PETSc.Sys.Print("  Done (z is at MAP with UFL expressions)")

    PETSc.Sys.Print("Recording forward at MAP for Hessian...")
    reset_manager()
    start_manager()
    J_gn = forward_gn(theta_map, phi_map)
    stop_manager()
    PETSc.Sys.Print(f"  J_GN at MAP = {float(J_gn):.6e} (should be ~0)")

    H_gn = CachedHessian(J_gn)

    def Hv_pair(v_theta, v_phi):
        _, _, ddJ = H_gn.action([theta_map, phi_map], [v_theta, v_phi])
        return ddJ[0].riesz_representation("L2"), ddJ[1].riesz_representation("L2")

    # ── Prior A^{-1} M (block-diagonal, fenics_ice convention) ──
    # A = delta*M + gamma*K where delta = 1/A, gamma = GAMMA*ELL^2/A
    delta_theta = Constant(1.0 / area_val)
    gamma_theta_eff = Constant(GAMMA_THETA * ELL**2 / area_val)
    delta_phi = Constant(1.0 / area_val)
    gamma_phi_eff = Constant(GAMMA_PHI * ELL**2 / area_val)
    PETSc.Sys.Print("Prior (fenics_ice convention):")
    PETSc.Sys.Print(
        f"  delta_theta={float(delta_theta):.6e}, gamma_theta={float(gamma_theta_eff):.6e}"
    )
    PETSc.Sys.Print(
        f"  delta_phi  ={float(delta_phi):.6e}, gamma_phi  ={float(gamma_phi_eff):.6e}"
    )

    PETSc.Sys.Print("Building prior A^{-1} M solvers...")

    def make_prior_solver(delta_c, gamma_c):
        r"""Build solver for A y = M z, i.e. y = A^{-1} M z."""
        y = Function(Q)
        rhs_f = Function(Q)
        trial_q, test_q = fd.TrialFunction(Q), fd.TestFunction(Q)
        a = (
            delta_c * inner(trial_q, test_q)
            + gamma_c * inner(grad(trial_q), grad(test_q))
        ) * dx
        L_form = inner(rhs_f, test_q) * dx
        prob = fd.LinearVariationalProblem(a, L_form, y)
        solver = fd.LinearVariationalSolver(
            prob,
            solver_parameters={
                "ksp_type": "cg",
                "ksp_rtol": 1e-10,
                "pc_type": "hypre",
                "pc_hypre_type": "boomeramg",
            },
        )

        def apply(v):
            rhs_f.assign(v)
            solver.solve()
            return y.copy(deepcopy=True)

        return apply

    Ainv_theta = make_prior_solver(delta_theta, gamma_theta_eff)
    Ainv_phi = make_prior_solver(delta_phi, gamma_phi_eff)

    # ── ARPACK eigensolve on K^{-1} H_GN ──
    PETSc.Sys.Print(f"ARPACK eigensolve for {K_LEADING} modes...")
    ndof = Q.dof_dset.size

    def matvec(x):
        v_theta = Function(Q)
        v_theta.dat.data[:] = x[:ndof]
        v_phi = Function(Q)
        v_phi.dat.data[:] = x[ndof:]
        h_theta, h_phi = Hv_pair(v_theta, v_phi)
        y_theta = Ainv_theta(h_theta)
        y_phi = Ainv_phi(h_phi)
        return np.concatenate([y_theta.dat.data_ro.copy(), y_phi.dat.data_ro.copy()])

    A_op = LinearOperator((2 * ndof, 2 * ndof), matvec=matvec, dtype=float)
    vals, vecs = eigsh(A_op, k=min(K_LEADING, 2 * ndof - 2), which="LM")

    idx = np.argsort(-vals)
    vals = vals[idx]
    vecs = vecs[:, idx]

    PETSc.Sys.Print(f"  Leading eigenvalues: {vals[:10]}")

    # ── Save modes as MixedFunctions ──
    QQ = MixedFunctionSpace([Q, Q])
    eig_fn = os.path.join(MESH_DIR, f"eigendec_{lc}.h5")
    with fd.CheckpointFile(eig_fn, "w") as chk:
        chk.save_mesh(mesh)
        for i in range(vals.size):
            mode = Function(QQ, name=f"mode_{i:04d}")
            mode.sub(0).dat.data[:] = vecs[:ndof, i]
            mode.sub(1).dat.data[:] = vecs[ndof:, i]
            chk.save_function(mode)
    PETSc.Sys.Print(f"Saved: {eig_fn}")

    eigval_fn = os.path.join(RESULTS_DIR, f"eigenvalues_{lc}.txt")
    with open(eigval_fn, "w") as f:
        for i, lam in enumerate(vals):
            f.write(f"{i} {lam:.16e}\n")
    PETSc.Sys.Print(f"Saved: {eigval_fn}")

    # Save prior hyperparameters for uq_sensitivity.py
    prior_fn = os.path.join(RESULTS_DIR, f"prior_params_{lc}.json")
    with open(prior_fn, "w") as f:
        json.dump(
            {
                "delta_theta": float(delta_theta),
                "gamma_theta": float(gamma_theta_eff),
                "delta_phi": float(delta_phi),
                "gamma_phi": float(gamma_phi_eff),
                "area": area_val,
                "GAMMA_THETA": GAMMA_THETA,
                "GAMMA_PHI": GAMMA_PHI,
                "ELL": ELL,
            },
            f,
            indent=2,
        )
    PETSc.Sys.Print(f"Saved: {prior_fn}")


def make_prior_covariance_dual(
    mesh, Q, *, delta_theta, gamma_theta, delta_phi, gamma_phi
):
    r"""Build prior covariance action Gamma = A^{-1} M A^{-1}.

    Following Isaac et al. (2015) / Recinos et al. (2023).
    Each returned callable applies the covariance operator for one control.
    """

    v = fd.TrialFunction(Q)
    w = fd.TestFunction(Q)

    def _build_cov(delta_c, gamma_c):
        a_form = (delta_c * inner(v, w) + gamma_c * inner(grad(v), grad(w))) * dx

        # Solver 1: A^{-1} (direct solve on A)
        A_mat = assemble(a_form)
        ksp = PETSc.KSP().create()
        ksp.setOperators(A_mat.petscmat)
        ksp.setType("preonly")
        pc = ksp.getPC()
        pc.setType("lu")
        pc.setFactorSolverType("mumps")
        ksp.setUp()

        # Solver 2: A^{-1} M (variational solve: A y = M z)
        y = Function(Q)
        rhs_f = Function(Q)
        L_form = inner(rhs_f, w) * dx
        prob = fd.LinearVariationalProblem(a_form, L_form, y)
        ainv_m_solver = fd.LinearVariationalSolver(
            prob,
            solver_parameters={
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
            },
        )

        def cov_action(x):
            r"""Compute Gamma x = A^{-1} M A^{-1} x."""
            # Step 1: tmp = A^{-1} x
            tmp = Function(Q)
            with x.dat.vec_ro as xv, tmp.dat.vec as tv:
                ksp.solve(xv, tv)
            # Step 2: y = A^{-1} M tmp
            rhs_f.assign(tmp)
            ainv_m_solver.solve()
            return y.copy(deepcopy=True)

        return cov_action

    cov_theta = _build_cov(Constant(delta_theta), Constant(gamma_theta))
    cov_phi = _build_cov(Constant(delta_phi), Constant(gamma_phi))
    return cov_theta, cov_phi


if __name__ == "__main__":
    main()
