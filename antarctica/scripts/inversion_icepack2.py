#!/usr/bin/env python3
"""
Joint inversion for basal friction and fluidity using icepack2 + tlm_adjoint.

Follows the Kangerd demo pattern from Shapero's dual-problems repo:
- 3-field Z = V * Sigma * T (no DG thickness in mixed space)
- Regularized Jacobian: J = J_r + α * J_1
- snes_divergence_tolerance = -1
- Sliding coefficient includes exp(m*theta)

Controls are log-deviations from PHYSICAL prior means - theta = log(C/C_w0)
on the balance-friction anchor, phi = log(A/A_prior) on a thermomechanical
fluidity prior - regularized with the fenics_ice Whittle-Matern prior
(icepack2_tools/prior.py). See antarctica/N3_FRAMEWORK.md.

Usage:
    python scripts/inversion_icepack2.py
    mpiexec -n 16 python scripts/inversion_icepack2.py
"""

import numpy as np
import os, sys, glob, json
from time import perf_counter

import firedrake as fd
from firedrake import (
    Constant,
    Function,
    max_value,
    sqrt,
    inner,
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
    COMM_WORLD,
    exp,
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

import rasterio, icepack
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
)
import colorcet as cc
import matplotlib.pyplot as plt

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")

# Repo root on the path so we can import the shared dual-friction operator.
sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.dual_friction import build_rc_residual, weertman_anchor
from icepack2_tools.prior import (
    regularization_form as reg_form,
    regularization_gradient_form as reg_grad_form,
)
from icepack2_tools.thermo_model import compute_fluidity_prior
from icepack2_tools.forcing import (load_racmo_smb_climatology,
                                    load_mean_annual_surface_temperature)
from mesh_naming import get_buffer_m, mesh_filename, bndids_filename

lc = int(os.environ.get("ISMIP7_LC", "8000"))
lc_coarse = int(os.environ.get("ISMIP7_LC_COARSE", str(lc * 10)))
buffer_m = get_buffer_m()

# Flow-law exponent default -- KEEP IN SYNC with simulation.py (the forward
# that loads this MAP must use the same rheology). THIS BRANCH: n=3 standard
# Glen, so no prefactor rescale (A4_FACTOR=1). See COMPOSITE_RHEOLOGY.md.
N_FLOW_DEFAULT = "3.0"
A4_FACTOR_DEFAULT = "1.0"


def map_n_tag():
    r"""`_n<N>` filename tag so MAPs at different flow exponents coexist;
    n=4 keeps the legacy untagged name. Mirrors simulation.map_n_tag."""
    n = float(os.environ.get("ISMIP7_N_FLOW", N_FLOW_DEFAULT))
    return "" if abs(n - 4.0) < 1e-9 else f"_n{int(round(n))}"

# Regularization -- fenics_ice Whittle-Matern prior (icepack2_tools/prior.py)
# on the log-deviation controls theta=log(C/C_w0), phi=log(A/A_prior). The
# controls sit on PHYSICAL prior means (balance friction, thermomechanical
# fluidity), so gamma is a PHYSICAL strength: gamma ~ 1e4 gives a prior std
# ~0.2-0.3 in log-units (the physical deviation scale; Recinos anchors
# sigma_C~0.22, sigma_A~0.26). The old gradient-only gamma=1 was "no prior"
# and let phi/theta blow up at n=3. Env-overridable.
GAMMA_THETA = float(os.environ.get("ISMIP7_GAMMA_THETA", "1e4"))
GAMMA_PHI = float(os.environ.get("ISMIP7_GAMMA_PHI", "1e4"))
L_REG = float(os.environ.get("ISMIP7_L_REG", "7.5e3"))

# Friction law: "budd" (power-law dual, default) or "regularized_coulomb"
# (Joughin/Schoof RC residual: grounded-only inference, exact-zero shelves).
FRICTION = os.environ.get("ISMIP7_FRICTION", "budd")
# Exact-zero-shelf residual laws share the C_w0/He/composite structure.
USE_RESIDUAL = FRICTION in ("regularized_coulomb", "budd")
USE_RC = USE_RESIDUAL  # geometry/anchor handling is shared
C0_RC = float(os.environ.get("ISMIP7_RC_C0", "0.5"))
# Buffer-node (h_clamp=0) coercivity controls; see dual_friction.build_rc_residual.
# h_visc_floor (membrane-only thickness floor) is the primary, bias-free cure;
# c_w0_floor is off by default (unnecessary once h_visc_floor is on).
RC_HVISC_FLOOR = float(os.environ.get("ISMIP7_RC_HVISC_FLOOR", "10.0"))
RC_CW0_FLOOR = float(os.environ.get("ISMIP7_RC_CW0_FLOOR", "0.0"))
# Budd N_hat knobs (fric_law="budd"): at the reference/inversion geometry
# N_hat=1 (with the PISM-delta grounded floor), so this inverts the exact-zero
# shelf He-gated law; the effective-pressure feedback is purely prognostic.
BUDD_DELTA = float(os.environ.get("ISMIP7_BUDD_DELTA", "0.02"))
BUDD_NHAT_CAP = float(os.environ.get("ISMIP7_BUDD_NHAT_CAP", "3.0"))
ALPHA_GL = (float(os.environ.get("ISMIP7_ALPHA_GL", "0.5"))
            if FRICTION == "budd" else 0.0)


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

    bndids_fn = os.environ.get("ISMIP7_BNDIDS", bndids_filename(lc_coarse, lc, buffer_m))
    with open(bndids_fn) as f:
        bnd_ids = json.load(f)

    use_calving_terminus = os.environ.get("ISMIP7_NO_CALVING_TERMINUS") is None

    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    # ── Load Data ──
    PETSc.Sys.Print("Loading data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:bed"), Q)
    # h_clamp default 0.0: invert against the *true* BedMachine geometry,
    # including h=0 over the buffered ocean region. Composite rheology
    # (added below) keeps the SNES nonsingular where h=0.
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "0.0"))
    H = Function(Q).interpolate(
        max_value(
            icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q),
            Constant(h_clamp),
        )
    )
    PETSc.Sys.Print(f"  H clamp: {h_clamp} m  "
                    f"(nodes h<=1m: {int((H.dat.data_ro <= 1.0).sum())} / "
                    f"{len(H.dat.data_ro)})")
    rho_ratio = Constant(917.0 / 1024.0)

    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q).interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))

    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V,
        fillvalue=0.0,
    )

    # Velocity-observation mask. icepack.interpolate zero-fills NODATA, so
    # without a mask unobserved regions read as "observed stationary" and
    # the inversion prescribes friction/stiff ice where there is no
    # constraint (0.7% of grounded area for MEaSUREs 450m v2). MEaSUREs
    # reports errors only where a velocity was measured, so ERR > 0 marks
    # real observations. ISMIP7_OBS_MASK=0 reverts to the unmasked misfit.
    obs_mask = Function(Q, name="obs_mask").assign(1.0)
    if os.environ.get("ISMIP7_OBS_MASK", "1") != "0":
        err = icepack.interpolate(
            (rasterio.open(f"netcdf:{vel_fn}:ERRX"),
             rasterio.open(f"netcdf:{vel_fn}:ERRY")),
            V,
            fillvalue=0.0,
        )
        emag = Function(Q).interpolate(sqrt(err[0] ** 2 + err[1] ** 2))
        obs_mask.dat.data[:] = (emag.dat.data_ro > 0.0).astype(float)
        n_no = COMM_WORLD.allreduce(int((obs_mask.dat.data_ro == 0.0).sum()))
        n_all = COMM_WORLD.allreduce(len(obs_mask.dat.data_ro))
        PETSc.Sys.Print(
            f"  Obs mask: {n_no}/{n_all} nodes without velocity obs "
            f"excluded from the misfit (regularization fills them)"
        )

    calving_ids = tuple(bnd_ids["calving"])

    # ── Rheology ──
    # Composite viscous rheology: flow exponent n_flow (this branch: n=3
    # standard Glen) main term + linear (n=1) regularization. Sliding stays
    # at Weertman m_slide=3.
    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    n_flow_val = float(os.environ.get("ISMIP7_N_FLOW", N_FLOW_DEFAULT))
    m_slide_val = float(os.environ.get("ISMIP7_M_SLIDE", "3.0"))
    a4_factor = float(os.environ.get("ISMIP7_A4_FACTOR", A4_FACTOR_DEFAULT))
    n_flow = Constant(n_flow_val)
    m_slide = Constant(m_slide_val)
    tau_c = Constant(0.1)
    PETSc.Sys.Print(
        f"  Rheology: flow n={n_flow_val:.1f}, sliding m={m_slide_val:.1f}, "
        f"A4_factor={a4_factor:.1f}"
    )

    u_speed = Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    )
    u_c = Constant(float(u_speed.dat.data_ro.mean()))
    PETSc.Sys.Print(f"  tau_c={float(tau_c):.3f} MPa, u_c={float(u_c):.1f} m/yr")

    sparams = {
        "snes_type": "newtonls",
        "snes_max_it": 200,
        "snes_linesearch_type": "nleqerr",
        "snes_divergence_tolerance": -1,
        "snes_stol": 0.0,
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 400,  # working memory increase
        "mat_mumps_icntl_24": 1,  # detect null pivots
        "mat_mumps_cntl_3": 1e-12,  # null pivot threshold
    }
    fc_params = {"quadrature_degree": 4}

    # ── Build form (Kangerd pattern: controls baked into sliding coefficient) ──
    z = Function(Z)
    z.sub(0).interpolate(Constant(0.1) * u_obs)

    theta = Function(Q, name="theta")  # log friction adjustment
    phi = Function(Q, name="phi")  # log fluidity adjustment

    # Warm-start theta/phi from a previous checkpoint if available
    warm_chk = os.environ.get("ISMIP7_WARM_START")
    if warm_chk:
        PETSc.Sys.Print(f"  Loading warm start from {warm_chk}")
        with fd.CheckpointFile(warm_chk, "r") as chk:
            chk_mesh = chk.load_mesh()
            theta_ws = chk.load_function(chk_mesh, name="log_friction")
            phi_ws = chk.load_function(chk_mesh, name="log_fluidity")
        theta.dat.data[:] = theta_ws.dat.data_ro
        phi.dat.data[:] = phi_ws.dat.data_ro
        PETSc.Sys.Print(
            f"    theta: [{theta.dat.data_ro.min():.3f}, {theta.dat.data_ro.max():.3f}]"
        )
        PETSc.Sys.Print(
            f"    phi:   [{phi.dat.data_ro.min():.3f}, {phi.dat.data_ro.max():.3f}]"
        )

    u_s, M_s, tau_s = split(z)

    fields = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": H,
        "surface": s,
    }

    # Sliding coefficient: Weertman with smooth Heaviside grounding mask
    # K = He * u_c / tau_c^n  (He = 0 floating, 1 grounded)
    # In the dual form friction_power: Ψ*(τ) = |τ|^{1+n} / ((1+n)*K)
    # He → 0 makes K → 0 which makes Ψ* → ∞, penalizing any basal
    # stress on floating ice. Add a small floor to keep K > 0.
    # Combined Ua approach: He multiplies τ in momentum_balance AND
    # K scales with 1/He so friction_power makes τ cheap on floating ice.
    # K large → τ unconstrained (floating), K small → τ penalized (grounded)
    # Together with floating=He in momentum_balance, this gives zero
    # effective friction on floating ice.
    # Smooth Heaviside sliding coefficient.
    # In the dual form: friction_power = K/(m+1)*|τ|^{m+1}
    #   Large K → τ penalized → zero friction (floating)
    #   Small K → τ allowed → normal friction (grounded)
    # He ≈ 1 grounded, He ≈ 0 floating. We use 1/(He + floor) to make
    # K large on floating ice. Floor = 0.0001 → K is 10000x larger on
    # shelves than grounded → effectively zero shelf friction.
    phi_eff = Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rho_W
            * g
            * max_value(Constant(0.0), -b)
            / (rho_I * g * max_value(H, Constant(1.0))),
            Constant(0.01),
        )
    )
    # Baseline sliding coefficient uses the Weertman m_slide=3 exponent.
    K_base = u_c / (phi_eff * tau_c) ** m_slide

    # Composite rheology following Goldsby-Kohlstedt 2001:
    #   ψ_visc  =     2·h·A /(n+1) · |M_dev|^(n+1)      (dislocation creep, n=n_flow)
    #          + α · 2·H_ref·A_1/2 · |M_dev|^2          (diffusion regularizer)
    #   ψ_fric  =     K_3 /4 · |τ|^4                    (Weertman, m=3)
    #          + α · K_1 /2 · |τ|^2                     (linear regularizer)
    # The α·(n=1, m=1) terms keep M and τ pinned where h→0 (calving front).
    # H_ref is constant so the viscous regularizer stays positive-definite
    # at h=0. α default 1e-2 because the inversion needs many forward
    # solves and SNES robustness is worth the small bias in θ, φ.
    # A_4 = a4_factor · A_3 chosen so the dislocation creep matches Glen
    # at τ = τ_c. Inverted φ absorbs any residual offset.
    alpha_reg = Constant(float(os.environ.get("ISMIP7_COMPOSITE_ALPHA", "1e-2")))
    H_ref = Constant(float(os.environ.get("ISMIP7_H_REF", "100.0")))
    PETSc.Sys.Print(
        f"  Composite rheology: alpha={float(alpha_reg):.1e}, "
        f"H_ref={float(H_ref):.0f} m"
    )

    # A4_base is the fluidity PRIOR MEAN, set below from the thermomechanical
    # A_prior (a Function) once C_w0 is available, so the control is
    # phi = log(A / A_prior). (Legacy constant baseline: A0 * a4_factor.)
    A4_base = None

    def _rheo_glen(theta_c, phi_c):
        # Main dislocation-creep flow law (n=n_flow, this branch n=3) + Weertman sliding (m=3)
        return {
            "flow_law_exponent": n_flow,
            "flow_law_coefficient": A4_base * exp(phi_c),
            "sliding_exponent": m_slide,
            "sliding_coefficient": K_base * exp(-m_slide * theta_c),
        }

    def _rheo_linear(theta_c, phi_c):
        # Linearize about τ = τ_c so the n=1 / m=1 powers agree with the
        # n_flow / m_slide powers at the calibration stress; this gives a
        # smooth crossover rather than a kink near τ_c.
        A_lin = A4_base * exp(phi_c) * tau_c ** (n_flow_val - 1)
        K_lin = u_c / (phi_eff * tau_c) * exp(-theta_c)
        return {
            "flow_law_exponent": Constant(1.0),
            "flow_law_coefficient": A_lin,
            "sliding_exponent": Constant(1.0),
            "sliding_coefficient": K_lin,
        }

    def _build_action(theta_c, phi_c, fields_):
        fields_reg = dict(fields_)
        fields_reg["thickness"] = H_ref
        rheo_glen = _rheo_glen(theta_c, phi_c)
        rheo_lin = _rheo_linear(theta_c, phi_c)
        L_ = (
            model.minimization.viscous_power(**fields_, **rheo_glen)
            + alpha_reg * model.minimization.viscous_power(
                **fields_reg, **rheo_lin)
            + model.minimization.friction_power(**fields_, **rheo_glen)
            + alpha_reg * model.minimization.friction_power(
                **fields_, **rheo_lin)
            + model.minimization.momentum_balance(**fields_)
        )
        if use_calving_terminus:
            L_ += model.minimization.calving_terminus(
                **fields_, outflow_ids=calving_ids
            )
        return L_

    # The residual laws need a fixed Weertman anchor C_w0 (driving-stress
    # balance); theta then inverts as an O(1) log-adjustment on top of it.
    C_w0 = weertman_anchor(H, s, u_obs, m_slide_val, Q)
    if USE_RESIDUAL:
        law_name = ("Budd N_hat (exact-zero shelf, delta="
                    f"{BUDD_DELTA:.3f}, alpha_gl={ALPHA_GL:.2f})"
                    if FRICTION == "budd"
                    else f"regularized Coulomb (c0={C0_RC})")
        PETSc.Sys.Print(
            f"  Friction: {law_name}; h_visc_floor={RC_HVISC_FLOOR:.0f}m; "
            f"C_w0 in [{float(C_w0.dat.data_ro.min()):.2e}, "
            f"{float(C_w0.dat.data_ro.max()):.2e}]"
        )
    else:
        PETSc.Sys.Print("  Friction: Budd power-law dual (legacy action)")

    # Physical FLUIDITY PRIOR MEAN A_prior(x): a fixed-velocity thermomechanical
    # (Stefan enthalpy) solve at the observed geometry/velocity, so the control
    # phi = log(A / A_prior) is a small deviation from a physically-motivated
    # fluidity rather than log(A / const). This is the Recinos/fenics_ice fix
    # for the n=3 blow-up (a constant A0 baseline forced phi to carry all the
    # spatial fluidity structure). Frictional heating uses the balance C_w0.
    if os.environ.get("ISMIP7_FLUIDITY_PRIOR", "thermo") == "thermo":
        acc_prior = load_racmo_smb_climatology(Q)
        T_srf = load_mean_annual_surface_temperature(Q)
        A_prior = compute_fluidity_prior(u_obs, H, s, b, C_w0, acc_prior, T_srf)
        A_prior.rename("fluidity_prior")
        PETSc.Sys.Print(
            f"  Fluidity prior A_prior in [{float(A_prior.dat.data_ro.min()):.2f}, "
            f"{float(A_prior.dat.data_ro.max()):.2f}] (thermomechanical)"
        )
    else:
        A_prior = Function(Q, name="fluidity_prior").interpolate(A0 * Constant(a4_factor))
        PETSc.Sys.Print("  Fluidity prior: constant A0*a4_factor (legacy)")
    A4_base = A_prior
    map_tag = ({"regularized_coulomb": "_rc", "budd": "_budd"}.get(FRICTION, "")
               + map_n_tag())

    def build_F(theta_c, phi_c):
        # Residual closure (tau linear, grounded-only theta via exp(theta*He),
        # exact-zero shelves): budd -> N_hat=1 at the reference geometry;
        # regularized_coulomb -> Coulomb cap. Legacy budd -> action derivative.
        if USE_RESIDUAL:
            return build_rc_residual(
                z, theta_c, phi_c, H=H, s=s, b=b, C_w0=C_w0,
                A4_base=A4_base, n_flow=n_flow, n_flow_val=n_flow_val,
                m_slide=m_slide_val, tau_c=tau_c, alpha=alpha_reg, H_ref=H_ref,
                fric_law=FRICTION, N_ref=None,
                nhat_floor=BUDD_DELTA, nhat_cap=BUDD_NHAT_CAP, alpha_gl=ALPHA_GL,
                c0=C0_RC, c_w0_floor=RC_CW0_FLOOR, h_visc_floor=RC_HVISC_FLOOR,
                calving_ids=calving_ids if use_calving_terminus else None,
            )
        return derivative(_build_action(theta_c, phi_c, fields), z)

    if use_calving_terminus:
        PETSc.Sys.Print("  Using calving_terminus BC")
    else:
        PETSc.Sys.Print("  NO calving_terminus BC (buffered mesh, h=0 at front)")
    F = build_F(theta, phi)

    # ── Warm start ──
    stop_manager()
    prob = NonlinearVariationalProblem(F, z, form_compiler_parameters=fc_params)
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)
    # Always use continuation — single solve at full exponents can fail
    # with checkpoint parameters that create ill-conditioned systems.
    # Ramp n_flow (1 → n_flow_val) and m_slide (1 → m_slide_val) together.
    PETSc.Sys.Print(
        f"Warm start (continuation n_flow 1→{n_flow_val:.1f}, "
        f"m_slide 1→{m_slide_val:.1f})..."
    )
    for t in np.linspace(0.0, 1.0, 5):
        n_flow.assign(1.0 + t * (n_flow_val - 1.0))
        m_slide.assign(1.0 + t * (m_slide_val - 1.0))
        slvr.solve()
    PETSc.Sys.Print("  Done")

    u_init = z.subfunctions[0]
    u_mag = Function(Q).interpolate(sqrt(inner(u_init, u_init)))
    PETSc.Sys.Print(f"  u_max = {float(u_mag.dat.data_ro.max()):.0f} m/yr")

    # ── Forward function for tlm_adjoint ──
    # Normalize by the OBSERVED area so the misfit magnitude stays
    # comparable between masked and unmasked runs.
    area_val = assemble(obs_mask * dx(mesh))

    def forward(theta_ctrl, phi_ctrl):
        clear_caches()
        F_ctrl = build_F(theta_ctrl, phi_ctrl)
        # Continuation inside annotation for robustness — ramp both
        # n_flow and m_slide on the same [0,1] parameter.
        for t in np.linspace(0.0, 1.0, 5):
            n_flow.assign(1.0 + t * (n_flow_val - 1.0))
            m_slide.assign(1.0 + t * (m_slide_val - 1.0))
            EquationSolver(
                F_ctrl == 0,
                z,
                solver_parameters=sparams,
                form_compiler_parameters=fc_params,
            ).solve()

        u_sol, _, _ = split(z)
        J = Functional(name="J")
        J.assign(
            0.5
            / area_val
            * obs_mask
            * ((u_sol[0] - u_obs[0]) ** 2 + (u_sol[1] - u_obs[1]) ** 2)
            * dx
        )
        return J

    # ── MPI helpers ──
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

    # ── L-BFGS-B Inversion ──
    max_iter = int(os.environ.get("ISMIP7_MAXITER", "500"))
    PETSc.Sys.Print(f"\nStarting L-BFGS-B inversion (theta + phi)...")
    PETSc.Sys.Print(f"  maxiter={max_iter}, nranks={COMM_WORLD.size}")

    global_ndof = len(func_to_global(theta))
    z_backup = z.copy(deepcopy=True)
    last_good_obj = [np.inf]
    iteration_count = [0]

    def objective_and_gradient(x_vec):
        t_iter = perf_counter()
        global_to_func(x_vec[:global_ndof], theta)
        global_to_func(x_vec[global_ndof:], phi)

        t_fwd = perf_counter()
        reset_manager()
        start_manager()
        try:
            J = forward(theta, phi)
        except fd.ConvergenceError:
            stop_manager()
            z.assign(z_backup)
            if iteration_count[0] == 0:
                # No successful evaluation yet: a large-obj + zero-grad return
                # would make L-BFGS-B declare convergence at iteration 0 and
                # write a theta=phi=0 garbage MAP. Fail loudly instead.
                raise RuntimeError(
                    "First forward solve failed - the inversion cannot start. "
                    "For a small mesh (e.g. 32 km) try fewer MPI ranks (MUMPS "
                    "is fragile at <~100 vertices/rank); also check the "
                    "fluidity prior."
                )
            PETSc.Sys.Print("  [!] Forward solve failed, returning large objective")
            return last_good_obj[0] * 10, np.zeros(2 * global_ndof)
        stop_manager()
        J_val = float(J)
        t_fwd = perf_counter() - t_fwd

        z_backup.assign(z)
        last_good_obj[0] = J_val

        t_adj = perf_counter()
        try:
            dJ_dtheta, dJ_dphi = compute_gradient(J, [theta, phi])
        except fd.ConvergenceError:
            # The adjoint jacobian solve can hit the SNES cap at rough
            # mid-optimization controls just like the forward (killed the
            # Jul 18 2500m Budd run at iter 60, ~11 hr in). Same recovery
            # as the forward guard: inflated objective + zero gradient
            # makes L-BFGS-B backtrack its line search.
            z.assign(z_backup)
            if iteration_count[0] == 0:
                raise RuntimeError(
                    "First adjoint solve failed - the inversion cannot start. "
                    "For a small mesh (e.g. 32 km) try fewer MPI ranks (MUMPS "
                    "is fragile at <~100 vertices/rank)."
                )
            PETSc.Sys.Print("  [!] Adjoint solve failed, returning large objective")
            return last_good_obj[0] * 10, np.zeros(2 * global_ndof)
        t_adj = perf_counter() - t_adj

        # Whittle-Matern prior energy + gradient (icepack2_tools/prior.py):
        # 0.5/area * gamma * (theta^2 + L^2 |grad theta|^2). The theta^2 mass
        # term (absent from the old gradient-only form) removes the null space
        # that let phi/theta drift to +-36 at n=3; combined with the physical
        # prior means it constrains the DEVIATION, not the amplitude.
        _psi = fd.TestFunction(Q)
        reg_theta = float(assemble(reg_form(theta, GAMMA_THETA, area_val, L_reg=L_REG)))
        reg_phi = float(assemble(reg_form(phi, GAMMA_PHI, area_val, L_reg=L_REG)))
        dR_theta = assemble(reg_grad_form(theta, _psi, GAMMA_THETA, area_val, L_reg=L_REG))
        dR_phi = assemble(reg_grad_form(phi, _psi, GAMMA_PHI, area_val, L_reg=L_REG))

        g_theta = func_to_global(dJ_dtheta) + func_to_global(dR_theta)
        g_phi = func_to_global(dJ_dphi) + func_to_global(dR_phi)

        total = J_val + reg_theta + reg_phi
        total_grad = np.concatenate([g_theta, g_phi])

        t_iter = perf_counter() - t_iter
        iteration_count[0] += 1
        PETSc.Sys.Print(
            f"  iter {iteration_count[0]:3d}: "
            f"misfit={J_val:.6e} reg_θ={reg_theta:.4e} reg_φ={reg_phi:.4e} "
            f"total={total:.6e} |grad|={np.linalg.norm(total_grad):.4e} "
            f"[fwd={t_fwd:.1f}s adj={t_adj:.1f}s total={t_iter:.1f}s]"
        )

        # Periodic checkpoint every 20 iterations
        if iteration_count[0] % 20 == 0:
            chk_fn = os.path.join(MESH_DIR, f"inversion_icepack2{map_tag}_{lc}.h5")
            with fd.CheckpointFile(chk_fn, "w") as chk:
                chk.save_mesh(mesh)
                chk.save_function(theta, name="log_friction")
                chk.save_function(phi, name="log_fluidity")
                chk.save_function(u_obs, name="velocity_obs")
                chk.save_function(obs_mask, name="obs_mask")
                chk.save_function(H, name="thickness")
                chk.save_function(b, name="bed")
                chk.save_function(s, name="surface")
                chk.save_function(A_prior, name="fluidity_prior")
                # Mesh provenance: the forward run derives the boundary_ids
                # sidecar path from these, never from its own environment.
                chk.set_attr("/", "lc_coarse", int(lc_coarse))
                chk.set_attr("/", "buffer_m", float(buffer_m))
            PETSc.Sys.Print(f"    [checkpoint saved: iter {iteration_count[0]}]")

        return total, total_grad

    x0 = np.concatenate([func_to_global(theta), func_to_global(phi)])
    result = scipy_minimize(
        objective_and_gradient,
        x0,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": max_iter, "ftol": 0, "gtol": 0, "disp": False},
    )

    PETSc.Sys.Print(f"\nOptimization finished: {result.message}")
    PETSc.Sys.Print(f"  {result.nit} iterations, {result.nfev} function evaluations")
    global_to_func(result.x[:global_ndof], theta)
    global_to_func(result.x[global_ndof:], phi)
    PETSc.Sys.Print(
        f"  theta range: [{float(theta.dat.data_ro.min()):.3f}, {float(theta.dat.data_ro.max()):.3f}]"
    )
    PETSc.Sys.Print(
        f"  phi range:   [{float(phi.dat.data_ro.min()):.3f}, {float(phi.dat.data_ro.max()):.3f}]"
    )

    # ── Save MAP immediately ──
    chk_fn = os.path.join(MESH_DIR, f"inversion_icepack2{map_tag}_{lc}.h5")
    with fd.CheckpointFile(chk_fn, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(theta, name="log_friction")
        chk.save_function(phi, name="log_fluidity")
        chk.save_function(u_obs, name="velocity_obs")
        chk.save_function(obs_mask, name="obs_mask")
        chk.save_function(H, name="thickness")
        chk.save_function(b, name="bed")
        chk.save_function(s, name="surface")
        chk.save_function(A_prior, name="fluidity_prior")
        # Mesh provenance: the forward run derives the boundary_ids sidecar
        # path from these, never from its own environment.
        chk.set_attr("/", "lc_coarse", int(lc_coarse))
        chk.set_attr("/", "buffer_m", float(buffer_m))
    PETSc.Sys.Print(f"Saved MAP: {chk_fn}")

    # ── Final forward solve ──
    PETSc.Sys.Print("\nFinal forward solve...")
    stop_manager()
    try:
        slvr.solve()
    except fd.ConvergenceError:
        PETSc.Sys.Print("  Final solve failed, using last optimization state")

    u_sol = z.subfunctions[0]
    u_sol_mag = Function(Q).interpolate(sqrt(u_sol[0] ** 2 + u_sol[1] ** 2))
    misfit = float(
        assemble(
            0.5
            / area_val
            * obs_mask
            * ((u_sol[0] - u_obs[0]) ** 2 + (u_sol[1] - u_obs[1]) ** 2)
            * dx
        )
    )
    PETSc.Sys.Print(f"  Final misfit (masked): {misfit:.6e}")

    # Update checkpoint with velocity
    with fd.CheckpointFile(chk_fn, "a") as chk:
        chk.save_function(u_sol, name="velocity")
    PETSc.Sys.Print(f"Saved velocity: {chk_fn}")

    # ── Plot ──
    PETSc.Sys.Print("Plotting...")
    u_speed = Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    )
    coords = mesh.coordinates.dat.data_ro
    xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
    ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
    pad = 0.02 * max(xmax - xmin, ymax - ymin)

    pw, ph, bot = 0.25, 0.75, 0.10
    gap1, gap2 = 0.03, 0.10
    cb_w, cb_gap = 0.012, 0.015
    cb_h = ph * 0.5
    cb_bot = bot + (ph - cb_h) / 2
    x0 = 0.04
    x1 = x0 + pw + gap1
    x_cb1 = x1 + pw + cb_gap
    x2 = x_cb1 + cb_w + gap2
    x_cb2 = x2 + pw + cb_gap

    fig = plt.figure(figsize=(22, 7))
    ax0 = fig.add_axes([x0, bot, pw, ph])
    ax1 = fig.add_axes([x1, bot, pw, ph])
    cax1 = fig.add_axes([x_cb1, cb_bot, cb_w, cb_h])
    ax2 = fig.add_axes([x2, bot, pw, ph])
    cax2 = fig.add_axes([x_cb2, cb_bot, cb_w, cb_h])
    for ax in [ax0, ax1, ax2]:
        ax.set_aspect("equal")
        ax.set_xlim(xmin - pad, xmax + pad)
        ax.set_ylim(ymin - pad, ymax + pad)
        fd.triplot(
            mesh,
            axes=ax,
            interior_kw={"linewidth": 0.1, "alpha": 0.3, "color": "k"},
            boundary_kw={"linewidth": 1.0, "color": "k"},
        )
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])

    fd.tripcolor(u_speed, vmin=0, vmax=1000, axes=ax0, cmap=cc.cm.CET_L19)
    ax0.set_title("Observed speed")
    cs = fd.tripcolor(u_sol_mag, vmin=0, vmax=1000, axes=ax1, cmap=cc.cm.CET_L19)
    ax1.set_title("Inverted speed (icepack2)")
    diff = Function(Q).interpolate(u_speed - u_sol_mag)
    cd = fd.tripcolor(diff, vmin=-500, vmax=500, axes=ax2, cmap=cc.cm.CET_CBTD1)
    ax2.set_title("Observed - Modeled")
    fig.colorbar(cs, cax=cax1, label="m/yr")
    fig.colorbar(cd, cax=cax2, label="m/yr")

    out_fn = os.path.join(FIG_DIR, f"inversion_icepack2_{lc}.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    PETSc.Sys.Print(f"Saved: {out_fn}")


if __name__ == "__main__":
    main()
