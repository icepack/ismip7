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
import os, sys, glob
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
    TestFunction,
    FacetNormal,
    dot,
    jump,
    ds,
    dS,
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
# colorcet/matplotlib are imported lazily at the plotting call below. They are
# needed only for the summary figure, so a missing optional plotting dependency
# must not abort a multi-hour inversion at import time (it did: colorcet is
# absent from venv-firedrake-2026 and every rank died before loading data).

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")

# Repo root on the path so we can import the shared dual-friction operator.
sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.dual_friction import build_rc_residual, weertman_anchor
from icepack2_tools.geometry import cg1_lift, sample_to_geometry
from icepack2_tools.grounding import height_above_flotation
from icepack2_tools.mpi_stats import (global_mean, global_range,
                                      global_max, global_size, global_count)
from icepack2_tools.naming import map_basename
from icepack2_tools.runconfig import (
    friction as _friction, geometry_space as _geometry_space,
    lc as _lc, lc_coarse as _lc_coarse, n_flow as _n_flow,
)
from icepack2_tools.prior import (
    regularization_form as reg_form,
    regularization_gradient_form as reg_grad_form,
)
from icepack2_tools.thermo_model import compute_fluidity_prior
from icepack2_tools.forcing import (load_racmo_smb_climatology,
                                    load_mean_annual_surface_temperature)
from mesh_naming import get_buffer_m, mesh_filename

lc = _lc()
lc_coarse = _lc_coarse()
buffer_m = get_buffer_m()

# Flow-law exponent: owned by icepack2_tools.runconfig, which the forward that
# loads this MAP reads too, so the two cannot drift apart. THIS BRANCH: n=3
# standard Glen, so no prefactor rescale (A4_FACTOR=1). See
# COMPOSITE_RHEOLOGY.md.
A4_FACTOR_DEFAULT = "1.0"
A4_FACTOR_N4 = "10.0"


def a4_factor_default():
    r"""Prefactor default derived from the flow exponent (10 at n=4, 1
    otherwise) so ISMIP7_N_FLOW=4 alone still gets the n=4 rescale.
    ISMIP7_A4_FACTOR still overrides. Mirrors simulation.a4_factor_default."""
    n = _n_flow()
    return A4_FACTOR_N4 if abs(n - 4.0) < 1e-9 else A4_FACTOR_DEFAULT

# Misfit normalization: "sigma" (default) divides each residual by the squared
# observational error of its own datum, so the misfit is a dimensionless chi^2
# and terms with different units can be traded off; "none" is the legacy
# dimensional functional. Read here (not only where the misfit is built)
# because the regularization defaults below are tied to it.
MISFIT_NORM = os.environ.get("ISMIP7_MISFIT_NORM", "sigma").lower()
if MISFIT_NORM not in ("sigma", "none"):
    raise ValueError(
        f"ISMIP7_MISFIT_NORM must be 'sigma' or 'none', got {MISFIT_NORM!r}"
    )

# Regularization -- fenics_ice Whittle-Matern prior (icepack2_tools/prior.py)
# on the log-deviation controls theta=log(C/C_w0), phi=log(A/A_prior). The
# controls sit on PHYSICAL prior means (balance friction, thermomechanical
# fluidity), so gamma is a PHYSICAL strength: gamma ~ 1e4 gives a prior std
# ~0.2-0.3 in log-units (the physical deviation scale; Recinos anchors
# sigma_C~0.22, sigma_A~0.26). The old gradient-only gamma=1 was "no prior"
# and let phi/theta blow up at n=3. Env-overridable.
#
# THE DEFAULT IS COUPLED TO ISMIP7_MISFIT_NORM. Normalizing by sigma divides
# the misfit by ~sigma^2 (MEaSUREs per-component sigma is median 2.6 m/yr, p90
# 6.6 m/yr, so roughly one order of magnitude), which weakens the prior by the
# same factor if gamma is held fixed. Carrying the unnormalized 1e4 into the
# chi^2 functional gave control ranges of +/-5 to 6 log-units; 1e5 under sigma
# recovered the plausible +/-3 seen on two converged runs (32 km and 2500 m).
# So the default tracks the norm rather than leaving the shipped configuration
# knowingly mistuned. Setting ISMIP7_GAMMA_* explicitly overrides either way.
GAMMA_DEFAULT = "1e5" if MISFIT_NORM == "sigma" else "1e4"
GAMMA_THETA = float(os.environ.get("ISMIP7_GAMMA_THETA", GAMMA_DEFAULT))
GAMMA_PHI = float(os.environ.get("ISMIP7_GAMMA_PHI", GAMMA_DEFAULT))
L_REG = float(os.environ.get("ISMIP7_L_REG", "7.5e3"))

# Friction law: "budd" (power-law dual, default) or "regularized_coulomb"
# (Joughin/Schoof RC residual: grounded-only inference, exact-zero shelves).
FRICTION = _friction()
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

    use_calving_terminus = os.environ.get("ISMIP7_NO_CALVING_TERMINUS") is None
    # Per-mesh sidecar, hard-checked against this mesh: a stale sidecar leaves
    # most of the front with no terminus back-pressure, and the inversion would
    # absorb that into theta/phi where the t=0 misfit cannot reveal it.
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

    # Geometry space, matching simulation.py. This MUST agree with the forward
    # that loads the MAP: the inversion absorbs whatever the front treatment
    # gets wrong into theta/phi, so a MAP inverted under CG1 geometry silently
    # carries the lumped lift's inflated calving-front thickness into the
    # friction field, and the t=0 velocity misfit cannot reveal it.
    geometry_space = _geometry_space()
    geom_dg = geometry_space == "dg0"
    Q_g = FunctionSpace(mesh, "DG", 0) if geom_dg else Q
    PETSc.Sys.Print(f"  Geometry space: {geometry_space.upper()}")

    # ── Load Data ──
    PETSc.Sys.Print("Loading data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    # Cell average onto the geometry space, NOT a centroid point sample --
    # see geometry.sample_to_geometry for the measurements behind that.
    b = sample_to_geometry(
        lambda sp: icepack.interpolate(
            rasterio.open(f"netcdf:{bm_fn}:bed"), sp), Q_g, Q)
    # h_clamp default 0.0: invert against the *true* BedMachine geometry,
    # including h=0 over the buffered ocean region. Composite rheology
    # (added below) keeps the SNES nonsingular where h=0.
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "0.0"))
    H = sample_to_geometry(
        lambda sp: icepack.interpolate(
            rasterio.open(f"netcdf:{bm_fn}:thickness"), sp),
        Q_g, Q, floor=h_clamp)
    PETSc.Sys.Print(f"  H clamp: {h_clamp} m  "
                    f"(nodes h<=1m: "
                    f"{global_count(H.dat.data_ro <= 1.0, mesh.comm)} / "
                    f"{global_size(H)})")
    rho_ratio = Constant(917.0 / 1024.0)

    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q_g).interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))

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
    err = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:ERRX"),
         rasterio.open(f"netcdf:{vel_fn}:ERRY")),
        V,
        fillvalue=0.0,
    )
    # ONE definition of "this node carries a real velocity observation",
    # consumed by both the mask below and the sigma normalization further down.
    # The two used to derive it independently from `err` and could drift apart.
    emag = Function(Q, name="obs_err_mag").interpolate(
        sqrt(err[0] ** 2 + err[1] ** 2))
    observed = np.asarray(emag.dat.data_ro > 0.0)

    obs_mask = Function(Q, name="obs_mask").assign(1.0)
    if os.environ.get("ISMIP7_OBS_MASK", "1") != "0":
        obs_mask.dat.data[:] = observed.astype(float)
        n_no = COMM_WORLD.allreduce(int((obs_mask.dat.data_ro == 0.0).sum()))
        n_all = COMM_WORLD.allreduce(len(obs_mask.dat.data_ro))
        PETSc.Sys.Print(
            f"  Obs mask: {n_no}/{n_all} nodes without velocity obs "
            f"excluded from the misfit (regularization fills them)"
        )

    # ── Observation-error normalization ──────────────────────────────────
    # Each misfit term is divided by the SQUARED observational error of its own
    # datum, so every term is dimensionless (a chi^2 density) and the relative
    # weight between them is a pure number rather than an accident of units.
    # Without this the velocity term carries (m/yr)^2 and cannot be traded off
    # against a dH/dt term in any principled way.
    #
    # Velocity: per-component MEaSUREs formal error, floored. The floor matters
    # -- ERR goes to ~0 in places, and an unfloored 1/sigma^2 would let a
    # handful of nodes dominate the whole functional.
    #
    # NOTE: normalizing changes the ABSOLUTE misfit scale, so a GAMMA tuned
    # against the other scale is NOT transferable. GAMMA_THETA/GAMMA_PHI
    # therefore DEFAULT differently per norm (see GAMMA_DEFAULT above);
    # ISMIP7_MISFIT_NORM=none recovers the legacy functional AND its 1e4
    # gamma. Measured at 32 km the MEaSUREs per-component sigma has median
    # 2.6 m/yr and p90 6.6 m/yr, so the shift is roughly one order of
    # magnitude, not two; the printout below reports the actual distribution.
    sigma_u_floor = float(os.environ.get("ISMIP7_SIGMA_U_FLOOR", "1.0"))  # m/yr
    if MISFIT_NORM == "sigma":
        # Where MEaSUREs reports no velocity, icepack.interpolate zero-fills
        # BOTH u_obs and ERR. Flooring such a node at sigma_u_floor would give
        # it the SMALLEST sigma in the field, i.e. the LARGEST 1/sigma^2
        # weight, on a fabricated zero velocity -- unobserved nodes would then
        # dominate the functional wherever the model flows fast. With the obs
        # mask on they are excluded anyway; the mask-off escape hatch has to be
        # made safe explicitly, by giving those nodes a sigma large enough that
        # they carry no weight rather than maximal weight.
        sigma_u_unobs = float(os.environ.get("ISMIP7_SIGMA_U_UNOBS", "1e4"))
        sig_ux = Function(Q, name="sigma_ux").interpolate(
            max_value(abs(err[0]), Constant(sigma_u_floor)))
        sig_uy = Function(Q, name="sigma_uy").interpolate(
            max_value(abs(err[1]), Constant(sigma_u_floor)))
        _unobs = ~observed
        sig_ux.dat.data[_unobs] = sigma_u_unobs
        sig_uy.dat.data[_unobs] = sigma_u_unobs
        _n_unobs = COMM_WORLD.allreduce(int(_unobs.sum()))
        # Order statistics have to be taken over the WHOLE field: .dat.data_ro
        # is this rank's owned dofs, so a bare np.median would print rank 0's
        # partition as if it were the ice-sheet-wide MEaSUREs distribution --
        # and this printout is the evidence cited for the gamma rescale.
        _parts = COMM_WORLD.gather(
            np.asarray(sig_ux.dat.data_ro, dtype="f8"), root=0)
        if COMM_WORLD.rank == 0:
            _all = np.concatenate(_parts)
            _stats = (float(np.median(_all)), float(np.percentile(_all, 90)))
        else:
            _stats = None
        _sig_med, _sig_p90 = COMM_WORLD.bcast(_stats, root=0)
        PETSc.Sys.Print(
            f"  Misfit normalization: chi^2 (per-datum sigma^2). "
            f"velocity sigma floor {sigma_u_floor:g} m/yr, "
            f"sigma_x median {_sig_med:.2f} p90 {_sig_p90:.2f} m/yr"
        )
        PETSc.Sys.Print(
            f"  Unobserved nodes (ERR==0): {_n_unobs} given "
            f"sigma={sigma_u_unobs:g} m/yr so the zero-filled velocity there "
            f"carries no weight"
        )
        PETSc.Sys.Print(
            f"  NOTE: absolute misfit is dimensionless; GAMMA_THETA="
            f"{GAMMA_THETA:g} GAMMA_PHI={GAMMA_PHI:g} default to the "
            f"sigma-normalized scale (see ISMIP7_MISFIT_NORM)."
        )
    else:
        sig_ux = Function(Q).assign(1.0)
        sig_uy = Function(Q).assign(1.0)
        PETSc.Sys.Print("  Misfit normalization: NONE (legacy, dimensional)")

    # ── Rheology ──
    # Composite viscous rheology: flow exponent n_flow (this branch: n=3
    # standard Glen) main term + linear (n=1) regularization. Sliding stays
    # at Weertman m_slide=3.
    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    n_flow_val = _n_flow()
    m_slide_val = float(os.environ.get("ISMIP7_M_SLIDE", "3.0"))
    a4_factor = float(os.environ.get("ISMIP7_A4_FACTOR", a4_factor_default()))
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
    # global_mean, not .dat.data_ro.mean(): the latter is the OWNED slice, so
    # under MPI each rank built the sliding coefficient from a different
    # reference speed and the same physical location got a different friction
    # depending on which rank owned it. Reaches only the legacy action path
    # (K_base/K_lin -> _rheo_glen/_rheo_linear); build_rc_residual anchors on
    # C_w0 and never sees u_c. See icepack2_tools/mpi_stats.
    u_c = Constant(global_mean(u_speed))
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
        _ws_t_lo, _ws_t_hi = global_range(theta)
        _ws_p_lo, _ws_p_hi = global_range(phi)
        PETSc.Sys.Print(f"    theta: [{_ws_t_lo:.3f}, {_ws_t_hi:.3f}]")
        PETSc.Sys.Print(f"    phi:   [{_ws_p_lo:.3f}, {_ws_p_hi:.3f}]")

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
    phi_eff = Function(Q_g).interpolate(
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
    # Under DG0 geometry the anchor's |grad s| comes from a CG1 reconstruction
    # inside weertman_anchor (a cell-wise surface has no cell gradient); the
    # anchor is a fixed reference scaling, not a force in the residual.
    C_w0 = weertman_anchor(H, s, u_obs, m_slide_val, Q_g)
    if USE_RESIDUAL:
        law_name = ("Budd N_hat (exact-zero shelf, delta="
                    f"{BUDD_DELTA:.3f}, alpha_gl={ALPHA_GL:.2f})"
                    if FRICTION == "budd"
                    else f"regularized Coulomb (c0={C0_RC})")
        C_w0_lo, C_w0_hi = global_range(C_w0)
        PETSc.Sys.Print(
            f"  Friction: {law_name}; h_visc_floor={RC_HVISC_FLOOR:.0f}m; "
            f"C_w0 in [{C_w0_lo:.2e}, {C_w0_hi:.2e}]"
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
        # The thermal prior is a smooth englacial calculation and it is the
        # prior MEAN of the CG1 control phi = log(A/A_prior), so it runs on CG1
        # geometry whatever the prognostic geometry space is. Two reasons not
        # to run it on DG0 and convert afterwards:
        #   * an L2 DG0->CG1 projection overshoots at the front and produced a
        #     NEGATIVE fluidity (A_prior min -9.88 against a DG0 range of
        #     [1.0, 446.7]), which is unphysical and poisons log(A/A_prior);
        #   * the Picard loop stalled on cell-wise inputs (dA plateaued at
        #     4.6e-2, never reaching rtol=1e-3).
        # cg1_lift is a convex combination of cell values, so it cannot
        # overshoot and keeps a positive field positive.
        if geom_dg:
            H_th = cg1_lift(H)
            b_th = cg1_lift(b)
            C_th = cg1_lift(C_w0)
            s_th = Function(Q).interpolate(
                max_value(b_th + H_th, (Constant(1.0) - rho_ratio) * H_th)
            )
        else:
            H_th, b_th, s_th, C_th = H, b, s, C_w0
        A_prior = compute_fluidity_prior(
            u_obs, H_th, s_th, b_th, C_th, acc_prior, T_srf
        )
        A_prior.rename("fluidity_prior")
        A_prior_lo, A_prior_hi = global_range(A_prior)
        PETSc.Sys.Print(
            f"  Fluidity prior A_prior in [{A_prior_lo:.2f}, "
            f"{A_prior_hi:.2f}] (thermomechanical)"
        )
    else:
        A_prior = Function(Q, name="fluidity_prior").interpolate(A0 * Constant(a4_factor))
        PETSc.Sys.Print("  Fluidity prior: constant A0*a4_factor (legacy)")
    A4_base = A_prior
    # The MAP name carries BOTH the flow exponent and the geometry space it was
    # inverted under, so n=3/n=4 and DG0/CG1 MAPs coexist on disk and a forward
    # cannot silently pair itself with a MAP whose front treatment differs.
    # Built by the shared helper the forward and the preflight gates use.
    # ISMIP7_MAP_OUT overrides the full output path: without it, EVERY run --
    # including a short smoke test -- writes to the production filename, and a
    # 3-iteration artifact silently replaces a converged MAP (this nearly
    # happened twice in Aug 2026 validation). Variant MAPs (e.g. the transient
    # dH/dt-constrained inversion) should also name themselves distinctly here
    # rather than shadow the velocity-only MAP the forwards auto-load.
    map_out = os.environ.get("ISMIP7_MAP_OUT")
    map_fn = os.path.basename(map_out) if map_out else map_basename(FRICTION, lc)
    # A bare filename (ISMIP7_MAP_OUT=map.h5) has no dirname; resolve it under
    # MESH_DIR like the non-override path rather than silently against the CWD.
    _map_dir = (os.path.dirname(map_out) or MESH_DIR) if map_out else MESH_DIR
    if map_out:
        PETSc.Sys.Print(f"  MAP output override: {map_out}")
    # Validate the destination NOW. The first write is the iteration-20
    # checkpoint, or for a short run the very end, so a typo'd or absent
    # directory would otherwise throw away hours of a multi-rank inversion --
    # the exact opposite of what this knob exists to protect against.
    if COMM_WORLD.rank == 0:
        _probe_err = None
        try:
            os.makedirs(_map_dir, exist_ok=True)
            _probe = os.path.join(_map_dir, f".map_write_probe.{os.getpid()}")
            with open(_probe, "w") as _fh:
                _fh.write("ok")
            os.remove(_probe)
        except OSError as exc:
            _probe_err = f"MAP output directory {_map_dir!r} is unusable: {exc}"
    else:
        _probe_err = None
    _probe_err = COMM_WORLD.bcast(_probe_err, root=0)
    if _probe_err:
        raise RuntimeError(_probe_err)
    PETSc.Sys.Print(f"  MAP output: {os.path.join(_map_dir, map_fn)}")

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
    PETSc.Sys.Print(f"  u_max = {global_max(u_mag):.0f} m/yr")

    # ── Forward function for tlm_adjoint ──
    # Normalize by the OBSERVED area so the misfit magnitude stays
    # comparable between masked and unmasked runs.
    area_val = assemble(obs_mask * dx(mesh))

    # ── Transient (dH/dt) constraint ─────────────────────────────────────
    # A velocity-only inversion never constrains div(h u), so the MAP can carry
    # a flux divergence wildly inconsistent with the observed geometry; the
    # forward then drifts or needs a large frozen apparent-MB term to stand in
    # for it. Adding one prognostic step and matching the resulting thickness
    # tendency to an observed mean dH/dt map closes that directly.
    #
    # Restricted to GROUNDED ice by default, for three reasons: floating dH/dt
    # from altimetry is noisy and firn/tide/ocean confounded (the observed field
    # integrates to +212 Gt/yr over shelves against -85.6 Gt/yr grounded, which
    # matches IMBIE-3), floating thickness change does not move VAF, and on
    # grounded ice ocean melt is identically zero.
    #
    # That last point is about the MISFIT, not the source term: the step still
    # applies the shelf melt (below) so the single prognostic solve matches the
    # forcing the forward experiment uses, but because the DG0 upwind operator
    # couples a cell only to its UPWIND neighbours and the flow runs
    # grounded->shelf, no grounded-cell equation ever sees a shelf value. The
    # melt source therefore cannot move the grounded-only misfit -- which is
    # why a missing per-basin K degrades to an SMB-only source rather than
    # aborting the inversion.
    dhdt_w = float(os.environ.get("ISMIP7_DHDT_WEIGHT", "0.0"))
    use_dhdt = dhdt_w > 0.0
    # Resolved (not merely requested) net-term provenance: ISMIP7_DHDT_NET_SIGMA
    # only reaches the objective when the dH/dt term itself is on, so these stay
    # at their off values unless the term is actually built below. save_map
    # stamps net_sigma_used, so a MAP can never claim a constraint it never saw.
    use_dhdt_net = False
    net_sigma_used = 0.0
    if use_dhdt:
        if not geom_dg:
            raise RuntimeError(
                "ISMIP7_DHDT_WEIGHT requires ISMIP7_GEOMETRY_SPACE=dg0: the "
                "prognostic step is the DG0 upwind FV operator and must act on "
                "the same field the forward transports."
            )
        from icepack2_tools.obs_dhdt import load_dhdt_obs

        dhdt_sigma = float(os.environ.get("ISMIP7_DHDT_SIGMA", "0.1"))   # m/yr
        dhdt_dt = float(os.environ.get("ISMIP7_DHDT_DT", "1.0"))          # yr
        clim0 = int(os.environ.get("ISMIP7_DHDT_CLIM_START", "2003"))
        clim1 = int(os.environ.get("ISMIP7_DHDT_CLIM_END", "2019"))

        dhdt_obs, dhdt_cov = load_dhdt_obs(Q_g)
        smb_obs = load_racmo_smb_climatology(Q_g, clim_start=clim0, clim_end=clim1)

        # Ocean melt for the prognostic step, from the SAME parameterisation
        # and forcing the experiments use: OI-climatology TF/so at draft depth
        # + per-basin calibrated K through the Burgard quadratic-mixed-slope
        # formula, via the shared make_climatology_ocean_callback -- not a
        # reimplementation. Evaluated ONCE at the frozen reference geometry
        # (the controls move; the reference geometry does not), so it is a
        # constant field like smb_obs. It is zero on grounded ice by
        # construction, so the grounded-only misfit is unchanged by it; what
        # it fixes is the source term the implicit step sees on shelf cells,
        # keeping the single prognostic step consistent with the forward
        # experiment. ISMIP7_DHDT_MELT=0 drops it (SMB-only source).
        melt_ref = Function(Q_g, name="ocean_melt_ref")
        if os.environ.get("ISMIP7_DHDT_MELT", "1") != "0":
            from icepack2_tools.forcing import (
                load_K_per_basin, make_climatology_ocean_callback)
            _k_lc = os.path.join(_ROOT, "results",
                                 f"calibrated_K_per_basin_{lc}.npz")
            _k_2500 = os.path.join(_ROOT, "results",
                                   "calibrated_K_per_basin_2500.npz")
            k_npz = os.environ.get(
                "ISMIP7_K_PER_BASIN_NPZ",
                _k_lc if os.path.exists(_k_lc) else _k_2500)
            if not os.path.exists(k_npz):
                # Warn, do not abort: the melt source only touches shelf cells
                # and the misfit is grounded-only, so an absent ocean
                # calibration cannot change this objective. Making it a hard
                # prerequisite would fail the whole inversion over an artifact
                # the term provably does not use.
                PETSc.Sys.Print(
                    f"  [!] dH/dt melt source: per-basin K not found at "
                    f"{k_npz}; falling back to an SMB-only source. The "
                    f"grounded-only misfit is unaffected (melt is zero on "
                    f"grounded ice and the DG0 upwind step never carries a "
                    f"shelf value upstream); shelf cells of the single "
                    f"prognostic step are no longer forcing-consistent with "
                    f"the forward. Calibrate K, or set ISMIP7_DHDT_MELT=0 to "
                    f"silence this."
                )
            else:
                _W2 = VectorFunctionSpace(mesh, "DG", 0)
                _xy = Function(_W2).interpolate(
                    fd.SpatialCoordinate(mesh)).dat.data_ro.reshape(-1, 2)
                geom_xy = (_xy[:, 0].copy(), _xy[:, 1].copy())
                K_field = load_K_per_basin(
                    k_npz, geom_xy[0], geom_xy[1], fill=0.0)
                K_field = K_field * float(
                    os.environ.get("ISMIP7_K_SCALE", "1.0"))
                _ctx = {"mesh": mesh, "Q": Q, "V": V, "Q_g": Q_g,
                        "geom_xy": geom_xy, "h": H, "b": b, "s": s,
                        "ocean_melt": melt_ref}
                make_climatology_ocean_callback(K_field)(_ctx, 0.0)
                _melt_gt = float(assemble(melt_ref * dx)) * 917.0 / 1e12
                PETSc.Sys.Print(
                    f"  dH/dt melt source: per-basin K "
                    f"({os.path.basename(k_npz)}), integrated "
                    f"{_melt_gt:.0f} Gt/yr at reference geometry")
        else:
            PETSc.Sys.Print("  dH/dt melt source: DISABLED (SMB-only)")

        # grounded indicator on the reference geometry, frozen (the control
        # fields move, the classification of where we trust dH/dt should not)
        haf_ref = Function(Q_g).interpolate(height_above_flotation(H, b))
        grounded = Function(Q_g)
        grounded.dat.data[:] = np.where(
            (haf_ref.dat.data_ro > 0.0) & (H.dat.data_ro > 1.0), 1.0, 0.0)
        dhdt_mask = Function(Q_g, name="dhdt_mask")
        dhdt_mask.dat.data[:] = dhdt_cov.dat.data_ro * grounded.dat.data_ro
        dhdt_area = assemble(dhdt_mask * dx(mesh))
        if dhdt_area <= 0.0:
            raise RuntimeError(
                "dH/dt constraint enabled but no cell is both grounded and "
                "observation-covered; check ISMIP7_OBS_KIT and the mesh."
            )
        # Global, not rank-local: this sits on the same line as dhdt_area,
        # which is a global assemble, and a rank-local count beside it once
        # read as the constraint having dropped 96% of its cells.
        _n_use = COMM_WORLD.allreduce(int((dhdt_mask.dat.data_ro > 0.5).sum()))
        PETSc.Sys.Print(
            f"  dH/dt constraint: weight={dhdt_w:g} sigma={dhdt_sigma:g} m/yr "
            f"dt={dhdt_dt:g} yr, {_n_use} grounded+observed cells, "
            f"{dhdt_area / 1e6 / 1e3:.0f} x10^3 km^2, "
            f"SMB clim {clim0}-{clim1}"
        )
        if variable_note := os.environ.get("ISMIP7_DHDT_VAR", "dhdt_smith"):
            if variable_note == "dhdt_smith":
                PETSc.Sys.Print(
                    "  dH/dt target: dhdt_smith 2003-2019 mean -- DECLARED "
                    "overlap with the post-2015 projection era (protocol "
                    "dates DA initial states before 2015; see obs_dhdt.py)"
                )
        h_next = Function(Q_g, name="h_next")

        # ── Net (integral) mass-balance term -- OFF BY DEFAULT ──────────
        # Explored Aug 2026 and mechanically validated (it moves the grounded
        # net tendency from +604 to -54 Gt/yr against an observed -73 at
        # 2500 m, at ~10% pointwise-RMS cost), then EXCLUDED from the
        # production objective by decision (Andrew, Aug 30):
        #   1. Assimilating the integrated mass trend forfeits it as
        #      INDEPENDENT validation -- the observed trend is the one number
        #      every downstream assessment (IMBIE/GRACE consistency) checks.
        #   2. proj - CTRL differencing under ISMIP7_APPARENT_MB removes the
        #      shared frozen drift by construction, so the bias this term
        #      fixes largely cancels in the reported signal anyway.
        #   3. The integral is bought with regionally structured adjustments
        #      to a residual that is div(h u) DATA NOISE (interior response
        #      times are 1e3-1e4 yr; a 16-yr window carries no dynamical
        #      interior signal), i.e. aggregated noise-fitting.
        # The machinery stays for experiments (set ISMIP7_DHDT_NET_SIGMA>0),
        # and the per-iteration net= DIAGNOSTIC prints regardless, so the
        # tendency bias is always visible without being penalised.
        dhdt_net_sigma = float(os.environ.get("ISMIP7_DHDT_NET_SIGMA", "0"))
        use_dhdt_net = dhdt_net_sigma > 0.0
        if use_dhdt_net:
            net_sigma_used = dhdt_net_sigma
            R_net = FunctionSpace(mesh, "R", 0)
            A_tot_net = float(assemble(Constant(1.0) * dx(mesh)))
            _rho_gt = 917.0 / 1e12
            PETSc.Sys.Print(
                f"  dH/dt NET constraint: sigma_net={dhdt_net_sigma:g} Gt/yr "
                f"(0.5*(net/sigma)^2 on the grounded+observed integral; "
                f"ISMIP7_DHDT_NET_SIGMA=0 disables)"
            )

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
        # chi^2 density: each residual divided by the squared error of its own
        # observation, so the term is dimensionless.
        integrand = (
            0.5
            / area_val
            * obs_mask
            * (
                (u_sol[0] - u_obs[0]) ** 2 / sig_ux ** 2
                + (u_sol[1] - u_obs[1]) ** 2 / sig_uy ** 2
            )
        )

        if use_dhdt:
            # ONE prognostic step, using the SAME DG0 upwind FV operator the
            # forward transports with -- that consistency is the point: the
            # inversion must be penalised for the divergence its own transport
            # scheme will produce, not an idealised one.
            w = TestFunction(Q_g)
            nrm = FacetNormal(mesh)
            un = dot(u_sol, nrm)
            unp = (un + abs(un)) / 2
            dt_c = Constant(dhdt_dt)
            F_h = (
                (h_next - H) / dt_c * w * dx
                + (unp("+") * h_next("+") - unp("-") * h_next("-")) * jump(w) * dS
                + unp * h_next * w * ds
                - (smb_obs - melt_ref) * w * dx
            )
            EquationSolver(
                F_h == 0,
                h_next,
                solver_parameters={
                    "snes_type": "ksponly",
                    "ksp_type": "gmres",
                    "pc_type": "bjacobi",
                    "sub_pc_type": "ilu",
                    "ksp_rtol": 1e-10,
                },
                form_compiler_parameters=fc_params,
            ).solve()
            dhdt_model = (h_next - H) / dt_c
            integrand = integrand + (
                0.5
                * Constant(dhdt_w)
                / dhdt_area
                * dhdt_mask
                * ((dhdt_model - dhdt_obs) / Constant(dhdt_sigma)) ** 2
            )

            if use_dhdt_net:
                # m_net = domain mean of the masked residual (R-space
                # projection; the test function is the constant 1, so the
                # 1x1 solve is exactly mean = integral/area).
                m_net = Function(R_net, name="dhdt_net_mean")
                r_net = TestFunction(R_net)
                EquationSolver(
                    (m_net - (dhdt_model - dhdt_obs) * dhdt_mask) * r_net * dx
                    == 0,
                    m_net,
                    solver_parameters={
                        "snes_type": "ksponly",
                        "ksp_type": "preonly",
                        "pc_type": "jacobi",
                    },
                    form_compiler_parameters=fc_params,
                    # The R-space (Real) 1x1 matrix assembles as a
                    # python-type PETSc Mat, which tlm_adjoint's linear-solver
                    # cache cannot copy (MatAXPY: "no method getrow for Mat of
                    # type python"). Caching a 1x1 solve buys nothing anyway.
                    cache_jacobian=False,
                    cache_adjoint_jacobian=False,
                ).solve()
                # net [Gt/yr] = mean * total area * rho; integrand constant
                # over the mesh, so integrating /A_tot recovers the scalar.
                _c_net = Constant(A_tot_net * _rho_gt / dhdt_net_sigma)
                integrand = integrand + (
                    Constant(0.5 / A_tot_net) * (m_net * _c_net) ** 2
                )

        J = Functional(name="J")
        J.assign(integrand * dx)
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

    # ── Per-term diagnostics ─────────────────────────────────────────────
    # forward() returns ONE functional, so with two terms summed the reported
    # misfit alone cannot say which is being fitted. These re-assemble each
    # chi^2 contribution from the state left behind by the last forward solve
    # (diagnostic only -- outside the annotated tape, so they do not enter the
    # gradient). Both are dimensionless and directly comparable.
    _u_now = z.subfunctions[0]
    _vel_chi2 = (
        0.5 / area_val * obs_mask
        * ((_u_now[0] - u_obs[0]) ** 2 / sig_ux ** 2
           + (_u_now[1] - u_obs[1]) ** 2 / sig_uy ** 2)
    ) * dx(mesh)
    if use_dhdt:
        _net_form = (((h_next - H) / Constant(dhdt_dt) - dhdt_obs)
                     * dhdt_mask * dx(mesh))
        _dhdt_chi2 = (
            0.5 / dhdt_area * dhdt_mask
            * (((h_next - H) / Constant(dhdt_dt) - dhdt_obs)
               / Constant(dhdt_sigma)) ** 2
        ) * dx(mesh)

    # Velocity chi^2 of the last accepted state, in ITS OWN accumulator. The
    # optimizer's J_val is the COMBINED objective under use_dhdt, so it cannot
    # serve as the reference for a velocity-only comparison: a large
    # ISMIP7_DHDT_WEIGHT inflates it arbitrarily and would slacken the
    # final-save guard below to the point of passing a corrupt velocity.
    last_good_vel_chi2 = [np.inf]

    def term_report():
        r"""``' vel=... dhdt=...'`` for the iteration line, or '' if disabled."""
        try:
            out = f" vel={last_good_vel_chi2[0]:.4e}"
            if use_dhdt:
                out += f" dhdt={float(assemble(_dhdt_chi2)):.4e}"
                out += (f" net={float(assemble(_net_form)) * 917.0 / 1e12:+.0f}Gt/yr")
            return out
        except Exception:
            return ""

    # ── MAP checkpoint writer ────────────────────────────────────────────
    # ONE writer for both the periodic checkpoint and the final save: the two
    # used to be copy-pasted, so a field or attribute added to one silently
    # missed the other.
    #
    # The attributes record which OBJECTIVE produced the controls. The filename
    # (icepack2_tools/naming.py) encodes only friction, lc, geometry space and
    # flow exponent, so a velocity-only MAP and a dH/dt-constrained one, or two
    # runs under different normalizations/priors, land on the SAME path and
    # overwrite each other with nothing on disk to tell them apart. This repo
    # has been bitten by exactly that class of look-alike MAP before (see
    # ../GEOMETRY_DISCRETIZATION.md on the geometry tag), and MAPs are
    # gitignored, so the checkpoint is the only place this provenance can live.
    def save_map(path):
        with fd.CheckpointFile(path, "w") as chk:
            chk.save_mesh(mesh)
            chk.save_function(theta, name="log_friction")
            chk.save_function(phi, name="log_fluidity")
            chk.save_function(u_obs, name="velocity_obs")
            chk.save_function(obs_mask, name="obs_mask")
            chk.save_function(H, name="thickness")
            chk.save_function(b, name="bed")
            chk.save_function(s, name="surface")
            chk.save_function(A_prior, name="fluidity_prior")
            # The .msh this MAP was inverted on. A CheckpointFile mesh is named
            # "firedrake_default", so this is how the forward names its own
            # mesh and picks the matching per-mesh boundary-id sidecar.
            chk.set_attr("/", "mesh_basename", os.path.basename(mesh_fn))
            # Mesh PARAMETERS as well as the basename (Dan/David's scheme,
            # merged from upstream/integration). The forward resolves its
            # boundary_ids sidecar from these rather than from its own
            # environment: ISMIP7_BUFFER_M / ISMIP7_LC_COARSE can drift, and a
            # mismatched sidecar puts the calving BC on the wrong facets, where
            # ds(absent_id) integrates to zero -- silently wrong physics with
            # no crash. The basename covers meshes outside the standard naming
            # pattern; the parameters let bndids_filename() rebuild the name.
            chk.set_attr("/", "lc", int(lc))
            chk.set_attr("/", "lc_coarse", int(lc_coarse))
            chk.set_attr("/", "buffer_m", float(buffer_m))
            chk.set_attr("/", "misfit_norm", MISFIT_NORM)
            chk.set_attr("/", "gamma_theta", float(GAMMA_THETA))
            chk.set_attr("/", "gamma_phi", float(GAMMA_PHI))
            chk.set_attr("/", "dhdt_weight", float(dhdt_w))
            chk.set_attr("/", "dhdt_net_sigma", net_sigma_used)

    # ── L-BFGS-B Inversion ──
    max_iter = int(os.environ.get("ISMIP7_MAXITER", "500"))
    PETSc.Sys.Print("\nStarting L-BFGS-B inversion (theta + phi)...")
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
        last_good_vel_chi2[0] = float(assemble(_vel_chi2))

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
            f"misfit={J_val:.6e}{term_report()} "
            f"reg_θ={reg_theta:.4e} reg_φ={reg_phi:.4e} "
            f"total={total:.6e} |grad|={np.linalg.norm(total_grad):.4e} "
            f"[fwd={t_fwd:.1f}s adj={t_adj:.1f}s total={t_iter:.1f}s]"
        )

        # Periodic checkpoint every 20 iterations
        if iteration_count[0] % 20 == 0:
            save_map(os.path.join(_map_dir, map_fn))
            PETSc.Sys.Print(f"    [checkpoint saved: iter {iteration_count[0]}]")

        return total, total_grad

    # ── Optimization metric (ISMIP7_GRAD_PRECOND) ────────────────────────
    # L-BFGS-B works in raw dof coordinates, i.e. the Euclidean l2 metric, and
    # that metric is MESH-DEPENDENT: a gradient entry scales with the dof's
    # cell area, so the fine grounding-line cells -- exactly where the controls
    # must move -- get the smallest gradient entries and converge slowest.
    # `mass` optimizes in u = M^(1/2) x (M = lumped mass), which is steepest
    # descent in L2 and makes the rate mesh-independent; peers compensate for
    # the same defect with brute iteration counts (ISSM/M1QN3 runs 1000+300 in
    # cycles). The natural upgrade is the prior metric (delta*M+gamma*K)^-1
    # (fenics_ice); the operator exists in icepack2_tools/prior.py. Default
    # off so in-flight runs stay comparable; flip after the current A/B.
    grad_precond = os.environ.get("ISMIP7_GRAD_PRECOND", "none").lower()
    if grad_precond not in ("none", "mass"):
        raise ValueError(f"ISMIP7_GRAD_PRECOND must be none|mass, got {grad_precond!r}")
    if grad_precond == "mass":
        _mv = assemble(TestFunction(Q) * dx).dat.data_ro
        _mloc = Function(Q); _mloc.dat.data[:] = _mv
        _sqrtm_one = np.sqrt(np.maximum(func_to_global(_mloc), 1e-30))
        _sqrtm = np.concatenate([_sqrtm_one, _sqrtm_one])
        PETSc.Sys.Print(
            f"  Optimization metric: lumped-mass Riesz (u = sqrt(M) x); "
            f"sqrt(m) range [{_sqrtm_one.min():.2e}, {_sqrtm_one.max():.2e}]"
        )
        _inner_og = objective_and_gradient
        def objective_and_gradient(u_vec):  # noqa: F811 (deliberate wrap)
            J, g_x = _inner_og(u_vec / _sqrtm)
            return J, g_x / _sqrtm
    x0 = np.concatenate([func_to_global(theta), func_to_global(phi)])
    if grad_precond == "mass":
        x0 = x0 * _sqrtm
    result = scipy_minimize(
        objective_and_gradient,
        x0,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": max_iter, "ftol": 0, "gtol": 0, "disp": False},
    )

    PETSc.Sys.Print(f"\nOptimization finished: {result.message}")
    PETSc.Sys.Print(f"  {result.nit} iterations, {result.nfev} function evaluations")
    _x_final = result.x / _sqrtm if grad_precond == "mass" else result.x
    global_to_func(_x_final[:global_ndof], theta)
    global_to_func(_x_final[global_ndof:], phi)
    _t_lo, _t_hi = global_range(theta)
    _p_lo, _p_hi = global_range(phi)
    PETSc.Sys.Print(f"  theta range: [{_t_lo:.3f}, {_t_hi:.3f}]")
    PETSc.Sys.Print(f"  phi range:   [{_p_lo:.3f}, {_p_hi:.3f}]")

    # ── Save MAP immediately ──
    chk_fn = os.path.join(_map_dir, map_fn)
    save_map(chk_fn)
    PETSc.Sys.Print(
        f"Saved MAP: {chk_fn} "
        f"(misfit_norm={MISFIT_NORM} gamma_theta={GAMMA_THETA:g} "
        f"gamma_phi={GAMMA_PHI:g} dhdt_weight={dhdt_w:g})"
    )

    # ── Final forward solve ──
    # The single-shot solve at full exponents can fail, and the old code then
    # SAVED the failed Newton state as "velocity" while claiming the last
    # optimization state was used (the message said so, but z was never
    # restored). The Aug 3 dg0 32 km MAP carries a 2745 m/yr-RMS velocity this
    # way -- discovered only when compare_dhdt.py scored it against MEaSUREs.
    # Now: restore the last good optimization state on failure, and refuse to
    # save any velocity whose misfit grossly disagrees with the optimizer's
    # (a forward re-solves the diagnostic from theta/phi anyway; a missing
    # velocity is an inconvenience, a silently wrong one poisons everything
    # downstream that trusts the checkpoint).
    PETSc.Sys.Print("\nFinal forward solve...")
    stop_manager()
    try:
        slvr.solve()
    except fd.ConvergenceError:
        PETSc.Sys.Print("  Final solve failed; restoring last good "
                        "optimization state")
        z.assign(z_backup)

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

    # Update checkpoint with velocity -- unless the state is inconsistent
    # with the optimization it claims to represent. BOTH sides of the test are
    # the SAME form, `_vel_chi2`, assembled on the final state and on the last
    # accepted optimization state: an earlier version compared the dimensional
    # misfit against a chi^2 and blocked a good velocity at 13.8x (a corrupt
    # one sits >1000x off). Comparing against last_good_obj would repeat that
    # error in the other direction, since under use_dhdt that is the combined
    # objective and its dH/dt part has nothing to do with velocity.
    _guard = float(assemble(_vel_chi2))
    _ref = max(float(last_good_vel_chi2[0]), 1e-30)
    if np.isfinite(_guard) and _guard <= 10.0 * _ref:
        with fd.CheckpointFile(chk_fn, "a") as chk:
            chk.save_function(u_sol, name="velocity")
        PETSc.Sys.Print(f"Saved velocity: {chk_fn}")
    else:
        PETSc.Sys.Print(
            f"WARNING: NOT saving velocity -- final-state velocity misfit "
            f"{_guard:.3e} is inconsistent with the last accepted "
            f"{_ref:.3e} (>10x, same metric). The MAP controls are saved and "
            f"valid; forwards re-solve the diagnostic from theta/phi and are "
            f"unaffected."
        )

    # ── Plot ──
    # Optional: the MAP is already written and the velocity saved above, so a
    # missing plotting dependency must not fail the run at this point.
    try:
        import colorcet as cc
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        PETSc.Sys.Print(f"Skipping summary figure ({exc}); MAP is saved.")
        return

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
