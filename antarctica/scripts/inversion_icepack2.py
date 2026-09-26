#!/usr/bin/env python3
"""
Joint inversion for basal friction and fluidity using icepack2 + tlm_adjoint.

Follows the Kangerd demo pattern from Shapero's dual-problems repo:
- 3-field Z = V * Sigma * T (no DG thickness in mixed space)
- Regularized Jacobian: J = J_r + α * J_1
- snes_divergence_tolerance = -3 (PETSC_UNLIMITED; the Kangerd demo's -1 is
  PETSC_DETERMINE, which reinstates the default 1e4 growth cutoff)
- Sliding coefficient includes exp(m*theta)

Controls are log-deviations from PHYSICAL prior means - theta = log(C/C_w0)
on the balance-friction anchor, phi = log(A/A_prior) on a thermomechanical
fluidity prior - regularized with the Whittle-Matern prior of Recinos et al. (2023)
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
    ln,
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
from types import SimpleNamespace

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
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")

# Repo root on the path so we can import the shared dual-friction operator.
sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.dual_friction import (
    build_rc_residual,
    effective_pressure,
    weertman_anchor,
)
from icepack2_tools.geometry import cg1_lift, sample_to_geometry
from icepack2_tools.transfer import interpolate_with_fill, meshes_match
from icepack2_tools.grounding import height_above_flotation
from icepack2_tools.mpi_stats import (global_mean, global_range,
                                      global_max, global_size, global_count)
from icepack2_tools.naming import map_basename
from icepack2_tools.runconfig import (
    deltat_per_basin_npz, k_per_basin_npz,
    obs_data_root,
    BUDD_SHELF_GATE,
    friction as _friction, geometry_space as _geometry_space,
    raster_sample as _raster_sample,
    lc as _lc, lc_coarse as _lc_coarse, n_flow as _n_flow,
)
DATA_DIR = obs_data_root()
from icepack2_tools.prior import (
    bilaplacian_aux_residual,
    bilaplacian_coeffs,
    bilaplacian_energy_form,
    prior_operator_coeffs,
    prior_operator_form,
)
from icepack2_tools.thermo_model import compute_fluidity_prior
from icepack2_tools.optimization import (FunctionalDecreaseStop,
                                         recorded_objective,
                                         resolve_log_vel_weight)
from icepack2_tools.forcing import (load_racmo_smb_climatology,
                                    load_mean_annual_surface_temperature)
from icepack2_tools.runconfig import (
    TARGET_MESH_GEOMETRY_METHOD,
    front_hmin,
    residual_stabilizers,
)
from icepack2_tools.front import facet_neighbours, ocean_drag_cells
from icepack2_tools.continuation import ladder, ramp_exponents
from icepack2_tools.solverconfig import (
    continuation_steps,
    diagnostic_solver_label,
    diagnostic_solver_mode,
    diagnostic_solver_parameters,
    final_solve_bounds,
    final_solve_parameters,
    linearization_state,
    nonlinear_solver_options,
    snes_atol_scale,
    snes_monitor_enabled,
)
from mesh_naming import get_buffer_m, mesh_filename
from timing_campaign import MATRIX_T_START, atomic_write_json

# petsc4py returns SNES converged reasons as plain ints on most builds.
_SNES_REASON_NAMES = {
    value: name
    for name, value in vars(PETSc.SNES.ConvergedReason).items()
    if isinstance(value, int) and not name.startswith("_")
}


def _snes_reason_name(reason):
    name = getattr(reason, "name", None)
    if name:
        return str(name)
    return _SNES_REASON_NAMES.get(int(reason), str(int(reason)))

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

# Regularization -- Whittle-Matern prior (icepack2_tools/prior.py)
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

# ── ISSM-convention logarithmic velocity misfit ──────────────────────────
# An absolute (or sigma-normalized) velocity misfit is dominated by wherever
# the observational error is smallest, which on MEaSUREs is the slow interior
# (median sigma 2.6 m/yr): a 10 m/yr error on a 20 m/yr datum costs more than
# a 200 m/yr error on a 600 m/yr tributary. Measured on the 2500 m MAP
# (antarctica/scripts/region_budget.py): the grounding-line flux carried by
# the inverted velocity is 140% of observed below 100 m/yr and 52-62% between
# 100 and 1500 m/yr, i.e. the tributaries that deliver most of the discharge
# are systematically slow, and the ice sheet then gains grounded mass.
#
# ISSM's remedy is standard practice: sum the absolute misfit with a
# LOGARITHMIC one (cost functions 101 + 103, `SurfaceAbsVelMisfit` +
# `SurfaceLogVelMisfit`), the latter being a RELATIVE error and therefore
# scale-free:
#
#     J_log = 1/2 int [ ln( (|u| + eps) / (|u_obs| + eps) ) ]^2 dA / A
#
# with `eps` = ISSM's `epsvel`, a floor that keeps the logarithm finite where
# the ice is not moving. Other inversions use relative errors to the same end.
#
# ISMIP7_LOG_VEL_WEIGHT: 0 (default, the pre-Sep-2026 objective), a number, or
# "auto", which scales the log term to equal the velocity chi^2 term at the
# STARTING state, the only weight that means anything before the first
# iteration. A chained link starts from its predecessor's checkpoint, so under
# "auto" a warm start that records a positive weight under the same misfit
# norm and eps supplies that weight, and every link minimises one objective
# (issue 68, optimization.resolve_log_vel_weight). The resolved value and its
# source are stamped in the MAP.
LOG_VEL_WEIGHT = os.environ.get("ISMIP7_LOG_VEL_WEIGHT", "0")
LOG_VEL_EPS = float(os.environ.get("ISMIP7_LOG_VEL_EPS", "1.0"))   # m/yr
GAMMA_THETA = float(os.environ.get("ISMIP7_GAMMA_THETA", GAMMA_DEFAULT))
GAMMA_PHI = float(os.environ.get("ISMIP7_GAMMA_PHI", GAMMA_DEFAULT))
L_REG = float(os.environ.get("ISMIP7_L_REG", "7.5e3"))

# ── Prior form (ISMIP7_PRIOR_FORM) ──────────────────────────────────────
# `laplacian` (default, everything inverted so far) uses A = delta*M + gamma*K
# as the prior precision itself. `bilaplacian` uses A M^-1 A, which is
# the squared-operator prior of Villa et al. (2021) ("LM^-1L")
# and what the Whittle-Matern SPDE needs in 2-D for the field to be a function
# rather than a distribution (alpha = 2 > d/2; the un-squared operator has no
# pointwise variance to speak of and does not converge under refinement).
#
# They are DIFFERENT PRIORS, not two spellings of one: gamma under `laplacian`
# and (sigma, rho) under `bilaplacian` are not convertible, the regularization
# magnitudes differ by orders of magnitude, and a MAP inverted under one is not
# comparable with a MAP inverted under the other. save_map stamps prior_form
# for that reason. Default stays `laplacian` so in-flight inversions and every
# existing MAP keep their meaning.
PRIOR_FORM = os.environ.get("ISMIP7_PRIOR_FORM", "laplacian").lower()
if PRIOR_FORM not in ("laplacian", "bilaplacian"):
    raise ValueError(
        f"ISMIP7_PRIOR_FORM must be laplacian|bilaplacian, got {PRIOR_FORM!r}"
    )
# Bi-Laplacian knobs are PHYSICAL: the log-deviation scale the control is
# expected to carry, and the distance over which it decorrelates. The defaults
# are the scales prior.py's docstring quotes for the un-squared form (prior std
# ~0.2-0.3 at the L_REG correlation length), so the two forms start from the
# same intent even though their gammas are incomparable.
PRIOR_SIGMA_THETA = float(os.environ.get("ISMIP7_PRIOR_SIGMA_THETA", "0.3"))
PRIOR_SIGMA_PHI = float(os.environ.get("ISMIP7_PRIOR_SIGMA_PHI", "0.3"))
PRIOR_RHO = float(os.environ.get("ISMIP7_PRIOR_RHO", str(L_REG)))

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
# shelf HAF-gated law; the effective-pressure feedback is purely prognostic.
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
    # num_vertices()/num_cells() count this rank's plex, halo included; the
    # coordinate dofs and the owned cell set are reduced to global totals.
    PETSc.Sys.Print(f"  {global_size(mesh.coordinates)} vertices, "
                    f"{mesh.comm.allreduce(mesh.cell_set.size)} cells")

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
    raster_sample = _raster_sample()
    Q_g = FunctionSpace(mesh, "DG", 0) if geom_dg else Q
    PETSc.Sys.Print(f"  Geometry space: {geometry_space.upper()}")

    # ── Load Data ──
    PETSc.Sys.Print("Loading data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    # Cell average onto the geometry space, NOT a centroid point sample --
    # see geometry.sample_to_geometry for the measurements behind that.
    PETSc.Sys.Print(f"  Raster sampling onto geometry cells: {raster_sample}")
    b = sample_to_geometry(
        rasterio.open(f"netcdf:{bm_fn}:bed"), Q_g, Q, method=raster_sample)
    # h_clamp default 0.0: invert against the *true* BedMachine geometry,
    # including h=0 over the buffered ocean region. Composite rheology
    # (added below) keeps the SNES nonsingular where h=0.
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "0.0"))
    H = sample_to_geometry(
        rasterio.open(f"netcdf:{bm_fn}:thickness"),
        Q_g, Q, floor=h_clamp, method=raster_sample)
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

    # Newton/line-search settings come from icepack2_tools.solverconfig so the
    # ISMIP7_SNES_* knobs the campaign exports mean the same thing here as in
    # the transient: newtonls, nleqerr, max_it 200, stol 0, and divergence
    # tolerance -3 (PETSC_UNLIMITED; the legacy -1 is PETSC_DETERMINE, which
    # restores PETSc's 1e4 growth cutoff and reports DIVERGED_DTOL on solves
    # that would otherwise reach their real result -- README "Timing
    # benchmark" §7). The linear solve is the full mixed-Jacobian MUMPS LU:
    # tlm_adjoint differentiates through it, so this is deliberately NOT the
    # transient's condensed scpc_mumps mode, and the MAP records that as
    # state_solver_mode beside the lane contract diagnostic_solver_mode.
    sparams = nonlinear_solver_options()
    sparams.update({
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 400,  # working memory increase
        "mat_mumps_icntl_24": 1,  # detect null pivots
        "mat_mumps_cntl_3": 1e-12,  # null pivot threshold
        # MUMPS prints its error return (INFOG(1), the workspace or pivot
        # code) instead of failing silently: a factorisation that fails
        # reaches SNES only as DIVERGED_LINEAR_SOLVE (job 1612624).
        "mat_mumps_icntl_4": 1,
    })
    state_solver_mode = "full_mumps"
    state_solver_parameters = json.dumps(sparams, sort_keys=True)
    # Optional SNES monitoring (ISMIP7_SNES_MONITOR=1, ISMIP7_SNES_LOG=file),
    # the transient runner's convention. Applies to every annotated forward
    # and, through tlm_adjoint, the adjoint linear solves.
    _solver_log = os.environ.get("ISMIP7_SNES_LOG") if snes_monitor_enabled() else None
    if _solver_log:
        os.makedirs(os.path.dirname(os.path.abspath(_solver_log)), exist_ok=True)
    _viewer = f"ascii:{_solver_log}::append" if _solver_log else None
    _monitor_options = {
        "snes_monitor": _viewer,
        "snes_converged_reason": _viewer,
        # Which way the linear solve went: a failed factorisation
        # (DIVERGED_PC_FAILED) or an inexact one (DIVERGED_ITS).
        "ksp_converged_reason": _viewer,
    } if snes_monitor_enabled() else {}
    sparams.update(_monitor_options)
    # The mode consumers of the published state must run (validate_cache_manifest
    # asserts it). Resolved now so an invalid environment fails here, not
    # inside the final save after hours of work.
    lane_solver_mode = diagnostic_solver_mode()
    fc_params = {"quadrature_degree": 4}

    # ── Build form (Kangerd pattern: controls baked into sliding coefficient) ──
    z = Function(Z)
    z.sub(0).interpolate(Constant(0.1) * u_obs)

    theta = Function(Q, name="theta")  # log friction adjustment
    phi = Function(Q, name="phi")  # log fluidity adjustment

    # Warm-start from a previous MAP or timing-cache checkpoint. Prefer
    # interpolate (not a raw .dat copy): the prepared timing caches are
    # published on 1 rank and the invert runs on many, so dof ownership differs.
    warm_chk = os.environ.get("ISMIP7_WARM_START", "").strip()
    skip_continuation = (
        os.environ.get("ISMIP7_SKIP_CONTINUATION", "0").strip() == "1"
    )
    warm_A_prior = None
    warm_loaded_z = False
    # Residual the warm start's writer reached under the shared F (stamped by
    # save_model_state and by save_map): the forwards' absolute tolerance.
    warm_recorded = None
    # The log-velocity term the warm start was minimised under, which an
    # "auto" weight is held to. None without a warm start.
    warm_objective = None

    # What a target dof outside the warm start's mesh takes: the prior for
    # the log controls (0), and for the fluidity prior mean the constant
    # baseline the forward fills the same ring with when it loads a MAP from
    # a smaller mesh (A = A_prior exp(phi) must stay positive there;
    # transfer.py has the measurement behind that), so the ring of a MAP
    # inverted on the production mesh is what production already runs.
    # Same-mesh warm starts miss nothing.
    warm_fill = {"fluidity_prior": float(A0) * a4_factor}

    def _warm_load(chk, source_mesh, name, space):
        source_field = chk.load_function(source_mesh, name=name)
        target = Function(space, name=name)
        n_missing, n_total, n_clamped = interpolate_with_fill(
            target, source_field, warm_fill.get(name, 0.0))
        if n_missing or n_clamped:
            PETSc.Sys.Print(
                f"    transfer {name}: {n_missing}/{n_total} target dofs "
                f"outside the warm start's mesh -> {warm_fill.get(name, 0.0)}; "
                f"{n_clamped} clamped to the source range")
        return target

    if warm_chk:
        PETSc.Sys.Print(f"  Loading warm start from {warm_chk}")
        with fd.CheckpointFile(warm_chk, "r") as chk:
            chk_mesh = chk.load_mesh()
            if chk.has_attr("/", "full_state_residual"):
                try:
                    warm_recorded = float(chk.get_attr("/", "full_state_residual"))
                except (TypeError, ValueError):
                    warm_recorded = None
            warm_objective = recorded_objective(chk)
            theta.assign(_warm_load(chk, chk_mesh, "log_friction", Q))
            phi.assign(_warm_load(chk, chk_mesh, "log_fluidity", Q))
            # A warm start on this mesh supplies its geometry, observations
            # and mixed state as well. One from another mesh (a 2 km MAP
            # warm-starting a 1 km inversion) supplies the controls and the
            # fluidity prior only: its cell-wise geometry would arrive
            # blocky, and this mesh's own BedMachine sample and raster
            # observations are what the forward that loads the MAP will
            # use. ISMIP7_WARM_START_GEOMETRY=0/1 overrides the default.
            same_mesh = meshes_match(chk_mesh, mesh)
            warm_geometry = os.environ.get(
                "ISMIP7_WARM_START_GEOMETRY", "1" if same_mesh else "0"
            ).strip() != "0"
            PETSc.Sys.Print(
                "    warm start is on " + ("this mesh" if same_mesh else "another mesh")
                + ("; taking its geometry, velocity_obs and state"
                   if warm_geometry else
                   "; taking its controls and fluidity prior only "
                   "(geometry and velocity_obs are this mesh's own)"))
            if not warm_geometry:
                raise_geometry = KeyError("warm start geometry not taken")
            try:
                if not warm_geometry:
                    raise raise_geometry
                H.assign(_warm_load(chk, chk_mesh, "thickness", Q_g))
                b.assign(_warm_load(chk, chk_mesh, "bed", Q_g))
                s.assign(_warm_load(chk, chk_mesh, "surface", Q_g))
                PETSc.Sys.Print(
                    "    geometry: thickness/bed/surface from warm start"
                )
            except (KeyError, RuntimeError, ValueError):
                PETSc.Sys.Print(
                    "    geometry: keeping BedMachine sample "
                    "(warm start has no thickness/bed/surface)"
                )
            try:
                if not warm_geometry:
                    raise raise_geometry
                u_obs.assign(_warm_load(chk, chk_mesh, "velocity_obs", V))
                PETSc.Sys.Print("    velocity_obs from warm start")
            except (KeyError, RuntimeError, ValueError):
                pass
            try:
                warm_A_prior = _warm_load(
                    chk, chk_mesh, "fluidity_prior", Q
                )
            except (KeyError, RuntimeError, ValueError):
                warm_A_prior = None
            try:
                if not warm_geometry:
                    raise raise_geometry
                u_ws = _warm_load(chk, chk_mesh, "velocity", V)
                M_ws = _warm_load(
                    chk, chk_mesh, "membrane_stress",
                    z.subfunctions[1].function_space(),
                )
                tau_ws = _warm_load(
                    chk, chk_mesh, "basal_stress",
                    z.subfunctions[2].function_space(),
                )
                z.subfunctions[0].assign(u_ws)
                z.subfunctions[1].assign(M_ws)
                z.subfunctions[2].assign(tau_ws)
                warm_loaded_z = True
                PETSc.Sys.Print(
                    "    mixed state: velocity/membrane/basal from warm start"
                )
            except (KeyError, RuntimeError, ValueError):
                warm_loaded_z = False
        if warm_loaded_z:
            # Full-n mixed state from prepare already sits at the physical
            # exponents; ramping 1→n would only discard that work.
            skip_continuation = True
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
    # Combined approach: He multiplies τ in momentum_balance AND
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
    # Budd pins N_hat=1 at the inversion geometry; freeze N_ref with the MAP /
    # timing-cache so forwards reproduce the inverted friction at t=0.
    N_ref = None
    if FRICTION == "budd":
        N_ref = Function(Q_g, name="N_ref").interpolate(
            max_value(effective_pressure(H, s), Constant(0.0))
        )
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
    # fluidity rather than log(A / const). This is the Recinos et al. (2023) fix
    # for the n=3 blow-up (a constant A0 baseline forced phi to carry all the
    # spatial fluidity structure). Frictional heating uses the balance C_w0.
    # When warm-starting from a prepare cache / MAP that already carries
    # fluidity_prior, reuse it: phi = log(A/A_prior) is meaningless against a
    # freshly recomputed prior.
    if warm_A_prior is not None:
        A_prior = warm_A_prior
        A_prior.rename("fluidity_prior")
        A_prior_lo, A_prior_hi = global_range(A_prior)
        PETSc.Sys.Print(
            f"  Fluidity prior A_prior in [{A_prior_lo:.2f}, "
            f"{A_prior_hi:.2f}] (from warm start)"
        )
    elif os.environ.get("ISMIP7_FLUIDITY_PRIOR", "thermo") == "thermo":
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

    # Floor-cell stabilizers shared with the forward (runconfig): without
    # them the inverted mixed state solved a DIFFERENT F from the one the
    # forward assembles at restart -- ||F||=1.3e1 here, 1.5e10 there, on the
    # same 2500/25000 state (2026-09-14) -- and the forward's fast path then
    # trusted that state. k_lim stays 0: the term is a rescue-only gain.
    stabilizers = residual_stabilizers()
    PETSc.Sys.Print(
        "  Residual stabilizers (shared with the forward): "
        + ", ".join(f"{key}={value:g}" for key, value in stabilizers.items())
    )
    # The forward's ocean-drag gate (front.ocean_drag_cells) on the geometry
    # this inversion fits, which is the forward's t=0 extent: the drag acts
    # only in open water a cell away from the ice, so the controls are fitted
    # to the same free front, with no drag under any floating ice, that the
    # forward then runs. Without it the MAP learned a front the drag held at
    # a tenth of its observed speed.
    drag_mask = Function(FunctionSpace(mesh, "DG", 0), name="drag_mask")
    _ice = Function(drag_mask.function_space()).project(H).dat.data_ro >= front_hmin()
    drag_mask.dat.data[:] = ocean_drag_cells(
        _ice, facet_neighbours(drag_mask.function_space()), _ice)
    PETSc.Sys.Print(
        f"  Ocean drag gate: {COMM_WORLD.allreduce(int(drag_mask.dat.data_ro.sum()))} "
        f"of {COMM_WORLD.allreduce(int(drag_mask.dat.data_ro.size))} cells, open water "
        f"a cell away from the ice (h < {front_hmin():g} m), as in the forward"
    )

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
                k_lim=0.0, **stabilizers, drag_mask=drag_mask,
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
    n_flow.assign(n_flow_val)
    m_slide.assign(m_slide_val)
    if skip_continuation:
        if warm_loaded_z:
            PETSc.Sys.Print(
                f"Warm start: accepting loaded mixed state at full "
                f"n_flow={n_flow_val:.1f}, m_slide={m_slide_val:.1f} "
                f"(no 1→n continuation)"
            )
        else:
            PETSc.Sys.Print(
                f"Warm start: single solve at full n_flow={n_flow_val:.1f}, "
                f"m_slide={m_slide_val:.1f} (ISMIP7_SKIP_CONTINUATION=1)"
            )
            slvr.solve()
    else:
        # Ramp n_flow (1 → n_flow_val) and m_slide (1 → m_slide_val) together:
        # a single solve at the full exponents can fail from a cold (u≈0)
        # guess. The ramp is not annotated, so it runs under the lane's
        # diagnostic solver (ISMIP7_DIAGNOSTIC_LINEAR_SOLVER), the one the
        # forward that loads this MAP cold-starts with, and climbs the
        # transient's step ladder. The full-Jacobian MUMPS LU that tlm_adjoint
        # differentiates through takes over at the converged state: on the
        # 1 km production mesh its factorisation failed mid-ramp
        # (DIVERGED_LINEAR_SOLVE, job 1612624) where the forward's condensed
        # GAMG had climbed the same ramp on the same mesh and MAP.
        ramp_params = diagnostic_solver_parameters(lane_solver_mode)
        ramp_params.update(_monitor_options)
        ramp_J, ramp_pre_jacobian = None, None
        if linearization_state(lane_solver_mode) == "frozen":
            from icepack2_tools.preconditioners import frozen_linearization
            ramp_J, ramp_pre_jacobian = frozen_linearization(F, z)
        ramp_solver = NonlinearVariationalSolver(
            NonlinearVariationalProblem(
                F, z, J=ramp_J, form_compiler_parameters=fc_params
            ),
            solver_parameters=ramp_params,
            options_prefix="ismip7_inversion_continuation_",
            pre_jacobian_callback=ramp_pre_jacobian,
        )
        ramp_ladder = ladder(continuation_steps())
        PETSc.Sys.Print(
            f"Warm start (continuation n_flow 1→{n_flow_val:.1f}, "
            f"m_slide 1→{m_slide_val:.1f}; "
            f"{diagnostic_solver_label(lane_solver_mode)}, "
            f"{linearization_state(lane_solver_mode)} linearization, "
            f"steps {'/'.join(str(s) for s in ramp_ladder)})..."
        )

        def _ramp_solve(attempt, step, steps, t):
            t0 = perf_counter()
            try:
                ramp_solver.solve()
            finally:
                snes = ramp_solver.snes
                PETSc.Sys.Print(
                    f"    rung {attempt} step {step}/{steps} "
                    f"n_flow={float(n_flow):.4g} m_slide={float(m_slide):.4g}: "
                    f"{_snes_reason_name(snes.getConvergedReason())} "
                    f"snes_its={snes.getIterationNumber()} "
                    f"linear_its={snes.getLinearSolveIterations()} "
                    f"fnorm={snes.getFunctionNorm():.3e} "
                    f"{perf_counter() - t0:.1f}s"
                )

        _, ramp_steps = ramp_exponents(
            _ramp_solve, z, n_flow, m_slide, n_flow_val, m_slide_val,
            ramp_ladder, report=PETSc.Sys.Print,
        )
        PETSc.Sys.Print(f"  Done ({ramp_steps} continuation steps)")

    # Forward-solve tolerance. Now that the stabilizers are shared, a
    # prepared warm start is already a converged solution of THIS F, so the
    # first annotated forward starts at the residual floor, where the
    # relative test can never pass and stol=0 disables the step exit: 200
    # silent MUMPS factorisations (the 2026-09-14 hang, moved one stage
    # earlier by the fix). Solve every forward to snes_atol_scale x the
    # residual the warm start's writer recorded -- the transient's own run
    # rule -- so a floor-level start confirms at iteration 0 and a moved
    # control vector gets an ordinary Newton solve to the same absolute
    # level the forward runs at. Without a record, keep the relative test
    # but enable the step-size exit so a floor-level start cannot grind.
    with assemble(F, form_compiler_parameters=fc_params).dat.vec_ro as _rv:
        f_warm = float(_rv.norm())
    if (
        warm_recorded is not None
        and np.isfinite(warm_recorded)
        and warm_recorded > 0.0
    ):
        forward_atol = snes_atol_scale() * warm_recorded
        sparams["snes_atol"] = forward_atol
        PETSc.Sys.Print(
            f"  Warm-start residual ||F||={f_warm:.3e} (writer recorded "
            f"{warm_recorded:.3e}); forward snes_atol={forward_atol:.3e} "
            f"({snes_atol_scale():g}x recorded)"
        )
    else:
        sparams["snes_stol"] = final_solve_bounds()["snes_stol"]
        PETSc.Sys.Print(
            f"  Warm-start residual ||F||={f_warm:.3e}; no recorded residual, "
            f"forwards use the relative test with snes_stol="
            f"{sparams['snes_stol']:g}"
        )
    state_solver_parameters = json.dumps(sparams, sort_keys=True)
    # The adjoint solves are LINEAR (one Newton step to rtol) and inherit the
    # forward's parameters by default. An absolute tolerance sized for the
    # forward residual lets them exit at iteration 0 whenever ||dJ/du|| is
    # small -- it is ~1e-3 here -- returning a zero adjoint, so the gradient
    # is the prior's alone and L-BFGS pulls the controls toward the prior
    # means while the misfit rises (job 10432790, 2026-09-14: adjoint time
    # 25 s -> 0.7 s, |grad| 47 -> 12, misfit +7% in four evaluations).
    adjoint_sparams = {
        key: value for key, value in sparams.items() if key != "snes_atol"
    }

    u_init = z.subfunctions[0]
    u_mag = Function(Q).interpolate(sqrt(inner(u_init, u_init)))
    PETSc.Sys.Print(f"  u_max = {global_max(u_mag):.0f} m/yr")

    # ── Forward function for tlm_adjoint ──
    # Normalize by the OBSERVED area so the misfit magnitude stays
    # comparable between masked and unmasked runs.
    area_val = assemble(obs_mask * dx(mesh))

    # ── Prior operators (ISMIP7_PRIOR_FORM) ──────────────────────────────
    # The single place that answers, for whichever form is active: the prior
    # energy R(theta), its dual gradient, and the metric
    # ISMIP7_GRAD_PRECOND=prior descends in. Everything downstream asks here
    # instead of rebuilding a form, so the cost the optimizer pays, the
    # gradient it follows and the metric it measures in cannot end up
    # describing three different priors -- the failure prior.py exists to
    # prevent, and the property Recinos et al. 2023 rely on when the UQ
    # eigendecomposes the misfit Hessian against the inversion's own prior.
    #
    #   laplacian   : B = A,        R = 0.5 theta' A theta        (one form)
    #   bilaplacian : B = A M^-1 A, R = 0.5 ||M^-1 A theta||^2_M  (one solve)
    #
    # A = delta*M + gamma*K either way; only the coefficients and the power
    # differ. The bi-Laplacian's (delta, gamma) come from a physical
    # (sigma, rho) through the closed forms of Villa et al. (2021), which the un-squared
    # operator does not possess.
    _prior_test = TestFunction(Q)
    if PRIOR_FORM == "bilaplacian":
        _prior_dg = {
            "theta": bilaplacian_coeffs(PRIOR_SIGMA_THETA, PRIOR_RHO),
            "phi": bilaplacian_coeffs(PRIOR_SIGMA_PHI, PRIOR_RHO),
        }
        PETSc.Sys.Print(
            f"  Prior: bi-Laplacian A M^-1 A; "
            f"sigma_theta={PRIOR_SIGMA_THETA:g} sigma_phi={PRIOR_SIGMA_PHI:g} "
            f"rho={PRIOR_RHO:g} m -> "
            f"delta={_prior_dg['theta'][0]:.4e} gamma={_prior_dg['theta'][1]:.4e}"
        )
    else:
        _prior_dg = {
            "theta": prior_operator_coeffs(GAMMA_THETA, area_val, L_REG),
            "phi": prior_operator_coeffs(GAMMA_PHI, area_val, L_REG),
        }
    _prior_aux = {k: Function(Q, name=f"prior_aux_{k}") for k in ("theta", "phi")}

    def _prior_energy_form(ctrl, which):
        """``R(ctrl)`` as something ``assemble`` or ``Functional.addto`` takes.

        Under `bilaplacian` this SOLVES ``M f = A ctrl`` into
        ``_prior_aux[which]`` as a side effect, so :func:`_prior_grad` -- which
        needs that ``f`` -- must be called after this and before the next
        control changes. EquationSolver annotates when a manager is running
        (the TAO path needs the solve on the tape for the gradient) and is an
        ordinary solve when one is not (the scipy path differentiates it by
        hand below).
        """
        d, g = _prior_dg[which]
        if PRIOR_FORM == "laplacian":
            return 0.5 * prior_operator_form(ctrl, ctrl, d, g)
        aux = _prior_aux[which]
        EquationSolver(
            bilaplacian_aux_residual(ctrl, aux, _prior_test, d, g) == 0,
            aux,
            form_compiler_parameters=fc_params,
        ).solve()
        return bilaplacian_energy_form(aux)

    def _prior_grad(ctrl, which):
        """``dR/dctrl``, assembled (a dual vector, as the misfit gradient is).

        ``A ctrl`` for the Laplacian; ``A M^-1 A ctrl = A f`` for the
        bi-Laplacian, with ``f`` the field :func:`_prior_energy_form` just
        solved for -- A is symmetric, so no second solve is needed.
        """
        d, g = _prior_dg[which]
        src = ctrl if PRIOR_FORM == "laplacian" else _prior_aux[which]
        return assemble(prior_operator_form(src, _prior_test, d, g))

    def _prior_metric_solvers(metric):
        """Per-control solvers for the prior COVARIANCE, the metric
        ISMIP7_GRAD_PRECOND=prior descends in: ``A^-1`` for the Laplacian,
        ``A^-1 M A^-1`` for the bi-Laplacian (the covariance
        action "L^-1 M L^-1"). Returns a callable taking the
        two assembled gradients and returning the two preconditioned
        directions."""
        _tr = fd.TrialFunction(Q)
        # A is symmetric positive definite by construction (delta, gamma > 0),
        # so Cholesky; it is constant, so this factors once and every
        # application below is a back-substitution.
        _fac = {"ksp_type": "preonly", "pc_type": "cholesky",
                "pc_factor_mat_solver_type": "mumps"}
        solvers = {
            which: fd.LinearSolver(
                assemble(prior_operator_form(_tr, _prior_test,
                                             *_prior_dg[which])),
                solver_parameters=_fac,
            )
            for which in ("theta", "phi")
        }
        # The consistent mass Riesz map, for metric == "mass": the same
        # mass-matrix preconditioner, reached through the same code path so
        # the two options differ only in the operator.
        _mass_solver = fd.LinearSolver(
            assemble(inner(_tr, _prior_test) * dx), solver_parameters=_fac
        )

        # First-step scaling, per control BLOCK. The two controls of a dual
        # inversion tend to have very different magnitudes, so scaling both
        # by their combined mean is a bad idea. theta (friction) and phi
        # (fluidity) are exactly that pair here, and a single combined scalar is what a 32 km probe took
        # theta to [-883, +14165] with -- it is a log deviation, so order 1.
        #
        # The scaling exists because L-BFGS's first step is -H_0 g at unit
        # length, with no curvature pair yet to rescale it, and on this path a
        # line search cannot recover: one evaluation outside the region where
        # the forward has a solution returns NaN and every later trial point
        # inherits it. A line search can cap the first step with amax; TAO's
        # lmvm gives no equivalent once H_0 is supplied, so the
        # bound goes on H_0 instead. ISMIP7_PRECOND_STEP0 is the largest change
        # the first step may make to a control, in that control's own units.
        _step0 = float(os.environ.get("ISMIP7_PRECOND_STEP0", "0.15"))
        _scale = {}

        def action(g_theta, g_phi):
            out = []
            for which, rhs in (("theta", g_theta), ("phi", g_phi)):
                x = Function(Q)
                if metric == "mass_consistent":
                    # The consistent mass Riesz map: the same mass-matrix
                    # preconditioner, the Riesz map of the L2 inner product. It removes
                    # the cell-size dependency and nothing else, which is why
                    # it is robust where the prior metric is delicate.
                    _mass_solver.solve(x, rhs)
                else:
                    solvers[which].solve(x, rhs)
                    if PRIOR_FORM == "bilaplacian":
                        # M x assembled against the test function IS the mass
                        # action, and it lands in the dual space the second
                        # A-solve wants, with no vector juggling.
                        y = Function(Q)
                        solvers[which].solve(
                            y, assemble(inner(x, _prior_test) * dx)
                        )
                        x = y
                if which not in _scale:
                    with x.dat.vec_ro as _x:
                        _n = _x.norm(PETSc.NormType.NORM_INFINITY)  # collective
                    _scale[which] = (_step0 / _n) if _n > 0.0 else 1.0
                    PETSc.Sys.Print(
                        f"  Metric scale [{which}]: alpha={_scale[which]:.4e} "
                        f"(unscaled first step |d{which}|_max={_n:.4e}, "
                        f"bounded to {_step0:g})"
                    )
                if _scale[which] != 1.0:
                    x.dat.data[:] *= _scale[which]
                out.append(x)
            return tuple(out)

        return action

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
    # The melt calibration the dH/dt melt source used, stamped into the MAP
    # under the attribute's original name, dhdt_melt_k_npz, so the artifact
    # carries the answer.
    k_npz_used = "none"
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
        # + the melt calibration through the Burgard quadratic-mixed-slope
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
            # The forward's melt calibration: per-basin deltaT at one K, the
            # tracked file unless another is named, else a legacy per-basin
            # K named with ISMIP7_K_PER_BASIN_NPZ. The tracked file is in
            # every checkout, so the source always carries melt.
            dT_npz = deltat_per_basin_npz()
            k_npz = None if dT_npz is not None else k_per_basin_npz()
            k_npz_used = dT_npz or k_npz
            _W2 = VectorFunctionSpace(mesh, "DG", 0)
            _xy = Function(_W2).interpolate(
                fd.SpatialCoordinate(mesh)).dat.data_ro.reshape(-1, 2)
            geom_xy = (_xy[:, 0].copy(), _xy[:, 1].copy())
            K_field = None
            if k_npz is not None:
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
                f"  dH/dt melt source: {os.path.basename(k_npz_used)}, "
                f"integrated {_melt_gt:.0f} Gt/yr at reference geometry")
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

    # Log-velocity misfit (ISSM convention).
    _eps_v = Constant(LOG_VEL_EPS)

    def _log_ratio(u_expr):
        sp = sqrt(u_expr[0] ** 2 + u_expr[1] ** 2 + Constant(1e-12))
        sp_obs = sqrt(u_obs[0] ** 2 + u_obs[1] ** 2 + Constant(1e-12))
        return ln((sp + _eps_v) / (sp_obs + _eps_v))

    _derived_w = None
    if LOG_VEL_WEIGHT.lower() == "auto":
        # Scale so the two velocity terms start out comparable: the ratio of
        # the chi^2 term to the unweighted log term at the state this run
        # starts from, used when the warm start holds no weight of its own.
        _u0 = z.subfunctions[0]
        _chi2_0 = float(assemble(
            (0.5 / area_val * obs_mask
             * ((_u0[0] - u_obs[0]) ** 2 / sig_ux ** 2
                + (_u0[1] - u_obs[1]) ** 2 / sig_uy ** 2)) * dx(mesh)))
        _log_0 = float(assemble(
            (0.5 / area_val * obs_mask * _log_ratio(_u0) ** 2) * dx(mesh)))
        _derived_w = (_chi2_0 / _log_0) if _log_0 > 0 else 0.0
    log_vel_w, log_vel_source, _log_vel_note = resolve_log_vel_weight(
        LOG_VEL_WEIGHT, _derived_w, warm_objective,
        misfit_norm=MISFIT_NORM, eps=LOG_VEL_EPS)
    if log_vel_source == "warm_start":
        PETSc.Sys.Print(
            f"  Log-velocity misfit (ISSM 103): auto weight {log_vel_w:.6g} "
            f"held from the warm start (chi2 {_chi2_0:.3e} / log "
            f"{_log_0:.3e} here would give {_derived_w:.4g}), "
            f"eps={LOG_VEL_EPS:g} m/yr")
    elif log_vel_source == "derived":
        PETSc.Sys.Print(
            f"  Log-velocity misfit (ISSM 103): auto weight {log_vel_w:.4g} "
            f"= chi2 {_chi2_0:.3e} / log {_log_0:.3e} at "
            f"the start, eps={LOG_VEL_EPS:g} m/yr")
    elif log_vel_w > 0.0 or _log_vel_note:
        PETSc.Sys.Print(
            f"  Log-velocity misfit (ISSM 103): weight {log_vel_w:g}, "
            f"eps={LOG_VEL_EPS:g} m/yr")
    if _log_vel_note:
        PETSc.Sys.Print(f"    {_log_vel_note}")

    def forward(theta_ctrl, phi_ctrl):
        clear_caches()
        F_ctrl = build_F(theta_ctrl, phi_ctrl)
        if skip_continuation:
            # Timing-matrix short invert starts from a full-n prepare cache;
            # stay at the physical exponents so each eval is one Newton solve.
            n_flow.assign(n_flow_val)
            m_slide.assign(m_slide_val)
            EquationSolver(
                F_ctrl == 0,
                z,
                solver_parameters=sparams,
                adjoint_solver_parameters=adjoint_sparams,
                form_compiler_parameters=fc_params,
            ).solve()
        else:
            # Continuation inside annotation for robustness — ramp both
            # n_flow and m_slide on the same [0,1] parameter.
            for t in np.linspace(0.0, 1.0, 5):
                n_flow.assign(1.0 + t * (n_flow_val - 1.0))
                m_slide.assign(1.0 + t * (m_slide_val - 1.0))
                EquationSolver(
                    F_ctrl == 0,
                    z,
                    solver_parameters=sparams,
                    adjoint_solver_parameters=adjoint_sparams,
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

        if log_vel_w > 0.0:
            # ISSM SurfaceLogVelMisfit: a relative (scale-free) error, so the
            # fit is not bought entirely in the slow interior.
            integrand = integrand + (
                Constant(0.5 * log_vel_w) / area_val * obs_mask
                * _log_ratio(u_sol) ** 2
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
    _log_chi2 = (0.5 / area_val * obs_mask * _log_ratio(_u_now) ** 2) * dx(mesh)
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

    def _residual_norm():
        """||F(z; theta, phi)|| at full n/m -- what every consumer recomputes."""
        with assemble(F, form_compiler_parameters=fc_params).dat.vec_ro as _rv:
            return float(_rv.norm())

    # Residual level the last accepted forward actually reached, and the
    # controls it was solved at. The publishing solve below is judged against
    # these: a state already at that level needs no Newton iterations.
    last_good_fnorm = [float("nan")]
    last_good_x = [None]

    # Provenance of the mixed state written by save_map(full_state=True):
    # its residual under the controls saved in the SAME file, and how the
    # publishing solve ended. The timing-cache restart fast path
    # (simulation.py) recomputes ||F|| itself; these attributes are the
    # audit trail for what it will find.
    full_state_solve = {
        "residual": float("nan"),
        "solve_reason": "",
        "solve_iterations": -1,
        "atol": float("nan"),
        "fnorm_ref": float("nan"),
    }

    def term_report():
        r"""``' vel=... dhdt=...'`` for the iteration line, or '' if disabled."""
        try:
            out = f" vel={last_good_vel_chi2[0]:.4e}"
            if log_vel_w > 0.0:
                out += f" log={float(assemble(_log_chi2)):.4e}"
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
    def save_map(path, *, full_state=False):
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
            if full_state:
                chk.save_function(z.subfunctions[0], name="velocity")
                chk.save_function(z.subfunctions[1], name="membrane_stress")
                chk.save_function(z.subfunctions[2], name="basal_stress")
                chk.save_function(H, name="H_init")
                chk.save_function(phi_eff, name="phi_eff")
                chk.save_function(C_w0, name="C_w0")
                if N_ref is not None:
                    chk.save_function(N_ref, name="N_ref")
            # The .msh this MAP was inverted on. A CheckpointFile mesh is named
            # "firedrake_default", so this is how the forward names its own
            # mesh and picks the matching per-mesh boundary-id sidecar.
            chk.set_attr("/", "mesh_basename", os.path.basename(mesh_fn))
            # Mesh PARAMETERS as well as the basename (the collaborators'
            # scheme, merged from upstream/integration). The forward resolves its
            # boundary_ids sidecar from these rather than from its own
            # environment: ISMIP7_BUFFER_M / ISMIP7_LC_COARSE can drift, and a
            # mismatched sidecar puts the calving BC on the wrong facets, where
            # ds(absent_id) integrates to zero -- silently wrong physics with
            # no crash. The basename covers meshes outside the standard naming
            # pattern; the parameters let bndids_filename() rebuild the name.
            chk.set_attr("/", "lc", int(lc))
            chk.set_attr("/", "lc_coarse", int(lc_coarse))
            chk.set_attr("/", "buffer_m", float(buffer_m))
            # The configuration theta/phi only mean anything under. The
            # derived MAP filename encodes all three, but ISMIP7_INVERSION
            # bypasses the name, so the forward needs them recorded to check
            # the MAP it was pointed at against the law it is about to run.
            chk.set_attr("/", "friction", str(FRICTION))
            chk.set_attr("/", "n_flow", float(n_flow_val))
            chk.set_attr("/", "geometry_space", str(geometry_space))
            chk.set_attr("/", "misfit_norm", MISFIT_NORM)
            chk.set_attr("/", "log_vel_weight", float(log_vel_w))
            # requested, derived (auto, at this run's start) or warm_start
            # (auto, held from the MAP this run warm-started from). A MAP
            # without it predates issue 68: if several chained links wrote it
            # under auto, each re-derived the weight, and log_vel_weight is
            # the last link's.
            chk.set_attr("/", "log_vel_weight_source", log_vel_source)
            chk.set_attr("/", "log_vel_eps", float(LOG_VEL_EPS))
            chk.set_attr("/", "gamma_theta", float(GAMMA_THETA))
            chk.set_attr("/", "gamma_phi", float(GAMMA_PHI))
            # Which prior this MAP was inverted under. laplacian and
            # bilaplacian are different priors with incomparable gammas, so a
            # MAP that does not say which one it paid cannot be interpreted.
            chk.set_attr("/", "prior_form", str(PRIOR_FORM))
            if PRIOR_FORM == "bilaplacian":
                chk.set_attr("/", "prior_sigma_theta", float(PRIOR_SIGMA_THETA))
                chk.set_attr("/", "prior_sigma_phi", float(PRIOR_SIGMA_PHI))
                chk.set_attr("/", "prior_rho", float(PRIOR_RHO))
            chk.set_attr("/", "dhdt_weight", float(dhdt_w))
            chk.set_attr("/", "dhdt_net_sigma", net_sigma_used)
            # Which per-basin K the dH/dt melt source used, or "none". The
            # fallback is deliberately non-fatal (the misfit is grounded-only
            # and melt is zero there), but a MAP that cannot say whether it
            # had the calibration cannot be told apart from one that did.
            chk.set_attr("/", "dhdt_melt_k_npz", str(k_npz_used))
            # How BedMachine was put onto the cells (runconfig.RASTER_SAMPLES).
            # theta/phi absorb the bed representation just as they absorb the
            # front treatment, so a forward must reproduce it.
            chk.set_attr("/", "raster_sample", raster_sample)
            if full_state:
                chk.set_attr("/", "t_yr", float(MATRIX_T_START))
                chk.set_attr("/", "friction", str(FRICTION))
                if str(FRICTION) == "budd":
                    chk.set_attr("/", "friction_gate", BUDD_SHELF_GATE)
                chk.set_attr("/", "geometry_space", str(geometry_space))
                chk.set_attr("/", "n_flow", float(n_flow_val))
                chk.set_attr("/", "a4_factor", float(a4_factor))
                # Two solver facts, kept apart: the mode the lanes consuming
                # this state must run (the cache contract), and the solver
                # that actually produced the state.
                chk.set_attr("/", "diagnostic_solver_mode", lane_solver_mode)
                chk.set_attr("/", "state_solver_mode", state_solver_mode)
                chk.set_attr(
                    "/", "state_solver_parameters", state_solver_parameters
                )
                chk.set_attr(
                    "/", "geometry_source", os.path.realpath(bm_fn)
                )
                chk.set_attr(
                    "/", "geometry_source_method", TARGET_MESH_GEOMETRY_METHOD
                )
                for key, value in full_state_solve.items():
                    chk.set_attr("/", f"full_state_{key}", value)

    # ── L-BFGS-B Inversion ──
    max_iter = int(os.environ.get("ISMIP7_MAXITER", "500"))
    # Relative-decrease stopping rule; 0 disables it and the budget above decides.
    ftol = float(os.environ.get("ISMIP7_FTOL", "1e-10"))
    min_iter = int(os.environ.get("ISMIP7_MIN_ITER", "3"))
    PETSc.Sys.Print("\nStarting L-BFGS-B inversion (theta + phi)...")
    PETSc.Sys.Print(f"  maxiter={max_iter} ftol={ftol:g} min_iter={min_iter}, "
                    f"nranks={COMM_WORLD.size}")

    global_ndof = len(func_to_global(theta))
    z_backup = z.copy(deepcopy=True)
    last_good_obj = [np.inf]
    last_x = [None]                      # controls of the last CONVERGED evaluation
    iteration_count = [0]
    timing_history = []
    timing_json = os.environ.get("ISMIP7_INVERSION_TIMING_JSON", "").strip()
    t_opt0 = perf_counter()

    def _eval_terms():
        """Assemble diagnostic term values for the JSON / print line."""
        terms = {"vel": float(last_good_vel_chi2[0])}
        if log_vel_w > 0.0:
            terms["log"] = float(assemble(_log_chi2))
        if use_dhdt:
            terms["dhdt"] = float(assemble(_dhdt_chi2))
            terms["net_gt_per_yr"] = (
                float(assemble(_net_form)) * 917.0 / 1e12
            )
        return terms

    def _write_timing_json(
        *, phase, message="", nit=None, nfev=None, final_solve=None
    ):
        if not timing_json or COMM_WORLD.rank != 0:
            return
        written = os.path.realpath(os.path.join(_map_dir, map_fn))
        published = os.environ.get("ISMIP7_MAP_OUT_FINAL", "").strip()
        payload = {
            "phase": phase,
            "map_path": os.path.realpath(published) if published else written,
            "map_path_written": written,
            "mesh_basename": os.path.basename(mesh_fn),
            "lc": int(lc),
            "lc_coarse": int(lc_coarse),
            "buffer_m": float(buffer_m),
            "ncores": int(COMM_WORLD.size),
            "maxiter": int(max_iter),
            "nit": int(nit if nit is not None else iteration_count[0]),
            "nfev": int(nfev if nfev is not None else iteration_count[0]),
            "message": str(message),
            "optimize_seconds": perf_counter() - t_opt0,
            "knobs": {
                "misfit_norm": MISFIT_NORM,
                "log_vel_weight_requested": LOG_VEL_WEIGHT,
                "log_vel_weight": float(log_vel_w),
                "log_vel_weight_source": log_vel_source,
                "log_vel_eps": float(LOG_VEL_EPS),
                "gamma_theta": float(GAMMA_THETA),
                "gamma_phi": float(GAMMA_PHI),
                "ftol": float(ftol),
                "min_iter": int(min_iter),
                "dhdt_weight": float(dhdt_w),
                "dhdt_net_sigma": float(net_sigma_used),
                "grad_precond": os.environ.get(
                    "ISMIP7_GRAD_PRECOND", "none"
                ).lower(),
                "warm_start": bool(warm_chk),
                "skip_continuation": bool(skip_continuation),
            },
            "evaluations": list(timing_history),
            # Set only by the "finished" record: how the publishing solve
            # ended and whether the full mixed state was written.
            "final_solve": final_solve,
        }
        atomic_write_json(timing_json, payload)

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
        last_x[0] = np.array(x_vec, copy=True)
        last_good_vel_chi2[0] = float(assemble(_vel_chi2))
        last_good_fnorm[0] = _residual_norm()
        last_good_x[0] = np.array(x_vec, copy=True)

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
        # Energy first, then gradient: under `bilaplacian` the energy solves
        # the M f = A theta the gradient reuses (see _prior_energy_form).
        reg_theta = float(assemble(_prior_energy_form(theta, "theta")))
        dR_theta = _prior_grad(theta, "theta")
        reg_phi = float(assemble(_prior_energy_form(phi, "phi")))
        dR_phi = _prior_grad(phi, "phi")

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

        if timing_json:
            try:
                terms = _eval_terms()
            except Exception:
                terms = {"vel": float(last_good_vel_chi2[0])}
            timing_history.append({
                "eval": iteration_count[0],
                "misfit": J_val,
                "reg_theta": reg_theta,
                "reg_phi": reg_phi,
                "total": total,
                "grad_norm": float(np.linalg.norm(total_grad)),
                "fwd_seconds": t_fwd,
                "adj_seconds": t_adj,
                "total_seconds": t_iter,
                "terms": terms,
            })
            _write_timing_json(phase="running", message="in progress")

        # Periodic checkpoint every 20 iterations
        if iteration_count[0] % 20 == 0:
            save_map(os.path.join(_map_dir, map_fn))
            PETSc.Sys.Print(f"    [checkpoint saved: iter {iteration_count[0]}]")

        return total, total_grad

    def _minimize_prior_metric():
        """L-BFGS in the prior metric, through TAO.

        The seed inverse Hessian is A^-1 with A the prior precision of the
        regularization this run pays (prior.prior_bilinear_form), applied per
        control because theta and phi carry different gamma. A is constant, so
        it is factored once here and every application is a back-substitution:
        milliseconds against a forward-plus-adjoint.

        Why TAO rather than the scipy path above: a metric can only reach
        L-BFGS-B as a change of variables, which needs a factor of A and not
        just its inverse -- the serial-only construction this replaces. TAO
        takes the inverse action directly, and tlm_adjoint wires H_0_action
        only for the lmvm and blmvm types, so the type is not free.

        What it costs, relative to the scipy path:

        * No line-search rescue. objective_and_gradient turns a failed forward
          or adjoint into an inflated objective with a zero gradient, which
          makes L-BFGS-B backtrack (it saved the Jul 18 2500 m Budd run at
          iteration 60). Through a ReducedFunctional the same failure raises
          and ends the run at the last periodic checkpoint. Use `none` or
          `mass` for a configuration whose solves are known to be marginal.
        * Per-iteration bookkeeping moves to the TAO monitor, so z_backup and
          the residual on record track the last EVALUATED point rather than
          the last accepted one. The publishing solve already tests that and
          reports `result.x != last accepted controls` when they differ.
        """
        from tlm_adjoint.firedrake import TAOSolver

        # The prior COVARIANCE action, from the same operator the objective
        # pays for: A^-1, or A^-1 M A^-1 under `bilaplacian`.
        _A_inv = _prior_metric_solvers(grad_precond)

        _nfev = [0]

        def forward_total(theta_ctrl, phi_ctrl):
            """The objective TAO differentiates: misfit plus BOTH prior terms.

            The scipy path adds the regularization outside the tape; here it
            has to be inside it, or TAO would descend on the misfit gradient
            in a metric built from a prior the objective never saw.
            """
            # Mirror this evaluation's controls and mixed state onto the
            # module-level Functions that save_map and the monitor read.
            # Raw dof writes: Function.assign under a running manager would
            # annotate, and these are bookkeeping, not part of the model.
            _nfev[0] += 1
            theta.dat.data[:] = theta_ctrl.dat.data_ro
            phi.dat.data[:] = phi_ctrl.dat.data_ro
            J = forward(theta_ctrl, phi_ctrl)
            J.addto(_prior_energy_form(theta_ctrl, "theta"))
            J.addto(_prior_energy_form(phi_ctrl, "phi"))
            for _zb, _z in zip(z_backup.subfunctions, z.subfunctions):
                _zb.dat.data[:] = _z.dat.data_ro
            return J

        gtol = float(os.environ.get("ISMIP7_GTOL", "0.0"))
        _step0_env = float(os.environ.get("ISMIP7_PRECOND_STEP0", "0.15"))
        if grad_precond == "mass_consistent":
            _desc = "consistent mass Riesz map (M^-1), after Recinos et al. (2023)"
        else:
            _desc = (
                f"{PRIOR_FORM} prior covariance "
                f"({'A^-1 M A^-1' if PRIOR_FORM == 'bilaplacian' else 'A^-1'})"
                " -- EXPERIMENTAL: the prior-preconditioned variant was not used in the reference work; "
                "its two prior-preconditioned H_0 attempts were left commented out as not working"
            )
        PETSc.Sys.Print(
            f"  Optimization metric: {_desc}; via TAO lmvm; "
            f"gatol={gtol:g} max_it={max_iter} step0={_step0_env:g}"
        )
        solver = TAOSolver(
            forward_total, [Q, Q],
            solver_parameters={
                "tao_type": "lmvm",
                "tao_max_it": max_iter,
                # The gradient norm TAO tests is the one M_inv_action defines,
                # i.e. sqrt(g' A^-1 g) -- mesh independent, unlike the raw l2
                # norm the scipy path prints. 0 leaves stopping to the ftol
                # rule in _monitor and to tao_max_it.
                "tao_gatol": gtol,
                "tao_grtol": 0.0,
                "tao_gttol": 0.0,
            },
            H_0_action=_A_inv, M_inv_action=_A_inv,
        )

        _t_last = [perf_counter()]
        _ftol_stop = FunctionalDecreaseStop(ftol, min_iter)

        def _monitor(tao):
            its, f_val, gnorm, _cnorm, _xdiff, _reason = tao.getSolutionStatus()
            iteration_count[0] = int(its)
            if _ftol_stop.update(iteration_count[0], f_val):
                tao.setConvergedReason(PETSc.TAO.ConvergedReason.CONVERGED_USER)
            now = perf_counter()
            t_iter = now - _t_last[0]
            _t_last[0] = now
            # z holds the last evaluated forward, mirrored into z_backup above.
            last_good_obj[0] = float(f_val)
            last_good_fnorm[0] = _residual_norm()
            last_good_vel_chi2[0] = float(assemble(_vel_chi2))
            _x = np.concatenate([func_to_global(theta), func_to_global(phi)])
            last_x[0] = _x
            last_good_x[0] = np.array(_x, copy=True)
            reg_theta = float(assemble(_prior_energy_form(theta, "theta")))
            reg_phi = float(assemble(_prior_energy_form(phi, "phi")))
            PETSc.Sys.Print(
                f"  iter {iteration_count[0]:3d}: "
                f"misfit={f_val - reg_theta - reg_phi:.6e} "
                f"reg_θ={reg_theta:.4e} reg_φ={reg_phi:.4e} "
                f"total={f_val:.6e} |grad|_A={gnorm:.4e} "
                f"dJ/J={_ftol_stop.criterion if _ftol_stop.criterion is not None else 0.0:.1e} "
                f"[total={t_iter:.1f}s]"
            )
            if timing_json:
                timing_history.append({
                    "eval": iteration_count[0],
                    "misfit": float(f_val - reg_theta - reg_phi),
                    "reg_theta": reg_theta,
                    "reg_phi": reg_phi,
                    "total": float(f_val),
                    "grad_norm": float(gnorm),
                    "total_seconds": t_iter,
                    "terms": {"vel": float(last_good_vel_chi2[0])},
                })
                _write_timing_json(phase="running", message="in progress")
            if iteration_count[0] > 0 and iteration_count[0] % 20 == 0:
                save_map(os.path.join(_map_dir, map_fn))
                PETSc.Sys.Print(f"    [checkpoint saved: iter {iteration_count[0]}]")

        solver.tao.setMonitor(_monitor)
        # TAO reports the iteration cap as a diverged reason and TAOSolver
        # turns any such reason into an exception -- but only AFTER writing the
        # solution back into (theta, phi). Reaching ISMIP7_MAXITER is how a
        # production inversion normally ends here, so the cap is read back from
        # TAO rather than treated as a failure.
        try:
            solver.solve([theta, phi])
        except RuntimeError:
            pass
        reason = int(solver.tao.getConvergedReason())
        message = (
            f"CONVERGED: relative functional decrease <= ftol={ftol:g}"
            if reason == int(PETSc.TAO.ConvergedReason.CONVERGED_USER)
            else "CONVERGED: gradient tolerance reached" if reason > 0
            else f"STOP: TAO reason {reason} (iteration limit is {max_iter})"
        )
        return SimpleNamespace(
            x=np.concatenate([func_to_global(theta), func_to_global(phi)]),
            nit=int(solver.tao.getIterationNumber()),
            # petsc4py exposes no evaluation count on TAO, and the number
            # that matters is the taped forwards this run actually paid for.
            nfev=_nfev[0],
            message=message,
        )

    # ── Optimization metric (ISMIP7_GRAD_PRECOND) ────────────────────────
    # L-BFGS-B works in raw dof coordinates, i.e. the Euclidean l2 metric, and
    # that metric is MESH-DEPENDENT: a gradient entry scales with the dof's
    # cell area, so the fine grounding-line cells -- exactly where the controls
    # must move -- get the smallest gradient entries and converge slowest.
    # `mass` optimizes in u = M^(1/2) x (M = lumped mass), which is steepest
    # descent in L2 and makes the rate mesh-independent; peers compensate for
    # the same defect with brute iteration counts (ISSM/M1QN3 runs 1000+300 in
    # cycles).
    #
    # `prior` goes further: it descends in the metric of the prior precision
    # A = delta_eff*M + gamma_eff*K (prior.prior_bilinear_form, the Hessian of
    # the regularization this run actually pays), so the reduced Hessian the
    # optimizer sees is A^-1 H_misfit + I -- clustered at 1 with one outlier
    # per data-informed mode. The iteration count then tracks the number of
    # those modes rather than the mesh.
    #
    # It is the MPI-parallel form of the whitening in the MISMIP TDDA
    # pipeline (scripts/run_inversion_weertman.py::PriorPreconditioner,
    # zeta = L^T theta with A = L L^T, validated there to 1e-16). That
    # construction is SERIAL-ONLY -- the PETSc native Cholesky factor's
    # solveForward/solveBackward and its full-local-vector indexing assume one
    # process -- so it cannot run at the 16-32 ranks these inversions use.
    # L-BFGS seeded with the initial inverse Hessian A^-1 is the same metric
    # and needs only a parallel solve against A, never a triangular factor.
    #
    # Scipy's L-BFGS-B takes no preconditioner, so `prior` runs through TAO
    # (tlm_adjoint's TAOSolver, whose H_0_action is exactly this seed). See
    # _minimize_prior_metric below for what that costs us.
    grad_precond = os.environ.get("ISMIP7_GRAD_PRECOND", "none").lower()
    _METRICS = ("none", "mass", "mass_consistent", "prior")
    if grad_precond not in _METRICS:
        raise ValueError(
            f"ISMIP7_GRAD_PRECOND must be one of {_METRICS}, "
            f"got {grad_precond!r}"
        )
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
    if grad_precond in ("mass_consistent", "prior"):
        result = _minimize_prior_metric()
    else:
        x0 = np.concatenate([func_to_global(theta), func_to_global(phi)])
        if grad_precond == "mass":
            x0 = x0 * _sqrtm
        result = scipy_minimize(
            objective_and_gradient,
            x0,
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": max_iter, "ftol": ftol, "gtol": 0},
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
        f"(misfit_norm={MISFIT_NORM} log_vel_weight={log_vel_w:g} "
        f"gamma_theta={GAMMA_THETA:g} gamma_phi={GAMMA_PHI:g} "
        f"dhdt_weight={dhdt_w:g})"
    )
    # The chain runner reads <ISMIP7_MAP_OUT>.done as "the MAP is on disk, do
    # not re-invert it". Write it here, the moment the checkpoint write returns:
    # the tail below (final solve, summary figure) runs for long enough that
    # the wall clock can kill the job inside it, and the runner's post-srun
    # rule would then never get to write the marker.
    if map_out and COMM_WORLD.rank == 0:
        with open(map_out + ".done", "w"):
            pass

    if timing_json:
        _write_timing_json(
            phase="final_solve",
            message=str(result.message),
            nit=result.nit,
            nfev=result.nfev,
        )

    # ── Final diagnostic ──
    # Publish the mixed state at the returned controls. When result.x is the
    # last evaluated point (the usual L-BFGS-B ending) z_backup already solves
    # F(z; result.x) = 0 to the level the annotated forwards reached, and a
    # fresh Newton solve from there sits at the rounding floor: the relative
    # test cannot pass, snes_stol=0 disables the step exit, and each
    # iteration is another full MUMPS factorisation. That ran the 2500/25000
    # timing invert silently to snes_max_it (200 iterations, ~30 min) and
    # then into a five-stage n-continuation retry -- which starts from a
    # full-n state and ramps m_slide, a float inside build_rc_residual, so it
    # only walked away from the solution. The transient runner solved the
    # same problem with a self-scaled absolute tolerance (simulation.py,
    # restart fast path); final_solve_parameters applies it here, bounds the
    # iteration counts, and always prints the converged reason, so this is
    # either a 0-iteration confirmation or a short, visible Newton solve from
    # a line-search neighbour. A plain NonlinearVariationalSolver (not
    # tlm_adjoint's EquationSolver, which discards its SNES) keeps the reason,
    # iteration count and function norm readable after a failure.
    PETSc.Sys.Print("\nFinal forward solve...")
    stop_manager()
    reset_manager()
    clear_caches()
    n_flow.assign(n_flow_val)
    m_slide.assign(m_slide_val)
    z.assign(z_backup)
    f_ref = float(last_good_fnorm[0])
    f0 = _residual_norm()
    x_same = bool(
        last_good_x[0] is not None
        and np.array_equal(last_good_x[0], _x_final)
    )
    final_sparams = final_solve_parameters(sparams, f_ref, viewer=_viewer)
    atol_final = float(final_sparams.get("snes_atol", float("nan")))
    PETSc.Sys.Print(
        f"  ||F(z_backup; result.x)|| = {f0:.3e}; last accepted forward "
        f"reached {f_ref:.3e}; result.x {'==' if x_same else '!='} last "
        f"accepted controls; snes_atol={atol_final:.3e} "
        f"snes_max_it={final_sparams['snes_max_it']}"
    )
    if "snes_atol" not in final_sparams:
        PETSc.Sys.Print(
            "  WARNING: no converged forward residual on record; the final "
            "solve falls back to the relative test alone"
        )
    final_solver = NonlinearVariationalSolver(
        NonlinearVariationalProblem(F, z, form_compiler_parameters=fc_params),
        solver_parameters=final_sparams,
    )
    final_solve_ok = False
    try:
        final_solver.solve()
        final_solve_ok = True
    except (fd.ConvergenceError, PETSc.Error) as exc:
        PETSc.Sys.Print(f"  Final solve raised {type(exc).__name__}: {exc}")
    final_reason = _snes_reason_name(final_solver.snes.getConvergedReason())
    final_its = int(final_solver.snes.getIterationNumber())
    if not final_solve_ok:
        # Never publish a half-converged Newton iterate.
        z.assign(z_backup)
    f_end = _residual_norm()
    PETSc.Sys.Print(
        f"  Final solve: {final_reason} after {final_its} Newton iterations, "
        f"||F|| {f0:.3e} -> {f_end:.3e}"
    )
    full_state_solve.update({
        "residual": f_end,
        "solve_reason": final_reason,
        "solve_iterations": final_its,
        "atol": atol_final,
        "fnorm_ref": f_ref,
    })
    # Whether the state on hand belongs to THIS MAP's controls; the figure
    # step below reads it under this name.
    final_state_ok = final_solve_ok

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

    # Rewrite the MAP with the full mixed state only when that state solves
    # the residual at the controls saved beside it. Timing caches are
    # published from this checkpoint without a second prepare.
    # The gate is whether the final solve CONVERGED at the saved controls,
    # not whether the state resembles the last optimizer eval. The old
    # guard compared assemble(_vel_chi2) against last_good_vel_chi2, but
    # the failure path had just done z.assign(z_backup) and
    # last_good_vel_chi2 is that same state's own metric -- so it compared
    # z_backup with itself and passed unconditionally. Every full-state MAP
    # written between then and 2026-09-14 carries controls from result.x
    # and a velocity solved for a different control vector; ||F|| at the
    # published state reached 1.5e10 and the forward blew up in 4 steps.
    _fnorm = f_end
    PETSc.Sys.Print(
        f"  Published-state residual ||F(z; theta, phi)|| = {_fnorm:.6e} "
        f"(final solve {'converged' if final_solve_ok else 'FAILED'}: "
        f"{final_reason})"
    )
    published = bool(final_solve_ok and np.isfinite(_fnorm))
    if published:
        save_map(chk_fn, full_state=True)
        PETSc.Sys.Print(f"Saved full mixed-state MAP: {chk_fn}")
    else:
        PETSc.Sys.Print(
            f"WARNING: NOT saving mixed state -- the final forward solve at "
            f"the saved controls did not converge ({final_reason} after "
            f"{final_its} iterations, ||F||={_fnorm:.3e}, result.x "
            f"{'==' if x_same else '!='} last accepted controls). A mixed "
            f"state from a different control vector would be accepted by "
            f"the restart fast path and never re-solved. Controls-only MAP "
            f"remains; timing-cache publish from this file will fail until "
            f"re-run."
        )

    if timing_json:
        _write_timing_json(
            phase="finished",
            message=str(result.message),
            nit=result.nit,
            nfev=result.nfev,
            final_solve={
                "reason": final_reason,
                "iterations": final_its,
                "fnorm_start": f0,
                "fnorm_end": f_end,
                "fnorm_ref": f_ref,
                "atol": atol_final,
                "result_x_is_last_accepted": x_same,
                "published": published,
            },
        )
        PETSc.Sys.Print(f"Inversion timing record -> {timing_json}")

    # ── Plot ──
    # Optional: the MAP is already written and the velocity saved above, so a
    # missing plotting dependency must not fail the run at this point.
    if not final_state_ok:
        PETSc.Sys.Print(
            "Skipping summary figure: the final solve did not converge, so the "
            "only state on hand is the last converged evaluation's, which "
            "belongs to different controls than this MAP. The MAP is saved."
        )
        return
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
