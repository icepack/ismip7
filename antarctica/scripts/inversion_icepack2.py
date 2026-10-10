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
    paused_manager,
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
    grounded_mask,
    floating_control_nodes,
    build_rc_residual,
    effective_pressure,
    rebase_log_friction,
    weertman_anchor,
)
from icepack2_tools.geometry import (
    cg1_lift, front_changed_nodes, raise_bed_to_lake_ice_base,
    sample_bed_thickness, sample_to_geometry,
)
from icepack2_tools.preconditioners import frozen_linearization, with_scpc_blocks
from icepack2_tools.taped_solve import StateSolverCache, taped_state_solve
from icepack2_tools.transfer import (
    harmonic_extension, interpolate_with_fill, load_checkpoint_mesh, meshes_match,
)
from icepack2_tools.grounding import height_above_flotation
from icepack2_tools.mpi_stats import (global_mean, global_range,
                                      global_max, global_size, global_count,
                                      global_rss_mib)
from icepack2_tools.naming import map_basename
from icepack2_tools.runconfig import (
    deltat_per_basin_npz, k_per_basin_npz,
    obs_data_root,
    BUDD_SHELF_GATE,
    friction as _friction, geometry_space as _geometry_space,
    raster_sample as _raster_sample,
    lc as _lc, lc_coarse as _lc_coarse, n_flow as _n_flow,
    eval_continuation, inversion_mesh_source, transfer_fill, drag_gate,
    DRAG_GATE_NONE, hvisc_floor, exact_front_version, warm_start_geometry,
    warm_start_state, ramp_slide_fixed, warm_start_fluidity,
    warm_start_state_fluidity, front_band_controls, front_band_extends,
)
DATA_DIR = obs_data_root()
from icepack2_tools.prior import (
    BilaplacianAuxSolver,
    bilaplacian_coeffs,
    bilaplacian_energy_form,
    prior_operator_coeffs,
    prior_operator_form,
)
from icepack2_tools.thermo_model import compute_fluidity_prior
from icepack2_tools.handoff import (
    OBJECTIVE_KEYS, OBJECTIVE_RECORD_KEYS, SUBELEMENT_SCHEME_VERSIONS,
    accepted_evaluation, frozen_in_control, handoff_gap, objective_mismatches)
from icepack2_tools.relaxation import (
    END_STATE_ATTR, RELAX_MAP_KEYS, anchor_ratio_counts, describe_relaxation,
    end_state_problems, inherited_geometry, is_relaxed)
from icepack2_tools.profiling import Spans
from icepack2_tools.optimization import (NOT_FINAL,
                                         FunctionalDecreaseStop,
                                         TrialFailures,
                                         recorded_objective,
                                         resolve_log_vel_weight)
from icepack2_tools.forcing import (load_racmo_smb_climatology,
                                    load_mean_annual_surface_temperature)
from icepack2_tools.runconfig import (
    TARGET_MESH_GEOMETRY_METHOD,
    anchor_length,
    file_sha256,
    front_hmin,
    lake_ice_base,
    residual_stabilizers,
)
from icepack2_tools.front import facet_neighbours, ocean_drag_cells, vertex_neighbours
from icepack2_tools.continuation import ladder, ramp_exponents
from icepack2_tools.solverconfig import (
    continuation_steps,
    diagnostic_solver_label,
    direct_forward_enabled,
    direct_forward_parameters,
    diagnostic_solver_mode,
    diagnostic_solver_parameters,
    mumps_analysis,
    final_solve_bounds,
    final_solve_parameters,
    inversion_adjoint_parameters,
    inversion_solver_mode,
    inversion_state_parameters,
    linearization_state,
    snes_atol_scale,
    snes_monitor_enabled,
    trial_rescue_rungs,
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

# ── Friction control (ISMIP7_FRICTION_CONTROL) ──────────────────────────
# `log` (default): theta = log(C / C_w0), a deviation from the balance anchor
# with prior mean zero, so the prior sets the friction's amplitude wherever
# the data are weak. `sqrt`: the control is alpha = sqrt(C) itself with a
# ZERO prior mean (Recinos et al. 2023 / fenics_ice: their mass term is
# negligible, the prior only smooths, and the data set the amplitude); the
# anchor is the initial guess only, as their driving-stress initial guess
# is. Residual laws and the bi-Laplacian prior only. The prior scale
# ISMIP7_PRIOR_SIGMA_ALPHA is in sqrt(MPa yr^(1/m) m^(-1/m)) (Recinos's
# m^(-1/6) yr^(1/6) Pa^(1/2) times 1e-3); "auto" takes the grounded median of
# the initial alpha, so the pointwise prior std is the typical friction and
# the amplitude is free in practice. ISMIP7_PRIOR_RHO_THETA is the friction
# correlation length (default ISMIP7_PRIOR_RHO; Recinos's l = sqrt(gamma/delta)
# is rho / sqrt(8)).
# `exp` (Rice, 27 Sep 2026: "instead of alpha squared we should use
# exp(alpha)"): C = C_ref exp(alpha), alpha a zero-mean control as under
# `sqrt`, but a log parameterisation, so a step in alpha is a multiplicative
# change of C everywhere: the weak beds under the ice streams, where C must
# fall by orders of magnitude, are as reachable as the stiff interior (under
# `sqrt` a step is additive in sqrt(C), small where C is small). C_ref is
# one scalar (ISMIP7_C_REF: "auto" = grounded median of the friction the
# start describes; 1 makes it exp(alpha) outright), so the prior mean is
# "typical friction" with no spatial structure, as Recinos's zero mean is.
# ISMIP7_PRIOR_SIGMA_ALPHA is then in log units ("auto" = 1, a factor e);
# the MAP carries alpha as log_friction and C_ref as a constant C_w0, so
# every consumer reads it as C_w0 exp(log_friction) with no new code.
FRICTION_CONTROL = os.environ.get("ISMIP7_FRICTION_CONTROL", "log").strip().lower()
if FRICTION_CONTROL not in ("log", "sqrt", "exp"):
    raise ValueError(
        f"ISMIP7_FRICTION_CONTROL must be log|sqrt|exp, got {FRICTION_CONTROL!r}")
if FRICTION_CONTROL in ("sqrt", "exp") and PRIOR_FORM != "bilaplacian":
    raise ValueError(
        f"ISMIP7_FRICTION_CONTROL={FRICTION_CONTROL} needs ISMIP7_PRIOR_FORM=bilaplacian")
# ISMIP7_GRAD_CHECK=1 runs its Taylor test on the TAO path only. Refused here
# on the scipy metrics, where it would run an ordinary inversion and write
# over ISMIP7_MAP_OUT.
if (os.environ.get("ISMIP7_GRAD_CHECK", "0").strip() == "1"
        and os.environ.get("ISMIP7_GRAD_PRECOND", "none").lower()
        not in ("mass_consistent", "prior")):
    raise ValueError(
        "ISMIP7_GRAD_CHECK=1 needs the TAO path "
        "(ISMIP7_GRAD_PRECOND=mass_consistent or prior)")
PRIOR_SIGMA_ALPHA = os.environ.get("ISMIP7_PRIOR_SIGMA_ALPHA", "auto").strip().lower()
C_REF = os.environ.get("ISMIP7_C_REF", "auto").strip().lower()
PRIOR_RHO_THETA = float(os.environ.get("ISMIP7_PRIOR_RHO_THETA", str(PRIOR_RHO)))

# ── Fluidity prior mean (ISMIP7_FLUIDITY_PRIOR) ─────────────────────────
# `pattyn` (default since 27 Sep 2026, Rice: "focus on the pattyn prior";
# its shelves start ten times closer to the observed speed than the thermal
# model's at 2 km): rate_factor of the depth-averaged Pattyn temperature raster
# (icepack2_tools/rheology_prior.py; ISMIP7_PATTYN_TEMP names the file,
# default <data root>/temp/Pattyn_2013.tif), the source Recinos et al. (2023)
# take their rheology prior mean from. Anything else: the constant A0.
FLUIDITY_PRIOR = os.environ.get("ISMIP7_FLUIDITY_PRIOR", "pattyn").strip().lower()

# ── Where the fluidity control acts (ISMIP7_FLUIDITY_CONTROL) ───────────
# `all` (default): phi = log(A / A_prior) everywhere. `floating`: phi acts
# on floating ice only, through the smooth grounded indicator (1 - He); on
# grounded ice the fluidity is the prior mean (a temperature field), and the
# velocity there is fitted by the friction alone. That is the split most
# ISMIP6/7 groups initialise with (a thermal rate factor on grounded ice,
# friction inverted; a shelf rheology inverted where there is no friction),
# and it removes the friction/rheology trade-off on grounded ice that let the
# data move the friction while the fluidity sat at its prior (Rice, 26 Sep
# 2026). phi is held at exactly zero on the nodes where the gate has shut it
# off (dual_friction.floating_control_nodes): the warm start's values there are
# zeroed and the gradient is projected out, so the MAP's log_fluidity carries
# no grounded variation for a reader to pair with this friction (issue #153).
# The L-BFGS-B path does the projection; the TAO path writes the controls
# itself, so it is refused.
FLUIDITY_CONTROL = os.environ.get("ISMIP7_FLUIDITY_CONTROL", "all").strip().lower()
if FLUIDITY_CONTROL not in ("all", "floating"):
    raise ValueError(
        f"ISMIP7_FLUIDITY_CONTROL must be all|floating, got {FLUIDITY_CONTROL!r}")
PHI_GROUNDED = "zero" if FLUIDITY_CONTROL == "floating" else "inverted"
if (FLUIDITY_CONTROL == "floating"
        and os.environ.get("ISMIP7_GRAD_PRECOND", "none").lower()
        in ("mass_consistent", "prior")):
    raise ValueError(
        "ISMIP7_FLUIDITY_CONTROL=floating holds grounded phi at zero on the "
        "L-BFGS-B path only (ISMIP7_GRAD_PRECOND=none or mass)")

# ── Which controls move (ISMIP7_INVERT) ──────────────────────────────────
# `both` (default): theta and phi descend together. `phi` / `theta`: only
# that control moves; the other stays at its current value (the prior mean,
# or the warm start's field) as a fixed coefficient of the forward. The
# OBJECTIVE is the same in every mode, so a staged inversion warm-starts each
# stage from the previous stage's MAP and the handoff check still applies.
# Staging is how the ISSM groups initialised ISMIP6 ("first invert for ice
# shelf viscosity and then for basal friction under grounded ice", Seroussi
# et al. 2020 Appendix C9): in the joint descent the friction gradient is
# ~50x the fluidity gradient per unit control, so L-BFGS fits the friction
# first and the fluidity sits at its prior for tens of iterations (the 26 Sep
# 32 km MAPs: phi rms 0.01 after 3 iterations). TAO path only.
INVERT = os.environ.get("ISMIP7_INVERT", "both").strip().lower()
if INVERT not in ("both", "phi", "theta"):
    raise ValueError(f"ISMIP7_INVERT must be both|phi|theta, got {INVERT!r}")

# ── Likelihood scale (ISMIP7_MISFIT_SCALE) ───────────────────────────────
# The data terms are normalised by the observed AREA, so the misfit is a mean
# chi^2 per node (1 at a perfect fit) while the prior energy is a sum over
# the whole domain; the prior therefore outweighs the data by about the
# number of observed nodes, and a MAP sits on its prior mean (91-94% of the
# 2 km MAP within 0.1 of it, Sep 2026). A Gaussian likelihood is the SUM of
# the per-datum chi^2, which is what fenics_ice assembles (Recinos et al.
# 2023) and why their zero-mean friction prior can leave the amplitude to the
# data. `nodes` multiplies every data term by the number of observed nodes
# (the sum over the control space's nodes); a number is a plain multiplier;
# 1 (default) keeps the mean. Recorded in the MAP as misfit_scale.
MISFIT_SCALE = os.environ.get("ISMIP7_MISFIT_SCALE", "1").strip().lower()

# ── Sub-element grounding scheme (ISMIP7_SUBELEMENT_FRICTION) ───────────
# ISSM's SEP2 through icepack_tools (icepack2_tools/subelement.py): the
# basal friction is integrated over the grounded part of each cell, so the
# grounding line is exact inside a cell instead of a cell-wise staircase.
# Budd then runs with N_hat = 1 on the grounded part (no N_ref, no delta
# floor). ISMIP7_EXACT_FRONT adds the exact depth-integrated push on a
# calving front inside the mesh: on by default under the scheme, since it
# is the shared residual's default, and off by default for cell-wise
# friction, where =1 adds dual_friction.front_cliff_correction (a grounded
# marine cliff otherwise gets 15 to 33 % too little push, issue #153) and =2
# its exposed-face form, which leaves no push against rock above the ice
# surface (issue #166; cell-wise friction only). The MAP records the version.
SUBELEMENT_FRICTION = os.environ.get("ISMIP7_SUBELEMENT_FRICTION", "0").strip() == "1"
# ISMIP7_SUBELEMENT_SCHEME: sep1 (ISSM's default: whole-cell quadrature,
# drag times the grounded fraction; the default here since 1 Oct 2026) or
# sep2 (grounded-part quadrature); recorded in the MAP and followed by the
# forward. On the same 2 km start SEP1 reached misfit 3.42e9 in 21 accepted
# iterations with no failed line-search trial, SEP2 4.21e9 with nine, at a
# quarter of SEP1's iterations per hour (NOTS scavenge chains, 1 Oct).
SUBELEMENT_SCHEME = os.environ.get("ISMIP7_SUBELEMENT_SCHEME", "sep1").strip().lower()
if SUBELEMENT_SCHEME not in ("sep2", "sep1"):
    raise ValueError(f"ISMIP7_SUBELEMENT_SCHEME must be sep2 or sep1, not {SUBELEMENT_SCHEME!r}")
# The form of that scheme this code builds, recorded with it in the MAP
# (icepack2_tools.handoff.SUBELEMENT_SCHEME_VERSIONS).
SUBELEMENT_SCHEME_VERSION = SUBELEMENT_SCHEME_VERSIONS[SUBELEMENT_SCHEME]
EXACT_FRONT = exact_front_version(os.environ.get(
    "ISMIP7_EXACT_FRONT", "1" if SUBELEMENT_FRICTION else "0"))
if SUBELEMENT_FRICTION and EXACT_FRONT == 2:
    raise ValueError(
        "exact_front version 2 (the exposed-face push, issue #166) is not in "
        "icepack_tools.momentum.front_cliff_correction; the sub-element scheme "
        "takes 0 or 1")
if MISFIT_SCALE != "nodes":
    try:
        float(MISFIT_SCALE)
    except ValueError:
        raise ValueError(
            f"ISMIP7_MISFIT_SCALE must be a number or 'nodes', got {MISFIT_SCALE!r}")

# Friction law: "budd" (power-law dual, default) or "regularized_coulomb"
# (Joughin/Schoof RC residual: grounded-only inference, exact-zero shelves).
FRICTION = _friction()
# Exact-zero-shelf residual laws share the C_w0/He/composite structure.
USE_RESIDUAL = FRICTION in ("regularized_coulomb", "budd")
USE_RC = USE_RESIDUAL  # geometry/anchor handling is shared
if SUBELEMENT_FRICTION and (not USE_RESIDUAL or FRICTION not in ("budd", "weertman")):
    raise ValueError("ISMIP7_SUBELEMENT_FRICTION=1 needs ISMIP7_FRICTION=budd (or "
                     "weertman): the scheme integrates a velocity-only stress")
if FRICTION_CONTROL in ("sqrt", "exp") and not USE_RESIDUAL:
    raise ValueError(f"ISMIP7_FRICTION_CONTROL={FRICTION_CONTROL} needs a residual friction "
                     "law (budd or regularized_coulomb)")
C0_RC = float(os.environ.get("ISMIP7_RC_C0", "0.5"))
# Buffer-node (h_clamp=0) coercivity controls; see dual_friction.build_rc_residual.
# h_visc_floor (membrane-only thickness floor) is the primary cure. The first
# water row carries it too and so couples the ice front to the ocean drag one
# cell out; runconfig.hvisc_floor's 2.5 m keeps a fifth of the 10 m coupling
# (GEOMETRY_DISCRETIZATION.md). c_w0_floor is off by default (unnecessary once
# h_visc_floor is on).
RC_HVISC_FLOOR = hvisc_floor()
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

    mesh_lc = lc
    mesh_lc_coarse = lc_coarse
    mesh_buffer_m = buffer_m
    mesh_fn, mesh_from_warm = inversion_mesh_source(
        mesh_filename(mesh_lc_coarse, mesh_lc, mesh_buffer_m))
    if mesh_from_warm:
        # The warm start's own mesh. mesh_fn becomes the .msh basename that
        # checkpoint recorded, which picks the boundary-id sidecar and is
        # stamped into this run's MAP, so the chain and the forward that loads
        # the result name the mesh the controls were inverted on.
        PETSc.Sys.Print(f"Loading mesh from the warm start: {mesh_fn}")
        (mesh, mesh_fn, mesh_lc, mesh_lc_coarse,
         mesh_buffer_m) = load_checkpoint_mesh(mesh_fn)
        PETSc.Sys.Print(
            f"  recorded mesh: {mesh_fn} "
            f"(lc={mesh_lc}, lc_coarse={mesh_lc_coarse}, "
            f"buffer_m={mesh_buffer_m:g})"
        )
    else:
        PETSc.Sys.Print(f"Loading mesh: {mesh_fn}")
        mesh = Mesh(mesh_fn)
    # num_vertices()/num_cells() count this rank's plex, halo included; the
    # coordinate dofs and the owned cell set are reduced to global totals.
    PETSc.Sys.Print(f"  {global_size(mesh.coordinates)} vertices, "
                    f"{mesh.comm.allreduce(mesh.cell_set.size)} cells")

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
    map_fn = (os.path.basename(map_out) if map_out
              else map_basename(FRICTION, mesh_lc))

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
    # h_clamp default 0.0: invert against the *true* BedMachine geometry,
    # including h=0 over the buffered ocean region. Composite rheology
    # (added below) keeps the SNES nonsingular where h=0.
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "0.0"))
    # The front sampling (vertex_front) rebuilds the marine front cells from
    # BedMachine's own mask (geometry.front_cells, issue #167).
    b, H, front_counts = sample_bed_thickness(
        bm_fn, Q_g, Q, floor=h_clamp, method=raster_sample)
    if front_counts is not None:
        PETSc.Sys.Print(f"  Front cells ({raster_sample}): {front_counts}")
    PETSc.Sys.Print(f"  H clamp: {h_clamp} m  "
                    f"(nodes h<=1m: "
                    f"{global_count(H.dat.data_ro <= 1.0, mesh.comm)} / "
                    f"{global_size(H)})")
    # Lake Vostok: BedMachine's bed there is the lake floor, so b + H would
    # sink the surface by the water column (runconfig.lake_ice_base).
    LAKE_ICE_BASE = lake_ice_base()
    ANCHOR_LENGTH = anchor_length()     # read here: the warm start checks it too
    if LAKE_ICE_BASE:
        n_lake = raise_bed_to_lake_ice_base(b, H, bm_fn, Q_g, Q, method=raster_sample)
        PETSc.Sys.Print(f"  Lake ice base: bed raised to BedMachine's ice base on "
                        f"{n_lake} geometry dofs over its subglacial lake")
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
    # benchmark" §7). The linear solve of every annotated forward, of the
    # adjoint tlm_adjoint solves against it, and of the publishing solve
    # follows ISMIP7_INVERSION_LINEAR_SOLVER: full_mumps is the full
    # mixed-Jacobian MUMPS LU, scpc_mumps and scpc_gamg (the default) the
    # transient's condensed modes (icepack2_tools/taped_solve.py). The MAP records it as
    # state_solver_mode beside the lane contract diagnostic_solver_mode.
    # The condensed factorizations (scpc_mumps, and schur_mumps's velocity
    # block) run MUMPS's own sequential analysis unless the environment names
    # one (solverconfig.mumps_analysis): under the distributed PT-Scotch
    # analysis the adjoint factorizations aborted the 2 km chains in
    # MUMPS_LOAD_RECV_MSGS. Set before the first solver options are resolved
    # below, and read now so an invalid value fails here.
    os.environ.setdefault("ISMIP7_MUMPS_ANALYSIS", "sequential")
    mumps_analysis()
    state_solver_mode = inversion_solver_mode()
    sparams = inversion_state_parameters(state_solver_mode)
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
    # A solver condensing with SCPC needs the (M, tau) structural-zero blocks
    # in its form (preconditioners.with_scpc_blocks), and an assembled one is
    # better without them: they are entries of the AIJ sparsity. So the taped
    # form F carries them under an scpc_* inversion solver and the startup
    # ramp's form under an scpc_* lane solver, each for its own solver.
    state_scpc = state_solver_mode.startswith("scpc_")
    lane_scpc = lane_solver_mode.startswith("scpc_")
    PETSc.Sys.Print(
        f"  Inversion linear solver: {diagnostic_solver_label(state_solver_mode)} "
        f"(ISMIP7_INVERSION_LINEAR_SOLVER={state_solver_mode}); startup ramp: "
        f"{diagnostic_solver_label(lane_solver_mode)}")
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
    warm_prior_origin = None
    # The warm start's fluidity prior under ISMIP7_WARM_START_PHI=physical:
    # phi is moved onto this run's prior so A = A_prior exp(phi) is kept.
    warm_A_prior_rebase = None
    warm_loaded_z = False
    warm_state_guess = False
    # the warm start's anchor length, when its theta is to be rebased onto
    # this run's anchor once that is built (ISMIP7_WARM_START_THETA=physical)
    warm_theta_anchor = None
    # alpha = sqrt(C) of a warm start inverted on the sqrt control
    warm_alpha = None
    # theta holds an exp-control warm start's alpha = ln(C / C_ref), a
    # deviation from the constant C_ref and not from any anchor
    warm_exp_alpha = False
    # Residual the warm start's writer reached under the shared F (stamped by
    # save_model_state and by save_map): the forwards' absolute tolerance.
    warm_recorded = None
    # Everything the warm start's writer recorded about ITS objective and the
    # objective value at its checkpointed iterate (icepack2_tools.handoff).
    warm_attrs = {}
    # Where the warm start's geometry came from (icepack2_tools.relaxation):
    # its record, whether it is a relaxation's end state, the anchor that
    # state's year ran with, and the record this run writes when it takes a
    # relaxed geometry (None: this run's own BedMachine sample).
    warm_geometry_attrs = {}
    warm_end_state = False
    warm_C_w0 = None
    run_geometry = None
    # The log-velocity term the warm start was minimised under, which an
    # "auto" weight is held to. None without a warm start.
    warm_objective = None

    # What a target dof outside the warm start's mesh takes, as the forward
    # does when it loads a MAP from a smaller mesh (ISMIP7_TRANSFER_FILL).
    # Under `extend` (the default) the controls continue harmonically from
    # the warm start's outline and the fluidity prior's logarithm does too
    # (transfer.harmonic_extension), so the bi-Laplacian prior pays nothing
    # for a step there; the constants below remain only where the extension
    # cannot reach. Under `constant` they are the fill: the prior for the log
    # controls (0), and for the fluidity prior mean the constant baseline,
    # which keeps A = A_prior exp(phi) positive (transfer.py has the
    # measurement behind that). Same-mesh warm starts miss nothing.
    warm_fill = {"fluidity_prior": float(A0) * a4_factor}
    warm_fill_mode = transfer_fill()
    # name -> extend its logarithm
    warm_extend = ({"log_friction": False, "log_fluidity": False,
                    "sqrt_friction": False, "fluidity_prior": True}
                   if warm_fill_mode == "extend" else {})

    def _warm_load(chk, source_mesh, name, space):
        source_field = chk.load_function(source_mesh, name=name)
        target = Function(space, name=name)
        fill = warm_fill.get(name, 0.0)
        n_missing, n_total, n_clamped = interpolate_with_fill(
            target, source_field, fill,
            extend=name in warm_extend, log=warm_extend.get(name, False))
        if n_missing or n_clamped:
            how = (f"a harmonic extension{' of its logarithm' if warm_extend[name] else ''} "
                   f"({fill:g} where it cannot reach)" if name in warm_extend else f"{fill:g}")
            PETSc.Sys.Print(
                f"    transfer {name}: {n_missing}/{n_total} target dofs "
                f"outside the warm start's mesh -> {how}; "
                f"{n_clamped} clamped to the source range")
        return target

    def _warm_state(chk, source_mesh):
        # a published MAP's state, else a periodic checkpoint's
        try:
            pre = ""
            u_ws = _warm_load(chk, source_mesh, "velocity", V)
        except (KeyError, RuntimeError, ValueError):
            pre = "ckpt_"
            u_ws = _warm_load(chk, source_mesh, "ckpt_velocity", V)
        M_ws = _warm_load(chk, source_mesh, f"{pre}membrane_stress",
                          z.subfunctions[1].function_space())
        tau_ws = _warm_load(chk, source_mesh, f"{pre}basal_stress",
                            z.subfunctions[2].function_space())
        return u_ws, M_ws, tau_ws

    if warm_chk:
        PETSc.Sys.Print(f"  Loading warm start from {warm_chk}")
        with fd.CheckpointFile(warm_chk, "r") as chk:
            # Under ISMIP7_MESH=checkpoint the compute mesh IS this file's
            # mesh, so its fields load onto it directly, with no transfer.
            chk_mesh = mesh if mesh_from_warm else chk.load_mesh()
            if chk.has_attr("/", "full_state_residual"):
                try:
                    warm_recorded = float(chk.get_attr("/", "full_state_residual"))
                except (TypeError, ValueError):
                    warm_recorded = None
            warm_objective = recorded_objective(chk)
            _warm_anchor = (float(chk.get_attr("/", "friction_anchor_length"))
                            if chk.has_attr("/", "friction_anchor_length") else 0.0)
            _warm_lake = (int(chk.get_attr("/", "lake_ice_base"))
                          if chk.has_attr("/", "lake_ice_base") else 0)
            _warm_rs = (str(chk.get_attr("/", "raster_sample")).lower()
                        if chk.has_attr("/", "raster_sample") else "vertex")
            for _key in OBJECTIVE_KEYS + OBJECTIVE_RECORD_KEYS:
                if chk.has_attr("/", _key):
                    warm_attrs[_key] = chk.get_attr("/", _key)
            for _key in RELAX_MAP_KEYS + (END_STATE_ATTR, "geometry_source",
                                          "geometry_source_method", "t_yr", "stalled"):
                if chk.has_attr("/", _key):
                    warm_geometry_attrs[_key] = chk.get_attr("/", _key)
            warm_end_state = bool(int(warm_geometry_attrs.get(END_STATE_ATTR, 0) or 0))
            if "objective_total" in warm_attrs:
                PETSc.Sys.Print(
                    f"    objective recorded at iteration "
                    f"{int(float(warm_attrs.get('objective_iteration', -1)))}: "
                    f"total {float(warm_attrs['objective_total']):.6e}")
            # theta is a log-deviation from the warm start's anchor. Under a
            # different anchor the same theta is a different friction:
            # ISMIP7_WARM_START_THETA=0 starts it at the new prior mean, and
            # =physical rebases it so the friction C_w0 exp(theta) is the warm
            # start's (the first solve then reproduces the warm start's).
            _theta_mode = os.environ.get("ISMIP7_WARM_START_THETA", "1").strip()
            if _theta_mode not in ("0", "1", "physical"):
                raise ValueError(
                    f"ISMIP7_WARM_START_THETA={_theta_mode!r}: use 1, 0 or physical")
            if _theta_mode == "0":
                PETSc.Sys.Print("    log_friction: prior mean (ISMIP7_WARM_START_THETA=0)")
            else:
                theta.assign(_warm_load(chk, chk_mesh, "log_friction", Q))
                if _theta_mode == "physical":
                    warm_theta_anchor = _warm_anchor
                elif _warm_anchor != ANCHOR_LENGTH:
                    PETSc.Sys.Print(
                        f"    WARNING: log_friction taken from a MAP whose anchor "
                        f"length is {_warm_anchor:g} m, not this run's "
                        f"{ANCHOR_LENGTH:g} m: the same theta means a different "
                        f"friction. ISMIP7_WARM_START_THETA=physical keeps the "
                        f"friction; =0 starts at the prior mean.")
            # =physical keeps the warm start's fluidity under this run's
            # own prior mean (ISMIP7_FLUIDITY_PRIOR): phi is rebased once the
            # prior is built, below.
            _phi_mode = os.environ.get("ISMIP7_WARM_START_PHI", "1").strip()
            if _phi_mode not in ("0", "1", "physical"):
                raise ValueError(
                    f"ISMIP7_WARM_START_PHI={_phi_mode!r}: use 1, 0 or physical")
            # ISMIP7_WARM_START_FLUIDITY: phi and the prior it deviates from
            # come from another MAP, so A = A_prior exp(phi) is that MAP's.
            _fl_chk = warm_start_fluidity()
            _fl_prior = None
            _fl_origin = None
            # the fluidity MAP's raster sampling (the front band's controls)
            _fl_rs = None
            # its mixed state, the first guess under ISMIP7_WARM_START_STATE=fluidity
            _fl_state = None
            if _phi_mode == "0":
                if _fl_chk:
                    raise ValueError(
                        "ISMIP7_WARM_START_FLUIDITY names a fluidity to start from "
                        "and ISMIP7_WARM_START_PHI=0 the prior mean: set one")
                PETSc.Sys.Print("    log_fluidity: prior mean (ISMIP7_WARM_START_PHI=0)")
            elif _fl_chk:
                with fd.CheckpointFile(_fl_chk, "r") as _fl:
                    _fl_fc = (str(_fl.get_attr("/", "fluidity_control")).lower()
                              if _fl.has_attr("/", "fluidity_control") else "all")
                    if _fl_fc != FLUIDITY_CONTROL:
                        raise ValueError(
                            f"ISMIP7_WARM_START_FLUIDITY={_fl_chk} was inverted "
                            f"under fluidity_control={_fl_fc}, this run under "
                            f"{FLUIDITY_CONTROL}")
                    _fl_rs = (str(_fl.get_attr("/", "raster_sample")).lower()
                              if _fl.has_attr("/", "raster_sample") else "vertex")
                    _fl_mesh = _fl.load_mesh()
                    phi.assign(_warm_load(_fl, _fl_mesh, "log_fluidity", Q))
                    _fl_prior = _warm_load(_fl, _fl_mesh, "fluidity_prior", Q)
                    _fl_origin = (str(_fl.get_attr("/", "fluidity_prior_origin"))
                                  if _fl.has_attr("/", "fluidity_prior_origin")
                                  else f"warm start {os.path.basename(_fl_chk)}")
                    if warm_start_state_fluidity(same_mesh=meshes_match(_fl_mesh, mesh)):
                        _fl_state = _warm_state(_fl, _fl_mesh)
                PETSc.Sys.Print(
                    f"    log_fluidity and fluidity_prior from {_fl_chk} "
                    f"(ISMIP7_WARM_START_FLUIDITY); log_friction from the warm start")
            else:
                phi.assign(_warm_load(chk, chk_mesh, "log_fluidity", Q))
            # A MAP inverted on the sqrt(C) control carries alpha itself; its
            # log_friction is zero by construction.
            _warm_fc = (str(chk.get_attr("/", "friction_control"))
                        if chk.has_attr("/", "friction_control") else "log")
            warm_exp_alpha = _warm_fc == "exp" and _theta_mode != "0"
            if _warm_fc == "log" and _theta_mode != "0":
                # The band a forward bounds a loaded MAP's log controls to
                # (simulation.py, ISMIP7_MAP_CLIP: 10 at n=3, 6 at n=4). The
                # 94-iteration sigma=3 MAP carries ~30 spikes up to theta=19
                # (C x 1e8) at the grounding line, and from them the first
                # n=1 solve of a sub-element sqrt run wandered at ||F|| 1e13-1e21.
                _clip = float(os.environ.get(
                    "ISMIP7_MAP_CLIP", "6.0" if float(n_flow_val) == 4.0 else "10.0"))
                if _clip > 0.0:
                    _d = theta.dat.data
                    _n_clip = COMM_WORLD.allreduce(int((np.abs(_d) > _clip).sum()))
                    np.clip(_d, -_clip, _clip, out=_d)
                    if _n_clip:
                        PETSc.Sys.Print(
                            f"    log_friction: {_n_clip} node(s) bounded to "
                            f"|theta|<={_clip:g} (ISMIP7_MAP_CLIP, as a forward does)")
            if _warm_fc == "sqrt" and _theta_mode != "0":
                warm_alpha = _warm_load(chk, chk_mesh, "sqrt_friction", Q)
                PETSc.Sys.Print("    sqrt_friction (alpha) from warm start")
            # A warm start on this mesh supplies its geometry, observations
            # and mixed state as well. One from another mesh (a 2 km MAP
            # warm-starting a 1 km inversion) supplies the controls and the
            # fluidity prior only: its cell-wise geometry would arrive
            # blocky, and this mesh's own BedMachine sample and raster
            # observations are what the forward that loads the MAP will
            # use. ISMIP7_WARM_START_GEOMETRY=0/1 overrides the default.
            same_mesh = meshes_match(chk_mesh, mesh)
            # A MAP inverted before the lake fix carries the lake bowl in its
            # geometry; taking it would undo the fix, so by default it is not.
            # Nor does one whose geometry was sampled another way: a front
            # sampling's cells differ at the front (issue #167), and a refit
            # from a vertex MAP is how a MAP gets them.
            warm_geometry = warm_start_geometry(
                same_mesh=same_mesh, same_lake=_warm_lake == int(LAKE_ICE_BASE),
                warm_sampling=_warm_rs, run_sampling=raster_sample)
            if same_mesh and _warm_rs != raster_sample:
                PETSc.Sys.Print(
                    f"    warm start records raster_sample={_warm_rs}, this run "
                    f"{raster_sample}: its geometry is not taken")
            # The warm start's controls were fitted against the front its own
            # sampling built. Under ISMIP7_WARM_START_FRONT_EXTEND=1 theta and
            # phi over the cells the front rule rebuilt or emptied become a
            # harmonic blend of the ice upstream and the nodes seaward of the
            # band, which keep the warm start's values, and the refit sets
            # them afresh (issue #167).
            # A phi taken from a MAP fitted under a front sampling
            # (ISMIP7_WARM_START_FLUIDITY) already belongs to this front and
            # is kept.
            if front_band_extends(same_mesh=same_mesh, geometry_taken=warm_geometry,
                                  run_sampling=raster_sample, warm_sampling=_warm_rs):
                _H_warm_rs = sample_to_geometry(
                    rasterio.open(f"netcdf:{bm_fn}:thickness"), Q_g, Q,
                    floor=h_clamp, method=_warm_rs)
                _front_nodes = front_changed_nodes(H, _H_warm_rs, Q)
                _band = front_band_controls(_fl_rs if _fl_chk else None)
                for _name in _band:
                    harmonic_extension({"theta": theta, "phi": phi}[_name],
                                       _front_nodes, 0.0, COMM_WORLD)
                PETSc.Sys.Print(
                    f"    front band: {' and '.join(_band)} continued harmonically on "
                    f"{global_count(_front_nodes, COMM_WORLD)} nodes of the cells "
                    f"the front rule changed, between the ice upstream and the "
                    f"warm start seaward (ISMIP7_WARM_START_FRONT_EXTEND=1)"
                    + ("; phi is the fluidity MAP's, fitted under "
                       f"{_fl_rs}" if "phi" not in _band else ""))
            if same_mesh and _warm_lake != int(LAKE_ICE_BASE):
                PETSc.Sys.Print(
                    f"    warm start records lake_ice_base={_warm_lake}, this run "
                    f"{int(LAKE_ICE_BASE)}: its geometry is "
                    + ("taken anyway (ISMIP7_WARM_START_GEOMETRY)" if warm_geometry
                       else "not taken"))
            PETSc.Sys.Print(
                "    warm start is on " + ("this mesh" if same_mesh else "another mesh")
                + ("; taking its geometry, velocity_obs and state"
                   if warm_geometry else
                   "; taking its controls and fluidity prior only "
                   "(geometry and velocity_obs are this mesh's own)"))
            if not warm_geometry:
                raise_geometry = KeyError("warm start geometry not taken")
            warm_geometry_loaded = False
            try:
                if not warm_geometry:
                    raise raise_geometry
                H.assign(_warm_load(chk, chk_mesh, "thickness", Q_g))
                b.assign(_warm_load(chk, chk_mesh, "bed", Q_g))
                s.assign(_warm_load(chk, chk_mesh, "surface", Q_g))
                warm_geometry_loaded = True
                PETSc.Sys.Print(
                    "    geometry: thickness/bed/surface from warm start"
                )
            except (KeyError, RuntimeError, ValueError):
                PETSc.Sys.Print(
                    "    geometry: keeping BedMachine sample "
                    "(warm start has no thickness/bed/surface)"
                )
            # A relaxed geometry (icepack2_tools.relaxation) comes with the
            # controls fitted to it, so it is taken or the run stops. From a
            # relaxation's end state the re-inversion keeps the MAP's theta
            # (IU, 6 Oct 2026, issue #162): the friction then follows the relaxed
            # driving stress through the anchor, and where basal drag carries
            # that stress the starting speed stays the MAP's.
            if warm_end_state:
                _why = end_state_problems(warm_geometry_attrs, same_mesh=same_mesh,
                                          geometry_taken=warm_geometry_loaded,
                                          map_out=map_fn)
                if _theta_mode == "physical":
                    _why.append(
                        "ISMIP7_WARM_START_THETA=physical rebases between anchor "
                        "lengths on this run's geometry; a re-inversion keeps the "
                        "MAP's theta, so leave it unset")
                if _why:
                    raise RuntimeError(
                        f"ISMIP7_WARM_START={warm_chk} is a relaxation's end state "
                        f"that cannot seed this re-inversion: " + "; ".join(_why))
                if _warm_fc == "log" and _theta_mode == "1":
                    try:
                        warm_C_w0 = _warm_load(chk, chk_mesh, "C_w0", Q_g)
                    except (KeyError, RuntimeError, ValueError):
                        warm_C_w0 = None
            if is_relaxed(warm_geometry_attrs):
                _sha = COMM_WORLD.bcast(
                    file_sha256(warm_chk) if (warm_end_state and COMM_WORLD.rank == 0)
                    else None, root=0)
                run_geometry = inherited_geometry(
                    warm_geometry_attrs, geometry_taken=warm_geometry_loaded,
                    warm_basename=os.path.basename(warm_chk), warm_sha256=_sha)
                PETSc.Sys.Print("    relaxed geometry: " + describe_relaxation(run_geometry))
            try:
                if not warm_geometry:
                    raise raise_geometry
                u_obs.assign(_warm_load(chk, chk_mesh, "velocity_obs", V))
                PETSc.Sys.Print("    velocity_obs from warm start")
            except (KeyError, RuntimeError, ValueError):
                pass
            # ISMIP7_WARM_START_PRIOR=0 recomputes the thermal prior instead,
            # e.g. after a change to its physics; phi is kept, as a deviation
            # from the new prior mean.
            if _phi_mode == "physical":
                warm_A_prior_rebase = (_fl_prior if _fl_prior is not None else
                                       _warm_load(chk, chk_mesh, "fluidity_prior", Q))
                warm_A_prior = None
                PETSc.Sys.Print("    fluidity_prior: recomputed below; log_fluidity "
                                "rebased onto it (ISMIP7_WARM_START_PHI=physical)")
            elif os.environ.get("ISMIP7_WARM_START_PRIOR", "1").strip() == "0":
                warm_A_prior = None
                PETSc.Sys.Print("    fluidity_prior: recomputed below "
                                "(ISMIP7_WARM_START_PRIOR=0); log_fluidity kept")
            elif _fl_prior is not None:
                warm_A_prior = _fl_prior
                warm_prior_origin = _fl_origin
            else:
                try:
                    warm_A_prior = _warm_load(
                        chk, chk_mesh, "fluidity_prior", Q
                    )
                    warm_prior_origin = (
                        str(chk.get_attr("/", "fluidity_prior_origin"))
                        if chk.has_attr("/", "fluidity_prior_origin")
                        else f"warm start {os.path.basename(warm_chk)}")
                except (KeyError, RuntimeError, ValueError):
                    warm_A_prior = None
            _state_load, warm_state_guess = warm_start_state(
                geometry_taken=warm_geometry, same_mesh=same_mesh)
            try:
                if _fl_state is not None:
                    u_ws, M_ws, tau_ws = _fl_state
                    warm_state_guess = True
                    _state_from = f"{_fl_chk} (ISMIP7_WARM_START_STATE=fluidity)"
                elif not _state_load:
                    raise raise_geometry
                else:
                    u_ws, M_ws, tau_ws = _warm_state(chk, chk_mesh)
                    _state_from = ("warm start (ISMIP7_WARM_START_STATE=1)"
                                   if warm_state_guess else "warm start")
                z.subfunctions[0].assign(u_ws)
                z.subfunctions[1].assign(M_ws)
                z.subfunctions[2].assign(tau_ws)
                warm_loaded_z = True
                PETSc.Sys.Print(
                    f"    mixed state: velocity/membrane/basal from {_state_from}"
                    + (", the first guess on this run's geometry"
                       if warm_state_guess else ""))
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
    C_w0 = weertman_anchor(H, s, u_obs, m_slide_val, Q_g, length=ANCHOR_LENGTH, b=b)
    if warm_C_w0 is not None:
        # How far the relaxed geometry moved the anchor from the one the
        # relaxation year ran with: R = C_w0 here / C_w0 there, on grounded
        # ice. Under the kept theta the first evaluation's friction is the
        # MAP's times R on these cells.
        _haf = Function(Q_g).interpolate(height_above_flotation(H, b))
        _n, _n1, _n3, _sum, _lo, _hi = anchor_ratio_counts(
            C_w0.dat.data_ro, warm_C_w0.dat.data_ro, _haf.dat.data_ro > 0.0)
        _n, _n1, _n3 = (COMM_WORLD.allreduce(v) for v in (_n, _n1, _n3))
        _sum = COMM_WORLD.allreduce(_sum)
        _lo, _hi = min(COMM_WORLD.allgather(_lo)), max(COMM_WORLD.allgather(_hi))
        PETSc.Sys.Print(
            f"  Anchor on the relaxed geometry: ln R over {_n} grounded dofs in "
            f"[{_lo:.3f}, {_hi:.3f}], mean |ln R| {_sum / max(_n, 1):.4f}, "
            f"|ln R| > 0.1 on {_n1} ({100.0 * _n1 / max(_n, 1):.2f} %), "
            f"> 0.3 on {_n3} ({100.0 * _n3 / max(_n, 1):.2f} %)")
    PETSc.Sys.Print(
        "  Friction anchor: " + ("local driving stress" if ANCHOR_LENGTH == 0.0 else
                                 f"grounded driving stress averaged over {ANCHOR_LENGTH / 1e3:g} km"))
    if warm_theta_anchor is not None and not (warm_exp_alpha and FRICTION_CONTROL == "exp"):
        # The warm start's anchor, rebuilt on this run's geometry, and its theta
        # moved onto this run's anchor: the friction the first solve sees is the
        # warm start's wherever the geometry is the same. An exp-control MAP's
        # anchor is its constant C_ref; an exp-control run keeps it as is.
        if warm_exp_alpha:
            C_prev = Function(Q_g).assign(Constant(float(warm_attrs["friction_c_ref"])))
            _from = "the exp control's C_ref"
        else:
            C_prev = weertman_anchor(H, s, u_obs, m_slide_val, Q_g,
                                     length=warm_theta_anchor, b=b)
            _from = f"a {warm_theta_anchor / 1e3:g} km anchor"
        theta_prev = theta.copy(deepcopy=True)
        theta.assign(rebase_log_friction(theta_prev, C_prev, C_w0, H, b))
        warm_exp_alpha = False
        _shift = Function(Q).interpolate(abs(theta - theta_prev))
        PETSc.Sys.Print(
            f"    log_friction rebased from {_from} onto "
            f"this run's: |shift| mean {global_mean(_shift):.3f}, max {global_max(_shift):.2f}")
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
            f"  Friction: {law_name}; h_visc_floor={RC_HVISC_FLOOR:g} m; "
            f"C_w0 in [{C_w0_lo:.2e}, {C_w0_hi:.2e}]"
        )
    else:
        PETSc.Sys.Print("  Friction: Budd power-law dual (legacy action)")

    # ── sqrt(C) friction control (ISMIP7_FRICTION_CONTROL=sqrt) ───────────
    # theta now IS alpha = sqrt(C). It starts from the friction the anchor
    # (and any warm-started log deviation, rebased above) describes, or from
    # the alpha of a warm start inverted on this control. The prior mean is
    # zero, so nothing but the start remembers the anchor.
    sigma_alpha_val = None
    if (FRICTION_CONTROL in ("sqrt", "exp") and PRIOR_SIGMA_ALPHA == "auto"
            and "prior_sigma_alpha" in warm_attrs
            and frozen_in_control(warm_attrs, "prior_sigma_alpha", FRICTION_CONTROL) is None):
        PETSc.Sys.Print(
            f"    sigma_alpha of the warm start ({warm_attrs['prior_sigma_alpha']}) is in the "
            f"{warm_attrs.get('friction_control', 'log')} control's units; not reused")
    if FRICTION_CONTROL == "sqrt":
        C_cg = cg1_lift(C_w0) if geom_dg else C_w0
        if warm_alpha is not None:
            theta.assign(warm_alpha)
        elif warm_exp_alpha:
            theta.interpolate(sqrt(Constant(float(warm_attrs["friction_c_ref"])) * exp(theta)))
        else:
            theta.interpolate(sqrt(max_value(C_cg * exp(theta), Constant(0.0))))
        theta.rename("alpha")
        H_cg = cg1_lift(H) if geom_dg else H
        b_cg = cg1_lift(b) if geom_dg else b
        _haf = Function(Q).interpolate(height_above_flotation(H_cg, b_cg)).dat.data_ro
        _grounded = (_haf > 0.0) & (H_cg.dat.data_ro > 10.0)
        _alpha_grounded = np.concatenate(COMM_WORLD.allgather(
            np.asarray(theta.dat.data_ro[_grounded], dtype=float)))
        _alpha_median = (float(np.median(_alpha_grounded))
                         if _alpha_grounded.size else float("nan"))
        # frozen within one control: re-deriving it from the warm-started
        # alpha would change the objective between links
        _warm_sigma = frozen_in_control(warm_attrs, "prior_sigma_alpha", FRICTION_CONTROL)
        if PRIOR_SIGMA_ALPHA == "auto" and _warm_sigma is not None:
            sigma_alpha_val = _warm_sigma
            PETSc.Sys.Print(f"    sigma_alpha frozen from the warm start: {sigma_alpha_val:.4e}")
        elif PRIOR_SIGMA_ALPHA == "auto":
            sigma_alpha_val = _alpha_median
        else:
            sigma_alpha_val = float(PRIOR_SIGMA_ALPHA)
        if not sigma_alpha_val > 0.0:
            raise ValueError(
                f"the prior sigma for alpha must be positive, got {sigma_alpha_val}")
        _a_lo, _a_hi = global_range(theta)
        PETSc.Sys.Print(
            f"  Friction control: alpha = sqrt(C) with a zero prior mean "
            f"(Recinos et al. 2023); initial alpha in [{_a_lo:.3e}, {_a_hi:.3e}], "
            f"grounded median {_alpha_median:.3e}; sigma_alpha={sigma_alpha_val:.3e} "
            f"({PRIOR_SIGMA_ALPHA}), rho_theta={PRIOR_RHO_THETA:g} m")
    elif FRICTION_CONTROL == "exp":
        # theta now IS alpha = ln(C / C_ref). It starts from the friction the
        # loaded state describes (the anchor with any log deviation, the
        # alpha^2 of a sqrt-control MAP, or the alpha of an exp-control MAP,
        # which is already this control) and the prior mean is zero.
        C_cg = cg1_lift(C_w0) if geom_dg else C_w0
        if warm_exp_alpha:
            C_start = Function(Q).interpolate(
                Constant(float(warm_attrs["friction_c_ref"])) * exp(theta))
        elif warm_alpha is not None:
            C_start = Function(Q).interpolate(warm_alpha ** 2)
        else:
            C_start = Function(Q).interpolate(C_cg * exp(theta))
        H_cg = cg1_lift(H) if geom_dg else H
        b_cg = cg1_lift(b) if geom_dg else b
        _haf = Function(Q).interpolate(height_above_flotation(H_cg, b_cg)).dat.data_ro
        _grounded = (_haf > 0.0) & (H_cg.dat.data_ro > 10.0)
        _C_grounded = np.concatenate(COMM_WORLD.allgather(
            np.asarray(C_start.dat.data_ro[_grounded], dtype=float)))
        _C_median = (float(np.median(_C_grounded)) if _C_grounded.size else float("nan"))
        _warm_c_ref = frozen_in_control(warm_attrs, "friction_c_ref", FRICTION_CONTROL)
        if _warm_c_ref is not None:
            # frozen: the control is a deviation from THIS reference
            c_ref_val = _warm_c_ref
            PETSc.Sys.Print(f"    C_ref frozen from the warm start: {c_ref_val:.4e}")
        elif C_REF == "auto":
            c_ref_val = _C_median
        else:
            c_ref_val = float(C_REF)
        if not c_ref_val > 0.0:
            raise ValueError(f"the friction reference C_ref must be positive, got {c_ref_val}")
        # ln floors at 1e-4 C_ref (alpha >= -9.2): the anchor is ~0 under
        # some grounded nodes, and -inf is not a start. Floating and ice-free
        # nodes, where the friction is inert, start at the prior mean. A chain
        # link from an exp-control warm start resumes its alpha unchanged.
        if not warm_exp_alpha:
            theta.interpolate(ln(max_value(C_start, Constant(1e-4 * c_ref_val)) / Constant(c_ref_val)))
            theta.dat.data[~_grounded] = 0.0
        theta.rename("alpha")
        _warm_sigma = frozen_in_control(warm_attrs, "prior_sigma_alpha", FRICTION_CONTROL)
        if PRIOR_SIGMA_ALPHA == "auto" and _warm_sigma is not None:
            sigma_alpha_val = _warm_sigma
            PETSc.Sys.Print(f"    sigma_alpha frozen from the warm start: {sigma_alpha_val:.4e}")
        elif PRIOR_SIGMA_ALPHA == "auto":
            sigma_alpha_val = 1.0
        else:
            sigma_alpha_val = float(PRIOR_SIGMA_ALPHA)
        if not sigma_alpha_val > 0.0:
            raise ValueError(
                f"the prior sigma for alpha must be positive, got {sigma_alpha_val}")
        _a_lo, _a_hi = global_range(theta)
        PETSc.Sys.Print(
            f"  Friction control: C = C_ref exp(alpha) with a zero prior mean on alpha; "
            f"C_ref={c_ref_val:.3e} ({C_REF}; grounded median of the start {_C_median:.3e}), "
            f"initial alpha in [{_a_lo:.3f}, {_a_hi:.3f}]; sigma_alpha={sigma_alpha_val:.3g} "
            f"log units ({PRIOR_SIGMA_ALPHA}), rho_theta={PRIOR_RHO_THETA:g} m")
    elif warm_exp_alpha:
        # an exp-control MAP warm-starting a log-control run: its alpha is a
        # deviation from C_ref; rebase it onto THIS anchor
        C_cg = cg1_lift(C_w0) if geom_dg else C_w0
        theta.interpolate(ln(max_value(Constant(float(warm_attrs["friction_c_ref"])) * exp(theta),
                                       Constant(1e-12))
                             / max_value(C_cg, Constant(1e-12))))
        PETSc.Sys.Print("    log_friction rebuilt from the warm start's exp control: "
                        "theta = ln(C_ref exp(alpha) / C_w0)")
    elif warm_alpha is not None:
        # a sqrt-control MAP warm-starting a log-control run: the deviation
        # from THIS anchor that reproduces its friction on grounded ice
        C_cg = cg1_lift(C_w0) if geom_dg else C_w0
        theta.interpolate(ln(max_value(warm_alpha ** 2, Constant(1e-12))
                             / max_value(C_cg, Constant(1e-12))))
        PETSc.Sys.Print("    log_friction rebuilt from the warm start's alpha: "
                        "theta = ln(alpha^2 / C_w0)")

    # Physical FLUIDITY PRIOR MEAN A_prior(x): a fixed-velocity thermomechanical
    # (Stefan enthalpy) solve at the observed geometry/velocity, so the control
    # phi = log(A / A_prior) is a small deviation from a physically-motivated
    # fluidity rather than log(A / const). This is the Recinos et al. (2023) fix
    # for the n=3 blow-up (a constant A0 baseline forced phi to carry all the
    # spatial fluidity structure). Frictional heating uses the balance C_w0.
    # When warm-starting from a prepare cache / MAP that already carries
    # fluidity_prior, reuse it: phi = log(A/A_prior) is meaningless against a
    # freshly recomputed prior.
    prior_origin = "constant A0*a4_factor"
    if warm_A_prior is not None:
        A_prior = warm_A_prior
        prior_origin = warm_prior_origin or "warm start"
        if FLUIDITY_PRIOR != "thermo":
            PETSc.Sys.Print(
                f"    WARNING: ISMIP7_FLUIDITY_PRIOR={FLUIDITY_PRIOR} but the warm "
                "start's fluidity_prior is kept (ISMIP7_WARM_START_PRIOR=0 recomputes)")
        A_prior.rename("fluidity_prior")
        A_prior_lo, A_prior_hi = global_range(A_prior)
        PETSc.Sys.Print(
            f"  Fluidity prior A_prior in [{A_prior_lo:.2f}, "
            f"{A_prior_hi:.2f}] (from warm start)"
        )
    elif FLUIDITY_PRIOR == "pattyn":
        from icepack2_tools.rheology_prior import (
            column_fill_temperature, fluidity_prior_from_temperature_raster)
        _pattyn_fn = os.environ.get(
            "ISMIP7_PATTYN_TEMP", os.path.join(DATA_DIR, "temp", "Pattyn_2013.tif"))
        if not os.path.exists(_pattyn_fn):
            raise FileNotFoundError(
                f"ISMIP7_FLUIDITY_PRIOR=pattyn: temperature raster {_pattyn_fn} "
                "not found (ISMIP7_PATTYN_TEMP names it)")
        # The raster is defined only under BedMachine ice (outside it is an
        # extrapolation halo); nodes it does not define take a column
        # temperature from the surface-temperature climatology and the
        # flotation state, on the same CG1 geometry the thermal prior uses.
        _H_cg = cg1_lift(H) if geom_dg else H
        _b_cg = cg1_lift(b) if geom_dg else b
        _T_fill = column_fill_temperature(
            load_mean_annual_surface_temperature(Q), _H_cg, _b_cg, Q)
        A_prior, _pinfo = fluidity_prior_from_temperature_raster(
            _pattyn_fn, Q, bm_fn=bm_fn, fill=_T_fill)
        _ice_node = _H_cg.dat.data_ro > 10.0
        _n_fill = COMM_WORLD.allreduce(int(_pinfo["nodes_filled"]))
        _n_fill_ice = COMM_WORLD.allreduce(int((_pinfo["filled"] & _ice_node).sum()))
        _n_ice = COMM_WORLD.allreduce(int(_ice_node.sum()))
        _n_tot = COMM_WORLD.allreduce(int(_pinfo["nodes_total"]))
        prior_origin = (
            "rate_factor of the depth-averaged Pattyn temperature "
            f"({os.path.basename(_pattyn_fn)}) under BedMachine ice, column "
            "temperature from the surface climatology where it is undefined; "
            "no strain heating, geothermal flux or water content of our own")
        A_prior_lo, A_prior_hi = global_range(A_prior)
        PETSc.Sys.Print(
            f"  Fluidity prior A_prior in [{A_prior_lo:.2f}, {A_prior_hi:.2f}] "
            f"(Pattyn depth-averaged temperature {_pattyn_fn}: "
            f"{_pinfo['pixels_replaced']} of {_pinfo['pixels_total']} pixels are "
            "undefined - outside BedMachine ice, missing or outside [200, 273.15] K; "
            f"{_n_fill} of {_n_tot} nodes took the surface-climatology column "
            f"temperature, {_n_fill_ice} of the {_n_ice} nodes with more than 10 m of ice)")
    elif FLUIDITY_PRIOR == "thermo":
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
        prior_origin = ("thermomechanical; over water the base is at the ice-ocean melting point "
                        "and supplies nothing more (no friction, no geothermal flux, no water content)")
        A_prior.rename("fluidity_prior")
        A_prior_lo, A_prior_hi = global_range(A_prior)
        PETSc.Sys.Print(
            f"  Fluidity prior A_prior in [{A_prior_lo:.2f}, "
            f"{A_prior_hi:.2f}] (thermomechanical)"
        )
    else:
        A_prior = Function(Q, name="fluidity_prior").interpolate(A0 * Constant(a4_factor))
        PETSc.Sys.Print("  Fluidity prior: constant A0*a4_factor (legacy)")
    if warm_A_prior_rebase is not None:
        _phi_prev = phi.copy(deepcopy=True)
        phi.interpolate(phi + ln(max_value(warm_A_prior_rebase, Constant(1e-30))
                                 / max_value(A_prior, Constant(1e-30))))
        _shift = Function(Q).interpolate(phi - _phi_prev)
        _s_lo, _s_hi = global_range(_shift)
        PETSc.Sys.Print(
            f"    log_fluidity rebased onto this prior (A kept): shift in "
            f"[{_s_lo:.2f}, {_s_hi:.2f}], mean {global_mean(_shift):.3f}")
    A4_base = A_prior
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
    # ISMIP7_DRAG_GATE=vertex also keeps it off water cells that touch the
    # ice at a vertex, whose drag acts on the ice's own front nodes.
    DRAG_GATE = drag_gate()
    _neighbours = (vertex_neighbours if DRAG_GATE == "vertex" else facet_neighbours)
    drag_mask.dat.data[:] = ocean_drag_cells(
        _ice, _neighbours(drag_mask.function_space()), _ice)
    _n_drag = COMM_WORLD.allreduce(int(drag_mask.dat.data_ro.sum()))
    PETSc.Sys.Print(
        f"  Ocean drag gate: {_n_drag} "
        f"of {COMM_WORLD.allreduce(int(drag_mask.dat.data_ro.size))} cells, open water "
        f"a cell away from the ice (h < {front_hmin():g} m, ISMIP7_DRAG_GATE={DRAG_GATE})"
    )
    # The gate the controls absorbed, which a forward then runs
    # (runconfig.forward_drag_gate): none when no cell was dragged.
    DRAG_RECORD = (DRAG_GATE if _n_drag and stabilizers["ocean_drag"] > 0.0
                   else DRAG_GATE_NONE)

    subelement = None
    if SUBELEMENT_FRICTION:
        from icepack2_tools.subelement import ice_indicator, subelement_from_geometry
        subelement = subelement_from_geometry(
            mesh, H, b, ice=ice_indicator(H, front_hmin()))
        _fr = subelement.fraction.dat.data_ro
        _n_part = COMM_WORLD.allreduce(int(((_fr > 0.0) & (_fr < 1.0)).sum()))
        _n_full = COMM_WORLD.allreduce(int((_fr == 1.0).sum()))
        PETSc.Sys.Print(
            f"  Sub-element grounding (ISSM {SUBELEMENT_SCHEME.upper()} version "
            f"{SUBELEMENT_SCHEME_VERSION}, icepack_tools): {_n_full} cells fully "
            f"grounded, {_n_part} partly grounded; {FRICTION} runs with N_hat = 1 on the "
            f"grounded part (no N_ref, no delta floor); exact front push "
            f"{f'version {EXACT_FRONT}' if EXACT_FRONT else 'off'}")

    _zero_theta = Constant(0.0)
    # the smooth grounded indicator of THIS geometry, for a floating-only
    # fluidity control (UFL on the cell fields; the residual builds its own
    # copy for the friction gate)
    _He_phi = grounded_mask(H, b) if FLUIDITY_CONTROL == "floating" else None
    # The nodes phi may move on; the rest hold phi = 0 (see FLUIDITY_CONTROL).
    phi_free = None
    if _He_phi is not None:
        phi_free = floating_control_nodes(H, b, Q)
        phi.dat.data[:] = phi.dat.data_ro * phi_free.dat.data_ro
        _n_free = COMM_WORLD.allreduce(int(phi_free.dat.data_ro.sum()))
        _n_all = COMM_WORLD.allreduce(int(phi_free.dat.data_ro.size))
        PETSc.Sys.Print(
            f"  Fluidity control: phi acts on floating ice only and is held at 0 "
            f"on {_n_all - _n_free} of {_n_all} nodes (grounded ice keeps the "
            f"prior fluidity)")

    def build_F(theta_c, phi_c, *, scpc_blocks=state_scpc):
        F_c = _build_residual(theta_c, phi_c)
        return with_scpc_blocks(F_c, z) if scpc_blocks else F_c

    def _build_residual(theta_c, phi_c):
        # Residual closure (tau linear, grounded-only theta via exp(theta*He),
        # exact-zero shelves): budd -> N_hat=1 at the reference geometry;
        # regularized_coulomb -> Coulomb cap. Legacy budd -> action derivative.
        if USE_RESIDUAL:
            if FRICTION_CONTROL == "sqrt":
                # C = alpha^2 outright; the He-gated log deviation is zero
                theta_arg, C_arg = _zero_theta, theta_c ** 2
            elif FRICTION_CONTROL == "exp":
                # C = C_ref exp(alpha) outright, as the sqrt control passes
                # alpha^2: the residual's He gate on a log deviation would
                # otherwise pull the friction to C_ref across the grounding
                # band, which is the anchor's role under `log` and no one's here
                theta_arg, C_arg = _zero_theta, Constant(c_ref_val) * exp(theta_c)
            else:
                theta_arg, C_arg = theta_c, C_w0
            if _He_phi is not None:
                phi_c = phi_c * (Constant(1.0) - _He_phi)
            if subelement is not None:
                from icepack2_tools.subelement import build_subelement_residual
                return build_subelement_residual(
                    z, theta_arg, phi_c, H=H, s=s, b=b, C_w0=C_arg,
                    A4_base=A4_base, n_flow=n_flow, n_flow_val=n_flow_val,
                    m_slide=m_slide_val, tau_c=tau_c, alpha=alpha_reg, H_ref=H_ref,
                    subelement=subelement, scheme=SUBELEMENT_SCHEME,
                    fric_law=FRICTION, nhat_cap=BUDD_NHAT_CAP,
                    alpha_gl=ALPHA_GL, c_w0_floor=RC_CW0_FLOOR,
                    h_visc_floor=RC_HVISC_FLOOR, k_lim=0.0, **stabilizers,
                    drag_mask=drag_mask,
                    calving_ids=calving_ids if use_calving_terminus else None,
                    exact_front=EXACT_FRONT,
                )
            return build_rc_residual(
                z, theta_arg, phi_c, H=H, s=s, b=b, C_w0=C_arg,
                A4_base=A4_base, n_flow=n_flow, n_flow_val=n_flow_val,
                m_slide=m_slide_val, tau_c=tau_c, alpha=alpha_reg, H_ref=H_ref,
                fric_law=FRICTION, N_ref=None,
                nhat_floor=BUDD_DELTA, nhat_cap=BUDD_NHAT_CAP, alpha_gl=ALPHA_GL,
                c0=C0_RC, c_w0_floor=RC_CW0_FLOOR, h_visc_floor=RC_HVISC_FLOOR,
                k_lim=0.0, **stabilizers, drag_mask=drag_mask,
                calving_ids=calving_ids if use_calving_terminus else None,
                exact_front=EXACT_FRONT, front_hmin=front_hmin(),
            )
        return derivative(_build_action(theta_c, phi_c, fields), z)

    if EXACT_FRONT and not SUBELEMENT_FRICTION:
        PETSc.Sys.Print(
            f"  Exact cliff push on internal fronts (ISMIP7_EXACT_FRONT={EXACT_FRONT}, "
            "dual_friction.front_cliff_correction: "
            f"{'the free-cliff push' if EXACT_FRONT == 1 else 'the face above the neighbour bed'})")
    if use_calving_terminus:
        PETSc.Sys.Print("  Using calving_terminus BC")
    else:
        PETSc.Sys.Print("  NO calving_terminus BC (buffered mesh, h=0 at front)")
    F = build_F(theta, phi)
    # The startup ramp's residual, with the blocks its own solver wants.
    F_ramp = (F if lane_scpc == state_scpc
              else build_F(theta, phi, scpc_blocks=lane_scpc))

    def _untaped_state_solver(F_form, params):
        """An unannotated solve of F_form under the taped solve's mode: under
        scpc_* with the Jacobian frozen at the Newton iterate, as the taped
        solve and the transient run it (preconditioners.frozen_linearization;
        tlm_adjoint refuses its callback while annotating, so call this with
        the manager stopped)."""
        J_form, pre_jacobian = None, None
        if state_solver_mode.startswith("scpc_"):
            J_form, pre_jacobian = frozen_linearization(F_form, z)
        return NonlinearVariationalSolver(
            NonlinearVariationalProblem(
                F_form, z, J=J_form, form_compiler_parameters=fc_params),
            solver_parameters=params,
            pre_jacobian_callback=pre_jacobian,
        )

    # ── Warm start ──
    stop_manager()
    slvr = _untaped_state_solver(F, sparams)
    n_flow.assign(n_flow_val)
    m_slide.assign(m_slide_val)
    ramp_params = diagnostic_solver_parameters(lane_solver_mode)
    ramp_params.update(_monitor_options)
    ramp_J, ramp_pre_jacobian = None, None
    if linearization_state(lane_solver_mode) == "frozen":
        ramp_J, ramp_pre_jacobian = frozen_linearization(F_ramp, z)
    ramp_solver = NonlinearVariationalSolver(
        NonlinearVariationalProblem(
            F_ramp, z, J=ramp_J, form_compiler_parameters=fc_params
        ),
        solver_parameters=ramp_params,
        options_prefix="ismip7_inversion_continuation_",
        pre_jacobian_callback=ramp_pre_jacobian,
    )
    ramp_ladder = ladder(continuation_steps())
    # The sliding exponent a ramp starts from: 1, or its target under
    # ISMIP7_RAMP_SLIDE_FIXED=1 (issue #167).
    ramp_m_start = m_slide_val if ramp_slide_fixed() else 1.0
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


    def _reramp_at_current_controls():
        r"""Re-climb the exponent ladder, unannotated, at the controls the
        module-level (theta, phi) hold, from the state z holds: the rescue
        for a line-search trial point where the single Newton solve of the
        annotated forward diverges (the 32 km sum-likelihood run stopped at
        iteration 5 on three such trials). Leaves the exponents at their
        targets whatever happens."""
        try:
            _, steps = ramp_exponents(
                _ramp_solve, z, n_flow, m_slide, n_flow_val, m_slide_val,
                ramp_ladder[:trial_rescue_rungs()], report=PETSc.Sys.Print,
                m_start=ramp_m_start)
        finally:
            n_flow.assign(n_flow_val)
            m_slide.assign(m_slide_val)
        return steps

    if skip_continuation:
        if warm_loaded_z and warm_state_guess:
            PETSc.Sys.Print(
                f"Warm start: the loaded mixed state is the first guess at full "
                f"n_flow={n_flow_val:.1f}, m_slide={m_slide_val:.1f}; the first "
                f"evaluation's forward solves it on this geometry (no 1→n "
                f"continuation)"
            )
        elif warm_loaded_z:
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
        # transient's step ladder. The selected inversion solver takes over at
        # the converged state. On the 1 km production mesh the full-Jacobian
        # MUMPS factorisation failed mid-ramp
        # (DIVERGED_LINEAR_SOLVE, job 1612624) where the forward's condensed
        # GAMG had climbed the same ramp on the same mesh and MAP.
        PETSc.Sys.Print(
            f"Warm start (continuation n_flow 1→{n_flow_val:.1f}, "
            f"m_slide {ramp_m_start:g}→{m_slide_val:.1f}; "
            f"{diagnostic_solver_label(lane_solver_mode)}, "
            f"{linearization_state(lane_solver_mode)} linearization, "
            f"steps {'/'.join(str(s) for s in ramp_ladder)})..."
        )

        _, ramp_steps = ramp_exponents(
            _ramp_solve, z, n_flow, m_slide, n_flow_val, m_slide_val,
            ramp_ladder, report=PETSc.Sys.Print, m_start=ramp_m_start,
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
    # The same rule drops every SNES option: under scpc_* tlm_adjoint hands
    # these to a LinearVariationalSolver, which a forward's newtonls would
    # turn into a line-searched Newton solve (inversion_adjoint_parameters).
    adjoint_sparams = inversion_adjoint_parameters(sparams)

    u_init = z.subfunctions[0]
    u_mag = Function(Q).interpolate(sqrt(inner(u_init, u_init)))
    PETSc.Sys.Print(f"  u_max = {global_max(u_mag):.0f} m/yr")

    # ── Forward function for tlm_adjoint ──
    # Normalize by the OBSERVED area so the misfit magnitude stays
    # comparable between masked and unmasked runs.
    area_val = assemble(obs_mask * dx(mesh))
    # ISMIP7_MISFIT_SCALE: the multiplier that turns the mean chi^2 into the
    # likelihood the prior is weighed against (see the knob's comment).
    if MISFIT_SCALE == "nodes":
        misfit_scale = float(COMM_WORLD.allreduce(
            int((obs_mask.dat.data_ro > 0.5).sum())))
    else:
        misfit_scale = float(MISFIT_SCALE)
    if not misfit_scale > 0.0:
        raise ValueError(f"ISMIP7_MISFIT_SCALE resolved to {misfit_scale}, not positive")
    PETSc.Sys.Print(
        f"  Misfit scale: {misfit_scale:g} x the mean chi^2 "
        f"({'sum over observed nodes' if MISFIT_SCALE == 'nodes' else 'ISMIP7_MISFIT_SCALE'})")

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
        _sigma_theta = (sigma_alpha_val if FRICTION_CONTROL in ("sqrt", "exp")
                        else PRIOR_SIGMA_THETA)
        _prior_dg = {
            "theta": bilaplacian_coeffs(_sigma_theta, PRIOR_RHO_THETA),
            "phi": bilaplacian_coeffs(PRIOR_SIGMA_PHI, PRIOR_RHO),
        }
        PETSc.Sys.Print(
            f"  Prior: bi-Laplacian A M^-1 A; "
            f"sigma_{'alpha' if FRICTION_CONTROL == 'sqrt' else 'theta'}="
            f"{_sigma_theta:g} rho_theta={PRIOR_RHO_THETA:g} m, "
            f"sigma_phi={PRIOR_SIGMA_PHI:g} rho={PRIOR_RHO:g} m -> "
            f"theta delta={_prior_dg['theta'][0]:.4e} gamma={_prior_dg['theta'][1]:.4e}"
        )
    else:
        _prior_dg = {
            "theta": prior_operator_coeffs(GAMMA_THETA, area_val, L_REG),
            "phi": prior_operator_coeffs(GAMMA_PHI, area_val, L_REG),
        }
    _prior_aux = {k: Function(Q, name=f"prior_aux_{k}") for k in ("theta", "phi")}
    _prior_aux_solver = BilaplacianAuxSolver(Q, form_compiler_parameters=fc_params)

    def _prior_energy_form(ctrl, which):
        """``R(ctrl)`` as something ``assemble`` or ``Functional.addto`` takes.

        Under `bilaplacian` this SOLVES ``M f = A ctrl`` into
        ``_prior_aux[which]`` as a side effect, so :func:`_prior_grad` -- which
        needs that ``f`` -- must be called after this and before the next
        control changes. The solve back-substitutes with M factored once
        (prior.BilaplacianAuxSolver). It is recorded while a manager annotates
        (the TAO path needs it on the tape for the gradient) and is an
        ordinary solve when one does not (the scipy path differentiates it by
        hand below, and the TAO monitor only reads the energy).
        """
        d, g = _prior_dg[which]
        if PRIOR_FORM == "laplacian":
            return 0.5 * prior_operator_form(ctrl, ctrl, d, g)
        aux = _prior_aux[which]
        _prior_aux_solver(ctrl, aux, d, g)
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
            # a block passed as None is a frozen control (ISMIP7_INVERT): no
            # solve, no entry in the result
            out = []
            for which, rhs in (("theta", g_theta), ("phi", g_phi)):
                if rhs is None:
                    continue
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
    log_vel_w = (0.0 if LOG_VEL_WEIGHT.lower() == "auto"
                 else float(LOG_VEL_WEIGHT))
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

    # How each objective evaluation reaches its state (README, "Inversion
    # solver"):
    # - direct (ISMIP7_DIRECT_FORWARD, the default): one untaped Newton solve
    #   at the full exponents from the last converged state, then a taped
    #   solve that starts converged; a lost trial goes to the caller's rescue
    #   or backtrack. The taped 5-stage ladder restarted every evaluation at
    #   n=1 with ||F|| 1e11-1e13 and failed at the same trial points, at up to
    #   an hour a try on the 2 km mesh.
    # - single: one taped solve at the full exponents, for a timing lane that
    #   asks for it by name (ISMIP7_SKIP_CONTINUATION=1) and, with the direct
    #   forward off, for a warm start that loaded its mixed state or
    #   ISMIP7_EVAL_CONTINUATION=0.
    # - ladder: the taped 5-stage n,m: 1->n ladder in every evaluation.
    _skip_by_name = os.environ.get("ISMIP7_SKIP_CONTINUATION", "0").strip() == "1"
    if direct_forward_enabled() and not _skip_by_name:
        eval_mode = "direct"
    elif skip_continuation or not eval_continuation():
        eval_mode = "single"
    else:
        eval_mode = "ladder"
    eval_full_n = eval_mode != "ladder"
    PETSc.Sys.Print("  Evaluations: " + {
        "direct": "direct forward, an untaped Newton solve at full n_flow/m_slide "
                  "from the last converged state, then a taped confirmation "
                  "(ISMIP7_DIRECT_FORWARD)",
        "single": "one taped solve at full n_flow/m_slide",
        "ladder": "1->n continuation in 5 taped solves",
    }[eval_mode])
    direct_sparams = direct_forward_parameters(sparams)

    # The untaped work of the last evaluation's state solves (taped_state_solve);
    # written into the timing record per evaluation.
    state_work = []
    # The paused Newton solver of the last taped form, kept while the forward
    # keeps solving that form (StateSolverCache): under the direct forward for
    # every inversion solver, otherwise under scpc_*.
    state_solver_cache = StateSolverCache()

    def _taped_state_solve(F_ctrl, *, direct=False):
        work = taped_state_solve(
            F_ctrl, z, state_solver_mode,
            direct_sparams if direct else sparams, adjoint_sparams,
            form_compiler_parameters=fc_params, cache=state_solver_cache,
            direct=direct,
        )
        state_work.append(work)
        if direct:
            PETSc.Sys.Print(
                f"    direct forward: "
                f"{_snes_reason_name(work['converged_reason'])} "
                f"snes_its={work['snes_iterations']} fnorm={work['fnorm']:.3e} "
                f"{work['seconds']:.1f}s")

    # The taped residual is one form for the run, over Functions that live
    # as long, so the symbolic work Firedrake caches on a form and the solver
    # state_solver_cache keeps for it serve every evaluation. The scipy
    # path's controls are the module-level (theta, phi) F is built over. TAO
    # hands the forward fresh copies of its controls each evaluation;
    # assignments on the tape carry them into (theta_eval, phi_eval), and the
    # gradient reaches the copies through them. This holds because build_F
    # writes the controls into UFL and evaluates nothing from them: a
    # Function interpolated from theta inside build_F would keep the first
    # evaluation's values.
    theta_eval = Function(theta.function_space(), name="theta_eval")
    phi_eval = Function(phi.function_space(), name="phi_eval")
    F_eval = []

    def _residual_at(theta_ctrl, phi_ctrl):
        if theta_ctrl is theta and phi_ctrl is phi:
            return F
        theta_eval.assign(theta_ctrl)
        phi_eval.assign(phi_ctrl)
        if not F_eval:
            F_eval.append(build_F(theta_eval, phi_eval))
        return F_eval[0]

    def forward(theta_ctrl, phi_ctrl):
        clear_caches()
        state_work.clear()
        F_ctrl = _residual_at(theta_ctrl, phi_ctrl)
        if eval_mode == "direct":
            n_flow.assign(n_flow_val)
            m_slide.assign(m_slide_val)
            # A lost trial raises ConvergenceError with z back at the last
            # converged state, for the caller's rescue (one untaped rung of
            # the ladder at these controls, then this solve) or backtrack.
            _taped_state_solve(F_ctrl, direct=True)
        elif eval_mode == "single":
            # Stay at the physical exponents so each eval is one Newton solve.
            n_flow.assign(n_flow_val)
            m_slide.assign(m_slide_val)
            _taped_state_solve(F_ctrl)
        else:
            # Continuation inside annotation for robustness — ramp both
            # n_flow and m_slide on the same [0,1] parameter.
            for t in np.linspace(0.0, 1.0, 5):
                n_flow.assign(1.0 + t * (n_flow_val - 1.0))
                m_slide.assign(1.0 + t * (m_slide_val - 1.0))
                _taped_state_solve(F_ctrl)

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
        J.assign(Constant(misfit_scale) * integrand * dx)
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

    # Under a floating-only fluidity control the optimiser sees phi on the
    # free nodes only: phi is written with the others at zero and its
    # gradient is zero there, so L-BFGS-B never moves them.
    _phi_free_global = func_to_global(phi_free) if phi_free is not None else None

    def set_phi(arr):
        global_to_func(arr, phi)
        if phi_free is not None:
            phi.dat.data[:] = phi.dat.data_ro * phi_free.dat.data_ro

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
    # The objective at the last ACCEPTED iterate (written into every
    # checkpoint), the resolved run settings (idem), and a handoff failure
    # raised after TAO returns (its monitor cannot raise through petsc4py).
    last_accepted = {}
    run_settings = {}
    _handoff_failed = [None]

    def _friction_reference():
        # What a consumer multiplies exp(log_friction) by: the anchor, or
        # under the sqrt control alpha^2 itself (log_friction is saved as 0).
        if FRICTION_CONTROL == "sqrt":
            return Function(Q_g, name="C_w0").project(theta ** 2)
        if FRICTION_CONTROL == "exp":
            # the constant reference: C = C_w0 exp(log_friction) reads exactly
            return Function(Q_g, name="C_w0").interpolate(Constant(c_ref_val))
        return C_w0

    def save_map(path, *, full_state=False):
        """Write the MAP to ``path`` atomically: into a sibling temporary file,
        renamed over ``path`` once every rank has closed it, so a link killed
        mid-write (a wall limit, a preempted scavenge job) leaves the previous
        checkpoint intact for its successor rather than a truncated file."""
        tmp = f"{path}.tmp"
        _write_map(tmp, full_state=full_state)
        COMM_WORLD.Barrier()
        if COMM_WORLD.rank == 0:
            os.replace(tmp, path)
        COMM_WORLD.Barrier()

    def _write_map(path, *, full_state=False):
        with fd.CheckpointFile(path, "w") as chk:
            chk.save_mesh(mesh)
            if FRICTION_CONTROL == "sqrt":
                # the control is alpha = sqrt(C): C_w0 (below) is alpha^2 and
                # the log deviation a consumer applies to it is zero
                chk.save_function(Function(Q, name="log_friction"),
                                  name="log_friction")
                chk.save_function(theta, name="sqrt_friction")
            else:
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
                chk.save_function(_friction_reference(), name="C_w0")
                if N_ref is not None:
                    chk.save_function(N_ref, name="N_ref")
            else:
                # A periodic checkpoint carries the accepted mixed state for
                # the next chain link only, under names no forward reads: the
                # link starts from it instead of re-ramping n=1->3 from rest
                # (~20 min of a one-hour 2 km link, 30 Sep 2026).
                chk.save_function(z.subfunctions[0], name="ckpt_velocity")
                chk.save_function(z.subfunctions[1], name="ckpt_membrane_stress")
                chk.save_function(z.subfunctions[2], name="ckpt_basal_stress")
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
            chk.set_attr("/", "lc", int(mesh_lc))
            chk.set_attr("/", "lc_coarse", int(mesh_lc_coarse))
            chk.set_attr("/", "buffer_m", float(mesh_buffer_m))
            # The configuration theta/phi only mean anything under. The
            # derived MAP filename encodes all three, but ISMIP7_INVERSION
            # bypasses the name, so the forward needs them recorded to check
            # the MAP it was pointed at against the law it is about to run.
            chk.set_attr("/", "friction", str(FRICTION))
            chk.set_attr("/", "n_flow", float(n_flow_val))
            chk.set_attr("/", "geometry_space", str(geometry_space))
            # theta is a log-deviation from THIS anchor, on THIS geometry: a
            # forward rebuilds C_w0 from them, so it takes both from here.
            chk.set_attr("/", "friction_anchor_length", float(ANCHOR_LENGTH))
            chk.set_attr("/", "lake_ice_base", int(LAKE_ICE_BASE))
            # How a warm start from another mesh filled the dofs beyond it
            # (provenance: ISMIP7_TRANSFER_FILL; same-mesh starts fill none).
            chk.set_attr("/", "warm_start_fill", warm_fill_mode)
            # The ocean-drag gate and the membrane floor the controls
            # absorbed: a forward runs both (runconfig.forward_drag_gate,
            # forward_hvisc_floor).
            chk.set_attr("/", "drag_gate", DRAG_RECORD)
            chk.set_attr("/", "h_visc_floor", float(RC_HVISC_FLOOR))
            # Which field the friction is: C_w0 exp(log_friction) on the
            # anchor (log), or sqrt_friction^2 with no anchor (sqrt).
            chk.set_attr("/", "friction_control", FRICTION_CONTROL)
            chk.set_attr("/", "fluidity_control", FLUIDITY_CONTROL)
            # Whether log_fluidity on grounded ice is held at zero (floating)
            # or inverted (all).
            chk.set_attr("/", "phi_grounded", PHI_GROUNDED)
            # Which controls this stage moved (provenance only: the objective
            # is the same in every ISMIP7_INVERT mode)
            chk.set_attr("/", "invert_controls", INVERT)
            # The grounding scheme and the front push the MAP was inverted
            # under: a forward follows them (icepack2_tools.subelement).
            chk.set_attr("/", "subelement_friction", int(SUBELEMENT_FRICTION))
            chk.set_attr("/", "subelement_scheme", SUBELEMENT_SCHEME)
            chk.set_attr("/", "subelement_scheme_version", int(SUBELEMENT_SCHEME_VERSION))
            chk.set_attr("/", "exact_front", int(EXACT_FRONT))
            if FRICTION_CONTROL in ("sqrt", "exp"):
                chk.set_attr("/", "prior_sigma_alpha", float(sigma_alpha_val))
                chk.set_attr("/", "prior_rho_theta", float(PRIOR_RHO_THETA))
            if FRICTION_CONTROL == "exp":
                chk.set_attr("/", "friction_c_ref", float(c_ref_val))
            chk.set_attr("/", "fluidity_prior_origin", str(prior_origin))
            chk.set_attr("/", "misfit_norm", MISFIT_NORM)
            chk.set_attr("/", "misfit_scale", float(misfit_scale))
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
            # The optimisation metric and the objective at the checkpointed
            # iterate, so the next link can verify it continues THIS
            # minimisation from THIS point (icepack2_tools.handoff).
            if run_settings.get("grad_precond") is not None:
                chk.set_attr("/", "grad_precond", str(run_settings["grad_precond"]))
            for _key, _val in last_accepted.items():
                chk.set_attr("/", f"objective_{_key}", float(_val))
            # How BedMachine was put onto the cells (runconfig.RASTER_SAMPLES).
            # theta/phi absorb the bed representation just as they absorb the
            # front treatment, so a forward must reproduce it.
            chk.set_attr("/", "raster_sample", raster_sample)
            # How each evaluation reached the full exponents: 1 for the
            # five-solve continuation, 0 for one solve. A solver fact, outside
            # the objective (handoff.OBJECTIVE_KEYS).
            chk.set_attr("/", "eval_continuation", int(not eval_full_n))
            chk.set_attr("/", "eval_mode", eval_mode)
            # Where the geometry came from, on every checkpoint, so the next
            # link of a chain and the forward know a relaxed geometry from
            # this mesh's own BedMachine sample (icepack2_tools.relaxation).
            if run_geometry is None:
                chk.set_attr("/", "geometry_source", os.path.realpath(bm_fn))
                chk.set_attr("/", "geometry_source_method", TARGET_MESH_GEOMETRY_METHOD)
            else:
                for _key, _val in run_geometry.items():
                    chk.set_attr("/", _key, _val)
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
                for key, value in full_state_solve.items():
                    chk.set_attr("/", f"full_state_{key}", value)

    # ── L-BFGS-B Inversion ──
    max_iter = int(os.environ.get("ISMIP7_MAXITER", "500"))
    # Relative-decrease stopping rule; 0 disables it and the budget above decides.
    ftol = float(os.environ.get("ISMIP7_FTOL", "1e-10"))
    min_iter = int(os.environ.get("ISMIP7_MIN_ITER", "3"))
    PETSc.Sys.Print("\nStarting L-BFGS-B inversion ("
                    + {"both": "theta + phi", "phi": "phi only, theta frozen",
                       "theta": "theta only, phi frozen"}[INVERT] + ")...")
    PETSc.Sys.Print(f"  maxiter={max_iter} ftol={ftol:g} min_iter={min_iter}, "
                    f"nranks={COMM_WORLD.size}")

    global_ndof = len(func_to_global(theta))
    z_backup = z.copy(deepcopy=True)
    last_good_obj = [np.inf]
    last_good_total = [np.inf]
    last_x = [None]                      # controls of the last CONVERGED evaluation
    trial_failures = TrialFailures()
    iteration_count = [0]
    timing_history = []
    timing_json = os.environ.get("ISMIP7_INVERSION_TIMING_JSON", "").strip()
    t_opt0 = perf_counter()
    # What an evaluation spends outside the forward and the adjoint (issue
    # 156: 27 s of a 43 s 1 km evaluation on Quartz), by span, each the
    # slowest rank's time. `spans` covers the work inside total_seconds;
    # `gap_spans` the work between two evaluations: the previous one's report,
    # term assembly, timing write and checkpoint, and the optimizer's own step.
    spans = Spans(COMM_WORLD)
    gap_spans = Spans(COMM_WORLD)
    t_body_end = [None]

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
        if not timing_json:
            return
        optimize_seconds = global_max(
            np.array([perf_counter() - t_opt0]), comm=COMM_WORLD
        )
        if COMM_WORLD.rank != 0:
            return
        written = os.path.realpath(os.path.join(_map_dir, map_fn))
        published = os.environ.get("ISMIP7_MAP_OUT_FINAL", "").strip()
        payload = {
            "phase": phase,
            "map_path": os.path.realpath(published) if published else written,
            "map_path_written": written,
            "mesh_basename": os.path.basename(mesh_fn),
            "lc": int(mesh_lc),
            "lc_coarse": int(mesh_lc_coarse),
            "buffer_m": float(mesh_buffer_m),
            "ncores": int(COMM_WORLD.size),
            "maxiter": int(max_iter),
            "nit": int(nit if nit is not None else iteration_count[0]),
            "nfev": int(nfev if nfev is not None else iteration_count[0]),
            "message": str(message),
            "optimize_seconds": optimize_seconds,
            "knobs": {
                "misfit_norm": MISFIT_NORM,
                "misfit_scale": float(misfit_scale),
                "misfit_scale_requested": MISFIT_SCALE,
                "log_vel_weight_requested": LOG_VEL_WEIGHT,
                "log_vel_weight": float(log_vel_w),
                "log_vel_weight_source": log_vel_source,
                "log_vel_eps": float(LOG_VEL_EPS),
                "gamma_theta": float(GAMMA_THETA),
                "friction_control": FRICTION_CONTROL,
                "fluidity_prior": FLUIDITY_PRIOR,
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
                "eval_continuation": bool(not eval_full_n),
                "eval_mode": eval_mode,
                "inversion_linear_solver": state_solver_mode,
                "inversion_snes_linesearch": sparams.get("snes_linesearch_type"),
                # None under full_mumps, whose LU has no Krylov tolerance
                "inversion_ksp_rtol": sparams.get("ksp_rtol"),
                "diagnostic_linear_solver": lane_solver_mode,
            },
            "evaluations": list(timing_history),
            # Set only by the "finished" record: how the publishing solve
            # ended and whether the full mixed state was written.
            "final_solve": final_solve,
        }
        atomic_write_json(timing_json, payload)

    # Periodic checkpoint interval: accepted iterates on the TAO path,
    # evaluations on the scipy path. A wall-clocked link resumes from the
    # last one, so every iteration past it is repeated: 2 km link 1643735
    # lost iterations 21 to 30 (about five hours) to the old interval of 20,
    # and at five, with iterations of 3-6 h, links 1656863 and 1662734 each
    # lost four. A 2 km MAP write costs a minute.
    _ckpt_every = max(1, int(os.environ.get("ISMIP7_CHECKPOINT_EVERY_IT", "1")))

    def objective_and_gradient(x_vec):
        t_iter = perf_counter()
        spans.clear()
        if t_body_end[0] is not None:
            gap_spans.add("gap", t_iter - t_body_end[0])
        with spans("set_controls"):
            global_to_func(x_vec[:global_ndof], theta)
            set_phi(x_vec[global_ndof:])

        t_fwd = perf_counter()
        reset_manager()
        start_manager()
        try:
            J = forward(theta, phi)
        except fd.ConvergenceError as exc:
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
            PETSc.Sys.Print(
                f"  [!] Forward solve failed ({exc}), returning large objective")
            t_body_end[0] = perf_counter()
            return trial_failures.failed(last_good_total[0], 2 * global_ndof)
        stop_manager()
        J_val = float(J)
        t_fwd = perf_counter() - t_fwd

        with spans("record_state"):
            z_backup.assign(z)
            last_good_obj[0] = J_val
            last_x[0] = np.array(x_vec, copy=True)
        with spans("vel_chi2"):
            last_good_vel_chi2[0] = float(assemble(_vel_chi2))
        with spans("residual_norm"):
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
            t_body_end[0] = perf_counter()
            return trial_failures.failed(last_good_total[0], 2 * global_ndof)
        t_adj = perf_counter() - t_adj
        trial_failures.succeeded()

        # Whittle-Matern prior energy + gradient (icepack2_tools/prior.py):
        # 0.5/area * gamma * (theta^2 + L^2 |grad theta|^2). The theta^2 mass
        # term (absent from the old gradient-only form) removes the null space
        # that let phi/theta drift to +-36 at n=3; combined with the physical
        # prior means it constrains the DEVIATION, not the amplitude.
        # Energy first, then gradient: under `bilaplacian` the energy solves
        # the M f = A theta the gradient reuses (see _prior_energy_form).
        with spans("prior_solve"):
            _reg_form_theta = _prior_energy_form(theta, "theta")
        with spans("prior_energy"):
            reg_theta = float(assemble(_reg_form_theta))
        with spans("prior_grad"):
            dR_theta = _prior_grad(theta, "theta")
        with spans("prior_solve"):
            _reg_form_phi = _prior_energy_form(phi, "phi")
        with spans("prior_energy"):
            reg_phi = float(assemble(_reg_form_phi))
        with spans("prior_grad"):
            dR_phi = _prior_grad(phi, "phi")

        with spans("gather_gradient"):
            g_theta = func_to_global(dJ_dtheta) + func_to_global(dR_theta)
            g_phi = func_to_global(dJ_dphi) + func_to_global(dR_phi)
            if _phi_free_global is not None:
                g_phi = g_phi * _phi_free_global

            total = J_val + reg_theta + reg_phi
            last_good_total[0] = total
            total_grad = np.concatenate([g_theta, g_phi])

        t_body_end[0] = perf_counter()
        t_iter = t_body_end[0] - t_iter
        _other = t_iter - t_fwd - t_adj
        # Collective, outside total_seconds: every rank reaches this line once
        # per successful evaluation.
        eval_spans, eval_durations = spans.reduce(
            durations={
                "fwd_seconds": t_fwd,
                "adj_seconds": t_adj,
                "total_seconds": t_iter,
                "other_seconds": _other,
            },
            unspanned_total=_other,
        )
        t_fwd = eval_durations["fwd_seconds"]
        t_adj = eval_durations["adj_seconds"]
        t_iter = eval_durations["total_seconds"]
        _other = eval_durations["other_seconds"]
        _unspanned = eval_durations["unspanned"]
        prev_gap = gap_spans.reduce()
        iteration_count[0] += 1
        # NB: this path records the last EVALUATED point (L-BFGS-B trial
        # points included); the TAO path below records the accepted iterate.
        last_accepted.update(iteration=iteration_count[0], total=total, misfit=J_val,
                             reg_theta=reg_theta, reg_phi=reg_phi)
        if iteration_count[0] == 1:
            _check_handoff(total)
        with gap_spans("report"):
            PETSc.Sys.Print(
                f"  iter {iteration_count[0]:3d}: "
                f"misfit={J_val:.6e}{term_report()} "
                f"reg_θ={reg_theta:.4e} reg_φ={reg_phi:.4e} "
                f"total={total:.6e} |grad|={np.linalg.norm(total_grad):.4e} "
                f"[fwd={t_fwd:.1f}s adj={t_adj:.1f}s total={t_iter:.1f}s]"
            )
            PETSc.Sys.Print(
                f"    other={_other:.2f}s: "
                + " ".join(f"{k}={v:.2f}" for k, v in eval_spans.items())
                + (" | before: " + " ".join(f"{k}={v:.2f}" for k, v in prev_gap.items())
                   if prev_gap else "")
            )

        if timing_json:
            with gap_spans("eval_terms"):
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
                "other_spans": {**eval_spans, "unspanned": _unspanned},
                # since the previous evaluation's body ended: its report,
                # terms, timing write and checkpoint, then the optimizer's
                # step ("gap" is the whole interval)
                "before_spans": prev_gap,
                "terms": terms,
                "state_solves": [w for w in state_work if w],
                # after the adjoint: its growth across evaluations is how a
                # solver leak shows (issue #159)
                "rss_mib": global_rss_mib(COMM_WORLD),
            })
            with gap_spans("timing_json"):
                _write_timing_json(phase="running", message="in progress")

        if iteration_count[0] % _ckpt_every == 0:
            with gap_spans("checkpoint"):
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
        _A_inv_untimed = _prior_metric_solvers(grad_precond)

        def _A_inv(g_theta, g_phi):
            with spans("metric_action"):
                return _A_inv_untimed(g_theta, g_phi)

        _nfev = [0]
        _ring = []

        # A forward that "converged" by the relative test from a wild initial
        # residual can sit at ||F|| 1e9 (stage 2 of the 26 Sep staged 32 km
        # run: accepted iterate at 3.4e9, published as converged). Such a
        # state is a failed solve for the objective's purposes: it takes the
        # trial-point rescue path like a diverged one. The ceiling is a
        # multiple of the residual the last accepted forward reached.
        _fnorm_ceiling_factor = float(os.environ.get("ISMIP7_FNORM_CEILING", "1e4"))

        def _forward_checked(theta_ctrl, phi_ctrl):
            with spans("forward"):
                J = forward(theta_ctrl, phi_ctrl)
            # Never below the warm start's converged residual: an accepted
            # forward that began at the converged state ends at the rounding
            # floor (5e-5 on the 2 km mesh), and 1e4x that rejected an
            # ordinary relative-test solve at ||F|| 1.8 (NOTS 1691937).
            f_ref = float(np.nanmax([float(last_good_fnorm[0]), f_warm]))
            if np.isfinite(f_ref) and f_ref > 0.0:
                # A check, not part of the objective, so off the tape:
                # recorded, this one assembly took 4 to 8 s a call at 32 km
                # and at 4 km, against 0.1 s unrecorded.
                with spans("residual_norm"), paused_manager():
                    f_now = _residual_norm()
                if not np.isfinite(f_now) or f_now > _fnorm_ceiling_factor * f_ref:
                    raise fd.ConvergenceError(
                        f"forward reported convergence at ||F||={f_now:.3e}, above "
                        f"{_fnorm_ceiling_factor:g}x the last accepted {f_ref:.3e}")
            return J

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
            try:
                J = _forward_checked(theta_ctrl, phi_ctrl)
            except fd.ConvergenceError as err:
                PETSc.Sys.Print(f"    ({str(err).splitlines()[0][:200]})")
                # A line-search trial point where the single Newton solve
                # of the forward diverges (the trial is far from the last
                # converged state). First rescue: re-climb the exponent
                # ladder at the trial controls, unannotated, from the last
                # converged state, then take the annotated solve from there.
                # Failing that, the scipy path's rescue: restore the state
                # and hand TAO an inflated objective with no control
                # dependence, so its line search backtracks instead of the
                # run ending here.
                z.assign(z_backup)
                if not np.isfinite(last_good_obj[0]):
                    raise RuntimeError(
                        "First forward solve failed - the inversion cannot start "
                        "(fewer MPI ranks for a small mesh; check the fluidity prior).")
                stop_manager()
                try:
                    if trial_rescue_rungs() == 0:
                        raise fd.ConvergenceError("trial rescue disabled")
                    PETSc.Sys.Print(
                        "  [!] Forward solve failed at a trial point; re-climbing "
                        "the continuation there")
                    with spans("reramp"):
                        _reramp_at_current_controls()
                    reset_manager()
                    start_manager()
                    J = _forward_checked(theta_ctrl, phi_ctrl)
                except fd.ConvergenceError as err:
                    PETSc.Sys.Print(f"    ({str(err).splitlines()[0][:200]})")
                    z.assign(z_backup)
                    reset_manager()
                    start_manager()
                    PETSc.Sys.Print(
                        "  [!] The continuation failed there too; returning an "
                        "inflated objective so the line search backtracks")
                    J = Functional(name="J_failed")
                    J.assign(float(10.0 * last_good_obj[0]))
                    return J
            with spans("prior_taped"):
                J.addto(_prior_energy_form(theta_ctrl, "theta"))
                J.addto(_prior_energy_form(phi_ctrl, "phi"))
            # The last few evaluations with their controls: the monitor picks
            # the one TAO accepted, so a checkpoint never holds a rejected
            # line-search trial point.
            # The mixed state goes in with the controls: the monitor puts
            # the ACCEPTED evaluation's state back into z, so the state on
            # record, its residual and the published velocity belong to the
            # controls the MAP saves. Without it z_backup held the last
            # EVALUATED trial's state, and after a failed line search the
            # publishing solve started from a state of other controls,
            # "converged" in 0 iterations at an atol scaled from its own
            # residual (32 km joint Pattyn run, 26 Sep: ||F|| 3e35 published).
            with spans("record_state"):
                _ring.append((float(J), theta_ctrl.dat.data_ro.copy(),
                              phi_ctrl.dat.data_ro.copy(),
                              [_z.dat.data_ro.copy() for _z in z.subfunctions]))
                del _ring[:-6]
                for _zb, _z in zip(z_backup.subfunctions, z.subfunctions):
                    _zb.dat.data[:] = _z.dat.data_ro
            return J

        _accepted_entry = [None]

        def _restore_accepted(f_val):
            """Put the evaluation TAO accepted at objective ``f_val`` into
            (theta, phi, z, z_backup); returns the entry or None. The last
            accepted entry stays pinned so a long run of rejected trials
            cannot push it out of the ring."""
            _hit = accepted_evaluation(
                _ring + ([_accepted_entry[0]] if _accepted_entry[0] is not None else []),
                float(f_val))
            if _hit is None:
                return None
            _accepted_entry[0] = _hit
            theta.dat.data[:] = _hit[1]
            phi.dat.data[:] = _hit[2]
            for _z, _zb, _d in zip(z.subfunctions, z_backup.subfunctions, _hit[3]):
                _z.dat.data[:] = _d
                _zb.dat.data[:] = _d
            return _hit

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
        if os.environ.get("ISMIP7_GRAD_CHECK", "0").strip() == "1":
            # ISMIP7_GRAD_CHECK=1: verify the taped gradient of the FULL
            # objective by a Taylor test at the start, print the spread of
            # the metric-scaled first step over the nodes, and stop. Written
            # for the exp-control question (Rice, 28 Sep 2026): is the
            # objective or its gradient wrong under C = C_ref exp(alpha)?
            from tlm_adjoint import taylor_test
            import time as _time
            _theta0 = theta.copy(deepcopy=True)
            _phi0 = phi.copy(deepcopy=True)
            reset_manager()
            start_manager()
            _J0 = forward_total(theta, phi)
            stop_manager()
            _dJ = compute_gradient(_J0, [theta, phi])
            _xt, _xp = _A_inv(_dJ[0], _dJ[1])
            _Hc = cg1_lift(H) if geom_dg else H
            _bc = cg1_lift(b) if geom_dg else b
            _hafc = Function(Q).interpolate(height_above_flotation(_Hc, _bc)).dat.data_ro
            _gr = (_hafc > 0.0) & (_Hc.dat.data_ro > 10.0)
            _fl = (_hafc <= 0.0) & (_Hc.dat.data_ro > 10.0)

            def _pct(arr, msk):
                a = np.concatenate(COMM_WORLD.allgather(
                    np.abs(np.asarray(arr[msk], dtype=float))))
                if a.size == 0:
                    return "n/a"
                q = np.percentile(a, [50, 90, 99, 100])
                return (f"median {q[0]:.3e} p90 {q[1]:.3e} p99 {q[2]:.3e} max {q[3]:.3e} "
                        f"(max/median {q[3] / max(q[0], 1e-300):.1e})")
            PETSc.Sys.Print(f"  Gradient check: objective {float(_J0):.6e}")
            PETSc.Sys.Print(f"    raw dJ/dtheta, grounded nodes: {_pct(_dJ[0].dat.data_ro, _gr)}")
            PETSc.Sys.Print(f"    raw dJ/dphi, floating nodes:   {_pct(_dJ[1].dat.data_ro, _fl)}")
            PETSc.Sys.Print(f"    first step H0 g, theta grounded: {_pct(_xt.dat.data_ro, _gr)}")
            PETSc.Sys.Print(f"    first step H0 g, phi floating:   {_pct(_xp.dat.data_ro, _fl)}")
            _seed = float(os.environ.get("ISMIP7_GRAD_CHECK_SEED", "1e-3"))
            _t0 = _time.perf_counter()
            _order = taylor_test(forward_total, [theta, phi], J_val=float(_J0), dJ=_dJ,
                                 seed=_seed, size=4)
            PETSc.Sys.Print(
                f"  Taylor test (seed {_seed:g}, 4 sizes): minimum order {_order:.3f} "
                f"(2 = the gradient is consistent with the objective; 1 = it is not) "
                f"[{_time.perf_counter() - _t0:.0f}s]")
            theta.dat.data[:] = _theta0.dat.data_ro
            phi.dat.data[:] = _phi0.dat.data_ro
            PETSc.Sys.Print("  GRAD_CHECK: stopped after the Taylor test; no MAP written")
            COMM_WORLD.barrier()
            sys.exit(0)
        if INVERT == "both":
            _tao_forward, _tao_spaces, _tao_x, _tao_action = (
                forward_total, [Q, Q], [theta, phi], _A_inv)
        else:
            # One control moves. The other enters the forward as a fixed
            # coefficient (not a control: no gradient, its prior term a
            # constant), so the objective TAO sees is forward_total's with
            # that field held where it is.
            _frozen = Function(Q, name="frozen_control")
            _frozen.dat.data[:] = (theta if INVERT == "phi" else phi).dat.data_ro
            if INVERT == "phi":
                def _tao_forward(phi_ctrl):
                    return forward_total(_frozen, phi_ctrl)

                def _tao_action(g_phi):
                    return _A_inv(None, g_phi)
                _tao_x = [phi]
            else:
                def _tao_forward(theta_ctrl):
                    return forward_total(theta_ctrl, _frozen)

                def _tao_action(g_theta):
                    return _A_inv(g_theta, None)
                _tao_x = [theta]
            _tao_spaces = [Q]
            PETSc.Sys.Print(
                f"  Controls: {INVERT} only (ISMIP7_INVERT); "
                f"{'phi' if INVERT == 'theta' else 'theta'} frozen at its current field")
        solver = TAOSolver(
            _tao_forward, _tao_spaces,
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
            H_0_action=_tao_action, M_inv_action=_tao_action,
        )

        _t_last = [perf_counter()]
        _ftol_stop = FunctionalDecreaseStop(ftol, min_iter)
        def _monitor(tao):
            its, f_val, gnorm, _cnorm, _xdiff, _reason = tao.getSolutionStatus()
            iteration_count[0] = int(its)
            # theta/phi/z mirror the LAST EVALUATION; put the ACCEPTED
            # iterate (controls AND mixed state) there before anything below
            # reads, measures or saves them.
            with spans("monitor_restore"):
                _hit = _restore_accepted(f_val)
            if _hit is None:
                PETSc.Sys.Print(
                    f"    WARNING: no recent evaluation matches the accepted "
                    f"objective {float(f_val):.6e}; the controls on record are "
                    "the last evaluated point")
            if its == 0:
                _check_handoff(float(f_val), tao=tao)
            # An accepted point identical to the last one is a line search
            # that found no new point (every trial diverged), not a
            # functional decrease of zero: it must not read as convergence.
            with spans("monitor_gather"):
                _x_acc = np.concatenate([func_to_global(theta), func_to_global(phi)])
            _unchanged = (its > 0 and last_good_x[0] is not None
                          and np.array_equal(_x_acc, last_good_x[0]))
            if _unchanged:
                PETSc.Sys.Print(
                    "    (the line search accepted no new point this iteration; "
                    "not a functional-decrease stop)")
            elif _ftol_stop.update(iteration_count[0], f_val):
                tao.setConvergedReason(PETSc.TAO.ConvergedReason.CONVERGED_USER)
            # z holds the last evaluated forward, mirrored into z_backup above.
            last_good_obj[0] = float(f_val)
            with spans("residual_norm"):
                last_good_fnorm[0] = _residual_norm()
            with spans("vel_chi2"):
                last_good_vel_chi2[0] = float(assemble(_vel_chi2))
            with spans("monitor_gather"):
                _x = np.concatenate([func_to_global(theta), func_to_global(phi)])
                last_x[0] = _x
                last_good_x[0] = np.array(_x, copy=True)
            with spans("prior_solve"):
                _reg_form_theta = _prior_energy_form(theta, "theta")
            with spans("prior_energy"):
                reg_theta = float(assemble(_reg_form_theta))
            with spans("prior_solve"):
                _reg_form_phi = _prior_energy_form(phi, "phi")
            with spans("prior_energy"):
                reg_phi = float(assemble(_reg_form_phi))
            # One TAO iteration: every evaluation since the last monitor call
            # (the line search's trials included), the adjoint(s), TAO's own
            # work and the bookkeeping above. The report, timing write and
            # checkpoint below land in the next iteration's.
            now = perf_counter()
            t_iter = now - _t_last[0]
            _t_last[0] = now
            # Collective: TAO calls the monitor on every rank.
            iter_spans, iter_durations = spans.reduce(
                durations={"total_seconds": t_iter},
                unspanned_total=t_iter,
            )
            t_iter = iter_durations["total_seconds"]
            _unspanned = iter_durations["unspanned"]
            last_accepted.update(iteration=int(its), total=float(f_val),
                                 misfit=float(f_val) - reg_theta - reg_phi,
                                 reg_theta=reg_theta, reg_phi=reg_phi)
            PETSc.Sys.Print(
                f"  iter {iteration_count[0]:3d}: "
                f"misfit={f_val - reg_theta - reg_phi:.6e} "
                f"reg_θ={reg_theta:.4e} reg_φ={reg_phi:.4e} "
                f"total={f_val:.6e} |grad|_A={gnorm:.4e} "
                f"dJ/J={_ftol_stop.criterion if _ftol_stop.criterion is not None else 0.0:.1e} "
                f"[total={t_iter:.1f}s]"
            )
            # "unspanned" is the adjoint(s) and TAO's own work, which run
            # inside TAOSolver where no span reaches.
            PETSc.Sys.Print(
                "    " + " ".join(f"{k}={v:.2f}" for k, v in iter_spans.items())
                + f" unspanned={_unspanned:.2f}s")
            if timing_json:
                timing_history.append({
                    "eval": iteration_count[0],
                    "misfit": float(f_val - reg_theta - reg_phi),
                    "reg_theta": reg_theta,
                    "reg_phi": reg_phi,
                    "total": float(f_val),
                    "grad_norm": float(gnorm),
                    "total_seconds": t_iter,
                    "iteration_spans": {**iter_spans, "unspanned": _unspanned},
                    "terms": {"vel": float(last_good_vel_chi2[0])},
                    "state_solves": [w for w in state_work if w],
                    "rss_mib": global_rss_mib(COMM_WORLD),
                })
                with spans("timing_json"):
                    _write_timing_json(phase="running", message="in progress")
            if iteration_count[0] > 0 and iteration_count[0] % _ckpt_every == 0:
                with spans("checkpoint"):
                    save_map(os.path.join(_map_dir, map_fn))
                PETSc.Sys.Print(f"    [checkpoint saved: iter {iteration_count[0]}]")

        solver.tao.setMonitor(_monitor)
        # TAO reports the iteration cap as a diverged reason and TAOSolver
        # turns any such reason into an exception -- but only AFTER writing the
        # solution back into (theta, phi). Reaching ISMIP7_MAXITER is how a
        # production inversion normally ends here, so the cap is read back from
        # TAO rather than treated as a failure.
        try:
            solver.solve(_tao_x)
        except RuntimeError:
            pass
        if _handoff_failed[0]:
            raise RuntimeError(_handoff_failed[0])
        # TAOSolver writes its current point back into (theta, phi) on exit,
        # which after an aborted line search is a rejected trial point. The
        # MAP must hold the last ACCEPTED iterate, which the monitor kept.
        if last_good_x[0] is not None:
            _x_now = np.concatenate([func_to_global(theta), func_to_global(phi)])
            if not np.array_equal(_x_now, last_good_x[0]):
                global_to_func(last_good_x[0][:global_ndof], theta)
                global_to_func(last_good_x[0][global_ndof:], phi)
                PETSc.Sys.Print(
                    "  TAO exited on a point other than the last accepted iterate; "
                    "the accepted iterate is what the MAP records")
            # and its mixed state, evaluated at those controls, is what the
            # publishing solve starts from
            if _restore_accepted(last_good_obj[0]) is None:
                PETSc.Sys.Print(
                    "    WARNING: the accepted iterate's mixed state is not on "
                    "record; the publishing solve starts from the last evaluated state")
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

    # ── One objective across restarts (icepack2_tools.handoff) ──────────
    # What this run minimises, compared with what the warm start's writer
    # minimised. A chain link that would silently change the objective is a
    # different inversion; ISMIP7_WARM_START_STRICT=1 (the chain runner's
    # default) refuses it, otherwise the differences are printed. The first
    # evaluation must then reproduce the recorded objective to
    # ISMIP7_WARM_START_OBJ_TOL (relative, default 1e-3: the re-solved state
    # differs from the writer's by the forward tolerance only).
    run_settings.update({
        "misfit_norm": MISFIT_NORM, "misfit_scale": float(misfit_scale),
        "log_vel_weight": float(log_vel_w), "log_vel_eps": float(LOG_VEL_EPS),
        "dhdt_weight": float(dhdt_w), "dhdt_net_sigma": net_sigma_used,
        "prior_form": PRIOR_FORM, "gamma_theta": float(GAMMA_THETA),
        "gamma_phi": float(GAMMA_PHI), "friction_control": FRICTION_CONTROL,
        "friction": str(FRICTION), "n_flow": float(n_flow_val),
        "geometry_space": str(geometry_space),
        "friction_anchor_length": float(ANCHOR_LENGTH),
        "lake_ice_base": int(LAKE_ICE_BASE),
        "fluidity_prior_origin": str(prior_origin), "grad_precond": grad_precond,
        "subelement_friction": int(SUBELEMENT_FRICTION), "exact_front": int(EXACT_FRONT),
        "subelement_scheme": SUBELEMENT_SCHEME,
        "subelement_scheme_version": int(SUBELEMENT_SCHEME_VERSION),
        "fluidity_control": FLUIDITY_CONTROL,
        "drag_gate": DRAG_RECORD, "h_visc_floor": float(RC_HVISC_FLOOR),
        "phi_grounded": PHI_GROUNDED, "raster_sample": str(raster_sample),
    })
    if PRIOR_FORM == "bilaplacian":
        run_settings.update({"prior_sigma_theta": float(PRIOR_SIGMA_THETA),
                             "prior_sigma_phi": float(PRIOR_SIGMA_PHI),
                             "prior_rho": float(PRIOR_RHO)})
        if FRICTION_CONTROL in ("sqrt", "exp"):
            run_settings.update({"prior_sigma_alpha": float(sigma_alpha_val),
                                 "prior_rho_theta": float(PRIOR_RHO_THETA)})
        if FRICTION_CONTROL == "exp":
            run_settings["friction_c_ref"] = float(c_ref_val)
    _strict = os.environ.get("ISMIP7_WARM_START_STRICT", "0").strip() == "1"
    _obj_tol = float(os.environ.get("ISMIP7_WARM_START_OBJ_TOL", "1e-3"))
    if warm_attrs:
        _mism = objective_mismatches(warm_attrs, run_settings)
        if _mism:
            PETSc.Sys.Print("  Handoff: this run's objective differs from the warm start's:")
            for _m in _mism:
                PETSc.Sys.Print(f"    {_m}")
            if _strict:
                raise RuntimeError(
                    "ISMIP7_WARM_START_STRICT=1: refusing to continue a "
                    "different objective from the warm start's controls "
                    f"({len(_mism)} setting(s) differ, listed above)")
        else:
            PETSc.Sys.Print("  Handoff: objective settings match the warm start's")
    _handoff_checked = [False]

    def _check_handoff(first_total, tao=None):
        """The first evaluation against the recorded objective."""
        if _handoff_checked[0]:
            return
        _handoff_checked[0] = True
        if "objective_total" not in warm_attrs:
            # A relaxation's end state carries the MAP's objective settings
            # and leaves its value behind, since that value belongs to the
            # MAP's geometry: the change is the relaxation's, not a gap.
            if warm_end_state and "relax_source_objective_total" in warm_geometry_attrs:
                PETSc.Sys.Print(
                    f"  Handoff: first objective {float(first_total):.6e} on the "
                    f"relaxed geometry; the source MAP recorded "
                    f"{float(warm_geometry_attrs['relax_source_objective_total']):.6e} "
                    f"on its own")
            return
        _rec = float(warm_attrs["objective_total"])
        _gap = handoff_gap(_rec, first_total)
        _msg = (f"  Handoff: first objective {float(first_total):.6e} vs recorded "
                f"{_rec:.6e} (relative gap {_gap:.2e}, tolerance {_obj_tol:g})")
        PETSc.Sys.Print(_msg)
        if _gap > _obj_tol:
            _why = (f"the warm start's objective is not reproduced (gap {_gap:.2e} "
                    f"> {_obj_tol:g}): the optimiser would not be continuing the "
                    "same minimisation")
            if _strict:
                _handoff_failed[0] = _why
                if tao is not None:
                    tao.setConvergedReason(PETSc.TAO.ConvergedReason.DIVERGED_USER)
                else:
                    raise RuntimeError(_why)
            else:
                PETSc.Sys.Print(f"    WARNING: {_why}")
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
        if INVERT != "both":
            raise ValueError(
                f"ISMIP7_INVERT={INVERT} needs the TAO path "
                "(ISMIP7_GRAD_PRECOND=mass_consistent or prior)")
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
    set_phi(_x_final[global_ndof:])
    _t_lo, _t_hi = global_range(theta)
    _p_lo, _p_hi = global_range(phi)
    PETSc.Sys.Print(f"  theta range: [{_t_lo:.3f}, {_t_hi:.3f}]")
    PETSc.Sys.Print(f"  phi range:   [{_p_lo:.3f}, {_p_hi:.3f}]")
    # A stop right after failed trials is no convergence test (TrialFailures),
    # so the MAP gets no done marker and the chain's next link resumes from
    # it. Printed before the MAP is saved: the runner reads this line, and a
    # kill between the two must not leave a "Saved MAP:" line without it.
    unfinished = trial_failures.stopped_on_failures
    if unfinished:
        PETSc.Sys.Print(
            f"  {NOT_FINAL} the optimizer stopped right after "
            f"{trial_failures.trailing} failed trial evaluations in a row "
            f"({trial_failures.total} in this run). The MAP holds the last "
            "accepted iterate; no done marker, so the chain's next link "
            "resumes from it.")

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
    if map_out and COMM_WORLD.rank == 0 and not unfinished:
        with open(map_out + ".done", "w"):
            pass

    if timing_json:
        _write_timing_json(
            phase="final_solve",
            message=str(result.message) + (
                f" ({NOT_FINAL} after {trial_failures.trailing} failed trials)"
                if unfinished else ""),
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
    # The publishing solve builds its own solver; release the forwards'.
    state_solver_cache.clear()
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
    final_solver = _untaped_state_solver(F, final_sparams)
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

    out_fn = os.path.join(FIG_DIR, f"inversion_icepack2_{mesh_lc}.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    PETSc.Sys.Print(f"Saved: {out_fn}")


if __name__ == "__main__":
    main()
