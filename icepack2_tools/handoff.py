r"""handoff.py - keeping one objective across an inversion's restarts.

A chained inversion warm-starts each link from the previous link's periodic
checkpoint. The optimiser only continues the same minimisation if the link
minimises the same objective from the point the previous link had accepted:

* every setting the objective depends on is recorded in the checkpoint and
  compared here on the warm start (:func:`objective_mismatches`); a link
  that would change one refuses under ``ISMIP7_WARM_START_STRICT=1``;
* the checkpoint holds the controls of the last ACCEPTED iterate, not of the
  last evaluated trial point of the line search
  (:func:`accepted_evaluation`);
* the checkpoint records the objective at that iterate, and the next link's
  first evaluation must reproduce it (:func:`handoff_gap`).

Rice (26 Sep 2026): "between restarts I think we need to check to make
sure the optimizer doesn't change."
"""
import math

# Checkpoint attributes the objective depends on. Auto-derived weights
# (log_vel_weight, prior_sigma_alpha) are FROZEN from the warm start rather
# than re-derived, so they appear here as plain values.
OBJECTIVE_KEYS = (
    "misfit_norm", "misfit_scale", "log_vel_weight", "log_vel_eps",
    "dhdt_weight", "dhdt_net_sigma", "dhdt_filter_l",
    "dhdt_speed_ratio", "dhdt_speed_slow", "dhdt_speed_fast",
    "prior_form", "gamma_theta", "gamma_phi",
    "prior_sigma_theta", "prior_sigma_phi", "prior_rho",
    "prior_sigma_alpha", "prior_rho_theta", "friction_c_ref",
    "friction_control", "friction", "n_flow", "geometry_space",
    "friction_anchor_length", "lake_ice_base", "fluidity_prior_origin",
    "grad_precond", "subelement_friction", "subelement_scheme",
    "subelement_scheme_version", "exact_front", "fluidity_control",
    "drag_gate", "h_visc_floor", "phi_grounded", "raster_sample",
    "vel_sample", "u_obs_max", "dhdt_reweight", "dhdt_reweight_scale", "kill_bergs",
    "velocity_sample",
)

# The residual each sub-element scheme builds today (icepack2_tools.subelement),
# recorded in the MAP as subelement_scheme_version. A change to a scheme's
# residual at the same control fields takes a new version.
#   sep2 1: quadrature over the grounded part (26 Sep 2026).
#   sep1 1: whole-cell quadrature, drag times the grounded fraction, theta
#           through the shared Weertman closure as exp(He theta) (PR 158,
#           73bd05a).
#   sep1 2: the same with exp(theta) ungated, as under SEP2 (PR 158, 1c13b38).
# The two sep1 forms agree where theta is zero: the sqrt and exp friction
# controls, which carry the friction in C_w0.
SUBELEMENT_SCHEME_VERSIONS = {"sep2": 1, "sep1": 2}

# Recorded with the controls: the objective at the checkpointed iterate.
OBJECTIVE_RECORD_KEYS = (
    "objective_total", "objective_misfit", "objective_reg_theta",
    "objective_reg_phi", "objective_iteration",
)


def _same(a, b, rel_tol):
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)
    if math.isnan(fa) and math.isnan(fb):
        return True
    return math.isclose(fa, fb, rel_tol=rel_tol, abs_tol=0.0)


def _truthy(value):
    try:
        return bool(int(float(value)))
    except (TypeError, ValueError):
        return False


def recorded_scheme_version(scheme, version, friction_control):
    r"""The version of ``SUBELEMENT_SCHEME_VERSIONS`` a sub-element record
    stands for. A record from before the version existed carries none: SEP2
    has had one form, and under the sqrt and exp controls both SEP1 forms
    agree, so those stand for the current version; a log-control SEP1 record
    cannot say which form it was inverted under, and stands for none."""
    if version is not None:
        return int(version)
    scheme = str(scheme)
    if scheme == "sep2":
        return SUBELEMENT_SCHEME_VERSIONS["sep2"]
    if str(friction_control) in ("sqrt", "exp"):
        return SUBELEMENT_SCHEME_VERSIONS.get(scheme)
    return None


def check_subelement_record(scheme, version, friction_control, source="the MAP"):
    r"""Refuse a sub-element MAP this code cannot rebuild: one recorded under
    another version of its scheme, or a log-control SEP1 MAP from before the
    version was recorded. Returns the version."""
    scheme = str(scheme)
    if scheme not in SUBELEMENT_SCHEME_VERSIONS:
        raise RuntimeError(f"{source} records sub-element scheme {scheme!r}; "
                           f"this code builds {sorted(SUBELEMENT_SCHEME_VERSIONS)}")
    current = SUBELEMENT_SCHEME_VERSIONS[scheme]
    got = recorded_scheme_version(scheme, version, friction_control)
    if got is None:
        raise RuntimeError(
            f"{source} was inverted under {scheme.upper()} with the "
            f"{friction_control} friction control before the scheme's version "
            f"was recorded, so it cannot say whether theta entered as exp(He "
            f"theta) (version 1) or exp(theta) (version 2, this code)")
    if got != current:
        raise RuntimeError(
            f"{source} was inverted under {scheme.upper()} version {got}; this "
            f"code builds version {current} (icepack2_tools.handoff."
            f"SUBELEMENT_SCHEME_VERSIONS)")
    return got


def objective_mismatches(recorded, current, rel_tol=1e-9):
    r"""The keys of ``current`` whose recorded value differs, as
    ``"key: recorded -> current"`` strings. Keys the record lacks are not
    mismatches (an older checkpoint), keys ``current`` lacks are ignored.
    A sub-element record with no ``subelement_scheme`` was written under
    SEP2, the only scheme before the record existed, and one with no
    ``subelement_scheme_version`` stands for ``recorded_scheme_version``; a
    record that cannot name its version (log-control SEP1) is a mismatch.
    The scheme and its version are compared only when both sides use
    sub-element friction, since otherwise they select nothing."""
    out = []
    if _truthy(recorded.get("subelement_friction")):
        if "subelement_scheme" not in recorded:
            recorded = {**recorded, "subelement_scheme": "sep2"}
        if "subelement_scheme_version" not in recorded:
            version = recorded_scheme_version(
                recorded["subelement_scheme"], None,
                recorded.get("friction_control", "log"))
            recorded = {**recorded, "subelement_scheme_version":
                        "unrecorded" if version is None else version}
    both_sub = (_truthy(recorded.get("subelement_friction"))
                and _truthy(current.get("subelement_friction")))
    for key in OBJECTIVE_KEYS:
        if key in ("subelement_scheme", "subelement_scheme_version") and not both_sub:
            continue
        if key in recorded and key in current:
            if not _same(recorded[key], current[key], rel_tol):
                out.append(f"{key}: {recorded[key]!r} -> {current[key]!r}")
    return out


def accepted_evaluation(ring, f_val, rel_tol=1e-10):
    r"""The entry of ``ring`` (an iterable of ``(value, *payload)`` for the
    most recent evaluations) whose value is ``f_val`` to ``rel_tol``: the
    evaluation the optimiser accepted. None if no entry matches, e.g. after
    a restart of the ring."""
    best, best_gap = None, math.inf
    for entry in ring:
        value = float(entry[0])
        gap = abs(value - f_val) / max(abs(f_val), 1e-300)
        if gap < best_gap:
            best, best_gap = entry, gap
    return best if best_gap <= rel_tol else None


def handoff_gap(recorded_total, first_total):
    r"""Relative difference between the objective a checkpoint recorded and
    the value the warm-started link evaluates first: zero when the link
    continues the same minimisation from the same point."""
    return abs(float(first_total) - float(recorded_total)) / max(
        abs(float(recorded_total)), 1e-300)


def frozen_in_control(recorded, key, friction_control):
    r"""The positive value ``recorded`` holds for ``key`` (an auto weight in
    the friction control's units: ``prior_sigma_alpha``, ``friction_c_ref``)
    if the record was minimised on the same ``friction_control``, else None.
    A record without a ``friction_control`` was inverted on the log
    control."""
    if str(recorded.get("friction_control", "log")) != str(friction_control):
        return None
    try:
        value = float(recorded.get(key, 0.0) or 0.0)
    except (TypeError, ValueError):
        return None
    return value if value > 0.0 else None
