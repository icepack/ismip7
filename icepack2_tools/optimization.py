r"""Stopping rules for the inversions, shared by the TAO and scipy paths, and
the log-velocity weight a warm start carries into an inversion."""
import math


class FunctionalDecreaseStop:
    r"""Relative-decrease stopping rule, for a TAO monitor to apply.

    PETSc's TAO tests only gradient norms. This stops on the relative decrease
    of the functional between successive iterates,

        (J_old - J_new) / max(|J_old|, |J_new|, 1) <= ftol,

    and never before ``min_iter`` iterations. scipy's L-BFGS-B ``ftol`` is the
    same expression, so both optimizer paths read one knob. 1e-10 converges a
    production inversion; the L-surface sweeps pass 1e-4.
    """

    def __init__(self, ftol, min_iter=3):
        self.ftol = float(ftol)
        self.min_iter = int(min_iter)
        self.criterion = None
        self._last = None

    def update(self, it, f):
        r"""Record iterate ``it``'s functional ``f``; True when the rule says stop."""
        f = float(f)
        if self._last is None:
            self._last = f
            return False
        self.criterion = (self._last - f) / max(abs(self._last), abs(f), 1.0)
        self._last = f
        return self.ftol > 0 and it >= self.min_iter and self.criterion <= self.ftol


# The root attributes that say which log-velocity term a MAP was minimised
# under: the weight, and the two knobs that set the scale of what it weighs
# (misfit_norm the units of the chi^2 term, log_vel_eps the floor inside the
# logarithm). The inversion's save_map writes all three.
LOG_VEL_ATTRS = ("log_vel_weight", "log_vel_eps", "misfit_norm")


def recorded_objective(chk):
    r"""The log-velocity attributes a MAP or inversion checkpoint records.

    ``chk`` is an open ``firedrake.CheckpointFile``, or anything with its
    ``has_attr``/``get_attr``. An attribute the file lacks is left out: MAPs
    written before 5 September 2026 carry no weight, and an adapted-mesh
    checkpoint keeps its geometry attributes only.
    """
    out = {}
    for key in LOG_VEL_ATTRS:
        if chk.has_attr("/", key):
            value = chk.get_attr("/", key)
            out[key] = value.decode() if isinstance(value, bytes) else value
    return out


def _held_weight(recorded, misfit_norm, eps):
    r"""The warm start's weight, or None and why an inversion cannot hold it."""
    if recorded is None:
        return None, ""
    if "log_vel_weight" not in recorded:
        return None, "it records no log-velocity weight"
    weight = float(recorded["log_vel_weight"])
    if not (math.isfinite(weight) and weight > 0.0):
        return None, f"it records weight {weight:g}"
    norm = recorded.get("misfit_norm")
    if norm is None or str(norm).lower() != misfit_norm.lower():
        return None, f"misfit_norm {norm} there, {misfit_norm} here"
    rec_eps = recorded.get("log_vel_eps")
    if rec_eps is None or not math.isclose(float(rec_eps), eps, rel_tol=1e-9):
        return None, f"log_vel_eps {rec_eps} there, {eps:g} here"
    return weight, ""


def resolve_log_vel_weight(requested, derived, recorded, *, misfit_norm, eps):
    r"""``ISMIP7_LOG_VEL_WEIGHT`` resolved, as ``(weight, source, note)``.

    ``requested`` is the knob, ``auto`` or a number. ``derived`` is the
    velocity chi^2 term over the unweighted log term at the state this run
    starts from, and is read under ``auto`` only. ``recorded`` is
    ``recorded_objective`` of the warm start, None without one.

    Under ``auto``, a warm start that records a positive weight under the same
    misfit norm and eps supplies that weight (source ``warm_start``). A chained
    link resumes from the checkpoint the link before it wrote, so every link
    minimises the objective the first one set; re-deriving at the checkpoint
    gave 2495 where the first link of IU's 32 km chain had 17452 (issue 68).
    Otherwise the weight is ``derived`` (source ``derived``), and ``note``
    says why the warm start's weight was not taken. A number is used as given
    (source ``requested``), and ``note`` says so when it replaces a different
    weight the warm start was minimised under.
    """
    held, why = _held_weight(recorded, misfit_norm, eps)
    if str(requested).strip().lower() == "auto":
        if held is not None:
            return held, "warm_start", ""
        note = f"the warm start's weight is not used: {why}" if why else ""
        return float(derived), "derived", note
    weight = float(requested)
    if held is not None and not math.isclose(weight, held, rel_tol=1e-9):
        return weight, "requested", (
            f"it replaces the weight {held:.6g} the warm start was minimised under")
    return weight, "requested", ""
