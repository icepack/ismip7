r"""Stopping rules for the inversions, shared by the TAO and scipy paths, and
the log-velocity weight a warm start carries into an inversion."""
import math

import numpy as np

# What the driver prints, ahead of its "Saved MAP:" line, when the optimizer
# stopped right after failed trial evaluations. The chain runner reads it and
# leaves the MAP without a done marker, so the next link resumes from it.
NOT_FINAL = "MAP not final:"


class TrialFailures:
    r"""Failed trial evaluations of the scipy path since the last one that
    succeeded.

    A forward or adjoint that fails at a line-search trial returns
    ``failed(J_last, n)``: ten times the last good total objective (misfit
    and regularisation, the value scipy minimises) and a zero gradient, so
    L-BFGS-B backtracks. A misfit alone would sit below the current total
    wherever regularisation is more than nine times the misfit, and the line
    search would accept the failed trial. When every trial near the current
    iterate fails, scipy still ends with a CONVERGENCE message (the projected
    gradient below pgtol, read off the zero gradient, or the relative
    reduction below factr) at a point whose gradient is not small: RC's
    final-product link 11461566 (issue #153) stopped that way at evaluation
    107 with a gradient norm of 1.03, its objective still falling about 1 an
    evaluation. ``stopped_on_failures`` after the optimizer returns says its
    stop came that way, and the MAP is then an evaluated point of an
    unfinished minimisation (``AcceptedIterate``)."""

    def __init__(self):
        self.trailing = 0
        self.total = 0

    def failed(self, last_objective, size):
        r"""Count a failed trial; the (objective, gradient) to return for it."""
        self.trailing += 1
        self.total += 1
        return last_objective * 10, np.zeros(size)

    def succeeded(self):
        self.trailing = 0

    @property
    def stopped_on_failures(self):
        return self.trailing > 0


class AcceptedIterate:
    r"""The scipy path's record of the evaluations a MAP can be published at.

    The driver keeps one mixed state, the last successful evaluation's, and
    scipy's ``result.x`` is rarely that evaluation. Over 600 wall-limited
    quadratics (scipy 1.17, every trial past a wall failing) it was the last
    evaluation in 5 and the iterate the L-BFGS-B callback accepted in all
    600. At 2 km both failure stops of 16 Oct 2026 (the DG1 velocity-only
    re-inversions 12286864 and 12291649) returned the last FAILED trial: the
    publishing solve, started from the same last good state as that trial,
    diverged as it did (||F|| 1.841e27 at Newton iteration 11 in 12291649),
    and no velocity was published.

    ``evaluated(x, value, payload)`` records a successful evaluation at ``x``
    with objective ``value``; ``iterate(xk)`` is the optimiser callback, at
    which the most recent evaluation is the iterate just accepted. The first
    evaluation, at the starting point, is accepted until a callback says
    otherwise. ``publish(x)`` is ``(x, payload, how)`` for the evaluation to
    publish: the one at ``x`` (the accepted, lowest or last evaluation, how
    "returned"), else the lowest-objective evaluation (how "lowest"), where
    the publishing solve has that evaluation's converged state to start
    from; None before any evaluation. ``x`` is compared bitwise, and every
    rank holds the same global vector, so all ranks agree without a
    reduction.
    """

    def __init__(self):
        self.last = None
        self.accepted = None
        self.lowest = None

    def evaluated(self, x, value, payload):
        self.last = (x, float(value), payload)
        if self.accepted is None:
            self.accepted = self.last
        if self.lowest is None or self.last[1] < self.lowest[1]:
            self.lowest = self.last

    def iterate(self, xk):
        if self.last is not None and np.array_equal(xk, self.last[0]):
            self.accepted = self.last

    def publish(self, x):
        for entry in (self.accepted, self.lowest, self.last):
            if entry is not None and np.array_equal(entry[0], x):
                return entry[0], entry[2], "returned"
        if self.lowest is None:
            return None
        return self.lowest[0], self.lowest[2], "lowest"


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
    (source ``requested``), as Rice's 1 km inversion passes 85,380, and
    ``note`` says so when it differs by more than 0.1 percent from the weight
    the warm start was minimised under; a copy of that weight rounded by hand
    is the same objective.
    """
    held, why = _held_weight(recorded, misfit_norm, eps)
    if str(requested).strip().lower() == "auto":
        if held is not None:
            return held, "warm_start", ""
        note = f"the warm start's weight is not used: {why}" if why else ""
        return float(derived), "derived", note
    weight = float(requested)
    if held is not None and not math.isclose(weight, held, rel_tol=1e-3):
        return weight, "requested", (
            f"it replaces the weight {held:.6g} the warm start was minimised under")
    return weight, "requested", ""
