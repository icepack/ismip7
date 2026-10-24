r"""``optimization.TrialFailures`` against scipy's L-BFGS-B, which the
inversion's scipy path runs with ``gtol=0``.

A failed forward or adjoint at a line-search trial returns ten times the last
good objective and a zero gradient. Issue 153's RC link 11461566 hit three in
a row at evaluation 107 and L-BFGS-B stopped with a CONVERGENCE message at a
gradient norm of 1.03; the chain then marked the MAP done. The tracker has to
tell that stop apart from a minimisation that met failed trials along the way
and recovered.

``optimization.AcceptedIterate`` keeps the evaluations a MAP can be published
at. scipy's result.x is rarely the last evaluation, and after failed trials
it can be a failed trial itself; the publishing solve and the MAP's record
need a successful evaluation's state and objective.
"""
import numpy as np
from scipy.optimize import minimize

from icepack2_tools.optimization import AcceptedIterate, TrialFailures

HESSIAN = np.diag(np.linspace(1.0, 50.0, 20))


def _minimise(fails, regularisation=0.0, scale_from="total"):
    r"""L-BFGS-B on a quadratic misfit plus a constant ``regularisation``, as
    the driver calls it; ``fails(n, x, x_good)`` says whether evaluation ``n``
    at ``x`` fails, and a failed trial scales from the last good
    ``scale_from`` ("misfit" or "total")."""
    tracker = TrialFailures()
    good = {"misfit": None, "total": None, "x": None, "n": 0}

    def objective_and_gradient(x):
        good["n"] += 1
        if good["x"] is not None and fails(good["n"], x, good["x"]):
            return tracker.failed(good[scale_from], x.size)
        misfit = 0.5 * x @ HESSIAN @ x
        good.update(misfit=misfit, total=misfit + regularisation, x=x.copy())
        tracker.succeeded()
        return good["total"], HESSIAN @ x

    result = minimize(objective_and_gradient, np.ones(20), jac=True,
                      method="L-BFGS-B",
                      options={"maxiter": 300, "ftol": 1e-12, "gtol": 0})
    return result, tracker, float(np.linalg.norm(HESSIAN @ result.x)), good["x"]


def test_a_stop_after_failed_trials_is_flagged():
    r"""Every trial fails from evaluation 13 on, as when the direct forward
    stalls at the residual floor near the accepted point. scipy calls the
    stop a success; the gradient there is far from small."""
    result, tracker, grad_norm, _x_good = _minimise(lambda n, x, x_good: n > 12)
    assert result.success
    assert grad_norm > 0.1
    assert tracker.stopped_on_failures
    assert tracker.trailing == tracker.total > 0


def test_failed_trials_the_line_search_recovers_from_are_not_flagged():
    r"""Trials that step too far fail and the line search backtracks past
    them; the minimisation converges and the tracker leaves it final."""
    result, tracker, grad_norm, _x_good = _minimise(
        lambda n, x, x_good: np.linalg.norm(x - x_good) > 0.05)
    assert tracker.total > 0
    assert grad_norm < 1e-4
    assert not tracker.stopped_on_failures


def test_a_minimisation_without_failures_is_not_flagged():
    result, tracker, grad_norm, _x_good = _minimise(lambda n, x, x_good: False)
    assert grad_norm < 1e-4
    assert tracker.total == 0
    assert not tracker.stopped_on_failures


def test_a_stop_after_failed_trials_holds_the_last_good_point_under_heavy_regularisation():
    r"""Regularisation a hundred times the starting misfit and every trial
    failing from evaluation 13 on. Scaled from the misfit alone, a failed
    trial sits below the current total, L-BFGS-B accepts it and the result is
    a point whose forward failed; scaled from the total, every failed trial is
    rejected and the result is the last point whose forward succeeded."""
    regularisation = 100 * 0.5 * np.trace(HESSIAN)
    for scale_from, holds in (("misfit", False), ("total", True)):
        result, tracker, _grad_norm, x_good = _minimise(
            lambda n, x, x_good: n > 12, regularisation=regularisation,
            scale_from=scale_from)
        assert tracker.stopped_on_failures
        assert np.array_equal(result.x, x_good) == holds, scale_from


def test_a_failed_trial_returns_ten_times_the_last_objective_and_no_gradient():
    objective, gradient = TrialFailures().failed(3.0, 4)
    assert objective == 30.0
    assert np.array_equal(gradient, np.zeros(4))


def _minimise_recorded(fails, regularisation=0.0, scale_from="total"):
    r"""``_minimise`` with every good evaluation recorded in an
    ``AcceptedIterate`` (its total as the payload, where the driver records
    its mixed state) and the record's callback; ``fails(n, x)``."""
    tracker, record = TrialFailures(), AcceptedIterate()
    good = {"misfit": None, "total": None, "x": None, "n": 0, "totals": []}

    def objective_and_gradient(x):
        good["n"] += 1
        if good["x"] is not None and fails(good["n"], x):
            return tracker.failed(good[scale_from], x.size)
        misfit = 0.5 * x @ HESSIAN @ x
        good.update(misfit=misfit, total=misfit + regularisation, x=x.copy())
        good["totals"].append(good["total"])
        record.evaluated(good["x"], good["total"], good["total"])
        tracker.succeeded()
        return good["total"], HESSIAN @ x

    result = minimize(objective_and_gradient, np.ones(20), jac=True,
                      method="L-BFGS-B", callback=record.iterate,
                      options={"maxiter": 300, "ftol": 1e-12, "gtol": 0})
    return result, tracker, record, good


def test_the_returned_controls_are_an_evaluation_other_than_the_last():
    r"""Every trial past the wall ``x.sum() < wall`` fails, as a forward that
    stops converging past some controls (the minimiser at 0 lies past it).
    At some walls L-BFGS-B stops right after a failed trial, at others it
    stops ABNORMAL with none trailing. Which wall does which belongs to the
    installed scipy's L-BFGS-B: the test was written where 5 and 2 gave the
    two stops, and scipy 1.15.3 stops right after a failed trial at 6.5 but
    not at 5. So the test takes the first wall of each kind, trying 5 and 2
    first. Both times result.x is an earlier evaluation than the last good
    one, whose state the publishing solve used to start from, and the record
    publishes that evaluation; after the failure stop result.fun is the
    failed trial's inflated value."""
    walls = [5.0, 2.0] + [float(w) for w in np.arange(0.5, 10.01, 0.5)
                          if w not in (5.0, 2.0)]
    seen = {}
    for wall in walls:
        result, tracker, record, good = _minimise_recorded(
            lambda n, x: x.sum() < wall)
        stopped = tracker.stopped_on_failures
        if stopped in seen:
            continue
        seen[stopped] = wall
        assert not np.array_equal(result.x, good["x"]), wall
        x_pub, total, how = record.publish(result.x)
        assert how == "returned" and np.array_equal(x_pub, result.x), wall
        assert total == 0.5 * result.x @ HESSIAN @ result.x, wall
        if stopped:
            assert result.fun > total
        if len(seen) == 2:
            break
    assert set(seen) == {True, False}, seen


def test_a_failed_trial_returned_publishes_the_lowest_evaluation():
    r"""Under the misfit scaling with heavy regularisation L-BFGS-B accepts
    failed trials and returns one, as the 2 km re-inversions 12286864 and
    12291649 returned theirs: no evaluation is at result.x, so the record
    publishes the lowest-objective evaluation, which the publishing solve
    has a converged state for."""
    result, tracker, record, good = _minimise_recorded(
        lambda n, x: n > 12, regularisation=100 * 0.5 * np.trace(HESSIAN),
        scale_from="misfit")
    assert tracker.stopped_on_failures
    x_pub, total, how = record.publish(result.x)
    assert how == "lowest"
    assert not np.array_equal(x_pub, result.x)
    assert total == min(good["totals"])


def test_the_starting_point_is_accepted_until_a_callback_moves_it():
    record = AcceptedIterate()
    x0, x1, x2 = np.zeros(3), np.ones(3), np.full(3, 2.0)
    assert record.publish(x0) is None
    record.evaluated(x0, 3.0, "start")
    record.evaluated(x1, 5.0, "trial")
    assert record.publish(x0)[1:] == ("start", "returned")
    record.iterate(x2)
    assert record.publish(x0)[1:] == ("start", "returned")
    record.iterate(x1.copy())
    record.evaluated(x2, 4.0, "next")
    assert record.publish(x1)[1:] == ("trial", "returned")
    # x0 is still the lowest evaluation, and what a failed point publishes
    assert record.publish(x0)[1:] == ("start", "returned")
    assert record.publish(np.full(3, 9.0))[1:] == ("start", "lowest")
