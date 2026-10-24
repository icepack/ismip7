"""Adaptive substepping controller (icepack2_tools/substep.py)."""
import numpy as np
import pytest

from icepack2_tools.solverconfig import substep_settings
from icepack2_tools.substep import SubstepController


def _run(ctrl, h, rate, dt, nsteps):
    """Drive a macro-step loop the way run_simulation does: each macro step
    of dt is taken as m substeps; a rejection rewinds the step and repeats it
    with the refined count. ``rate(h)`` is the lagged tendency (the velocity
    solved at the previous geometry). Returns the trajectory and the m used."""
    traj, ms = [h.copy()], []
    for _ in range(nsteps):
        entry = h.copy()
        m = ctrl.m
        while True:
            h = entry.copy()
            ctrl.begin_attempt()
            ok = True
            for _j in range(m):
                before = h.copy()
                h = h + (dt / m) * rate(h)
                if ctrl.exceeds(ctrl.observe(before, h, dt / m)):
                    ok = False
                    err = ctrl.err_step
                    break
            if ok:
                break
            m = ctrl.refine(err)
            assert m is not None, "substeps exhausted"
        ms.append(m)
        ctrl.accept()
        traj.append(h.copy())
    return np.array(traj), ms


def test_smooth_evolution_coarsens_to_one_substep():
    ctrl = SubstepController(tol=1.0, m_init=8, quiet_steps=3)
    h = np.full(5, 1000.0)
    _, ms = _run(ctrl, h, lambda h: np.full_like(h, -2.0), 0.05, 20)
    assert ms[0] == 8 and ms[-1] == 1


def test_flip_flop_is_rejected_and_refined():
    ctrl = SubstepController(tol=1.0, m_init=1)
    h0 = np.full(4, 500.0)
    ctrl.begin_attempt()
    ctrl.observe(h0, h0 + 3.0, 0.05)
    ctrl.accept()
    ctrl.begin_attempt()
    err = ctrl.observe(h0 + 3.0, h0 - 1.0, 0.05)    # flips sign and grows
    assert err == pytest.approx(3.5) and ctrl.exceeds(err)
    assert ctrl.refine(err) == 4                     # sqrt(3.5)/0.9 = 2.1: x4
    assert ctrl.refine(30.0) == 32                   # sqrt(30)/0.9 = 6.1: x8 (the cap)


def test_rejection_rewinds_the_history():
    ctrl = SubstepController(tol=1.0)
    h = np.full(3, 100.0)
    ctrl.begin_attempt()
    ctrl.observe(h, h + 1.0, 0.1)
    ctrl.accept()
    ctrl.begin_attempt()
    ctrl.observe(h + 1.0, h - 50.0, 0.1)
    ctrl.refine(ctrl.err_step)
    ctrl.begin_attempt()
    # measured against the accepted increment (+1 over 0.1), not the rejected one
    err = ctrl.observe(h + 1.0, h - 1.0, 0.05)
    assert err == pytest.approx(0.05 / 0.15 * abs(-2.0 - 0.5 * 1.0))


def test_a_damped_reversal_is_accepted():
    """-1 < g < 0: the increment reverses but shrinks, which is stable."""
    ctrl = SubstepController(tol=0.1)
    h = np.full(2, 900.0)
    ctrl.begin_attempt()
    ctrl.observe(h, h + 4.0, 0.05)
    ctrl.accept()
    ctrl.begin_attempt()
    assert ctrl.observe(h + 4.0, h + 1.0, 0.05) == 0.0


def test_a_rate_jump_without_reversal_is_accepted():
    """A cell going afloat speeds its thinning by 2500 m/yr: the increment
    jumps but keeps its sign, which no step size makes smoother."""
    ctrl = SubstepController(tol=1.0)
    h = np.full(2, 1500.0)
    tau = 0.0125
    ctrl.begin_attempt()
    ctrl.observe(h, h - 5.0 * tau, tau)
    ctrl.accept()
    ctrl.begin_attempt()
    err = ctrl.observe(h - 5.0 * tau, h - 5.0 * tau - 2505.0 * tau, tau)
    assert err == 0.0 and not ctrl.exceeds(err)


def test_removed_cells_do_not_count():
    ctrl = SubstepController(tol=1.0, hmin=10.0)
    h = np.array([800.0, 300.0])
    ctrl.begin_attempt()
    ctrl.observe(h, h - 1.0, 0.05)
    ctrl.accept()
    ctrl.begin_attempt()
    calved = np.array([800.1, 0.0])                  # cell 1 removed by calving
    assert ctrl.observe(h - 1.0, calved, 0.05) == pytest.approx(1.05)


def test_skipped_cells_sit_out_this_comparison_and_the_next():
    """A cell an advected front crossed changes its transport with its side
    (only the cells inside export): it is left out of the comparison of the
    substep that crossed it and of the next one, which measures against it."""
    ctrl = SubstepController(tol=1.0)
    h = np.full(2, 400.0)
    ctrl.begin_attempt()
    ctrl.observe(h, h + 0.5, 0.05)
    ctrl.accept()
    ctrl.begin_attempt()
    crossed = np.array([False, True])
    flipped = h + 0.5 + np.array([0.6, -4.0])        # cell 1 stops importing
    err = ctrl.observe(h + 0.5, flipped, 0.05, skip=crossed)
    assert err == 0.0
    # the next comparison is against the skipped increment: still out
    assert ctrl.observe(flipped, flipped + np.array([0.6, 9.0]), 0.05) == 0.0
    # and the one after that counts it again
    later = flipped + np.array([0.6, 9.0])
    assert ctrl.observe(later, later + np.array([0.6, -20.0]), 0.05) == pytest.approx(14.5)


def test_ceiling_and_nonfinite():
    ctrl = SubstepController(tol=1.0, m_init=4, m_max=4)
    ctrl.begin_attempt()
    assert ctrl.refine(None) is None
    assert ctrl.exceeds(float("nan")) and ctrl.exceeds(float("inf"))
    ctrl = SubstepController(tol=1.0)
    h = np.full(2, 100.0)
    ctrl.begin_attempt()
    ctrl.observe(h, h + 1.0, 0.1)
    ctrl.begin_attempt()
    assert ctrl.exceeds(ctrl.observe(h + 1.0, np.array([101.0, np.nan]), 0.1))


def test_smooth_decay_never_rejects():
    """Backward-Euler decay of a thick column at a coarse step: no reversal."""
    ctrl = SubstepController(tol=1e-6)
    h = np.array([1000.0, 400.0])
    for _ in range(10):
        new = 50.0 + (h - 50.0) / (1 + 0.5)
        ctrl.begin_attempt()
        assert ctrl.observe(h, new, 0.5) == 0.0
        ctrl.accept()
        h = new


def test_lagged_coupling_blows_up_fixed_but_not_adaptive():
    """h relaxes to h* at rate lam with the tendency lagged one step (forward
    Euler): unstable for dt*lam > 2, the flip-flop the forward shows."""
    lam, dt, hstar = 100.0, 0.05, 1000.0               # dt*lam = 5
    rate = lambda h: -lam * (h - hstar)
    h0 = np.array([hstar + 5.0, hstar - 5.0, hstar])
    h = h0.copy()
    for _ in range(40):
        h = h + dt * rate(h)
    assert np.abs(h - hstar).max() > 1e6               # the fixed step diverges
    ctrl = SubstepController(tol=1.0, quiet_steps=1000)
    traj, ms = _run(ctrl, h0.copy(), rate, dt, 40)
    assert np.abs(traj[-1] - hstar).max() < 1.0
    # the run's first substep has no history to compare with; every later
    # macro step runs with tau*lam < 2, i.e. m > 2.5
    assert min(ms[1:]) >= 4


def test_settings_from_environment(monkeypatch):
    monkeypatch.delenv("ISMIP7_SUBSTEP_ADAPT", raising=False)
    assert substep_settings() is None
    monkeypatch.setenv("ISMIP7_SUBSTEP_ADAPT", "1")
    monkeypatch.setenv("ISMIP7_SUBSTEP_TOL", "0.5")
    monkeypatch.setenv("ISMIP7_SUBSTEP_INIT", "4")
    s = substep_settings()
    assert s["tol"] == 0.5 and s["m_init"] == 4 and s["m_max"] == 64
    assert s["ceiling"] == 20.0
    SubstepController(**s)


def test_a_finite_estimate_at_the_ceiling_is_accepted():
    # refining helps a little each time but never enough: at m_max a finite
    # estimate within the ceiling multiple is accepted instead of stalling
    ctrl = SubstepController(tol=1.0, m_init=2, m_max=4, ceiling=20.0)
    ctrl.begin_attempt()
    assert ctrl.exceeds(3.0)
    assert ctrl.refine(3.0) == 4
    ctrl.begin_attempt()
    # 1.5 m: refining halved the estimate (so not the ineffective rule) but it is still above tol
    assert not ctrl.exceeds(1.5) and ctrl.ceiling_now and ctrl.ceiling_accepted == 1
    ctrl.accept()
    assert not ctrl.ceiling_now
    # beyond the multiple, or non-finite, it still rejects and refine() is exhausted
    ctrl.begin_attempt()
    assert ctrl.exceeds(25.0) and ctrl.refine(25.0) is None
    assert ctrl.exceeds(float("inf"))
    # ceiling 0 is the old stall
    ctrl = SubstepController(tol=1.0, m_init=4, m_max=4, ceiling=0.0)
    ctrl.begin_attempt()
    assert ctrl.exceeds(2.0) and ctrl.refine(2.0) is None
    with pytest.raises(ValueError):
        SubstepController(tol=1.0, ceiling=-1.0)


def test_refinement_that_does_not_help_is_accepted():
    """A reversal that persists at any step (a grounding-zone cell flickering
    across flotation): once doubling m fails to cut the estimate by 30 %, the
    retry is accepted rather than doubled again; a blow-up still rejects."""
    ctrl = SubstepController(tol=1.0)
    h = np.full(1, 700.0)
    ctrl.begin_attempt()
    ctrl.observe(h, h + 2.0, 0.0125)
    ctrl.accept()
    ctrl.begin_attempt()
    err = ctrl.observe(h + 2.0, h - 0.8, 0.0125)
    assert ctrl.exceeds(err)
    ctrl.refine(err)
    ctrl.begin_attempt()
    assert not ctrl.exceeds(0.9 * err) and ctrl.tolerated == 1
    assert ctrl.exceeds(0.5 * err)                   # refining helped: keep going
    assert ctrl.exceeds(float("inf"))
    ctrl.accept()
    assert ctrl.m == 1                               # back to the count it rejected
    ctrl.begin_attempt()
    assert not ctrl.exceeds(err)                     # the same estimate is carried
    assert ctrl.exceeds(1.1 * err)                   # a larger one tests refining again
    ctrl.accept()
    ctrl.begin_attempt()
    ctrl.err_step = 0.5
    ctrl.accept()                                    # within the tolerance: the rule ends
    ctrl.begin_attempt()
    assert ctrl.exceeds(0.9 * err)


def test_a_persistent_flicker_keeps_m_bounded():
    """The 2004.5 flicker: an estimate that refining barely reduces must not
    ratchet m up one doubling per macro step until the ceiling stops the run."""
    estimate = lambda m: 1.1 + 1.2 / m               # 1.4, 1.25, 1.175 at 4, 8, 16
    ctrl = SubstepController(tol=1.0, m_init=4, m_max=64)
    ms = []
    for _ in range(30):
        m = ctrl.m
        while True:
            ctrl.begin_attempt()
            err = estimate(m)
            ctrl.err_step = err
            if not ctrl.exceeds(err):
                break
            m = ctrl.refine(err)
            assert m is not None, "substeps exhausted"
        ms.append(m)
        ctrl.accept()
    assert max(ms) == 8 and ms[1:] == [4] * 29 and ctrl.rejections == 1


def test_location_of_the_largest_estimate():
    ctrl = SubstepController(tol=1.0)
    h = np.full(3, 500.0)
    xy = np.array([[0.0, 0.0], [1e3, 2e3], [3e3, 4e3]])
    ctrl.begin_attempt()
    ctrl.observe(h, h + np.array([1.0, 1.0, 1.0]), 0.05, xy)
    ctrl.accept()
    ctrl.begin_attempt()
    ctrl.observe(h + 1.0, h + 1.0 + np.array([-1.5, -3.0, 0.5]), 0.05, xy)
    assert ctrl.err_xy == (1e3, 2e3)
