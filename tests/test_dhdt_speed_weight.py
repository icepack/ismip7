r"""``obs_dhdt.speed_weight``, the velocity terms' weight by observed speed
under ``ISMIP7_DHDT_SPEED_RATIO``, and its place in the objective a warm
start's handoff compares."""
import numpy as np
import pytest

from icepack2_tools.handoff import objective_mismatches
from icepack2_tools.obs_dhdt import speed_weight


def test_the_weight_is_one_on_fast_ice_and_one_over_the_ratio_on_slow_ice():
    w = speed_weight([0.0, 5.0, 10.0, 100.0, 300.0], 50.0, 10.0, 100.0)
    assert np.allclose(w, [1 / 50, 1 / 50, 1 / 50, 1.0, 1.0])


def test_the_weight_is_log_linear_in_speed_between_the_two_speeds():
    # at the geometric mean of the two speeds the weight is ratio^-1/2
    assert np.isclose(speed_weight([np.sqrt(10.0 * 100.0)], 100.0, 10.0, 100.0)[0], 0.1)
    log_w = np.log(speed_weight(np.geomspace(10.0, 100.0, 7), 30.0, 10.0, 100.0))
    assert np.all(np.diff(log_w) > 0)
    assert np.allclose(np.diff(log_w, 2), 0.0)


def test_a_ratio_of_one_weighs_every_speed_alike():
    assert np.array_equal(speed_weight([0.0, 3.0, 50.0, 1e4], 1.0, 10.0, 100.0), np.ones(4))


def test_the_weight_refuses_a_ratio_below_one_or_crossed_speeds():
    with pytest.raises(ValueError):
        speed_weight([1.0], 0.5, 10.0, 100.0)
    with pytest.raises(ValueError):
        speed_weight([1.0], 10.0, 100.0, 10.0)


def test_a_changed_speed_ratio_is_another_objective():
    assert objective_mismatches({"dhdt_speed_ratio": 1.0}, {"dhdt_speed_ratio": 50.0})
    assert objective_mismatches(
        {"dhdt_speed_ratio": 50.0, "dhdt_speed_fast": 100.0},
        {"dhdt_speed_ratio": 50.0, "dhdt_speed_fast": 300.0})
    # a MAP written before the knob records no ratio, so none is compared
    assert not objective_mismatches({}, {"dhdt_speed_ratio": 50.0})
