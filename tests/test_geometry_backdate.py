r"""A cold start before the 2015 geometry undoes the Smith et al. (2020) mean
thickness change on grounded ice (issue #117): the start years that ask for it,
and the cells it touches."""
import numpy as np
import pytest

from icepack2_tools.obs_dhdt import backdate_thickness
from icepack2_tools.runconfig import geometry_backdate_years

RHO = 917.0 / 1024.0


def test_a_2003_start_undoes_twelve_years(monkeypatch):
    monkeypatch.delenv("ISMIP7_GEOMETRY_BACKDATE", raising=False)
    assert geometry_backdate_years(2003.0) == 12.0
    assert geometry_backdate_years(2010.0) == 5.0
    assert geometry_backdate_years(2015.0) == 0.0     # the projections and the control
    assert geometry_backdate_years(2016.0) == 0.0


def test_a_start_before_the_dhdt_window_is_refused_unless_named(monkeypatch):
    monkeypatch.delenv("ISMIP7_GEOMETRY_BACKDATE", raising=False)
    with pytest.raises(ValueError, match="past the window"):
        geometry_backdate_years(1850.0)
    monkeypatch.setenv("ISMIP7_GEOMETRY_BACKDATE", "0")
    assert geometry_backdate_years(1850.0) == 0.0
    monkeypatch.setenv("ISMIP7_GEOMETRY_BACKDATE", "7.5")
    assert geometry_backdate_years(2003.0) == 7.5


def test_only_grounded_observed_cells_change():
    # cells: grounded thinning, grounded thickening, floating, grounded without
    # observations, grounded thickening that would go negative
    h = np.array([1000.0, 500.0, 400.0, 800.0, 5.0])
    bed = np.array([-200.0, 100.0, -600.0, -100.0, 50.0])
    dhdt = np.array([-2.0, 0.5, -3.0, -1.0, 1.0])
    observed = np.array([1.0, 1.0, 1.0, 0.0, 1.0])
    new, changed = backdate_thickness(h, bed, dhdt, observed, 12.0, RHO)
    assert list(changed) == [True, True, False, False, True]
    assert new[0] == pytest.approx(1024.0)      # thinning undone: thicker in 2003
    assert new[1] == pytest.approx(494.0)       # thickening undone: thinner
    assert new[2] == 400.0                      # floating: kept (noisy shelf altimetry)
    assert new[3] == 800.0                      # no observations: kept
    assert new[4] == 0.0                        # floored at zero
    assert np.array_equal(h, [1000.0, 500.0, 400.0, 800.0, 5.0])   # input untouched
