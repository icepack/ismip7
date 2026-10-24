r"""A retreat-only front (ISMIP7_FRONT_ADVANCE=none): the law retreats the
front and nothing advances past the t=0 extent."""
import numpy as np
import pytest

from icepack2_tools import runconfig
from icepack2_tools.front import front_removal_mask


def test_the_knob_parses_and_refuses_the_rest(monkeypatch):
    monkeypatch.delenv("ISMIP7_FRONT_ADVANCE", raising=False)
    assert runconfig.front_advance() == "free"
    monkeypatch.setenv("ISMIP7_FRONT_ADVANCE", "None")
    assert runconfig.front_advance() == "none"
    monkeypatch.setenv("ISMIP7_FRONT_ADVANCE", "retreat")
    with pytest.raises(ValueError):
        runconfig.front_advance()


def test_retreat_only_removes_beyond_the_initial_extent_as_well():
    law = np.array([False, True, False, False])     # the front passed cell 1
    init = np.array([False, False, True, True])     # cells 2, 3 were ice-free at t=0
    assert np.array_equal(front_removal_mask(law, init, False), law)
    assert np.array_equal(front_removal_mask(law, init, True),
                          np.array([False, True, True, True]))
    assert np.array_equal(front_removal_mask(law, None, True), law)


def test_retreat_only_holds_the_t0_extent_on_the_sea_only():
    r"""ISMIP7_FRONT_MARINE: a land-terminating front advances with the ice,
    so a cell beyond the t=0 extent on land keeps the ice that reaches it,
    and one in the sea is emptied as before."""
    law = np.array([False, True, False, False])     # the front passed cell 1
    init = np.array([False, False, True, True])     # cells 2, 3 were ice-free at t=0
    sea = np.array([True, True, True, False])       # cell 3 is land
    assert np.array_equal(front_removal_mask(law, init, True, sea=sea),
                          np.array([False, True, True, False]))
    # the law's own removals and the held cells are not the extent's
    calved = np.array([True, False, False, False])
    assert np.array_equal(front_removal_mask(law, init, True, calved, sea=sea),
                          np.array([True, True, True, False]))
    assert np.array_equal(front_removal_mask(law, init, False, sea=sea), law)


def test_land_fronts_do_not_calve_by_default(monkeypatch):
    r"""ISMIP7_FRONT_MARINE: on by default, every law made with marine_only
    unless its parameters say otherwise; 0 makes the law as it is."""
    calving = pytest.importorskip("icepack_tools.calving")
    if "marine_only" not in calving.make("hfb").p:
        pytest.skip("icepack_tools predates the land gate on every law")
    monkeypatch.delenv("ISMIP7_FRONT_MARINE", raising=False)
    monkeypatch.setenv("ISMIP7_CALVING", "hfb")
    monkeypatch.setenv("ISMIP7_CALVING_PARAMS", "sigma_max=0.15")
    assert runconfig.front_marine()
    law = runconfig.calving_law_object()
    assert law.p["marine_only"] is True
    assert "land-terminating fronts do not calve" in law.describe()
    monkeypatch.setenv("ISMIP7_CALVING_PARAMS", "sigma_max=0.15,marine_only=false")
    assert runconfig.calving_law_object().p["marine_only"] is False
    monkeypatch.setenv("ISMIP7_FRONT_MARINE", "0")
    monkeypatch.setenv("ISMIP7_CALVING_PARAMS", "sigma_max=0.15")
    assert not runconfig.front_marine()
    assert runconfig.calving_law_object().p["marine_only"] is False


def test_a_cell_the_front_emptied_stays_empty():
    from icepack2_tools.front import held_calved
    none = np.zeros(4, bool)
    h_old = np.array([300.0, 0.5, 0.0, 120.0])       # cell 1 was below the threshold already
    h_new = np.array([0.0, 0.0, 8.0, 119.0])         # 0 emptied; 2 was ice-free and took inflow
    beyond = np.array([True, True, False, False])    # the law passed cells 0 and 1
    calved = held_calved(none, h_old, h_new, 1.0, sliver=none, beyond=beyond,
                         calv_frac=np.zeros(4))
    assert np.array_equal(calved, [True, False, False, False])
    # once held, the mask only grows; a shed of the whole cell holds it too
    calved = held_calved(calved, np.array([0.0, 0.0, 8.0, 119.0]), np.array([5.0, 0.0, 0.0, 118.0]),
                         1.0, sliver=none, beyond=none, calv_frac=np.array([0, 0, 1.0, 0.2]))
    assert np.array_equal(calved, [True, False, True, False])
    law = np.array([False, False, False, False])
    assert np.array_equal(front_removal_mask(law, None, True, calved), calved)
    assert np.array_equal(front_removal_mask(law, None, False, calved), law)


def test_only_a_sliver_at_the_front_is_held():
    """A sliver is held where the law sheds (a front cell); the same sliver
    left by SMB or melt away from the front, or a cell melted to zero, is not."""
    from icepack2_tools.front import held_calved, retreat_slivers
    h_old = np.array([40.0, 1.02, 30.0, 40.0])
    h_after = np.array([0.6, 0.98, 0.5, 0.0])        # front sliver; SMB sliver; melt sliver; melted through
    sliver = retreat_slivers(h_after, h_old, 1.0)
    h_new = np.where(sliver, 0.0, h_after)
    calv_frac = np.array([0.3, 0.0, 0.0, 0.0])
    calved = held_calved(np.zeros(4, bool), h_old, h_new, 1.0, sliver=sliver,
                         beyond=np.zeros(4, bool), calv_frac=calv_frac)
    assert np.array_equal(calved, [True, False, False, False])
    # without a shedding law there is no front cell, so no sliver is held
    assert not held_calved(np.zeros(4, bool), h_old, h_new, 1.0, sliver=sliver,
                           beyond=np.zeros(4, bool)).any()


def test_the_front_anchor_defaults_to_extent_and_takes_advect(monkeypatch):
    r"""ISMIP7_FRONT_ANCHOR: extent by default (the front re-derived from the
    thickness), advect for a front carried by the ice velocity and the
    calving rate; anything else is refused at startup."""
    from icepack2_tools.runconfig import front_anchor
    monkeypatch.delenv("ISMIP7_FRONT_ANCHOR", raising=False)
    assert front_anchor() == "extent"
    monkeypatch.setenv("ISMIP7_FRONT_ANCHOR", " Advect ")
    assert front_anchor() == "advect"
    monkeypatch.setenv("ISMIP7_FRONT_ANCHOR", "lagrangian")
    import pytest
    with pytest.raises(ValueError, match="ISMIP7_FRONT_ANCHOR"):
        front_anchor()


def test_the_advected_front_reinitialises_on_model_time(monkeypatch):
    r"""ISMIP7_FRONT_REINIT_YR: 2.5 model years by default (every fifth
    advance stood the forward's front still), a positive number otherwise."""
    from icepack2_tools.runconfig import front_reinit_interval
    import pytest
    monkeypatch.delenv("ISMIP7_FRONT_REINIT_YR", raising=False)
    assert front_reinit_interval() == 2.5
    monkeypatch.setenv("ISMIP7_FRONT_REINIT_YR", " 0.5 ")
    assert front_reinit_interval() == 0.5
    for bad in ("0", "-1", "yearly", "nan"):
        monkeypatch.setenv("ISMIP7_FRONT_REINIT_YR", bad)
        with pytest.raises(ValueError, match="ISMIP7_FRONT_REINIT_YR"):
            front_reinit_interval()

