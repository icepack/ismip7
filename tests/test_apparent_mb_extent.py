r"""The frozen apparent-MB reference follows the LIVE ice extent.

Exercises :func:`icepack2_tools.front.clear_reference_where_ice_free`, the rule
that keeps ``a_ref`` off cells that hold no ice. Without it a free calving law
cannot retreat: ``a_ref`` at a t=0 front cell is the terminus outflow, large
and positive, so a cell the front calves is regrown by the transport source on
the very next step and re-booked to calving, oscillating at the t=0 front.

The reference is the DG0 flux divergence, so these are the per-cell arrays the
transport actually adds to its source term.

Serial, no data files, no run setup.
"""

import numpy as np

from icepack2_tools.front import clear_reference_where_ice_free


def test_calved_cell_inside_t0_extent_is_not_regrown():
    r"""The reported failure: a cell INSIDE the t=0 extent that the front has
    just calved keeps a large positive reference, so the next step's source
    regrows it above the extent threshold."""
    front_hmin, dt = 1.0, 0.1
    # Cell 1 sat at the t=0 front, so its reference is the terminus outflow.
    a_ref = np.array([0.5, 40.0, 0.3])
    h = np.array([500.0, 0.0, 480.0])          # cell 1 was just calved
    ice_free = np.array([False, True, False])

    # Before the fix a_ref was masked only at t=0, so cell 1 kept 40 m/yr.
    regrown = h + a_ref * dt
    assert regrown[1] > front_hmin, "guard: the un-masked source does regrow it"

    clear_reference_where_ice_free(a_ref, ice_free)
    assert a_ref[1] == 0.0
    assert (h + a_ref * dt)[1] == 0.0


def test_cells_that_still_hold_ice_keep_their_reference():
    r"""The rule is targeted: it must not disturb the balancing reference on
    the ice, which is what holds a control run steady."""
    a_ref = np.array([0.5, 40.0, -2.0])
    clear_reference_where_ice_free(
        a_ref, np.array([False, True, False]))
    assert a_ref.tolist() == [0.5, 0.0, -2.0]


def test_reference_stays_zero_when_the_front_readvances():
    r"""A cell that re-enters the ice does NOT get its reference back: the
    frozen reference was only ever defined on the t=0 ice."""
    a_ref = np.array([40.0])
    clear_reference_where_ice_free(a_ref, np.array([True]))
    # The front advances back over it: phi turns negative, so it is no longer
    # reported ice-free, and the rule simply does not touch it again.
    clear_reference_where_ice_free(a_ref, np.array([False]))
    assert a_ref[0] == 0.0


def test_nothing_ice_free_is_a_no_op():
    r"""An interior step with a stationary front leaves the field alone."""
    a_ref = np.array([0.5, 40.0, -2.0])
    before = a_ref.copy()
    clear_reference_where_ice_free(a_ref, np.zeros(3, dtype=bool))
    assert np.array_equal(a_ref, before)


# ── where the forcing acts ──────────────────────────────────────────────────
from icepack2_tools.front import unforced_cells  # noqa: E402


def test_open_ocean_is_not_forced_and_ice_and_land_are():
    r"""Ice on any bed is forced; an ice-free cell on a bed below sea level is
    open ocean and is not; an ice-free cell on land stays forced, so ice can
    still grow there."""
    h = np.array([300.0, 5.0, 0.0, 0.0, 0.0])
    bed = np.array([-800.0, 200.0, -400.0, 150.0, 0.0])
    assert unforced_cells(h, bed).tolist() == [False, False, True, False, False]


def test_cells_a_front_rule_holds_ice_free_are_not_forced():
    r"""Whatever a front rule empties, land or ocean, gets no forcing: SMB
    there would make ice the rule removes and books as calving."""
    h = np.array([300.0, 0.0, 0.0, 0.0])
    bed = np.array([-800.0, 150.0, 150.0, -400.0])
    beyond = np.array([False, True, False, False])
    ls_ice_free = np.array([False, False, False, True])
    assert unforced_cells(h, bed, beyond, None, ls_ice_free).tolist() == \
        [False, True, False, True]


from icepack2_tools.front import applied_forcing  # noqa: E402


def test_an_emptied_shelf_cell_gets_no_reference():
    r"""A floating cell the collapse mask empties becomes open ocean. Its t=0
    reference, about the shelf's t=0 melt, must not regrow it each advance;
    the neighbouring ice keeps its reference, and the stored field is kept
    so the reference applies again if ice returns."""
    dt = 0.1
    h = np.array([400.0, 0.0, 250.0])          # cell 1 was just collapsed
    bed = np.array([-600.0, -700.0, 100.0])
    accum = np.array([0.3, 0.2, 0.4])
    melt = np.array([5.0, 8.0, 0.0])
    a_ref = np.array([4.0, 7.9, -0.1])
    stored = a_ref.copy()

    forced = np.where(unforced_cells(h, bed), 0.0, 1.0)
    smb, mlt, ref = applied_forcing(forced, accum, melt, a_ref)
    src = smb - mlt + ref

    assert (h + a_ref * dt)[1] > 0.0, "guard: the unmasked reference regrows it"
    assert src[1] == 0.0 and ref[1] == 0.0
    assert ref[[0, 2]].tolist() == a_ref[[0, 2]].tolist()
    assert np.array_equal(a_ref, stored)

    # Ice flows back into the cell: the stored reference applies again.
    h[1] = 50.0
    forced = np.where(unforced_cells(h, bed), 0.0, 1.0)
    _, _, ref = applied_forcing(forced, accum, melt, a_ref)
    assert ref[1] == a_ref[1]


def test_no_reference_passes_through_as_none():
    smb, mlt, ref = applied_forcing(np.array([1.0, 0.0]), np.array([0.3, 0.3]),
                                    np.array([1.0, 1.0]))
    assert ref is None
    assert smb.tolist() == [0.3, 0.0] and mlt.tolist() == [1.0, 0.0]
