r"""Front bookkeeping the thickness transport does around the level set.

The level set decides where the front IS (see :mod:`icepack2_tools.levelset`).
These are the rules the DG0 transport applies to the cells it touches once the
front has moved. They live here, as pure functions over the cell arrays, so
they can be exercised without standing up a whole run.
"""

import numpy as np

__all__ = ["retreat_slivers", "clear_reference_where_ice_free", "clamp_thickness",
           "unforced_cells", "applied_forcing", "front_connected", "facet_neighbours",
           "ocean_drag_cells",
           "collapse_cell_counts",
           "collapse_banner", "collapse_csv_fields", "COLLAPSE_MARKER",
           "COLLAPSE_CSV_COLUMNS", "FRONT_OWNER_MARKER"]

# The prefix of the line a run prints naming the mechanism that owns the
# calving front (an ISMIP7_CALVING law, an external law with its parameters,
# or the legacy mask). core_report.py lifts it into the run record.
FRONT_OWNER_MARKER = "Calving front owner:"

# The prefix of every line a run prints about the collapse forcing: the mode
# at startup, then the cell counts. core_report.py lifts the lines carrying it
# into the run record, so keep it on one line with what it introduces.
COLLAPSE_MARKER = "Ice-shelf collapse forcing:"

# The timeseries columns collapse_cell_counts fills, in its return order.
COLLAPSE_CSV_COLUMNS = ("collapse_flagged_cells", "collapse_removed_cells",
                        "collapse_held_cells")


def unforced_cells(h, bed, *ice_free):
    r"""The cells the surface and ocean forcing must not act on.

    ``h`` is the per-cell thickness at the start of the advance, ``bed`` the
    per-cell bed elevation, and each ``ice_free`` argument a boolean mask (or
    ``None``) of cells a front rule holds ice-free. Returns a boolean mask:
    open ocean (no ice on a bed below sea level) and every such held cell.

    Forcing there has no ice to act on, and the model used to count it anyway.
    Measured on a 2 km control, the open ocean beyond a fixed front received
    57 Gt/yr of SMB and 87 to 90 Gt/yr of melt. The positivity limit withheld
    the melt that had no ice to melt and booked it as ``clamp`` (+64 Gt/yr),
    while the melt column and the ISMIP7 ``libmassbffl`` field still counted
    it; where the SMB won, it made ice on open water that the front removed
    on the next advance and booked as calving, 31 of the 48 Gt/yr the run
    called calving. None of it touched the ice.

    Ice-free land (a bed at or above sea level) stays forced unless a front
    rule holds it, so ice can still grow on it as before.
    """
    out = (np.asarray(h) <= 0.0) & (np.asarray(bed) < 0.0)
    for mask in ice_free:
        if mask is not None:
            out = out | np.asarray(mask, dtype=bool)
    return out


def applied_forcing(forced, accum, ocean_melt, a_ref=None):
    r"""The SMB, melt and apparent-MB reference an advance applies.

    ``forced`` is 1 where the forcing acts and 0 on the cells
    :func:`unforced_cells` names; the other arguments are the raw fields.
    Works on cell arrays and on UFL expressions alike. Returns
    ``(smb, melt, ref)``, each masked by ``forced``; ``ref`` is ``None`` when
    ``a_ref`` is.

    The reference is masked with the forcing. A cell emptied mid-run (a
    collapsed shelf, or one melted through to h=0) turns into open ocean, and
    its t=0 reference, about the shelf's t=0 melt, would regrow it every
    advance for the front to remove and book as calving. ``a_ref`` itself is
    left as it is, so it applies again if ice returns to the cell.
    """
    ref = None if a_ref is None else forced * a_ref
    return forced * accum, forced * ocean_melt, ref


def clamp_thickness(h, h_clamp, *ice_free):
    r"""Floor the cell thicknesses to ``h_clamp``, except where there is no ice.

    ``h`` is the per-cell thickness array, mutated in place and returned; each
    ``ice_free`` argument is a boolean mask (or ``None``) naming cells that
    hold no ice. The caller measures the mass the floor added by integrating
    before and after.

    The exemption is the whole point. A cell outside the ice domain floored to
    ``h_clamp`` is handed that much fresh ice out of nothing, and whatever rule
    empties it takes it away again on the very next advance, so the run
    fabricates a steady ``h_clamp / dt`` of both clamp mass and calving flux
    there forever. Every rule that empties a cell therefore has to name it
    here:

    * the level set's ice-free extent, not only the cells calved this step:
      under a free law the cells the front just passed are a handful, so
      flooring the rest of the buffer would fabricate ice across every
      never-glaciated cell;
    * the cells beyond a pinned or fixed front;
    * the FLOATING cells the ISMIP7 collapse mask empties, which is the same
      case: the shelf is gone there, and a floored cell would be re-calved
      into ``licalvf`` on every advance for the rest of the run.
    """
    floor = np.full_like(h, h_clamp)
    for mask in ice_free:
        if mask is not None:
            floor[mask] = 0.0
    np.maximum(h, floor, out=h)
    return h


def clear_reference_where_ice_free(a_ref, ice_free):
    r"""Zero the frozen apparent-MB reference in the cells that hold no ice.

    ``a_ref`` is the per-cell reference array, mutated in place; ``ice_free``
    is the level set's ice-free mask for the CURRENT extent.

    The reference is the t=0 flux divergence, so at a t=0 front cell it is the
    terminus outflow: large and positive. Applied against the live extent it
    closes both directions of the same rule - no balancing reference where
    there is no ice:

    * advance - a frozen SINK outside the extent re-empties the cells a free
      front advances into;
    * retreat - a frozen SOURCE inside the t=0 extent regrows the cells a free
      front has just calved, so the front cannot retreat and its calving tally
      counts the regrown ice again every step.

    Zeroing in place is deliberate: a cell that later re-enters the ice keeps
    ``a_ref = 0``, because the frozen reference was only ever defined on the
    t=0 ice.
    """
    a_ref[ice_free] = 0.0
    return a_ref


def retreat_slivers(h_new, h_old, front_hmin):
    r"""Mask of cells holding a retreat remainder the front should take.

    ``h_new`` and ``h_old`` are cell thicknesses after and before one transport
    advance, over the same owned cells. A cell qualifies when it HELD ice
    (``h_old > front_hmin``) and ended the step below the extent threshold
    (``0 < h_new <= front_hmin``). That remainder is what the sub-cell calving
    shed and the melt left behind, so it is the front's own loss and belongs
    to the calving tally.

    A cell that was ice-free is never touched, however little it holds. That
    sub-threshold inflow is how a free front ADVANCES: it accumulates over
    successive steps until it crosses ``front_hmin`` and joins the extent.
    Removing it every step would pin the front wherever the one-step influx
    is under the threshold.
    """
    return (h_new > 0.0) & (h_new <= front_hmin) & (h_old > front_hmin)


def front_connected(flagged, open_water, neighbours_of, any_rank=bool):
    r"""The ``flagged`` cells that open water can reach through flagged cells.

    The ISMIP7 collapse mask lights up near the grounding line of the Ross
    and Filchner-Ronne shelves long before it reaches their fronts. Emptying
    every flagged floating cell therefore opens holes hundreds of kilometres
    behind the front, and a hole is open ocean to the momentum balance: the
    DG0 facet term puts the full terminus water pressure on its rim and the
    glaciers feeding it lose their buttressing at once. Groups that apply the
    mask that way report 40 % to 100 % more sea level by 2300 from this
    alone, and the ice in such a hole has nowhere to go: it cannot calve into
    a shelf that is still there (discussion #30, September 2026). The other
    end-member the groups settled on is this one: a shelf collapses from its
    front, so a flagged cell goes only once the ocean has reached it.

    ``flagged`` and ``open_water`` are per-cell boolean arrays (the owned
    cells of this rank). ``neighbours_of(mask)`` returns the cells that share
    a facet with a cell of ``mask``, across rank boundaries; ``any_rank`` is
    the collective OR that ends the sweep on every rank at once. A flagged
    cell that is itself open water (one emptied on an earlier advance) is
    reached by definition, which is what keeps it named to the thickness
    floor. One sweep grows the set by one layer of cells, so the count is the
    width of the collapsing region in cells, and zero once it has collapsed.
    """
    reached = flagged & open_water
    water = open_water | reached
    while True:
        new = flagged & ~reached & neighbours_of(water)
        if not any_rank(new.any()):
            return reached
        reached |= new
        water |= new


_COLLAPSE_RULES = {
    "none": "no collapse mask is read and no cell is removed",
    "mask": "floating cells flagged by the mask are removed and booked as calving",
    "mask_front": "flagged floating cells are removed once open water reaches them, and booked as calving",
}


def collapse_banner(mode):
    r"""The line a run prints once at startup, under every ``ISMIP7_FRACTURE``
    mode. ``none`` prints too: a log with no such line is then a run that
    predates it, and the mode of any later run can be read off its log."""
    return f"{COLLAPSE_MARKER} ISMIP7_FRACTURE={mode} ({_COLLAPSE_RULES[mode]})"


def collapse_cell_counts(flagged, removed, count=np.count_nonzero):
    r"""``(flagged, removed, held)`` cell counts for one transport advance.

    ``flagged`` is the floating cells the year's mask names. ``removed`` is
    the ones the mode empties: all of them under ``mask``, the ones
    :func:`front_connected` returns under ``mask_front``. An emptied cell
    stays in both for as long as it is flagged and afloat (that is what names
    it to the thickness floor), so ``removed`` accumulates over a run. Held
    cells are the rest, flagged floating ice the advance leaves standing: zero
    under ``mask`` by construction, and under ``mask_front`` the size of the
    holes ``mask`` would have opened behind the front. Both arguments are
    ``None`` when no mask is read (``ISMIP7_FRACTURE=none``) and every count is
    zero.

    ``count`` reduces a boolean cell array to a number of cells. A run passes
    a global count (``icepack2_tools.mpi_stats.global_count``), which is
    collective: every rank calls this, or none does. The default counts the
    local array, for the tests.
    """
    if flagged is None:
        return 0, 0, 0
    return int(count(flagged)), int(count(removed)), int(count(flagged & ~removed))


def collapse_csv_fields(header, counts):
    r"""What a timeseries row appends for ``counts``: one field per entry of
    ``COLLAPSE_CSV_COLUMNS`` when ``header``, the file's own first line,
    carries them, and nothing otherwise. A series begun before the columns
    existed keeps its header on resume, so its rows keep their width and the
    counts reach the log alone."""
    if COLLAPSE_CSV_COLUMNS[-1] not in header:
        return ""
    return "," + ",".join(str(int(c)) for c in counts)


def ocean_drag_cells(ice, neighbours_of, extent0):
    r"""The cells the floor-cell ocean drag may act on: open water that holds
    no ice now, shares no facet with a cell that does, and lies outside the
    t=0 ice extent.

    ``ice`` is the per-cell boolean "holds ice now" (thickness at least the
    front threshold), ``extent0`` the same test on the t=0 thickness, and
    ``neighbours_of(mask)`` the cells sharing a facet with a cell of ``mask``
    (:func:`facet_neighbours`). Every other cell, floating or grounded, gets
    no drag.

    The drag exists to give the ice-free buffer some velocity coercivity. Left
    on everywhere below ``h_ocean`` it also acted on thin floating ice and on
    the water row the front's vertices share, and held the front back: at
    1875 in the 1 km CESM2-WACCM historical the floating front moved at 0.12
    of the observed speed (median; 0.09 summed), against 0.83 one cell further
    in, and only 19 Gt/yr crossed into the ice-free cells as calving. Keeping
    the t=0 extent drag-free as well means a front that retreats inside it
    meets no drag either, before a restart and after.
    """
    return ~ice & ~neighbours_of(ice) & ~extent0


def facet_neighbours(Q_dg):
    r"""``neighbours_of`` for :func:`front_connected` on the DG0 space
    ``Q_dg``: the cells sharing an interior facet with a cell of the mask.

    Assembled as a facet integral of the mask's indicator against the DG0
    test function on the far side, the same transfer the ISMIP7 output uses
    to book the grounding-line flux into the first floating cell. Writing the
    indicator through ``dat.data`` marks its halo stale, so the assembly
    exchanges it and a neighbour across a partition boundary counts.
    """
    import firedrake as fd

    indicator = fd.Function(Q_dg)
    phi = fd.TestFunction(Q_dg)
    form = (indicator("+") * phi("-") + indicator("-") * phi("+")) * fd.dS
    touched = fd.Cofunction(Q_dg.dual())

    def neighbours_of(mask):
        indicator.dat.data[:] = mask
        fd.assemble(form, tensor=touched)
        # a sum of facet lengths over the neighbours in the mask: exactly
        # zero with none, at least one facet length with any
        return touched.dat.data_ro > 0.0

    return neighbours_of
