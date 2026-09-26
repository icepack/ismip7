r"""ISMIP7 output, part one: what the forward records every year.

The data request (``isschecker/data/ISMIP7_variable_request.csv`` in
``ismip/ISM_SimulationChecker``, bundled here as
``icepack2_tools/ismip7_variable_request.csv``) asks for yearly 2D fields on
the 8 km AIS grid and yearly scalars. Regridding is a serial post-processing step
(``antarctica/scripts/write_ismip7_output.py``); this module is the part that
runs inside the parallel forward, gated by ``ISMIP7_OUTPUT=1``:

* state variables (``ST``) are snapshots at the end of each year, stamped
  1 January of the following year by the writer;
* flux variables (``FL``) are the year's means, accumulated every transport
  advance from the sources the transport APPLIED, stamped 1 July. Where the
  positivity limiter withholds part of a net sink that would draw a cell
  below ``h_clamp``, the withheld part comes off the SMB, melt and reference
  sinks in proportion to their sizes (``book_advance``), so a cell the
  limiter holds at the floor books only the melt it had ice for. The run's
  ``clamp`` budget column still carries the withheld total. Booking the
  requested sources instead put 1.43 % of the floating values of a
  full-length 32 km ssp585 below the request's -0.008 kg m-2 s-1, an
  isschecker 0.5.1 error, because by 2200 the limiter held back 99 % of the
  requested melt;
* the scalars are the integrals of the same fields over true area, written to
  a CSV as the run goes so an early stop loses nothing. A cell counts its
  map-plane area times af2 = (1/k)^2 of EPSG:3031 at its centroid, the factor
  ``ismip7-scalars`` weights every pixel by (``regrid.area_factor``). The
  cell means and the run's own mass budget stay map-plane, so at 32 km the
  scalars' mass sits about 2.6 % above the timeseries' ``mass_gt``.

Everything is kept on the model's own mesh in Firedrake checkpoints, ONE PER
YEAR: ``<results>/<experiment>_<lc>_ismip7_annual_<year>.h5``, each holding
that year's fields in the model's units (m, m/yr, MPa). The writer globs
them, sorts by the year in the name, and converts to the request's SI units
under the fill policies. ``base`` is not among them: it is the one submitted
variable the writer rebuilds on the 8 km grid (see ``GRID_DERIVED``).

One file per year rather than one growing file is what makes the series
durable. This is the only copy of what gets submitted, and six chained links
append ~285 years into it; a single file rewritten in place has no atomic
swap (HDF5 does not journal), so a job killed mid-write by a node failure or
a queue kill would take every earlier year with it, while re-copying a
multi-GB file once per year to get an atomic rename costs O(N^2) I/O. Each
year is instead written to ``<name>.tmp`` and renamed into place, the
``_save_state`` pattern, so a kill can only ever lose the year in flight.
The cost is one copy of the mesh per file, which is the price of never
being able to lose the years already banked.

The year in progress rides in the run's OWN checkpoint, not here: see
``state_fields`` below.

Conventions (from the request and discussions #16, #19, #22):

* ``acabf`` is the surface mass balance the transport applied: the forcing
  SMB (RACMO climatology plus the re-referenced anomaly) less its share of any
  sink the positivity limiter withheld. The apparent-mass-balance reference ``a_ref`` is NOT
  part of it: it cancels the discrete flux
  divergence spike by spike (up to ~1000 m/yr at the Pine Island grounding
  zone) and reported as SMB it would sit two orders of magnitude outside the
  request's range. It is recorded separately as ``acabf_correction`` (m/yr
  ice, not a request variable) so the grid budget can be closed by anyone
  who needs it, and the README states the convention.
* ``libmassbffl`` is the ocean melt the transport applied on floating cells
  (negative = loss), less than the parameterization's melt wherever the
  limiter holds a cell at the floor;
  ``libmassbfgr`` is zero (no grounded basal melt in the model).
* ``lifmassbf`` is the melt of ice that flowed into a marine cell holding no
  ice at either end of the year (issue #109, option 3 of the 25 September
  2026 meeting): grounded ice that goes afloat into an emptied shelf cell and
  melts there. The request's ``no_floating_ice`` fill would blank that melt
  from ``libmassbffl``; ``lifmassbf`` is never filled. ``year_end`` moves the
  inflow's share of such a cell's melt out of ``libmassbffl``
  (``split_front_melt``), and the share the frozen apparent-MB reference and
  the SMB supplied stays behind, so the two fields sum to the melt the
  transport applied. Each year file is stamped ``FRONT_MELT_ATTR``; the
  writer refuses a series that mixes this booking with the earlier one.
* ``licalvf`` is the ice removed at the front, negative = loss, booked in the
  cell it was removed from (whole-cell removal, sub-cell shed, retreat
  slivers), the same tallies as the ``calv`` budget column.
* ``ligroundf`` is the flux across the grounding line, booked as a specific
  mass flux into the first FLOATING cell (discussion #22), signed with the
  grounded sheet as the reference: positive for grounded ice going afloat,
  negative where floating ice flows onto grounded ice (a pinning point or an
  ice rumple), so the grid sum of ``ligroundf * area`` is the net
  grounding-line discharge. The group settled the reference on 22 September
  2026.
* ``dlithkdt`` is the change in thickness over the year divided by the year.
* the three area fractions are cell indicators here (0 or 1 per DG0 cell);
  the conservative regridding to 8 km turns them into fractions.
"""
import os

import numpy as np
from firedrake import Function, TestFunction, assemble, dS, dx
import firedrake as fd

from .mpi_stats import global_range
from .regrid import area_factor

# icepack's year (365.25 days): the model's own time unit, so every
# model-to-SI conversion the submission carries uses it. The time axis
# in the files is the standard calendar regardless.
SECONDS_PER_YEAR = 31557600.0
RHO_I = 917.0

#: request name -> (type, how this model produces it)
VARIABLES_2D = {
    "lithk": "ST", "orog": "ST", "topg": "ST", "base": "ST",
    "sftgif": "ST", "sftgrf": "ST", "sftflf": "ST",
    "xvelmean": "ST", "yvelmean": "ST", "xvelsurf": "ST", "yvelsurf": "ST",
    "xvelbase": "ST", "yvelbase": "ST", "strbasemag": "ST",
    "acabf": "FL", "libmassbfgr": "FL", "libmassbffl": "FL", "dlithkdt": "FL",
    "licalvf": "FL", "ligroundf": "FL", "lifmassbf": "FL",
}
SCALARS = ("lim", "limnsw", "iareagr", "iareafl", "tendacabf", "tendlibmassbfgr",
           "tendlibmassbffl", "tendlicalvf", "tendlifmassbf", "tendligroundf")

#: submitted variables the writer REBUILDS on the 8 km grid rather than
#: regridding from the model mesh, so the forward never banks them: the
#: checker requires orog == base + lithk pixel by pixel, and the elevations
#: are covered-part means while lithk is a whole-pixel mean, so only
#: base := orog - lithk taken on the grid satisfies it.
GRID_DERIVED = ("base",)
#: what year_end writes, and what the writer reads back per year
VARIABLES_BANKED = tuple(v for v in VARIABLES_2D if v not in GRID_DERIVED)

#: the attribute on every year file naming how its melt is booked, and the
#: value this module writes; a year file without it comes from a forward that
#: wrote ``lifmassbf`` as zero and kept all the melt in ``libmassbffl``
FRONT_MELT_ATTR = "front_melt"
FRONT_MELT = "inflow_share_of_empty_marine_cells"


def split_front_melt(libmassbffl, acabf, corr, dh, empty):
    r"""A year's booked melt split into ``(libmassbffl, lifmassbf)`` (issue #109).

    Every array is per cell: ``libmassbffl``, ``acabf`` and ``corr`` (the
    apparent-MB reference) are the year's sums of what the transport applied,
    ``dh`` is the change in thickness over the same span, all in metres of
    ice, and ``empty`` marks the cells that take part, marine cells holding
    no ice at either end of the year.

    The thickness budget closes on the booked sources, so the inflow a cell
    kept, net of what flowed on and of what the front removed as
    ``licalvf``, is the remainder ``I = dh - acabf - corr - libmassbffl``.
    The limiter nets the inflow, the SMB and the reference into one source,
    which leaves nothing to tell their ice apart, so the melt ``M`` is shared
    among what they supplied in proportion, the split ``leftout_melt.py``
    (issue #105) measured with, and ``lifmassbf`` is the inflow's share:

        lifmassbf = -M I+ / (I+ + acabf+ + corr+),    M = max(-libmassbffl, 0)

    and zero where nothing was supplied or outside ``empty``. Refreezing, a
    positive ``libmassbffl``, never moves. ``lifmassbf`` is the negated
    product of non-negative factors, so it never exceeds zero, the bound the
    request's range puts on it, and the two fields sum to ``libmassbffl``.
    The inflow that the front removes in the advance it arrives in never
    reaches the melt, and ``I`` leaves it out.

    The closure is exact under DG0 geometry, where every other change of a
    cell's thickness is booked and ``ISMIP7_H_CLAMP`` is 0; what escapes is
    round-off and the melt of a cell that grounds inside a subcycled step.
    Under CG1 geometry the forcing masks the melt per node while the booking
    counts it per cell, so melt goes unbooked near the grounding line and
    ``I`` comes out short there.
    """
    melt = np.maximum(-libmassbffl, 0.0)
    kept = np.maximum(dh - acabf - corr - libmassbffl, 0.0)
    supply = kept + np.maximum(acabf, 0.0) + np.maximum(corr, 0.0)
    share = np.divide(kept, supply, out=np.zeros_like(supply), where=supply > 0.0)
    lifmassbf = np.where(empty, -(melt * share), 0.0)
    return libmassbffl - lifmassbf, lifmassbf


class AnnualOutput:
    r"""Accumulates a year of a forward run and writes it at the year end.

    ``begin_step`` / ``commit_step`` bracket one time step so a rewound
    (subcycled) attempt does not double-count: ``_advance`` books into the
    step tallies, and only a completed step is added to the year.

    A chained run resumes into the same output: ``out_path`` names the stem
    ``<...>_ismip7_annual.h5`` and the years already on disk are the files
    ``<...>_ismip7_annual_<year>.h5`` beside it, so the written years are
    found by globbing rather than carried in an attribute. The partly
    accumulated year survives the link boundary too: ``state_fields`` /
    ``state_attrs`` hand it to the run's own checkpoint and ``resume`` takes
    it back, so a job that stops at 2021.4 goes on accumulating 2021 rather
    than losing four months of flux or relabelling them. That state is
    stamped with the series it belongs to, so a restart from another
    experiment's checkpoint (a projection starting from the historical
    endpoint carries the historical run's state) counts as a cold start here
    rather than as this series' resume.

    On this series' own resume, years at or after the resumed one are stale,
    left behind by an unclean kill that ran past the last state checkpoint,
    so they are discarded and re-simulated. Any other run refuses to touch a
    banked series: it names the files and stops, and removing them is the
    operator's own call. A year is never renamed, and a resume that would
    leave a HOLE in the series (further than one past the last kept year) is
    an error.

    A run that starts part-way through a year with no record of its earlier
    months does not bank that year at all: accumulation begins at the next
    1 January. The state carries the model's own year plus a skip flag, so a
    checkpoint taken inside that window resumes the skip rather than reading
    the year back as an accumulation in progress.
    """

    #: the per-cell year sums carried across a chained resume
    ACCUMULATORS = ("acabf", "acabf_correction", "libmassbffl", "licalvf",
                    "ligroundf")
    #: a cell holds ice when it is thicker than this (m): the forward's
    #: ``ice_cells`` for ``year_end``, and the same test on the thickness
    #: the year began with, which ``split_front_melt`` needs
    ICE_THICKNESS = 1.0
    #: dataset names of the in-progress year inside the run's checkpoint
    STATE_PREFIX = "ismip7_acc_"
    STATE_THICKNESS = "ismip7_year_start_thickness"
    STATE_YEAR = "ismip7_year"
    STATE_YEAR_TIME = "ismip7_year_time"
    #: which series the carried state belongs to: the annual stem's basename,
    #: which is <experiment>_<lc>_ismip7_annual.h5 and so names the run
    STATE_SERIES = "ismip7_series"
    #: set while the run is inside a partial year it will not bank, so a
    #: checkpoint taken in that window resumes the skip instead of reading
    #: the year back as an accumulation in progress
    STATE_SKIPPING = "ismip7_skipping"

    def __init__(self, mesh, Q_dg, out_path, scalars_path, first_year, rho_ratio,
                 comm=None, log=None, resume=None):
        self.mesh, self.Q_dg = mesh, Q_dg
        self.out_path, self.scalars_path = out_path, scalars_path
        self.year = int(np.floor(float(first_year) + 1e-6))   # the year the model time lies in
        self._skip_first_year_end = False
        self.rho_ratio = float(rho_ratio)
        self.comm = comm or mesh.comm
        self.log = log or (lambda s: None)
        # owned cells only: dat.data_ro is the owned slice, dof_count counts
        # the halo too (5650 vs 3772 on one of two ranks of the 32 km mesh)
        self.cell_area = assemble(TestFunction(Q_dg) * dx).dat.data_ro.copy()
        # The scalars integrate over true area (issue #97): each cell's
        # map-plane area times af2 at its centroid, the same owned cells in
        # the same order. cell_area stays map-plane, since book_advance
        # divides by it for the cell means.
        X = fd.SpatialCoordinate(mesh)
        x, y = (Function(Q_dg).interpolate(X[i]).dat.data_ro for i in (0, 1))
        self.area_factor = area_factor(x, y)
        self.true_area = self.cell_area * self.area_factor
        lo, hi = global_range(self.area_factor, comm=self.comm)
        self.log(f"  ISMIP7 output: scalars over true area, af2 from {lo:.4f} to "
                 f"{hi:.4f} on this mesh")
        n = len(self.cell_area)
        self.year_acc = {k: np.zeros(n) for k in self.ACCUMULATORS}
        self.step_acc = {k: np.zeros(n) for k in self.year_acc}
        self.year_time = 0.0
        self.h_year_start = None
        # A restart checkpoint carrying ISMIP7 attributes is only THIS
        # series' resume when it names this series. A projection restarting
        # from the historical endpoint carries the historical run's
        # accumulation state, which belongs to a different submitted series,
        # so for this output it is a cold start.
        if resume is not None and resume.get("series") != os.path.basename(out_path):
            self.log(f"  ISMIP7 output: the restart checkpoint carries the "
                     f"{resume.get('series') or 'unnamed'} series' accumulation "
                     f"state, not {os.path.basename(out_path)}; starting a new "
                     f"series here")
            resume = None
        if resume is not None and int(resume["year"]) != self.year:
            raise ValueError(
                f"the restart checkpoint was accumulating ISMIP7 year "
                f"{int(resume['year'])} but its timeline resumes at "
                f"t={float(first_year)!r}, which falls in year {self.year}; "
                f"the checkpoint's ISMIP7 state does not belong to it."
            )
        # A run part-way through a year it has no record of cannot bank that
        # year: a fraction of a year reported as the year's mean is a wrong
        # submitted value. Accumulation begins at the next 1 January instead
        # and the partial year is dropped. Both entries to that state land
        # here: a cold start at a non-integer time, and a resume from a
        # checkpoint taken while the skip was already in progress, which
        # carries the model's own year plus the skip flag so the window is
        # re-entered rather than read back as an accumulation.
        skipping = resume is not None and resume.get("skipping")
        if skipping or (resume is None and abs(float(first_year) - self.year) > 1e-6):
            self.year += 1
            self._skip_first_year_end = True
            resume = None
        if resume is not None:
            for k in self.ACCUMULATORS:
                self.year_acc[k][:] = resume["acc"][k]
            self.h_year_start = np.asarray(resume["h_year_start"]).copy()
            self.year_time = float(resume["year_time"])
        self._phi = TestFunction(Q_dg)
        self._n = fd.FacetNormal(mesh)
        self._gl_cof = fd.Cofunction(Q_dg.dual())
        # A chained job re-enters run_simulation and rebuilds this object, so
        # the years the earlier links banked are read back off disk. On a
        # RESUME, years at or after the one the checkpoint was accumulating
        # belong to a trajectory the restart abandons: an unclean kill banks
        # years past the last state checkpoint (the cadences differ, 5 yr vs
        # 1 yr by default), and keeping them would splice two runs into one
        # series, so they are discarded and re-simulated. Anything else is a
        # cold start for this output: nothing says the years on disk are
        # wrong, and this is the only copy of what gets submitted, so it
        # refuses rather than deleting a banked series. Removing them is the
        # operator's own deliberate act, which leaves an audit trail that an
        # in-process discard would not.
        self._written_years = self.years_on_disk(out_path)
        stale = [y for y in self._written_years if y >= self.year]
        if stale and resume is None:
            stem, ext = os.path.splitext(out_path)
            raise ValueError(
                f"{os.path.basename(out_path)} already holds ISMIP7 years "
                f"{stale[0]}-{stale[-1]}, but this run starts at {self.year} "
                f"without resuming THIS series' accumulation state, so it "
                f"would rewrite a banked submission series. Resume from a "
                f"checkpoint of this experiment, or move the series aside: "
                f"{stem}_{{{stale[0]}..{stale[-1]}}}{ext} and "
                f"{os.path.basename(scalars_path)}."
            )
        if stale:
            self._written_years = [y for y in self._written_years if y < self.year]
            if self.comm.rank == 0:
                for y in stale:
                    os.remove(self.year_path(out_path, y))
                self._trim_scalars(scalars_path, self.year)
            self.comm.barrier()
            self.log(f"  ISMIP7 output: discarded {len(stale)} year(s) "
                     f"({stale[0]}-{stale[-1]}) that the restart abandons; "
                     f"they will be re-simulated")
        if self._written_years:
            last = self._written_years[-1]
            if self.year > last + 1:
                raise ValueError(
                    f"{os.path.basename(out_path)} holds years "
                    f"{self._written_years[0]}-{last}, so the next year to "
                    f"accumulate is {last + 1}, but this run resumes inside "
                    f"year {self.year}: the submitted series would have a "
                    f"hole in it. Resume from a checkpoint inside {last + 1}."
                )
        if self._skip_first_year_end:
            self.log(f"  ISMIP7 output: inside {self.year - 1} at "
                     f"t={float(first_year)!r} with no record of its earlier "
                     f"months, so that year is NOT banked; accumulation "
                     f"begins at {self.year}-01-01.")
        if self._written_years:
            self.log(f"  ISMIP7 output: continuing {os.path.basename(out_path)} "
                     f"({len(self._written_years)} years through {last}; "
                     f"{self.year_time:.2f} yr of {self.year} carried over)")
        if self.comm.rank == 0:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            new = not os.path.exists(scalars_path)
            self._csv = open(scalars_path, "a")
            if new:
                self._csv.write("year," + ",".join(SCALARS) + "\n"); self._csv.flush()
        else:
            self._csv = None

    # ---- per-step bookkeeping -------------------------------------------
    def begin_step(self):
        for a in self.step_acc.values():
            a[:] = 0.0
        self.step_time = 0.0

    def book_advance(self, dt, accum, ocean_melt, a_ref, h_dg, u, grounded_cells,
                     withheld=None):
        r"""Called by ``_advance`` after the transport solve, BEFORE removal:
        books the sources this advance APPLIED (the SMB, the melt and the
        apparent-mass-balance reference) and the grounding-line flux with the
        velocity the transport used.

        ``withheld`` is the part of each cell's net sink that the positivity
        limiter held back in this advance (m/yr, never negative), the limited
        source minus the requested one. It comes off the three sinks, the
        negative SMB, the melt and the negative reference, in proportion to
        their sizes, so the booked sources sum to the source the transport
        applied. ``None`` books the requested sources."""
        smb = assemble(accum * self._phi * dx).dat.data_ro / self.cell_area        # m/yr, cell mean
        melt = assemble(ocean_melt * self._phi * dx).dat.data_ro / self.cell_area
        corr = (assemble(a_ref * self._phi * dx).dat.data_ro / self.cell_area
                if a_ref is not None else np.zeros_like(smb))
        if withheld is not None:
            sinks = np.maximum(-smb, 0.0) + np.maximum(melt, 0.0) + np.maximum(-corr, 0.0)
            keep = 1.0 - np.divide(np.minimum(withheld, sinks), sinks,
                                   out=np.zeros_like(sinks), where=sinks > 0.0)
            smb = np.where(smb < 0.0, smb * keep, smb)
            melt = np.where(melt > 0.0, melt * keep, melt)
            corr = np.where(corr < 0.0, corr * keep, corr)
        self.step_acc["acabf"] += smb * dt
        if a_ref is not None:
            self.step_acc["acabf_correction"] += corr * dt
        self.step_acc["libmassbffl"] += -melt * dt * (~grounded_cells)
        # grounding-line flux into the first floating cell: the upwind facet
        # flux across every facet whose two cells differ in grounding, the
        # same flux the DG0 transport moved, booked to the floating side (its
        # test function) and signed positive from grounded to floating. Ice
        # flowing from a shelf onto a pinning point books negative, so the
        # sum over the mesh is the net discharge and an ice rumple's
        # throughput cancels instead of counting as loss.
        g = Function(self.Q_dg); g.dat.data[:] = grounded_cells.astype(float)
        un = fd.dot(u, self._n); un_plus = (un + abs(un)) / 2
        flux = un_plus("+") * h_dg("+") - un_plus("-") * h_dg("-")   # upwind, "+" to "-"
        phi = self._phi
        form = (flux * g("+") * (1 - g("-")) * phi("-")
                - flux * g("-") * (1 - g("+")) * phi("+")) * dS
        assemble(form, tensor=self._gl_cof)
        self.step_acc["ligroundf"] += self._gl_cof.dat.data_ro / self.cell_area * dt   # m/yr equivalent
        self.step_time += dt

    def book_removal(self, cells, thickness_removed):
        r"""Ice removed at the front this advance (m of thickness per cell)."""
        self.step_acc["licalvf"][cells] -= thickness_removed

    def commit_step(self):
        for k in self.year_acc:
            self.year_acc[k] += self.step_acc[k]
        self.year_time += self.step_time

    @staticmethod
    def _trim_scalars(scalars_path, first_stale_year):
        r"""Drop the rows of the discarded years, so the CSV and the year
        files describe the same trajectory."""
        import csv as _csv
        if not os.path.exists(scalars_path):
            return
        with open(scalars_path) as f:
            reader = _csv.DictReader(f)
            header = reader.fieldnames
            kept = [r for r in reader if int(r["year"]) < first_stale_year]
        tmp = scalars_path + ".tmp"
        with open(tmp, "w") as f:
            writer = _csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(kept)
        os.replace(tmp, scalars_path)

    # ---- where a year lives ------------------------------------------------
    @staticmethod
    def year_path(out_path, year):
        r"""The checkpoint holding one year, derived from the stem."""
        stem, ext = os.path.splitext(out_path)
        return f"{stem}_{int(year)}{ext}"

    @staticmethod
    def years_on_disk(out_path):
        r"""The years already written beside ``out_path``, in order. Partial
        ``.tmp`` files are not matched, so a year killed mid-write is simply
        absent and gets rewritten."""
        import glob
        import re
        stem, ext = os.path.splitext(out_path)
        pattern = re.compile(re.escape(stem) + r"_(\d+)" + re.escape(ext) + r"$")
        years = [int(m.group(1)) for f in glob.glob(f"{stem}_*{ext}")
                 for m in [pattern.match(f)] if m]
        return sorted(years)

    # ---- carried across a chained resume ----------------------------------
    def state_fields(self):
        r"""The year in progress as DG0 Functions, for the run's checkpoint.

        Firedrake redistributes these on load, so a chained link may run on a
        different rank count than the one that wrote them."""
        fields = {}
        for k in self.ACCUMULATORS:
            f = Function(self.Q_dg, name=self.STATE_PREFIX + k)
            f.dat.data[:] = self.year_acc[k]
            fields[f.name()] = f
        h0 = Function(self.Q_dg, name=self.STATE_THICKNESS)
        if self.h_year_start is not None:
            h0.dat.data[:] = self.h_year_start
        fields[h0.name()] = h0
        return fields

    def state_attrs(self):
        r"""The year the model is IN, how much of it is in the sums, whether
        that year is being skipped, and which series it belongs to.

        The year recorded is the model's own, not the first year that will be
        banked: while a partial year is being skipped those differ by one, and
        a checkpoint stamped with the banked year would not match the timeline
        it was written alongside, so the resume would refuse its own state."""
        return {self.STATE_YEAR: int(self.year - 1 if self._skip_first_year_end
                                     else self.year),
                self.STATE_YEAR_TIME: float(self.year_time),
                self.STATE_SKIPPING: int(self._skip_first_year_end),
                self.STATE_SERIES: os.path.basename(self.out_path)}

    @classmethod
    def read_state(cls, chk, mesh):
        r"""The counterpart of ``state_fields``/``state_attrs``: the resume
        dict to hand back to the constructor, or None if the checkpoint was
        written by a run without ISMIP7 output."""
        if not chk.has_attr("/", cls.STATE_YEAR):
            return None
        return {
            "year": int(chk.get_attr("/", cls.STATE_YEAR)),
            "year_time": float(chk.get_attr("/", cls.STATE_YEAR_TIME)),
            "skipping": bool(int(chk.get_attr("/", cls.STATE_SKIPPING))
                             if chk.has_attr("/", cls.STATE_SKIPPING) else 0),
            "series": (str(chk.get_attr("/", cls.STATE_SERIES))
                       if chk.has_attr("/", cls.STATE_SERIES) else ""),
            "acc": {k: chk.load_function(mesh, name=cls.STATE_PREFIX + k).dat.data_ro.copy()
                    for k in cls.ACCUMULATORS},
            "h_year_start": chk.load_function(mesh, name=cls.STATE_THICKNESS).dat.data_ro.copy(),
        }

    # ---- year end ---------------------------------------------------------
    def start_year(self, h_dg):
        self.h_year_start = h_dg.dat.data_ro.copy()
        for a in self.year_acc.values():
            a[:] = 0.0
        self.year_time = 0.0

    def year_end(self, h_dg, s, b, u, tau, grounded_cells, ice_cells):
        r"""Write this year's state snapshot and flux means, then start the next."""
        if self._skip_first_year_end:
            self._skip_first_year_end = False
            self.start_year(h_dg)
            return
        yr = self.year
        T = self.year_time if self.year_time > 0 else 1.0
        fields = {}
        Q = self.Q_dg
        def dg(arr):
            f = Function(Q); f.dat.data[:] = arr; return f
        # thickness only where the mask says ice (h > 1 m): the checker
        # requires lithk == 0 wherever sftgif == 0, and sub-metre inflow in
        # buffer cells is the front's bookkeeping, not ice
        fields["lithk"] = dg(h_dg.dat.data_ro * ice_cells)
        fields["orog"] = Function(Q).interpolate(s)
        fields["topg"] = Function(Q).interpolate(b)
        fields["sftgif"] = dg(ice_cells.astype(float))
        fields["sftgrf"] = dg((ice_cells & grounded_cells).astype(float))
        fields["sftflf"] = dg((ice_cells & ~grounded_cells).astype(float))
        fields["strbasemag"] = Function(Q).interpolate(fd.sqrt(fd.inner(tau, tau)))   # MPa
        ux = Function(Q).interpolate(u[0]); uy = Function(Q).interpolate(u[1])
        for name in ("xvelmean", "xvelsurf", "xvelbase"):
            fields[name] = ux
        for name in ("yvelmean", "yvelsurf", "yvelbase"):
            fields[name] = uy
        for name in ("acabf", "acabf_correction", "licalvf", "ligroundf"):
            fields[name] = dg(self.year_acc[name] / T)                            # m/yr ice
        fields["libmassbfgr"] = dg(np.zeros_like(self.cell_area))
        h0 = self.h_year_start if self.h_year_start is not None else h_dg.dat.data_ro
        dh = h_dg.dat.data_ro - h0
        fields["dlithkdt"] = dg(dh / T)                                              # m/yr
        # front melt (issue #109): the inflow's share of the melt in marine
        # cells holding no ice at either end of the year
        empty = ((h0 <= self.ICE_THICKNESS) & ~ice_cells
                 & (fields["topg"].dat.data_ro < 0.0))
        melt, front = split_front_melt(
            self.year_acc["libmassbffl"], self.year_acc["acabf"],
            self.year_acc["acabf_correction"], dh, empty)
        fields["libmassbffl"] = dg(melt / T)
        fields["lifmassbf"] = dg(front / T)
        final_path = self.year_path(self.out_path, yr)
        tmp = final_path + ".tmp"
        with fd.CheckpointFile(tmp, "w") as chk:
            chk.save_mesh(self.mesh)
            for name, f in fields.items():
                chk.save_function(f, name=name)
            chk.set_attr("/", FRONT_MELT_ATTR, FRONT_MELT)                       # every rank
        self.comm.barrier()
        if self.comm.rank == 0:
            os.replace(tmp, final_path)
        self.comm.barrier()
        self._written_years.append(yr)
        # scalars, from the same fields over true area (kg, m2, kg/s)
        area = self.true_area
        def integ(arr):
            return self.comm.allreduce(float((arr * area).sum()))
        # limnsw is the mass of the ice ABOVE FLOTATION: the request defines it
        # as that volume times the ice density, so the integrand is the
        # thickness above flotation, h - h_f with h_f = max(-b, 0) / rho_ratio,
        # not the height above flotation s - s_float (which is rho_ratio times
        # smaller on marine beds and wrong outright where the bed is dry).
        topg = fields["topg"].dat.data_ro
        haf_thickness = np.maximum(
            h_dg.dat.data_ro - np.maximum(-topg, 0.0) / self.rho_ratio, 0.0)
        row = {
            "lim": integ(h_dg.dat.data_ro) * RHO_I,
            "limnsw": integ(haf_thickness * grounded_cells) * RHO_I,
            "iareagr": integ((ice_cells & grounded_cells).astype(float)),
            "iareafl": integ((ice_cells & ~grounded_cells).astype(float)),
            "tendacabf": integ(fields["acabf"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
            "tendlibmassbfgr": 0.0,
            "tendlibmassbffl": integ(fields["libmassbffl"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
            "tendlicalvf": integ(fields["licalvf"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
            "tendlifmassbf": integ(fields["lifmassbf"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
            "tendligroundf": integ(fields["ligroundf"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
        }
        # the melt that stays in libmassbffl on those cells: the share the
        # reference and the SMB supplied, which the fill drops wherever the
        # pixel holds no floating ice
        kept_gt = integ(fields["libmassbffl"].dat.data_ro * empty) * RHO_I / 1e12
        if self._csv is not None:
            self._csv.write(f"{yr}," + ",".join(f"{row[k]:.6e}" for k in SCALARS) + "\n"); self._csv.flush()
        self.log(f"  ISMIP7 output: year {yr} written ({len(fields)} fields; over true "
                 f"area, GL flux {row['tendligroundf'] * SECONDS_PER_YEAR / 1e12:+.0f} Gt/yr, "
                 f"calving {row['tendlicalvf'] * SECONDS_PER_YEAR / 1e12:+.0f} Gt/yr, "
                 f"front melt {row['tendlifmassbf'] * SECONDS_PER_YEAR / 1e12:+.0f} Gt/yr with "
                 f"{kept_gt:+.0f} Gt/yr more melt left in libmassbffl on the same cells)")
        self.year = yr + 1
        self.start_year(h_dg)

    def close(self):
        if self._csv is not None:
            self._csv.close()
