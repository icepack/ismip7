#!/usr/bin/env python3
r"""The organisers' scalar processing set against the model's own scalars
(issue #13).

    python antarctica/scripts/compare_scalars.py \
        --submission SUB/AIS/RICE/icepack2/CORE/C007 \
        --tool OUT/tool/nc/AIS/RICE/icepack2/CORE/C007 \
        --datapath ISMIP7/Output-Processing/Data/AIS \
        --params SUB/AIS/RICE/icepack2/params.nc \
        --refyear 2016 [--native-csv RUN_ismip7_scalars.csv] [--native-af2] \
        [--overlap RUN_ismip7_annual.h5.overlap.npz] \
        --out-csv OUT/compare/comparison.csv --out-md OUT/compare/comparison.md

``ismip7-scalars`` (ismip/ismip7-scalar-processing) recomputes the ten scalars
of the variable request from a submission's gridded files, and derives three
sea-level contributions from them. The model integrates the ten itself on its
own mesh (``icepack2_tools/ismip7_output.py``). This sets the two side by side
year by year and splits every difference into named parts:

    T - N = (T - R) + d_area + d_mm + d_fill + d_resid

- N, the model's value, read from the scalar files that sit beside the
  gridded ones.
- T, the tool's value.
- R, the tool's expressions replayed here on the same files as netCDF4 returns
  them. Masked arrays promote Python scalars, so nearly all of the tool's
  arithmetic is float64 and a faithful replay agrees to roundoff. T - R
  therefore checks the reading and the time alignment.
- d_area, R minus the same with af2 = 1. af2 = (1/k)^2 is the EPSG:3031 area
  factor, which map-plane native scalars do not carry. A forward with the
  area factor integrates its scalars over true area (issue #97), and the
  writer stamps which area a tree's scalars carry. With --native-af2, which
  the stamp must agree with, every replay below keeps af2, d_area is zero
  and the column R_noaf2 is R. The model takes af2 at its cell centroids and
  the tool at pixel centres, so the exact-sum gates then allow af2's largest
  change between neighbouring pixels (5.9e-4 at 8 km) of each sum's L1; a
  missing or doubled factor is 1 to 5 %.
- d_mm, the maximum-extent mask the tool applies to lim, limnsw and sea level.
- d_fill, the writer's fill convention undone. The tool sums every flux over
  whole pixels. A tree whose flux files carry ``flux_pixel_mean =
  whole_pixel`` holds whole-pixel means (issue #96), so nothing is undone and
  d_fill is zero. In a tree written before, acabf is a mean over the covered
  part of a pixel, so it is multiplied by the pixel coverage from the writer's
  overlap cache, and libmassbffl a mean over the part floating at year end, so
  it is multiplied by sftflf.
- d_resid, what is left against N. For tendlicalvf, tendlifmassbf,
  tendligroundf, the ice area iareagr + iareafl and tendacabf (whole-pixel
  means, or the overlap cache) it is zero up to float32 and the CSV's seven
  digits, because the writer's remap is conservative. A writer that writes
  floating cells within a centimetre of flotation as grounded (isschecker's
  elevation tolerance) moves area from iareafl to iareagr, so each of the two
  may differ while their sum holds; the moved area is reported. For lim
  d_resid is the mass in cells of 1 m or less, which lithk leaves out; for
  limnsw the pixel averaging of a nonlinear integrand; for tendlibmassbffl the
  melt booked where the grid holds none: in pixels with no floating ice at
  year end, or, in a tree written before, outside the year-end floating mask.
  Since the forward books the melt of ice flowing into cells holding no ice at
  either end of the year as lifmassbf (issue #109), that residual is the melt
  of shelf ice gone within the year and the share the frozen apparent-MB
  reference supplied; a tree from before carries the front melt in it too,
  and its lifmassbf is zero.

The sea-level contributions get native counterparts from the model's lim and
limnsw, with the tool's ocean area A_O = 3.625e14 m2:
slvaf = -d(limnsw)/(rho_fw A_O), and with a fixed bed
slg20 = -[d(limnsw)/rho_sw + d(lim)(1/rho_fw - 1/rho_sw)]/A_O, which is also
sla20, since the tool's slc_A2020 reduces to slc_G2020 cell by cell when the
bed does not move.

The reference is either the stamped year the tool was given with --refyear,
for a run that is its own reference (``--hist ssp585 --refyear 2016``), or the
last state of the historical named by --hist-submission, which is the tool's
default when it pairs a projection with its historical run.

Gates are identities only; the named differences are reported, and only
--strict fails on them. Exit status 0 when every gate holds, 1 when one fails,
2 when an input is missing or the files do not line up.
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

# The model's densities (RHO_I in icepack2_tools/ismip7_output.py, and
# simulation.py's rho_ratio = 917/1024) and its year, restated here so that the
# comparison runs without Firedrake.
RHO_I = 917.0
RHO_SW = 1024.0
RHO_FW = 1000.0
SECONDS_PER_YEAR = 31557600.0
# The tool normalises every submission by this ocean area (Gregory et al. 2019).
OCEAN_AREA = 3.625e14

ST_SCALARS = ("lim", "limnsw", "iareagr", "iareafl")
FL_SCALARS = (("tendacabf", "acabf"), ("tendlibmassbfgr", "libmassbfgr"),
              ("tendlibmassbffl", "libmassbffl"), ("tendlicalvf", "licalvf"),
              ("tendlifmassbf", "lifmassbf"), ("tendligroundf", "ligroundf"))
SEA_LEVEL = ("slvaf", "slg20", "sla20")
SCALARS = ST_SCALARS + tuple(s for s, _ in FL_SCALARS) + SEA_LEVEL
UNITS = dict({s: "kg" for s in ("lim", "limnsw")}, iareagr="m2", iareafl="m2",
             **{s: "kg s-1" for s, _ in FL_SCALARS}, **{s: "m" for s in SEA_LEVEL})
# forbidden-policy fields: whole-pixel means, so a grid sum is the mesh sum
# (tendacabf holds too in a tree of whole-pixel means, under its own gate).
# lifmassbf carries front melt since issue #109 and is zero in a tree from
# before, which the same gate passes.
EXACT = ("tendlicalvf", "tendlifmassbf", "tendligroundf")
ZERO = ("tendlibmassbfgr",)
GRID_FILES = {"af2": "af2_AIS_{res}000m_v1.nc", "maxmask1": "maxmask1_AIS_{res}000m_v0.nc"}
MARKERS = (2015, 2050, 2100, 2200, 2300)
# The writer's global attribute on its gridded flux files (FLUX_MEAN_ATTR and
# FLUX_MEAN in write_ismip7_output.py, restated so that this runs without
# Firedrake). A tree written before it carries none.
FLUX_MEAN_ATTR = "flux_pixel_mean"
WHOLE_PIXEL = "whole_pixel"
COVERED_PART = "covered_part"
# The writer's stamp on the ten scalar files (SCALAR_AREA_ATTR in
# write_ismip7_output.py): the area the run's scalars integrate over.
SCALAR_AREA_ATTR = "scalar_area"
TRUE_AREA, MAP_PLANE = "true_area", "map_plane"

# Tolerances. F32 is twice float32's half ulp: every gridded value and every
# submitted scalar was rounded to float32 once. CSV is twice the half digit of
# the scalars CSV's seven significant figures.
F32 = 1.2e-7
CSV_DIGITS = 1.0e-6
REPLAY = 1.0e-9          # T - R, relative to the L1 of the integrand
SEA_LEVEL_ABS = 1.0e-9   # m, T - R and the slg20 - slvaf identity
# sla20 against slg20: slc_A2020 rounds the flotation thickness in another
# order than slc_G2020, so the two agree to the algebra and not bit for bit;
# a micrometre of sea level is far below any bed or mask error they catch.
A2020_ABS = 1.0e-6
RESID_SHARE = 0.02       # --strict: limnsw and sea-level residuals


class InputError(RuntimeError):
    r"""A missing or inconsistent input: exit status 2."""


def netcdf():
    import netCDF4
    return netCDF4


def one_file(folder, var):
    r"""The single ``<var>_*.nc`` in ``folder``; ours and the tool's share names."""
    hits = sorted(glob.glob(os.path.join(folder, f"{var}_*.nc")))
    if len(hits) != 1:
        raise InputError(f"{len(hits)} files match {var}_*.nc in {folder}; one is needed")
    return hits[0]


def nominal_years(ds, flux):
    r"""Time values and nominal years of an ISMIP7 file. A state is stamped
    1 January of the year after its nominal year, a flux 1 July of the year
    itself (ismip/ismip7-time-encoding); anything else is refused."""
    nc4 = netcdf()
    t = ds.variables["time"]
    values = np.asarray(t[:], dtype=np.float64)
    dates = nc4.num2date(values, t.units, calendar=getattr(t, "calendar", "standard"))
    want = (7, 1) if flux else (1, 1)
    for d in dates:
        if (d.month, d.day) != want:
            raise InputError(f"{ds.filepath()}: time stamp {d} is not "
                             f"{'1 July' if flux else '1 January'}")
    return values, [d.year if flux else d.year - 1 for d in dates]


def read_scalar(folder, var, flux):
    r"""``{nominal year: (time value, value)}`` from a scalar file."""
    nc4 = netcdf()
    path = one_file(folder, var)
    with nc4.Dataset(path) as ds:
        values, years = nominal_years(ds, flux)
        data = np.ma.filled(np.ma.asarray(ds.variables[var][:], dtype=np.float64), np.nan)
    return path, {y: (float(t), float(v)) for y, t, v in zip(years, values, data)}


class Gridded:
    r"""One gridded variable of a submission, read a year at a time."""

    def __init__(self, folder, var, flux):
        nc4 = netcdf()
        self.path = one_file(folder, var)
        self.ds = nc4.Dataset(self.path)
        self.var = var
        self.values, years = nominal_years(self.ds, flux)
        self.index = {y: k for k, y in enumerate(years)}

    def year(self, y):
        if y not in self.index:
            raise InputError(f"{self.path} has no year {y}")
        return self.ds.variables[self.var][self.index[y], :, :]

    def last(self):
        return self.ds.variables[self.var][len(self.values) - 1, :, :]

    def close(self):
        self.ds.close()


def read_grid(datapath, name, res, x, y):
    r"""One auxiliary grid as the tool reads it (a masked array), checked
    against the submission's axes: the tool reads it by array position."""
    nc4 = netcdf()
    path = os.path.join(datapath, GRID_FILES[name].format(res=res))
    if not os.path.exists(path):
        raise InputError(f"missing {path}")
    with nc4.Dataset(path) as ds:
        arr = ds.variables[name][:, :]
        gx, gy = np.asarray(ds.variables["x"][:]), np.asarray(ds.variables["y"][:])
    if gx.shape != x.shape or gy.shape != y.shape or not (
            np.allclose(gx, x, rtol=0, atol=1e-3) and np.allclose(gy, y, rtol=0, atol=1e-3)):
        raise InputError(f"{path}: its x/y are not the submission's grid "
                         f"(the tool pairs pixels by array position)")
    if np.ma.count_masked(arr):
        raise InputError(f"{path}: {np.ma.count_masked(arr)} masked cells")
    return path, arr


def read_params(path):
    r"""The densities the tool used, from the params.nc it was given."""
    nc4 = netcdf()
    if not os.path.exists(path):
        raise InputError(f"missing {path}")
    with nc4.Dataset(path) as ds:
        return {k: float(ds.variables[k][()]) for k in ("rhoi", "rhow", "rhof")}


def pixel_coverage(npz, shape, pixel_area):
    r"""Covered share of every pixel, from the writer's cached overlap
    operator (pixels by model cells, exact overlap areas, CSR arrays)."""
    z = np.load(npz)
    npix = int(z["shape"][0])
    if npix != shape[0] * shape[1]:
        raise InputError(f"{npz}: {npix} rows for a {shape[0]} x {shape[1]} grid")
    rows = np.repeat(np.arange(npix), np.diff(z["indptr"]))
    covered = np.bincount(rows, weights=z["data"], minlength=npix)
    return (covered / pixel_area).reshape(shape)


def flux_means(grids):
    r"""How the writer took the pixel means of acabf and libmassbffl:
    WHOLE_PIXEL when their files say so, COVERED_PART for a tree written
    before the attribute (acabf over the covered part of a pixel, libmassbffl
    over the part floating at year end)."""
    seen = {}
    for v in ("acabf", "libmassbffl"):
        ds = grids[v].ds
        seen[v] = ds.getncattr(FLUX_MEAN_ATTR) if FLUX_MEAN_ATTR in ds.ncattrs() else None
    if seen["acabf"] != seen["libmassbffl"]:
        raise InputError(f"acabf and libmassbffl disagree on {FLUX_MEAN_ATTR}: "
                         f"{seen['acabf']!r} and {seen['libmassbffl']!r}")
    if seen["acabf"] not in (None, WHOLE_PIXEL):
        raise InputError(f"{FLUX_MEAN_ATTR} = {seen['acabf']!r} is no convention this "
                         f"comparison knows")
    return WHOLE_PIXEL if seen["acabf"] else COVERED_PART


def scalar_area_stamp(folder):
    r"""The area the writer found the native scalars of ``folder`` to
    integrate over, from its stamp on the lim file; None for a tree written
    before the stamp."""
    nc4 = netcdf()
    with nc4.Dataset(one_file(folder, "lim")) as ds:
        return ds.getncattr(SCALAR_AREA_ATTR) if SCALAR_AREA_ATTR in ds.ncattrs() else None


def af2_step(af2):
    r"""The largest relative change of af2 between neighbouring pixels.

    The model takes af2 at its cell centroids and the tool at pixel centres.
    The tool's factor for a cell is the overlap-weighted mean over the pixels
    it covers, which is af2 within about half a pixel of the centroid, plus
    the curvature of af2 over the cell (4e-5 for a 180 km cell). So over any
    sum the two differ by less than this share of its L1. Over the 8 km grid
    it is 5.9e-4."""
    a = np.ma.filled(af2, np.nan).astype(np.float64)
    steps = [np.abs(a[1:, :] - a[:-1, :]) / np.minimum(a[1:, :], a[:-1, :]),
             np.abs(a[:, 1:] - a[:, :-1]) / np.minimum(a[:, 1:], a[:, :-1])]
    return float(max(np.nanmax(s) for s in steps))


class Replay:
    r"""The tool's integrals (ismip7_scalars 0.1.0, ``scalars.py`` and
    ``slc/``), written as the same expressions on the same masked arrays.

    ``af2`` and ``maxmask1`` are swapped for ones to take the area factor and
    the maximum-extent mask out one at a time (af2 stays under --native-af2);
    the arithmetic keeps its types.
    """

    def __init__(self, af2, maxmask1, area_m2, rho):
        self.mm = maxmask1
        self.rhoi, self.rhow, self.rhof = rho["rhoi"], rho["rhow"], rho["rhof"]
        region = maxmask1 * 0 + 1
        self.A = region * af2 * area_m2            # compute_st_series, compute_slc_series
        self.W = region * (af2 * area_m2)          # run_fl_scalars: region_mask * weight_base
        self.ref = None

    def state(self, lithk, topg, sftgrf, sftflf):
        H = lithk * self.mm
        hf = np.maximum(-topg, 0) * self.rhow / self.rhoi
        return {"lim": float(np.sum(H * self.A) * self.rhoi),
                "limnsw": float(np.sum(np.maximum(H - hf, 0) * self.A) * self.rhoi),
                "iareagr": float(np.sum(sftgrf * self.A)),
                "iareafl": float(np.sum(sftflf * self.A))}

    def flux(self, field):
        return float(np.einsum("yx,yx->", np.ma.filled(field, 0.0), self.W))

    def _vaf(self, H, B, S):                      # slc_vaf.get_vaf
        hf = np.maximum(S - B, 0.0) * self.rhow / self.rhoi
        return np.sum(np.maximum(H - hf, 0.0) * self.A)

    def _vaf_g20(self, H, B):                     # slc_G2020.get_vaf_G2020
        hf = np.minimum(B, 0.0) * self.rhow / self.rhoi
        return np.sum(np.maximum(H + hf, 0.0) * self.A)

    def _pov(self, B):
        return np.sum(np.maximum(-B, 0.0) * self.A)

    def _den(self, H):
        return np.sum(H * (self.rhoi / self.rhof - self.rhoi / self.rhow) * self.A)

    def reference(self, lithk_ref, topg_ref):
        H0 = lithk_ref * self.mm
        S0 = topg_ref * 0.0                       # sea level fixed at zero
        self.S0 = S0
        self.ref = {"vaf": self._vaf(H0, topg_ref, S0), "vaf_g20": self._vaf_g20(H0, topg_ref),
                    "pov": self._pov(topg_ref), "den": self._den(H0)}
        hf0 = np.maximum(-topg_ref, 0) * self.rhow / self.rhoi
        self.lim0 = float(np.sum(H0 * self.A) * self.rhoi)
        self.limnsw0 = float(np.sum(np.maximum(H0 - hf0, 0) * self.A) * self.rhoi)

    def sea_level(self, lithk, topg):
        r"""(slvaf, slg20) against the reference state, in m."""
        H, B, AO, r = lithk * self.mm, topg, OCEAN_AREA, self.ref
        slvaf = -(self._vaf(H, B, self.S0) / AO * self.rhoi / self.rhof
                  - r["vaf"] / AO * self.rhoi / self.rhof)
        af = -(self._vaf_g20(H, B) / AO * self.rhoi / self.rhow
               - r["vaf_g20"] / AO * self.rhoi / self.rhow)
        pov = -(self._pov(B) / AO - r["pov"] / AO)
        den = -(self._den(H) / AO - r["den"] / AO)
        return float(slvaf), float(af + pov + den)


def same_field(a, b):
    r"""Bitwise equality of two masked planes, masks included."""
    ma, mb = np.ma.getmaskarray(a), np.ma.getmaskarray(b)
    return np.array_equal(ma, mb) and np.array_equal(np.ma.getdata(a)[~ma], np.ma.getdata(b)[~mb])


class Gate:
    r"""An identity checked every year: pass while |excess| <= tolerance."""

    def __init__(self, name):
        self.name, self.worst, self.checked = name, None, 0

    def check(self, year, excess, tol):
        self.checked += 1
        ratio = abs(excess) / tol if tol > 0 else (0.0 if excess == 0 else np.inf)
        if self.worst is None or ratio > self.worst[0]:
            self.worst = (ratio, year, excess, tol)

    @property
    def ok(self):
        return self.worst is None or self.worst[0] <= 1.0


def native_sea_level(lim, limnsw, lim0, limnsw0):
    r"""slvaf and slg20 (= sla20 with a fixed bed) from the model's own lim
    and limnsw, in the tool's conventions."""
    slvaf = -(limnsw - limnsw0) / (RHO_FW * OCEAN_AREA)
    slg20 = -((limnsw - limnsw0) / RHO_SW + (lim - lim0) * (1.0 / RHO_FW - 1.0 / RHO_SW)) / OCEAN_AREA
    return slvaf, slg20


def compare(a):
    nc4 = netcdf()
    if os.path.exists(a.tool) and os.path.samefile(a.submission, a.tool):
        raise InputError("--tool is the submission folder; the tool's files carry our names "
                         "and must be read from its own output tree")
    rho = read_params(a.params)
    fl_names = [v for _, v in FL_SCALARS] + ["dlithkdt"]
    grids = {v: Gridded(a.submission, v, flux=False) for v in ("lithk", "topg", "sftgrf", "sftflf")}
    grids.update({v: Gridded(a.submission, v, flux=True) for v in fl_names})
    means = flux_means(grids)
    whole = means == WHOLE_PIXEL
    # the area the native scalars integrate over: the switch, held to the
    # writer's stamp wherever the stamp exists
    told = TRUE_AREA if a.native_af2 else MAP_PLANE
    stamps = {}
    for folder in [a.submission] + ([a.hist_submission] if a.hist_submission else []):
        stamps[folder] = scalar_area_stamp(folder)
        if stamps[folder] not in (None, told):
            raise InputError(f"{folder}: the writer stamped the native scalars {stamps[folder]}, "
                             f"and the comparison was told {told} "
                             f"({'--native-af2' if a.native_af2 else 'no --native-af2'})")
    lith = grids["lithk"]
    x = np.asarray(lith.ds.variables["x"][:], dtype=np.float64)
    y = np.asarray(lith.ds.variables["y"][:], dtype=np.float64)
    dx = abs(float(x[1] - x[0]))
    res = f"{round(dx / 1000):02d}"                  # naming.resolution_string
    area_m2 = (float(res) * 1000.0) ** 2              # scalars.run
    af2_path, af2 = read_grid(a.datapath, "af2", res, x, y)
    mm_path, mm = read_grid(a.datapath, "maxmask1", res, x, y)
    cov = pixel_coverage(a.overlap, af2.shape, dx * dx) if a.overlap else None

    years = sorted(lith.index)
    native, tool, paths = {}, {}, {}
    for s in ST_SCALARS + tuple(s for s, _ in FL_SCALARS):
        flux = s not in ST_SCALARS
        paths[s], native[s] = read_scalar(a.submission, s, flux)
        _, tool[s] = read_scalar(a.tool, s, flux)
    for s in SEA_LEVEL:
        _, tool[s] = read_scalar(a.tool, s, False)
    for s in SCALARS:
        have = native[s] if s in native else native["lim"]
        for yr in years:
            if yr not in tool[s] or yr not in have:
                raise InputError(f"{s} has no value for {yr} in "
                                 f"{'the tool output' if yr not in tool[s] else 'the submission'}")
            if tool[s][yr][0] != have[yr][0]:
                raise InputError(f"{s} {yr}: the tool's time {tool[s][yr][0]} is not ours "
                                 f"{have[yr][0]}")

    # the reference state and the model's own lim and limnsw there
    if a.hist_submission:
        hist = {v: Gridded(a.hist_submission, v, flux=False) for v in ("lithk", "topg")}
        lithk_ref, topg_ref = hist["lithk"].last(), hist["topg"].last()
        hlast = max(hist["lithk"].index)
        _, h_lim = read_scalar(a.hist_submission, "lim", False)
        _, h_limnsw = read_scalar(a.hist_submission, "limnsw", False)
        lim0, limnsw0 = h_lim[hlast][1], h_limnsw[hlast][1]
        reference = f"the last state of {a.hist_submission} (nominal {hlast})"
        for g in hist.values():
            g.close()
    else:
        yref = a.refyear - 1
        if yref not in lith.index:
            raise InputError(f"--refyear {a.refyear}: no state stamped {a.refyear} in {lith.path}")
        lithk_ref, topg_ref = lith.year(yref), grids["topg"].year(yref)
        lim0, limnsw0 = native["lim"][yref][1], native["limnsw"][yref][1]
        reference = f"the state stamped {a.refyear} (nominal {yref}), the run's own"

    # With --native-af2 both sides carry the area factor, so the replays that
    # take the named differences out one at a time keep it (d_area is zero)
    # and the exact sums allow for where each side samples it.
    af2_n = af2 if a.native_af2 else np.ones_like(af2)
    w64 = np.ma.filled(af2_n, 1.0).astype(np.float64)   # 1.0 exactly without the switch
    allow = af2_step(af2) if a.native_af2 else 0.0
    ones_mm = np.ones_like(mm)
    replays = {"R": Replay(af2, mm, area_m2, rho), "Rn": Replay(af2_n, mm, area_m2, rho),
               "P": Replay(af2_n, ones_mm, area_m2, rho)}
    for r in replays.values():
        r.reference(lithk_ref, topg_ref)

    csv_rows = None
    if a.native_csv:
        with open(a.native_csv) as fh:
            csv_rows = {int(r["year"]): r for r in csv.DictReader(fh)}

    gates = {name: Gate(name) for name in (
        "densities", "scalar files against the CSV", "T against the replay R",
        "forbidden-policy sums against N", "tendacabf against N",
        "zero fluxes", "sla20 = slg20 (fixed bed)", "slg20 - slvaf identity",
        "topg constant in time", "lim change against dlithkdt", "lim residual sign",
        "T differs from N", "af2 > 0 under ice and flux")}
    for k, want in (("rhoi", RHO_I), ("rhow", RHO_SW), ("rhof", RHO_FW)):
        gates["densities"].check(k, rho[k] - want, 0.0)

    rows, carries, lim_prev = [], np.zeros(af2.shape, dtype=bool), None
    for yr in years:
        f = {v: g.year(yr) for v, g in grids.items()}
        carries |= np.ma.filled(f["lithk"] > 0, False)
        for v in fl_names[:-1]:
            carries |= np.ma.filled(f[v] != 0, False)
        gates["topg constant in time"].check(yr, 0.0 if same_field(f["topg"], topg_ref) else 1.0, 0.0)

        out = {}
        for key, rp in replays.items():
            out[key] = rp.state(f["lithk"], f["topg"], f["sftgrf"], f["sftflf"])
            out[key].update({s: rp.flux(f[v]) for s, v in FL_SCALARS})
            out[key]["slvaf"], out[key]["slg20"] = rp.sea_level(f["lithk"], f["topg"])
            out[key]["sla20"] = out[key]["slg20"]    # identity, fixed bed
        C = dict(out["P"])
        if not whole:
            if cov is not None:
                C["tendacabf"] = float(np.sum(np.ma.filled(f["acabf"], 0.0).astype(np.float64)
                                              * cov * w64) * area_m2)
            C["tendlibmassbffl"] = float(np.sum(np.ma.filled(f["libmassbffl"], 0.0).astype(np.float64)
                                                * np.ma.filled(f["sftflf"], 0.0).astype(np.float64)
                                                * w64) * area_m2)
        lithk64 = np.ma.filled(f["lithk"], 0.0).astype(np.float64)
        L1 = {"lim": RHO_I * np.sum(np.abs(lithk64) * w64) * area_m2, "limnsw": abs(out["P"]["limnsw"]),
              "iareagr": abs(out["P"]["iareagr"]), "iareafl": abs(out["P"]["iareafl"])}
        for s, v in FL_SCALARS:
            L1[s] = float(np.sum(np.abs(np.ma.filled(f[v], 0.0).astype(np.float64)) * w64) * area_m2)

        N = {s: native[s][yr][1] for s in native}
        N["slvaf"], N["slg20"] = native_sea_level(N["lim"], N["limnsw"], lim0, limnsw0)
        N["sla20"] = N["slg20"]
        for s in SCALARS:
            T, R, Rn, P = tool[s][yr][1], out["R"][s], out["Rn"][s], out["P"][s]
            rows.append({"scalar": s, "year": yr, "N": N[s], "T": T, "R": R, "R_noaf2": Rn,
                         "P": P, "C": C[s], "T_minus_N": T - N[s], "T_minus_R": T - R,
                         "d_area": R - Rn, "d_mm": Rn - P, "d_fill": P - C[s],
                         "d_resid": C[s] - N[s], "L1": L1.get(s, np.nan)})
            if s == "sla20":
                pass                                 # its own identity gate, below
            elif s in SEA_LEVEL:
                gates["T against the replay R"].check(f"{s} {yr}", T - R, SEA_LEVEL_ABS)
            else:
                gates["T against the replay R"].check(f"{s} {yr}", T - R, REPLAY * max(L1[s], abs(R)))
                gates["T differs from N"].check(f"{s} {yr}", float(T == N[s] and N[s] != 0.0), 0.0)
        # float32 and the CSV's digits, plus where the two sides sample af2
        # when both carry it (allow is zero otherwise)
        rel = F32 + allow
        for s in EXACT:
            gates["forbidden-policy sums against N"].check(
                f"{s} {yr}", C[s] - N[s], CSV_DIGITS * abs(N[s]) + rel * L1[s])
        # the two areas may trade near flotation, their sum may not, and area
        # only ever moves to the grounded side
        area_tol = CSV_DIGITS * (abs(N["iareagr"]) + abs(N["iareafl"])) + rel * (L1["iareagr"] + L1["iareafl"])
        gates["forbidden-policy sums against N"].check(
            f"ice area {yr}", (C["iareagr"] + C["iareafl"]) - (N["iareagr"] + N["iareafl"]), area_tol)
        gates["forbidden-policy sums against N"].check(
            f"grounded area gained {yr}", min(C["iareagr"] - N["iareagr"], 0.0), area_tol)
        if whole or cov is not None:
            gates["tendacabf against N"].check(
                yr, C["tendacabf"] - N["tendacabf"],
                CSV_DIGITS * abs(N["tendacabf"]) + rel * L1["tendacabf"])
        for s in ZERO:
            gates["zero fluxes"].check(f"{s} {yr}", max(abs(N[s]), abs(tool[s][yr][1]), L1[s]), 0.0)
        gates["lim residual sign"].check(yr, max(C["lim"] - N["lim"], 0.0),
                                         CSV_DIGITS * abs(N["lim"]) + rel * L1["lim"])
        slg, slv = tool["slg20"][yr][1], tool["slvaf"][yr][1]
        gates["sla20 = slg20 (fixed bed)"].check(yr, tool["sla20"][yr][1] - slg, A2020_ABS)
        R0 = replays["R"]
        gates["slg20 - slvaf identity"].check(
            yr, (slg - slv) - (1.0 / rho["rhof"] - 1.0 / rho["rhow"])
            * ((tool["limnsw"][yr][1] - R0.limnsw0) - (tool["lim"][yr][1] - R0.lim0)) / OCEAN_AREA,
            SEA_LEVEL_ABS)
        if lim_prev is not None:
            dh64 = np.ma.filled(f["dlithkdt"], 0.0).astype(np.float64)
            d_lim = RHO_I * SECONDS_PER_YEAR * float(np.sum(dh64 * w64)) * area_m2
            l1 = RHO_I * SECONDS_PER_YEAR * float(np.sum(np.abs(dh64) * w64)) * area_m2
            gates["lim change against dlithkdt"].check(
                yr, (N["lim"] - lim_prev) - d_lim, CSV_DIGITS * abs(N["lim"]) + rel * l1)
        lim_prev = N["lim"]
        if csv_rows is not None:
            if yr not in csv_rows:
                raise InputError(f"{a.native_csv} has no row for {yr}")
            for s in native:
                c = float(csv_rows[yr][s])
                gates["scalar files against the CSV"].check(f"{s} {yr}", N[s] - c, F32 * abs(c))
    bad = int(np.sum(carries & ~np.ma.filled(af2 > 0, False)))
    gates["af2 > 0 under ice and flux"].check("all years", float(bad), 0.0)
    for g in grids.values():
        g.close()
    return {"rows": rows, "gates": gates, "rho": rho, "reference": reference,
            "years": years, "paths": paths,
            "grids": {"af2": (af2_path, str(af2.dtype)), "maxmask1": (mm_path, str(mm.dtype))},
            "maxmask_pixels_dropped": int(np.sum(carries & ~np.ma.filled(mm > 0, False))),
            "coverage": a.overlap, "flux_means": means,
            "native_af2": a.native_af2, "af2_allow": allow, "stamps": stamps}


def warnings_for(res, strict_share=RESID_SHARE):
    r"""The named differences worth reading, one line each. They are the
    findings, and they only fail the run under --strict."""
    by = {}
    for r in res["rows"]:
        by.setdefault(r["scalar"], []).append(r)
    scale = {s: max(abs(r["N"]) for r in rows) for s, rows in by.items()}
    out = []

    def share(s, key):
        m = max(by[s], key=lambda r: abs(r[key]))
        return m[key] / scale[s] if scale[s] else 0.0, m["year"]

    whole = res["flux_means"] == WHOLE_PIXEL
    if res["native_af2"]:
        # the exact sums, as a share of their L1: how closely the model's af2
        # at its cell centroids reproduces the tool's at pixel centres
        exact = list(EXACT) + (["tendacabf"] if whole or res["coverage"] else [])
        worst = max([(abs(r["d_resid"]) / r["L1"], r["scalar"], r["year"])
                     for s in exact for r in by[s] if r["L1"] > 0]
                    + [(abs(g["d_resid"] + f["d_resid"]) / (g["L1"] + f["L1"]), "ice area", g["year"])
                       for g, f in zip(by["iareagr"], by["iareafl"]) if g["L1"] + f["L1"] > 0])
        out.append(f"area factor: the native scalars carry it too (--native-af2), so the area "
                   f"term is zero; the exact sums' worst |C - N| is {worst[0]:.2e} of their L1 "
                   f"({worst[1]} {worst[2]}), against an allowance of {res['af2_allow']:.2e}")
    else:
        area = [f"{s} {100 * share(s, 'd_area')[0]:+.2f}%" for s in by if scale[s]]
        out.append("area factor (af2), largest share of max |N|: " + ", ".join(area))
    mm = [f"{s} {100 * share(s, 'd_mm')[0]:+.3f}%" for s in ("lim", "limnsw", "slvaf", "slg20")
          if scale.get(s)]
    out.append(f"maximum-extent mask: {res['maxmask_pixels_dropped']} pixels carry ice or "
               f"flux outside maxmask1; " + ", ".join(mm))
    if whole:
        out.append(f"fill convention: acabf and libmassbffl are whole-pixel means "
                   f"({FLUX_MEAN_ATTR}), the tool's own, so nothing is undone and the "
                   f"fill term is zero")
    else:
        for s in ("tendacabf", "tendlibmassbffl"):
            v, yr = share(s, "d_fill")
            out.append(f"fill convention, {s}: the tool's whole-pixel sum minus the "
                       f"coverage-weighted one is {100 * v:+.2f}% of max |N| ({yr})")
    moved = [r for r in by["iareagr"]
             if r["d_resid"] > CSV_DIGITS * abs(r["N"]) + (F32 + res["af2_allow"]) * r["L1"]]
    if moved:
        m = max(moved, key=lambda r: r["d_resid"])
        out.append(f"near flotation: {len(moved)} years write floating area as grounded, at most "
                   f"{m['d_resid']:.3e} m2 ({m['year']})")
    if not res["coverage"] and not whole:
        out.append("no --overlap: tendacabf's fill convention is not undone, so its "
                   "residual carries it")
    m = max(by["tendlibmassbffl"], key=lambda r: abs(r["d_resid"]))
    where = ("in pixels with no floating ice at year end, which the fill leaves out" if whole
             else "outside the writer's year-end floating mask")
    out.append(f"tendlibmassbffl: the model's value carries "
               f"{-m['d_resid'] * SECONDS_PER_YEAR / 1e12:+.1f} Gt/yr {where} ({m['year']})")
    front = max(by["tendlifmassbf"], key=lambda r: abs(r["N"]))
    out.append(f"tendlifmassbf: the front melt booked in cells holding no ice at either end "
               f"of the year (issue #109) is at most "
               f"{front['N'] * SECONDS_PER_YEAR / 1e12:+.1f} Gt/yr ({front['year']})"
               if front["N"] else
               "tendlifmassbf: zero in every year (a tree from before the front-melt booking "
               "of issue #109, or no cell emptied)")
    for s in ("limnsw",) + SEA_LEVEL:
        v, yr = share(s, "d_resid")
        if abs(v) > strict_share:
            out.append(f"{s}: residual {100 * v:+.2f}% of max |N| ({yr}), above "
                       f"{100 * strict_share:.0f}%")
    return out


def write_csv(path, rows):
    keys = list(rows[0])
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.9e}" if isinstance(v, float) else v) for k, v in r.items()})


def write_markdown(path, res, a, warns, status):
    nc4 = netcdf()
    by = {}
    for r in res["rows"]:
        by.setdefault(r["scalar"], []).append(r)
    years = res["years"]
    lines = [f"# Scalar comparison: {a.label or os.path.basename(os.path.normpath(a.submission))}", "",
             f"- submission: `{a.submission}`",
             f"- tool output: `{a.tool}`",
             f"- grids: `{res['grids']['af2'][0]}` ({res['grids']['af2'][1]}), "
             f"`{res['grids']['maxmask1'][0]}` ({res['grids']['maxmask1'][1]})",
             f"- densities in params.nc: {res['rho']['rhoi']:g} / {res['rho']['rhow']:g} / "
             f"{res['rho']['rhof']:g}",
             f"- reference: {res['reference']}",
             f"- years: {years[0]} to {years[-1]} ({len(years)})",
             f"- flux means: whole pixel (`{FLUX_MEAN_ATTR}`)" if res["flux_means"] == WHOLE_PIXEL
             else f"- flux means: acabf over the covered part of a pixel, libmassbffl over the part "
                  f"floating at year end (no `{FLUX_MEAN_ATTR}`)",
             f"- pixel coverage: `{res['coverage']}`" if res["coverage"] else "- pixel coverage: none given",
             f"- native scalars: over {'true' if res['native_af2'] else 'map-plane'} area; the "
             f"writer's stamp: " + ", ".join(f"{v or 'none'} (`{k}`)" for k, v in res["stamps"].items()),
             f"- comparison: numpy {np.__version__}, netCDF4 {nc4.__version__}"]
    lines += [f"- {n}" for n in (a.note or [])]
    lines += ["", f"Exit status {status}.", "", "## Gates", "", "| gate | result | worst |", "|---|---|---|"]
    for g in res["gates"].values():
        if g.worst is None:
            lines.append(f"| {g.name} | not checked | |")
            continue
        ratio, where, excess, tol = g.worst
        lines.append(f"| {g.name} | {'pass' if g.ok else 'FAIL'} | {excess:.3e} against "
                     f"{tol:.3e} at {where} |")
    lines += ["", "## Named differences", ""] + [f"- {w}" for w in warns]
    lines += ["", "## By scalar", "",
              "T - N and its parts, as shares of max |N| over the run; sea level in mm.", ""]
    marks = sorted(set(y for y in MARKERS if y in years) | {years[0], years[-1]})
    for s in SCALARS:
        rows = {r["year"]: r for r in by[s]}
        scale = max(abs(r["N"]) for r in by[s])
        sea = s in SEA_LEVEL
        worst = max(by[s], key=lambda r: abs(r["T_minus_N"]))
        lines += [f"### {s} ({'mm' if sea else UNITS[s]})", ""]
        if scale == 0.0 and not sea:
            lines += [f"zero on both sides, {len(by[s])} years", ""]
            continue
        if sea:
            def f(v):
                return f"{1e3 * v:+.3f}"
            lines.append("| year | N | T | T - N | area | max mask | fill | residual |")
        else:
            def f(v):
                return f"{100 * v / scale:+.3f}%"
            lines.append(f"max |N| = {scale:.4e}")
            lines.append("")
            lines.append("| year | N | T | T - N | area | max mask | fill | residual |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for yr in sorted(set(marks) | {worst["year"]}):
            r = rows[yr]
            n, t = (f"{1e3 * r['N']:.3f}", f"{1e3 * r['T']:.3f}") if sea else \
                (f"{r['N']:.4e}", f"{r['T']:.4e}")
            tag = " (worst)" if yr == worst["year"] else ""
            lines.append(f"| {yr}{tag} | {n} | {t} | {f(r['T_minus_N'])} | {f(r['d_area'])} | "
                         f"{f(r['d_mm'])} | {f(r['d_fill'])} | {f(r['d_resid'])} |")
        lines.append("")
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--submission", required=True, help="our experiment folder, .../CORE/<counter>")
    ap.add_argument("--tool", required=True, help="the tool's folder for the same experiment, "
                    "<outpath>/nc/AIS/<source>/<ism>/CORE/<counter>")
    ap.add_argument("--datapath", required=True, help="the auxiliary grids (af2, maxmask1)")
    ap.add_argument("--params", required=True, help="the params.nc the tool read, the upload's "
                    "AIS/<source_id>/<ism_id>/params.nc")
    ref = ap.add_mutually_exclusive_group(required=True)
    ref.add_argument("--refyear", type=int, help="the stamped year given to the tool")
    ref.add_argument("--hist-submission", help="the historical's experiment folder, when the "
                     "tool paired the run with it")
    ap.add_argument("--native-csv", help="the run's _ismip7_scalars.csv, checked against the "
                    "scalar files")
    ap.add_argument("--native-af2", action="store_true",
                    help="the native scalars integrate over true area, map-plane cell area times "
                         "af2 at the centroid, as a forward with the area factor writes them; "
                         "the writer's scalar_area stamp must agree")
    ap.add_argument("--overlap", help="the writer's <annual>.overlap.npz, for pixel coverage; "
                    "only a tree without whole-pixel flux means needs it, to undo acabf's "
                    "covered-part means")
    ap.add_argument("--label", help="title of the summary")
    ap.add_argument("--note", action="append", help="a provenance line for the summary")
    ap.add_argument("--strict", action="store_true", help="fail on the named differences too")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-md", required=True)
    a = ap.parse_args(argv)
    try:
        res = compare(a)
    except InputError as e:
        print(f"compare_scalars: {e}", file=sys.stderr)
        return 2
    warns = warnings_for(res)
    failed = [g.name for g in res["gates"].values() if not g.ok]
    status = 1 if failed or (a.strict and warns) else 0
    write_csv(a.out_csv, res["rows"])
    write_markdown(a.out_md, res, a, warns, status)
    for g in res["gates"].values():
        print(f"  {'pass' if g.ok else 'FAIL'}  {g.name}")
    for w in warns:
        print(f"  note  {w}")
    print(f"compare_scalars exit status: {status} ({a.out_md})")
    return status


if __name__ == "__main__":
    sys.exit(main())
