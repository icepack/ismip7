#!/usr/bin/env python
r"""ISMIP7 output, part two: the model's yearly fields onto the 8 km AIS grid
in the request's files.

    python antarctica/scripts/write_ismip7_output.py ANNUAL.h5 --out-dir DIR
        --esm CESM2-WACCM --scenario ssp585 [--exp C007]
        [--source-id RICE] [--ism-id icepack2] [--set-id CORE] [--scalars CSV]

Serial. ANNUAL.h5 is the STEM ``<experiment>_<lc>_ismip7_annual.h5`` that
``icepack2_tools.ismip7_output`` names its output after; the years themselves
are the files ``<experiment>_<lc>_ismip7_annual_<year>.h5`` beside it (one per
year, each written atomically, so an interrupted run cannot damage the years
already banked) and this globs them in year order. It writes one NetCDF per
variable under
``DIR/AIS/<source_id>/<ism_id>/<set_id>/<exp>/``, named
``<var>_AIS_<source_id>_<ism_id>_m001_<ESM>_f001_<scenario>_<exp>_<y0>-<y1>.nc``.

Regridding is conservative: a supermesh mixed mass matrix between the model's
DG0 cells and a triangulated copy of the 8 km grid gives the exact area of
every (cell, pixel) overlap; a pixel's value is the area-weighted mean of the
cells under it. A flux (``FL``) is a mean over the whole pixel of everything
the cells under it booked, so its value times the pixel area sums to the
model's integral, which is how the organisers' scalar tool sums a flux (forum
thread 50). The flux files say so in the global attribute
``flux_pixel_mean``. For a flux the request's fill policy
(``isschecker/data/ISMIP7_variable_request.csv``) decides only where the
value is fill: ``outside_domain`` (``acabf``) where the model covers no part
of the pixel, ``no_floating_ice`` (``libmassbffl``) where no ice floats at
year end. Melt booked in such a pixel leaves ``libmassbffl``. The forward
books the melt of ice that flowed into a marine cell holding no ice at either
end of the year as ``lifmassbf`` instead (issue #109), whose ``forbidden``
policy never fills, so what leaves is the melt of shelf ice gone within the
year and the share the frozen apparent-MB reference supplied. The summary
lines report both, at their largest and at 2100, 2200 and 2300. Every year
file names its booking in the ``front_melt`` attribute, and a series that
mixes two bookings is refused. A state variable (``ST``) takes its mean over the area
the policy names: ``forbidden`` (thickness, fractions) over the whole pixel
with the uncovered part counting as zero, so sums over the grid are the
model's sums; ``outside_domain`` (elevations) over the covered part, filling
pixels the model does not cover; ``no_ice`` and friends over the ice part,
filling pixels without it. The overlap operator is cached next to the input
(``<annual>.overlap.npz``) because it depends on the mesh only.

The ten scalars are the run's CSV as it stands, after one check of the area
they integrate over. A current forward weights each cell by af2 = (1/k)^2 of
EPSG:3031 at its centroid, as ``ismip7-scalars`` weights its pixels, and a
forward from before summed map-plane area. Every year's ``iareagr`` and
``iareafl`` must match the annual file's sums under one of the two, the same
one in every year, and the scalar files name it in the global attribute
``scalar_area``.

Model-to-SI conversions use icepack's year, 365.25 days (31557600 s), which
is the model's own time unit; the time axis in the files is the standard
calendar regardless.

``acabf`` is written as the forcing surface mass balance, always. The
apparent-mass-balance reference stays where the forward put it, as
``acabf_correction`` in the annual file: it is not a request variable, and
folded into the SMB it would sit two orders of magnitude outside the
request's range. A reader who wants a grid budget that closes adds the two
from the annual file; the submission files never carry the sum.

Time follows ismip/ismip7-time-encoding: ``days since 1850-01-01`` on the
standard calendar; state variables are stamped 1 January of the following
year, fluxes 1 July with bounds over the year; the initial state is not
written.
"""
import argparse
import csv
import os
import re
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

import icepack2_tools.dual_friction  # noqa: F401,E402  (icepack2 -> irksome import order)
from icepack2_tools.ismip7_output import (AnnualOutput, FRONT_MELT_ATTR, RHO_I, SCALARS,  # noqa: E402
                                          SECONDS_PER_YEAR, VARIABLES_2D, VARIABLES_BANKED)
from icepack2_tools.regrid import (ISMIP7_DX, ISMIP7_NX, ISMIP7_NY, ISMIP7_X0, ISMIP7_Y0,  # noqa: E402
                                   area_factor)
import firedrake as fd  # noqa: E402

import netCDF4 as _nc4
FILL = _nc4.default_fillvals["f4"]          # the checker wants the netCDF4 default fill for the dtype
TIME_UNITS = "days since 1850-01-01"        # the checker's exact spelling
PIXEL_AREA = ISMIP7_DX * ISMIP7_DX
REQUEST = os.path.join(os.path.dirname(os.path.dirname(_ROOT)), "icepack2_tools", "ismip7_variable_request.csv")   # not under a data/ dir: .gitignore ignores those
# The global attribute on every gridded flux file, saying its pixel means are
# whole-pixel means (issue #96). compare_scalars.py reads it; a tree written
# before it carries none, and its acabf and libmassbffl are means over the
# covered part and over the floating part of a pixel.
FLUX_MEAN_ATTR = "flux_pixel_mean"
FLUX_MEAN = "whole_pixel"
# The global attribute on the ten scalar files, saying which area the run's
# scalars CSV integrates over (issue #97): ``true_area`` from a forward that
# weights each cell by af2, ``map_plane`` from one before. The writer reads
# which it is off the annual files, and compare_scalars.py holds its
# --native-af2 to it.
SCALAR_AREA_ATTR = "scalar_area"
TRUE_AREA, MAP_PLANE = "true_area", "map_plane"
# What a year file without the forward's front_melt stamp holds: a forward
# from before issue #109 wrote lifmassbf as zero and kept all the melt in
# libmassbffl.
UNSTAMPED_MELT = "none (lifmassbf zero, all melt in libmassbffl)"
# The years at which the melt summary reports its values, besides its
# largest year and the series' last.
MELT_MARKERS = (2100, 2200, 2300)


def standard_name(meta):
    r"""The CF standard name to write, or "" to omit the attribute.

    A blank standard_name column is the request stating that no standard name
    is assigned to that variable, so the attribute is omitted rather than
    written empty (an empty standard_name is not valid CF). The Comment column
    proposes names for some of them ("Should get standard name: ..."), but a
    proposal is not a registered CF name and a CF checker rejects it, so it is
    not used.
    """
    return (meta.get("standard_name") or "").strip()


def request_table():
    with open(REQUEST) as f:
        rows = list(csv.DictReader(f))
    return {r["Variable Name"]: r for r in rows}


# model units -> request units
CONVERT = {
    "m": 1.0, "1": 1.0,
    "m s-1": 1.0 / SECONDS_PER_YEAR,                      # m/yr -> m/s
    "kg m-2 s-1": RHO_I / SECONDS_PER_YEAR,               # m ice/yr -> kg m-2 s-1
    "Pa": 1.0e6,                                          # MPa -> Pa
}

# isschecker's ELEVATION_TOLERANCE (checker.py at 0.5.1): a wholly grounded
# pixel may sit this far from topg, and a wholly floating one must sit further
# above it than this.
FLOTATION_TOLERANCE_M = 1.0e-2


def ground_near_flotation(cells, tol=FLOTATION_TOLERANCE_M):
    r"""Write as grounded, in place, the floating cells whose base lies within
    ``tol`` of the bed, and return how many there were.

    The checker reads a wholly floating pixel less than 1 cm above ``topg`` as
    ice resting on the bed, an error at any share, and allows a wholly
    grounded one the same 1 cm. A DG0 cell just past flotation has its base
    millimetres above the bed: 1.9 to 9.2 mm in the full-length 32 km ssp585
    runs of September 2026, 13 to 24 pixel-years each. Only the two masks
    change; the geometry stays as the model had it. The melt booked on such a
    cell stays in ``libmassbffl`` while its pixel keeps other floating ice,
    and leaves it where the rule empties the pixel of floating ice; the
    summary line reports how much: at most 2.0 Gt/yr in the full-length 32 km
    control of September 2026 and 0.6 Gt/yr in its ssp585.
    """
    floating = cells["sftflf"] > 0.5
    near = floating & (cells["orog"] - cells["lithk"] - cells["topg"] <= tol)
    cells["sftflf"][near] = 0.0
    cells["sftgrf"][near] = 1.0
    return int(near.sum())


def scalar_area(row, sums, rtol=1.0e-6):
    r"""Which area one year's row of the scalars CSV integrates over.

    ``sums`` maps each convention, ``TRUE_AREA`` and ``MAP_PLANE``, to the
    grounded and floating ice areas the year's annual file gives under it.
    The CSV carries seven digits and the two conventions sit about 2 % apart
    over the ice sheet, so exactly one of them matches a row of the same run.
    Neither matching means the CSV belongs to another run."""
    def close(a, b):
        return abs(a - b) <= rtol * max(abs(a), abs(b))
    hits = [k for k, (gr, fl) in sums.items()
            if close(float(row["iareagr"]), gr) and close(float(row["iareafl"]), fl)]
    if len(hits) != 1:
        raise ValueError(
            f"year {row['year']}: the scalars CSV's iareagr {row['iareagr']} and iareafl "
            f"{row['iareafl']} match {'none' if not hits else 'more than one'} of the annual "
            f"file's area sums, "
            + ", ".join(f"{k} {gr:.6e} and {fl:.6e}" for k, (gr, fl) in sums.items())
            + ": the CSV and the annual files do not come from one run")
    return hits[0]


def series_area(area_of):
    r"""The one area convention of a series, from ``{year: convention}``.

    A chained run whose links straddled the change to true-area scalars
    wrote both conventions into one CSV, and the submission would carry the
    mix; it is refused."""
    kinds = sorted(set(area_of.values()))
    if len(kinds) == 1:
        return kinds[0]
    spans = {k: [y for y, v in sorted(area_of.items()) if v == k] for k in kinds}
    raise ValueError(
        "the scalars CSV mixes two areas: "
        + "; ".join(f"{k} in {len(ys)} years, {ys[0]} to {ys[-1]}" for k, ys in spans.items())
        + ". The run's links straddled the change to true-area scalars; run the "
          "series again on one version of the code.")


def series_front_melt(booking_of):
    r"""The one melt booking of a series, from ``{year: front_melt stamp}``.

    A year file written before the forward booked front melt carries no
    stamp and counts as ``UNSTAMPED_MELT``. A chained run whose links
    straddled that change would submit ``lifmassbf`` as zero in some years
    and the melt of the same cells in others, so any mix is refused."""
    kinds = sorted(set(booking_of.values()))
    if len(kinds) == 1:
        return kinds[0]
    spans = {k: [y for y, v in sorted(booking_of.items()) if v == k] for k in kinds}
    raise ValueError(
        "the annual files mix melt bookings: "
        + "; ".join(f"{k} in {len(ys)} years, {ys[0]} to {ys[-1]}" for k, ys in spans.items())
        + ". The run's links straddled the change that books front melt as "
          "lifmassbf (issue #109); run the series again on one version of the code.")


def melt_booking(W, cells, afloat_before, afloat):
    r"""One year's melt as the gridded fields carry it, in Gt/yr over the
    cells' map-plane area.

    ``left out`` is the melt booked in pixels with no floating ice at year
    end, which the ``no_floating_ice`` fill leaves out of ``libmassbffl``,
    and ``near flotation`` its part in pixels only the near-flotation rule
    emptied (``afloat_before`` and ``afloat`` are each pixel's floating area
    before and after that rule). ``front melt`` is the grid sum of
    ``lifmassbf``, which no fill touches."""
    melt_px = W @ cells["libmassbffl"]                     # m3/yr of ice per pixel
    gt = RHO_I / 1e12
    return {
        "left out": float(melt_px[afloat <= 0.0].sum()) * gt,
        "near flotation": float(melt_px[(afloat_before > 0.0) & (afloat <= 0.0)].sum()) * gt,
        "front melt": float((W @ cells["lifmassbf"]).sum()) * gt,
    }


def melt_summary(per_year, markers=MELT_MARKERS):
    r"""The summary lines for ``{year: melt_booking(...)}``: each quantity at
    its largest, then at the marker years the series holds and at its last."""
    what = {
        "left out": "libmassbffl leaves out the melt booked in pixels with no "
                    "floating ice at year end",
        "near flotation": "of it, in pixels the near-flotation rule left with no "
                          "floating ice",
        "front melt": "lifmassbf carries the front melt of cells holding no ice "
                      "at either end of the year",
    }
    years = sorted(per_year)
    shown = [y for y in markers if y in per_year]
    if years and years[-1] not in shown:
        shown.append(years[-1])
    lines = []
    for key, text in what.items():
        if not years:
            lines.append(f"  {text}: no years")
            continue
        top = max(years, key=lambda y: abs(per_year[y][key]))
        at = ", ".join(f"{per_year[y][key]:+.1f} in {y}" for y in shown)
        lines.append(f"  {text}: largest {per_year[top][key]:+.1f} Gt/yr ({top}); {at}")
    return lines


def grid_mesh():
    r"""A triangulated copy of the 8 km grid whose pixel centres are the
    ISMIP7 points, and the pixel index of every triangle."""
    half = ISMIP7_DX / 2
    Lx = ISMIP7_X0 + (ISMIP7_NX - 1) * ISMIP7_DX + half
    Ly = ISMIP7_Y0 + (ISMIP7_NY - 1) * ISMIP7_DX + half
    mesh = fd.RectangleMesh(ISMIP7_NX, ISMIP7_NY, Lx, Ly, originX=ISMIP7_X0 - half, originY=ISMIP7_Y0 - half,
                            diagonal="left", reorder=False)
    Q = fd.FunctionSpace(mesh, "DG", 0)
    xc = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro
    i = np.floor((xc[:, 0] - (ISMIP7_X0 - half)) / ISMIP7_DX).astype(int)
    j = np.floor((xc[:, 1] - (ISMIP7_Y0 - half)) / ISMIP7_DX).astype(int)
    return mesh, Q, j * ISMIP7_NX + i                     # pixel = row-major (y, x)


def overlap_operator(mesh_src, cache):
    r"""Sparse (npixels x ncells) matrix of overlap areas, cached."""
    import scipy.sparse as sp
    if cache and os.path.exists(cache):
        z = np.load(cache)
        return sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    from firedrake.supermeshing import assemble_mixed_mass_matrix
    Q_src = fd.FunctionSpace(mesh_src, "DG", 0)
    gmesh, Q_g, pix = grid_mesh()
    M = assemble_mixed_mass_matrix(Q_src, Q_g)            # rows: grid triangles, cols: model cells
    indptr, indices, data = M.getValuesCSR()
    Mt = sp.csr_matrix((data, indices, indptr), shape=(Q_g.dim(), Q_src.dim()))
    P = sp.csr_matrix((np.ones(len(pix)), (pix, np.arange(len(pix)))), shape=(ISMIP7_NX * ISMIP7_NY, len(pix)))
    W = (P @ Mt).tocsr()
    if cache:
        np.savez(cache, data=W.data, indices=W.indices, indptr=W.indptr, shape=np.array(W.shape))
    return W


def regrid(W, values, policy, mask=None, whole_pixel=False):
    r"""Pixel values under the request's fill policy; NaN where filled.

    ``whole_pixel`` divides everything the cells under a pixel carry by the
    pixel's area and keeps the policy for the fill alone: a pixel is filled
    where the policy's own mean would have no area to divide by."""
    if policy == "forbidden":
        return (W @ values) / PIXEL_AREA
    if policy == "outside_domain":
        # mean over the covered part of the pixel; any coverage counts, so
        # every pixel that carries ice (sftgif > 0) also carries elevations
        num, den = W @ values, W @ np.ones_like(values)
    else:
        # no_ice / no_grounded_ice / no_floating_ice: mean over the masked part
        m = mask.astype(float)
        num, den = (W @ values if whole_pixel else W @ (values * m)), W @ m
    out = np.full(num.shape, np.nan)
    ok = den > 0.0
    out[ok] = num[ok] / (PIXEL_AREA if whole_pixel else den[ok])
    return out


def pixel_values(W, values, meta, masks):
    r"""One variable on the 8 km grid, flat, in the request's units.

    A flux is a whole-pixel mean under any fill policy, so its value times the
    pixel area sums to the model's integral. A state variable is the mean its
    policy names: ``acabf`` and ``orog`` share ``outside_domain``, and the
    elevations stay covered-part means, which ``base := orog - lithk`` in
    ``main`` rests on."""
    policy = meta["fill_policy"]
    return (regrid(W, values, policy, masks.get(policy), whole_pixel=meta["Type"] == "FL")
            * CONVERT[meta["units"]])


GLOBAL = {}


def _global_attrs(ds, meta):
    ds.Conventions = "CF-1.8"
    for k, v in GLOBAL.items():
        setattr(ds, k, v)


def days_since_1850(year, month, day):
    import cftime
    return cftime.date2num(cftime.datetime(year, month, day, calendar="standard"),
                           "days since 1850-01-01 00:00:00", calendar="standard")


def create_2d(path, var, meta, years, is_flux):
    r"""Create the file with its axes, crs and attributes, and return the
    open dataset and the empty variable to write time slices into. The years
    are written one at a time so a 286-year submission never holds more than
    one year of grid in memory."""
    import netCDF4
    x = ISMIP7_X0 + ISMIP7_DX * np.arange(ISMIP7_NX)
    y = ISMIP7_Y0 + ISMIP7_DX * np.arange(ISMIP7_NY)
    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("time", None); ds.createDimension("y", ISMIP7_NY); ds.createDimension("x", ISMIP7_NX)
    if is_flux:
        ds.createDimension("nv", 2)
    vx = ds.createVariable("x", "f8", ("x",)); vx[:] = x
    vx.units = "m"; vx.standard_name = "projection_x_coordinate"; vx.long_name = "x coordinate of projection"; vx.axis = "X"
    vy = ds.createVariable("y", "f8", ("y",)); vy[:] = y
    vy.units = "m"; vy.standard_name = "projection_y_coordinate"; vy.long_name = "y coordinate of projection"; vy.axis = "Y"
    vt = ds.createVariable("time", "f4", ("time",))
    vt.units = TIME_UNITS; vt.calendar = "standard"; vt.standard_name = "time"; vt.long_name = "time"; vt.axis = "T"
    if is_flux:
        vt.bounds = "time_bnds"
        vb = ds.createVariable("time_bnds", "f4", ("time", "nv"))
        vt[:] = [days_since_1850(yr, 7, 1) for yr in years]
        vb[:] = [[days_since_1850(yr, 1, 1), days_since_1850(yr + 1, 1, 1)] for yr in years]
    else:
        vt[:] = [days_since_1850(yr + 1, 1, 1) for yr in years]
    crs = ds.createVariable("crs", "i4")
    crs.grid_mapping_name = "polar_stereographic"; crs.latitude_of_projection_origin = -90.0
    crs.standard_parallel = -71.0; crs.straight_vertical_longitude_from_pole = 0.0
    crs.false_easting = 0.0; crs.false_northing = 0.0; crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257223563; crs.epsg_code = "EPSG:3031"
    v = ds.createVariable(var, "f4", ("time", "y", "x"), zlib=True, complevel=4, fill_value=np.float32(FILL))
    _sn = standard_name(meta)
    if _sn:
        v.standard_name = _sn
    v.long_name = meta["long_name"]; v.units = meta["units"]
    v.grid_mapping = "crs"; v.coordinates = "y x"
    if is_flux:
        v.cell_methods = "time: mean"
    else:
        v.cell_methods = "time: point"
    _global_attrs(ds, meta)
    if is_flux:
        ds.setncattr(FLUX_MEAN_ATTR, FLUX_MEAN)
    return ds, v


def write_slice(v, k, plane):
    r"""One year of one gridded variable, filled where the policy left NaN."""
    arr = np.asarray(plane, dtype="f4").copy()
    arr[~np.isfinite(arr)] = FILL
    v[k] = arr


def write_scalar(path, var, meta, years, values, is_flux, attrs=None):
    import netCDF4
    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        for k, v in (attrs or {}).items():
            ds.setncattr(k, v)
        ds.createDimension("time", None)
        vt = ds.createVariable("time", "f4", ("time",))
        vt.units = TIME_UNITS; vt.calendar = "standard"; vt.standard_name = "time"; vt.long_name = "time"; vt.axis = "T"
        if is_flux:
            ds.createDimension("nv", 2); vt.bounds = "time_bnds"
            vb = ds.createVariable("time_bnds", "f4", ("time", "nv"))
            vt[:] = [days_since_1850(yr, 7, 1) for yr in years]
            vb[:] = [[days_since_1850(yr, 1, 1), days_since_1850(yr + 1, 1, 1)] for yr in years]
        else:
            vt[:] = [days_since_1850(yr + 1, 1, 1) for yr in years]
        v = ds.createVariable(var, "f4", ("time",), fill_value=FILL)
        arr = np.array(values, dtype="f4"); arr[~np.isfinite(arr)] = FILL
        v[:] = arr
        _sn = standard_name(meta)
        if _sn:
            v.standard_name = _sn
        v.long_name = meta["long_name"]; v.units = meta["units"]
        v.cell_methods = "time: mean" if is_flux else "time: point"
        _global_attrs(ds, meta)


# The core set's counters (the ISMIP7 filenames and conventions document;
# discussion #17). The counter names the experiment in every filename and in
# the directory, it was typed by hand, and C007 for a run that is C008 is a
# well-formed submission of the wrong experiment that no checker can catch.
CORE_COUNTER = {
    ("CESM2-WACCM", "historical"): "C001", ("MRI-ESM2-0", "historical"): "C002",
    ("CESM2-WACCM", "ssp370"): "C003", ("MRI-ESM2-0", "ssp370"): "C004",
    ("CESM2-WACCM", "ssp126"): "C005", ("MRI-ESM2-0", "ssp126"): "C006",
    ("CESM2-WACCM", "ssp585"): "C007", ("MRI-ESM2-0", "ssp585"): "C008",
    ("CESM2-WACCM", "ctrl"): "C009", ("MRI-ESM2-0", "ctrl"): "C010",
}
# OCX has no ESM, and what belongs in the filename's forcing field for it is
# not settled upstream: isschecker 0.5.1 checks that field against a list of
# CMIP models and its experiment table has no ocx row. So the counter follows
# from the scenario alone.
OCX_COUNTER = "C011"


def set_counter(esm, scenario, set_id, given=None):
    r"""The ``<set_counter>`` of a submission, e.g. ``C007``.

    For the core set it follows from the forcing, so ``given`` is optional
    and refused when it names another core experiment. Outside the core set
    (``ESM``, ``PPE``) nothing here knows the numbering, and ``given`` is
    required.
    """
    if given is not None and not re.fullmatch(r"[CEP]\d{3}", given):
        raise ValueError(f"--exp must look like C007, E041 or P132, got {given!r}")
    if set_id != "CORE":
        if given is None:
            raise ValueError(f"--exp is required for the {set_id} set: only the core set's "
                             f"counters follow from the forcing")
        return given
    expected = OCX_COUNTER if scenario.lower() == "ocx" else CORE_COUNTER.get((esm, scenario))
    if expected is None:
        raise ValueError(f"{esm} {scenario} is not a core experiment "
                         f"({', '.join(f'{e} {sc}' for e, sc in CORE_COUNTER)}, ocx); "
                         f"pass --set-id ESM or PPE with its --exp")
    if given is not None and given != expected:
        raise ValueError(f"--exp {given} names another experiment: {esm} {scenario} "
                         f"is {expected} in the core set")
    return expected


def submission_id(name, value):
    r"""``source_id`` and ``ism_id`` go into the filename between underscores,
    and the conventions allow neither an underscore nor a dot in them
    (discussion #17: ``ISSMv2026p2``, not ``ISSM_2026.2``). An underscore
    shifts every field after it, so the file parses as another run."""
    if not re.fullmatch(r"[A-Za-z0-9-]+", value):
        raise ValueError(f"--{name} {value!r}: letters, digits and hyphens only, "
                         f"no underscore, dot or space (ISMIP7 filename conventions)")
    return value


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("annual"); ap.add_argument("--out-dir", required=True)
    ap.add_argument("--esm", required=True); ap.add_argument("--scenario", required=True)
    ap.add_argument("--exp", default=None,
                    help="the set counter, e.g. C007. The core set's follows from --esm and "
                         "--scenario, so it is optional there and checked when given")
    ap.add_argument("--source-id", default=os.environ.get("ISMIP7_SOURCE_ID", "RICE"))
    ap.add_argument("--ism-id", default=os.environ.get("ISMIP7_ISM_ID", "icepack2"))
    ap.add_argument("--set-id", default="CORE")
    ap.add_argument("--contact-name", default=os.environ.get("ISMIP7_CONTACT_NAME", "Andrew Hoffman"))
    ap.add_argument("--contact-email", default=os.environ.get("ISMIP7_CONTACT_EMAIL", "ah301@rice.edu"))
    ap.add_argument("--scalars", default=None, help="the *_ismip7_scalars.csv (default: next to the annual file)")
    a = ap.parse_args()
    try:
        a.exp = set_counter(a.esm, a.scenario, a.set_id, a.exp)
        submission_id("source-id", a.source_id)
        submission_id("ism-id", a.ism_id)
    except ValueError as e:
        ap.error(str(e))
    req = request_table()
    years = AnnualOutput.years_on_disk(a.annual)
    if not years:
        hint = ("; pass the stem the run was named after."
                if not os.path.exists(a.annual) else
                ". A current forward writes each year to its own "
                "<stem>_<year>.h5 and carries the partial year in the run's "
                "own checkpoint, so the stem path stays a name. A stem file "
                "with no yearly files beside it comes from a forward from "
                "before 14 September 2026, which banked the years INSIDE that "
                "file, a layout this writer does not read: re-run the forward "
                "with ISMIP7_OUTPUT=1 on current code.")
        raise FileNotFoundError(
            f"no yearly checkpoints {os.path.basename(a.annual)[:-3]}_<year>.h5 "
            f"beside {a.annual}{hint}"
        )
    if years != list(range(years[0], years[-1] + 1)):
        missing = sorted(set(range(years[0], years[-1] + 1)) - set(years))
        raise FileNotFoundError(
            f"the yearly checkpoints beside {a.annual} run {years[0]}-{years[-1]} "
            f"with {len(missing)} missing ({missing[0]}..{missing[-1]}); a "
            f"submission needs a contiguous series."
        )
    with fd.CheckpointFile(AnnualOutput.year_path(a.annual, years[0]), "r") as chk:
        mesh = chk.load_mesh()
    print(f"{os.path.basename(a.annual)}: years {years[0]}-{years[-1]}, "
          f"{len(VARIABLES_2D)} variables", flush=True)
    W = overlap_operator(mesh, a.annual + ".overlap.npz")
    print(f"  overlap operator {W.shape}, {W.nnz} entries, pixels covered {int((W.sum(axis=1) > 0).sum())}", flush=True)
    # each cell's area and af2 at its centroid, in the operator's column
    # order, to tell which area the scalars CSV integrates over
    cell_area = np.asarray(W.sum(axis=0)).ravel()
    X = fd.SpatialCoordinate(mesh)
    cell_af2 = area_factor(*(fd.Function(fd.FunctionSpace(mesh, "DG", 0)).interpolate(X[i]).dat.data_ro
                             for i in (0, 1)))

    scal = a.scalars or a.annual.replace("_ismip7_annual.h5", "_ismip7_scalars.csv")
    if not os.path.exists(scal):
        raise FileNotFoundError(
            f"no scalars csv at {scal}; all ten scalar variables are "
            f"mandatory, so the submission cannot be written without it. "
            f"Pass --scalars if the run's csv is named differently."
        )
    with open(scal) as f:
        rows = {int(r["year"]): r for r in csv.DictReader(f)}
    missing_rows = [yr for yr in years if yr not in rows]
    if missing_rows:
        raise FileNotFoundError(
            f"{scal} has no row for year(s) {missing_rows[0]}"
            + (f"..{missing_rows[-1]}" if len(missing_rows) > 1 else "")
            + " although the gridded output has them; the scalar series "
              "would be submitted with a hole in it."
        )

    outdir = os.path.join(a.out_dir, "AIS", a.source_id, a.ism_id, a.set_id, a.exp)
    os.makedirs(outdir, exist_ok=True)
    tag = f"AIS_{a.source_id}_{a.ism_id}_m001_{a.esm}_f001_{a.scenario}_{a.exp}_{years[0]}-{years[-1]}"
    GLOBAL.update({"model": a.ism_id, "group": a.source_id, "crs": "EPSG:3031",
                   "contact_name": a.contact_name,
                   "contact_email": a.contact_email,
                   "source_id": a.source_id, "ism_id": a.ism_id, "experiment_id": a.exp,
                   "forcing": f"{a.esm} {a.scenario}",
                   "title": f"ISMIP7 AIS {a.exp} {a.esm} {a.scenario}, {a.source_id} {a.ism_id}"})

    # One open file per gridded variable, one year of grid in memory at a
    # time: the cubes for a 286-year submission on the production mesh would
    # otherwise be tens of GB in a serial script. Everything is written to
    # <name>.nc.tmp and renamed into place only once the last year and the
    # scalars are through, so a failure part-way leaves no half-filled file
    # that looks like a finished submission, and any previous files stand.
    paths = {var: os.path.join(outdir, f"{var}_{tag}.nc") for var in VARIABLES_2D}
    for var in SCALARS:
        paths[var] = os.path.join(outdir, f"{var}_{tag}.nc")
    tmps = {var: path + ".tmp" for var, path in paths.items()}
    handles = {}
    try:
        for var in VARIABLES_2D:
            handles[var] = create_2d(tmps[var], var, req[var], years,
                                     req[var]["Type"] == "FL")
        stats = {var: [np.inf, -np.inf, 0] for var in VARIABLES_2D}
        grounded_near = 0
        # each year's melt as the gridded fields carry it (melt_booking)
        melt_of = {}
        # each year's front_melt stamp, held to one booking per series
        booking_of = {}
        area_of = {}
        for k, yr in enumerate(years):
            with fd.CheckpointFile(AnnualOutput.year_path(a.annual, yr), "r") as chk:
                ymesh = chk.load_mesh()
                cells = {var: chk.load_function(ymesh, name=var).dat.data_ro.copy()
                         for var in VARIABLES_BANKED}
                booking_of[yr] = (str(chk.get_attr("/", FRONT_MELT_ATTR))
                                  if chk.has_attr("/", FRONT_MELT_ATTR) else UNSTAMPED_MELT)
            series_front_melt(booking_of)
            # the overlap operator is built once from the first year's mesh;
            # a series whose links remeshed cannot be regridded through it
            if len(cells["lithk"]) != W.shape[1]:
                raise ValueError(
                    f"year {yr} has {len(cells['lithk'])} cells but the overlap "
                    f"operator was built from year {years[0]}'s {W.shape[1]}: "
                    f"the mesh changed inside the series, so one conservative "
                    f"operator cannot cover it."
                )
            # the forward's own masks, before the near-flotation rule below
            gr, fl = cells["sftgrf"] * cell_area, cells["sftflf"] * cell_area
            area_of[yr] = scalar_area(rows[yr], {
                TRUE_AREA: (float((gr * cell_af2).sum()), float((fl * cell_af2).sum())),
                MAP_PLANE: (float(gr.sum()), float(fl.sum()))})
            afloat_before = W @ (cells["sftflf"] > 0.5).astype(float)
            grounded_near += ground_near_flotation(cells)
            masks = {"no_ice": cells["sftgif"] > 0.5,
                     "no_grounded_ice": cells["sftgrf"] > 0.5,
                     "no_floating_ice": cells["sftflf"] > 0.5}
            afloat = W @ masks["no_floating_ice"].astype(float)
            melt_of[yr] = melt_booking(W, cells, afloat_before, afloat)

            def plane(var):
                return pixel_values(W, cells[var], req[var], masks).reshape(
                    ISMIP7_NY, ISMIP7_NX).astype("f4")

            # the checker requires orog == base + lithk pixel by pixel and
            # orog >= 0; the elevations are covered-part means while lithk is a
            # whole-pixel mean (diluted where the pixel is partly covered), so
            # the base absorbs the dilution: base := orog - lithk on the grid.
            # Rebuilding orog instead put negative surfaces on partly covered
            # floating pixels (base < 0 plus a diluted thickness).
            done = {"orog": plane("orog"), "lithk": plane("lithk")}
            done["base"] = done["orog"] - done["lithk"]
            for var in VARIABLES_2D:
                p_yr = done.get(var)
                if p_yr is None:
                    p_yr = plane(var)
                write_slice(handles[var][1], k, p_yr)
                finite = np.isfinite(p_yr)
                st = stats[var]
                if finite.any():
                    st[0] = min(st[0], float(np.min(p_yr[finite])))
                    st[1] = max(st[1], float(np.max(p_yr[finite])))
                st[2] += int(finite.sum())
        print(f"  {grounded_near} cell-years within {FLOTATION_TOLERANCE_M:g} m of "
              f"flotation written as grounded", flush=True)
        print(f"  melt booking: {series_front_melt(booking_of)} in every year", flush=True)
        for line in melt_summary(melt_of):
            print(line, flush=True)
        for var in VARIABLES_2D:
            handles.pop(var)[0].close()
            lo, hi, nfin = stats[var]
            pct = 100.0 * nfin / (len(years) * ISMIP7_NY * ISMIP7_NX)
            rng = f"[{lo:.3e}, {hi:.3e}]" if nfin else "[all fill]"
            print(f"  {var:12s} {req[var]['units']:11s} {rng} "
                  f"{pct:5.1f}% valid  -> {os.path.basename(paths[var])}", flush=True)

        area = series_area(area_of)
        for var in SCALARS:
            meta = req[var]
            values = [float(rows[yr][var]) for yr in years]
            write_scalar(tmps[var], var, meta, years, values, meta["Type"] == "FL",
                         {SCALAR_AREA_ATTR: area})
        print(f"  {len(SCALARS)} scalars from {os.path.basename(scal)}, over {area.replace('_', ' ')} "
              f"in every year (iareagr and iareafl match the annual files' sums)", flush=True)

        for var, tmp in tmps.items():
            os.replace(tmp, paths[var])
    except BaseException:
        for ds, _ in handles.values():
            ds.close()
        for tmp in tmps.values():
            if os.path.exists(tmp):
                os.remove(tmp)
        raise
    print(f"wrote {outdir}")


if __name__ == "__main__":
    main()
