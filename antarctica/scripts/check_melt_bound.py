#!/usr/bin/env python
r"""Does the forward apply the melt its calibration was fitted to, and how
far past the variable request's bound does that melt reach?

    python antarctica/scripts/check_melt_bound.py [--npz <calibration.npz>]
        [--match-tol 1e-3]
        [--ocx [main|cold|warm|vary]] [--ocx-years 2000,2015,2025] [--ocx-tol 0.25]

Without ``--npz`` the calibration is the one a forward reads
(``runconfig.deltat_per_basin_npz``): the tracked
``antarctica/calibration/deltaT_per_basin_1000_K6.500e-05.npz`` unless
another is named. ``--npz`` takes an offsets file or a legacy per-basin K.

The forward half is the forward's own melt at the reference geometry: bed and
thickness sampled onto the DG0 cells, the surface from flotation,
s = max(b + H, (1 - 917/1024) H), as ``simulation.setup_model`` builds them on
a cold start and on a target mesh, melted by the forward's own
``forcing.make_climatology_ocean_callback`` with the OI climatology. With an
offsets file its melt is summed per IMBIE2 basin and set against the totals
the offsets were fitted to (the file's ``M_obs``), and the check exits 1 when
any basin is off by more than ``--match-tol``. That is the measurement issue
30 asked for: a calibration fitted through the forward's melt path, and the
forward's melt on the same cells. On a mesh other than the one the offsets
were fitted on, the table measures how far they carry over: at 32 km the
K50 offsets put the basins at 0.33 to 1.74 times their totals (job 10644430),
and the fit on IU's build of the production mesh misses basins 1 and 7 on
Rice's build by 0.18 and 0.17 percent (job 10649416). The table also gives what the forward booked on ice-free cells
that pass the flotation test, open ocean and bare land apart, before it
melted only cells holding ice (``forcing.melt_receiving``).

The ISMIP7 variable request gives ``libmassbffl`` an AIS minimum of
-0.008 kg m-2 s-1 with severity ``error``. In ice-equivalent thickness that is
275.3 m/yr, and a 10-year adaptive-mesh ssp585 reached -0.0117 (402.6 m/yr) on
grounding-zone cells. Two readings fit that: the parameterisation is too strong
somewhere, or one hot cell sets the value of its whole 8 km pixel. The writer's
flux means are whole-pixel means (issue #96), which weight a cell by its share
of the pixel, so a hot cell smaller than its pixel reaches the grid diluted.
Each row below reports the maximum, the 99th percentile, the area mean, the
integrated total against the total the calibration was fitted to, and how
much floating AREA sits past the bound.

With a legacy per-basin K file the check also rebuilds the calibration half,
calibrate_melt's CG1 fit: BedMachine interpolated with its raster ``surface``
and ``mask``, and grad(draft) projected onto CG1.

* capped: sin_alpha capped at the npz's ``sin_alpha_cap`` when it is finite,
  otherwise at calibrate_melt's CG1 default 5e-3, the slope K was fitted
  against;
* uncapped: the same slope with no cap.

and, under ``local``, the forward half with its cell slope capped at the same
value, which is where a cap inside ``forcing.compute_sin_alpha`` would act.

Both halves take their slope from ``ISMIP7_MELT_SLOPE`` as the forward and
calibrate_melt do. The capped and uncapped rows above belong to ``local``.
Under ``ant``, the default, the forward half melts with the one constant
``ISMIP7_SIN_ALPHA_ANT`` and the calibration half with the constant the K file
records, warning when the two differ by more than 1 percent; no cap applies,
and each half gives one row. A K file
whose recorded ``melt_slope`` differs from the run's is reported, since its K
does not transfer.

The two halves use different floating masks and different quadrature, nodal
area weights against cell areas, so their totals compare in magnitude and
differ in detail.

The bound is ``min_value_ais`` for ``libmassbffl`` in the same bundled request
table the writer reads, converted with the writer's year and ice density, so it
is the bound the compliance checker applies.

Scope: the reference state, with the OI thermal-forcing climatology and the
BedMachine geometry. A projection's thermal forcing warms above the
climatology and its shelves thin.

``--ocx`` adds the check core 11 needs before it runs on the ISMIP7 OCX
product. Every K here is fitted to the OI climatology, and the OCX ocean is a
different field: discussion #48 (17 September 2026) reports the OCX ``main``
thermal forcing so far from the Zhou climatology around Mertz that the
region's melt halves against a calibration made on the climatology, possibly
because OCX was built from an older extrapolated climatology. So the forward
half is melted a second time, uncapped as the forward runs, with thermal
forcing and salinity from the OCX ocean through the forward's own reader, for
each of ``--ocx-years``, and the two melts are set side by side per IMBIE2
basin and per 256 km block, the blocks because a basin total dilutes one
shelf. A basin or block whose OCX melt is off the climatology's by more than
``--ocx-tol`` is flagged and the exit status is 1, so a chain can stop on it.

Serial. Reuses calibrate_melt's loaders, so it needs the same inputs: a MAP for
the mesh (ISMIP7_INV_H5), the OI climatology, the IMBIE2 basins and BedMachine.

The measured rows below are the legacy per-basin K files that issue 26
settled. They predate the ``h > 0`` test and the seawater flotation test in
the forward half (icepack/ismip7#66): they counted ice-free cells as floating
and grounded the deep-draft shelf.

Measured on the adaptive 2 km mesh, September 2026, with
calibrated_K_per_basin_2000.npz calibrated against the re-released observation
table:

* calibration half, capped at 5e-3, over 1 512 899 km2 and 47 288 nodes:
  maximum 71.1 m/yr, 99th percentile 22.2, area mean 0.77, 1067 Gt/yr, which
  is the 1067.4 Gt/yr target K was fitted to, and nothing past the bound;
* calibration half, uncapped: maximum 1804.9 m/yr, 99th percentile 256.3, area
  mean 4.18, 5803 Gt/yr, and 421 nodes past the bound over 2920.6 km2 with a
  median node area of 6.33 km2;
* forward half, uncapped, over 1 631 466 km2 and 85 820 cells: maximum
  1144.1 m/yr, 99th percentile 59.0, area mean 1.16, 1732 Gt/yr, and 97 cells
  past the bound over 388.2 km2 with a median cell area of 3.61 km2;
* forward half, capped at 5e-3: maximum 57.8 m/yr, 99th percentile 13.0, area
  mean 0.43, 646 Gt/yr, and nothing past the bound.

With calibrated_K_per_basin_2500.npz, the coefficient file the 10-year ssp585
run read, everything else unchanged:

* calibration half, capped: maximum 54.0 m/yr, area mean 0.61, 841 Gt/yr, and
  nothing past the bound;
* calibration half, uncapped: maximum 1522.7 m/yr, area mean 3.34, 4634 Gt/yr,
  and 320 nodes past the bound over 2193.7 km2;
* forward half, uncapped: maximum 869.8 m/yr, area mean 0.92, 1380 Gt/yr, and
  63 cells past the bound over 248.6 km2 with a median cell area of 3.54 km2;
* forward half, capped: maximum 43.9 m/yr, area mean 0.34, 510 Gt/yr, and
  nothing past the bound.

The calibration half reproduces the target its own K was fitted to, 841
against 865 Gt/yr for the 2500 file, which is the internal consistency check.
The ratio of forward to calibration is the durable result, stable across both
coefficient files: 1.62 with the 2000 file and 1.64 with the 2500 file. The
10-year run booked 1860 Gt/yr with the 2500 file, above the 1380 Gt/yr its
forward half applies at the reference state, since that run carries warmer
ssp585 thermal forcing over evolving geometry; the gap is a consistent
residual. An earlier claim that the forward half lands within 7% of that run
compared two different coefficient files and does not hold.

Capping the forward's slope undershoots the target by 39% with the 2000 file,
so neither convention on its own reconciles the halves, which also differ in
floating area, mask and quadrature. Recalibrating through the forward's own
cell by cell melt path reconciles them by construction, which is how the
tracked calibration was fitted (issue 26).

An earlier form of this script lifted the forward's slope onto CG1 nodes and
melted it with CG1 forcing and the raster mask, so its forward rows,
4293 Gt/yr uncapped and 1028 Gt/yr capped, reproduced neither half and are
superseded.
"""
import argparse
import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))
sys.path.insert(0, _ROOT)

import numpy as np                                                    # noqa: E402
import firedrake as fd                                                # noqa: E402
from firedrake import assemble, dx                                    # noqa: E402
from firedrake.petsc import PETSc                                     # noqa: E402

import calibrate_melt as cm                                           # noqa: E402
from icepack2_tools.forcing import (quadratic_mixed_slope,            # noqa: E402
                                    compute_sin_alpha, is_floating,
                                    melt_receiving, melt_slope,
                                    sin_alpha_ant, ISMIP7Ocean, OCX,
                                    OCX_OCEAN_VARIANTS,
                                    describe_forcing_provenance,
                                    describe_melt_calibration,
                                    imbie2_basin_path, load_K_per_basin,
                                    load_deltaT_per_basin,
                                    make_climatology_ocean_callback,
                                    _basin_on_mesh, _RHO_I)
from icepack2_tools.mpi_stats import global_size                      # noqa: E402
from icepack2_tools.runconfig import (deltat_per_basin_npz,           # noqa: E402
                                      k_per_basin_npz,
                                      melt_calibration_contract,
                                      raster_sample)
# The same year and density the writer converts with, so the bound compared
# here is the one the checker applies.
from icepack2_tools.ismip7_output import RHO_I, SECONDS_PER_YEAR      # noqa: E402
from icepack2_tools.regrid import ISMIP7_DX                           # noqa: E402
from write_ismip7_output import request_table                         # noqa: E402

# The ice to seawater density ratio simulation.py builds the surface with.
RHO_RATIO = 917.0 / 1024.0


def m_per_yr(kg_m2_s):
    r"""kg m-2 s-1 of ice to m/yr of ice thickness."""
    return kg_m2_s * SECONDS_PER_YEAR / RHO_I


# Side of the blocks the OCX comparison also sums over: 32 pixels of the 8 km
# forcing grid, a few shelves wide, so one shelf is not lost in a basin total.
BLOCK_M = 256.0e3


def block_ids(x, y, size=BLOCK_M):
    r"""An integer id per point naming the ``size`` square it falls in, and
    the ``{id: (x centre, y centre)}`` of the squares that occur."""
    ix, iy = np.floor(x / size).astype(int), np.floor(y / size).astype(int)
    ids = ix * 100000 + iy
    centres = {int(i): ((a + 0.5) * size, (b + 0.5) * size)
               for i, a, b in zip(ids, ix, iy)}
    return ids, centres


def melt_by_group(groups, melt, area, floating):
    r"""``{group: Gt/yr}`` of ``melt`` (m/yr of ice) over the floating dofs."""
    gt = np.where(floating, melt * area, 0.0) * RHO_I / 1e12
    return {int(g): float(gt[groups == g].sum()) for g in np.unique(groups[floating])}


def off_by_more_than(reference, other, tol, floor_gt):
    r"""The groups whose ``other`` melt is off ``reference`` by more than
    ``tol`` (a fraction). Groups melting less than ``floor_gt`` either way are
    left out: a ratio of two near-zero totals flags nothing worth reading."""
    flagged = []
    for g, ref in reference.items():
        new = other.get(g, 0.0)
        if max(ref, new) < floor_gt:
            continue
        if ref <= 0.0 or abs(new / ref - 1.0) > tol:
            flagged.append(g)
    return flagged


def load_calibration(npz_arg):
    r"""The melt calibration to check and its kind, ``"deltaT"`` or ``"K"``.

    ``--npz`` names a file of either kind, told apart by its keys; without it
    this is the calibration a forward reads (runconfig.deltat_per_basin_npz:
    the tracked file unless another is named, else a legacy per-basin K named
    with ISMIP7_K_PER_BASIN_NPZ). A named file is put in the variable a
    forward reads it from, so the forward's own callback melts with it."""
    if npz_arg:
        with np.load(npz_arg, allow_pickle=True) as data:
            knob = ("ISMIP7_DELTAT_PER_BASIN_NPZ" if "deltaT_basin" in data.files
                    else "ISMIP7_K_PER_BASIN_NPZ")
        for name in ("ISMIP7_DELTAT_PER_BASIN_NPZ", "ISMIP7_K_PER_BASIN_NPZ"):
            os.environ.pop(name, None)
        os.environ[knob] = npz_arg
    path = deltat_per_basin_npz()
    if path is not None:
        return path, "deltaT"
    return k_per_basin_npz(), "K"


def forward_half(mesh, npz_path, kind):
    r"""The forward's melt at the reference geometry, cell by cell.

    The geometry is ``calibrate_melt.forward_cells``, the cells the
    calibrations fit on, built as ``simulation.setup_model`` builds a cold
    start and a target mesh: bed and thickness sampled onto the DG0 cells
    (``ISMIP7_RASTER_SAMPLE``), the surface from flotation. The melt is the
    forward's own, ``forcing.make_climatology_ocean_callback`` on that
    geometry, so it is the field a forward writes into ``ctx["ocean_melt"]``
    with the OI climatology the calibrations are fitted against.

    ``melt_ice_free`` is what the callback booked on ice-free cells that pass
    the flotation test before it melted only cells holding ice: the same
    callback run with a millimetre of ice on those cells, whose draft stays at
    the surface. ``bed`` tells their open ocean from their bare land."""
    c = cm.forward_cells(mesh)
    Q, Q_g, V, x, y = c["Q"], c["Q_g"], c["V"], c["x"], c["y"]
    b_dg, h_dg, s_dg = c["b"], c["h"], c["s"]
    K_field = load_K_per_basin(npz_path, x, y, fill=0.0) if kind == "K" else None
    # One callback, called twice: it keeps its climatology and offsets.
    callback = make_climatology_ocean_callback(K_field)

    def melt_on(h, s):
        ctx = {"mesh": mesh, "Q": Q, "V": V, "Q_g": Q_g, "geom_xy": (x, y),
               "h": h, "b": b_dg, "s": s, "ocean_melt": fd.Function(Q_g),
               "raster_sample": raster_sample()}
        callback(ctx, 0.0)
        return ctx["ocean_melt"].dat.data_ro.copy()

    PETSc.Sys.Print("  the forward's callback on the reference geometry:")
    melt = melt_on(h_dg, s_dg)
    b_np, h_np, s_np = b_dg.dat.data_ro, h_dg.dat.data_ro, s_dg.dat.data_ro
    ice_free = is_floating(s_np, b_np) & ~(h_np > 0)
    h_trace = fd.Function(Q_g).assign(h_dg)
    h_trace.dat.data[ice_free] = 1e-3
    s_trace = fd.Function(Q_g).interpolate(
        fd.max_value(b_dg + h_trace, (1.0 - RHO_RATIO) * h_trace))
    melt_ice_free = np.where(ice_free, melt_on(h_trace, s_trace), 0.0)

    draft = np.minimum(s_np - h_np, 0.0)
    if kind == "deltaT":
        dT, K = load_deltaT_per_basin(npz_path, x, y, fill=0.0)
        K = np.full(len(x), K)
    else:
        dT, K = np.zeros(len(x)), K_field
    return {"x": x, "y": y, "draft": draft, "bed": b_np.copy(),
            "sin_a": compute_sin_alpha({"Q": Q, "V": V, "Q_g": Q_g,
                                        "h": h_dg, "s": s_dg}),
            "floating": melt_receiving(s_np, b_np, h_np),
            "area": assemble(fd.TestFunction(Q_g) * dx).dat.data_ro,
            "tf": cm._grid_interp(cm.CLIM_TF, "tf", x, y, draft=draft),
            "sal": cm._grid_interp(cm.CLIM_SO, "so", x, y, draft=draft),
            "K": K, "dT": dT, "melt": melt,
            "ice_free": ice_free, "melt_ice_free": melt_ice_free}


def match_calibration(g, npz_path, tol):
    r"""The forward's melt per IMBIE2 basin against the totals the offsets
    were fitted to (the file's ``M_obs``), stamped through the basins the
    forward stamps the offsets through. Prints the table and returns the
    basins off by more than ``tol``, a fraction."""
    with np.load(npz_path) as data:
        bids = np.asarray(data["basin_ids"]).astype(int)
        fitted = np.asarray(data["M_obs"]).astype(float)
        recorded = str(data["imbie2_nc"]) if "imbie2_nc" in data else None
    basin = _basin_on_mesh(g["x"], g["y"], imbie2_basin_path(recorded))
    gt = g["melt"] * g["area"] * _RHO_I / 1e12
    free = g["melt_ice_free"] * g["area"] * _RHO_I / 1e12
    PETSc.Sys.Print(
        f"\n=== the forward's melt against its calibration, per IMBIE2 basin "
        f"(tolerance {100 * tol:g} percent) ===\n"
        f"  basin   fitted Gt/yr   forward Gt/yr    ratio")
    bad = []
    for bid, want in zip(bids, fitted):
        got = float(gt[basin == bid].sum())
        ratio = got / want if want > 0 else float("inf")
        off = not abs(ratio - 1.0) <= tol
        if off:
            bad.append(int(bid))
        PETSc.Sys.Print(f"  {bid:5d}   {want:12.3f}   {got:13.3f}   {ratio:7.4f}"
                        + ("   <-- OFF" if off else ""))
    inside = np.isin(basin, bids)
    PETSc.Sys.Print(
        f"  total   {fitted.sum():12.3f}   {float(gt[inside].sum()):13.3f}\n"
        f"  outside the fitted basins: {float(gt[~inside].sum()):.3f} Gt/yr\n"
        f"  ice-free cells that pass the flotation test, which the forward "
        f"melted before it melted only cells holding ice, none of it booked "
        f"now:")
    for name, where in (("open ocean (bed below sea level)", g["bed"] < 0.0),
                        ("bare land (bed at or above sea level)", g["bed"] >= 0.0)):
        cells = g["ice_free"] & where
        part = np.where(cells, free, 0.0)
        refreezing = float(-part[part < 0].sum())
        PETSc.Sys.Print(
            f"    {name}: {int(cells.sum())} cells, melt "
            f"{float(part[part > 0].sum()):.1f} Gt/yr, refreezing "
            f"{refreezing if refreezing > 0.0 else 0.0:.1f} Gt/yr")
    return bad


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", default=None,
                    help="a melt calibration, offsets (calibrate_deltaT.py, "
                         "select_melt_parameters.py) or a legacy per-basin K "
                         "(calibrate_melt.py); default: the one a forward reads")
    ap.add_argument("--match-tol", type=float, default=1e-3,
                    help="flag a basin whose forward melt is off the total its "
                         "offsets were fitted to by more than this fraction "
                         "(default %(default)s); exit 1 if any is")
    ap.add_argument("--ocx", nargs="?", const="main", default=None,
                    choices=OCX_OCEAN_VARIANTS,
                    help="also melt the forward half with this OCX ocean and "
                         "compare it with the climatology K was fitted to "
                         "(discussion #48); exit 1 if any region is off")
    ap.add_argument("--ocx-years", default="2000,2015,2025",
                    help="comma list of OCX years to compare (default %(default)s)")
    ap.add_argument("--ocx-tol", type=float, default=0.25,
                    help="flag a region whose OCX melt is off by more than "
                         "this fraction (default %(default)s)")
    ap.add_argument("--ocx-floor", type=float, default=2.0,
                    help="ignore regions melting less than this many Gt/yr "
                         "either way (default %(default)s)")
    a = ap.parse_args()

    bound = float(request_table()["libmassbffl"]["min_value_ais"])
    bound_m_yr = abs(m_per_yr(bound))
    PETSc.Sys.Print(f"=== melt against the request bound {bound} kg m-2 s-1 "
                    f"({bound_m_yr:.1f} m/yr ice) ===")

    npz_path, kind = load_calibration(a.npz)
    d = np.load(npz_path, allow_pickle=(kind == "K"))
    for line in describe_melt_calibration(npz_path if kind == "deltaT" else None,
                                          npz_path if kind == "K" else None):
        PETSc.Sys.Print(f"  {line}")
    if "obs_csv" in d:
        PETSc.Sys.Print(f"  calibrated against {d['obs_csv']}")
    if kind == "deltaT" and "M_obs" in d:
        target = (f" (the offsets were fitted to "
                  f"{float(np.sum(d['M_obs'])):.1f})")
    elif "obs_total_gtyr" in d:
        target = f" (K was fitted to {float(d['obs_total_gtyr']):.0f})"
    else:
        target = ""
        PETSc.Sys.Print(f"  {os.path.basename(npz_path)} records no observed "
                        f"total, so the rows print no target")

    mesh = cm._load_mesh()
    if mesh.comm.size > 1:
        raise SystemExit("check_melt_bound.py is serial: its tables sum the "
                         "dofs of one process")
    # num_vertices()/num_cells() count this rank's plex, halo included; the
    # coordinate dofs and the owned cell set are reduced to global totals.
    PETSc.Sys.Print(f"  Mesh: {global_size(mesh.coordinates)} vertices, "
                    f"{mesh.comm.allreduce(mesh.cell_set.size)} cells")

    # The forward half: the forward's own melt, cell by cell under DG0.
    forward = forward_half(mesh, npz_path, kind)
    running = melt_slope()
    as_run = (f"forward half, the constant {sin_alpha_ant():g} on DG0 cells, "
              f"as the forward runs" if running == "ant" else
              "forward half, local, uncapped, as the forward runs")
    cases = [(as_run, forward, forward["sin_a"], "cells", forward["melt"])]

    if kind == "K":
        # A legacy per-basin K was fitted on CG1 nodes: the calibration half
        # rebuilds that fit, BedMachine interpolated onto CG1 nodes with its
        # raster surface and mask, and the slope under this run's convention.
        def k_at(xs, ys):
            r"""The per-basin K the forward stamps onto the mesh, rebuilt from
            K_basin so this runs against any mesh."""
            basin = np.round(
                cm._grid_interp(cm.IMBIE2_NC, "basinNumber", xs, ys)).astype(int)
            K = np.zeros(len(xs))
            for bid, kb in zip(d["basin_ids"], d["K_basin"]):
                if np.isfinite(kb):
                    K[basin == int(bid)] = kb
            return K

        c = cm.calibration_geometry(mesh)
        calibration = {"x": c["x"], "y": c["y"], "draft": c["draft"],
                       "sin_a": c["sin_a"], "floating": c["floating"],
                       "area": c["area"], "K": k_at(c["x"], c["y"]),
                       "tf": cm._grid_interp(cm.CLIM_TF, "tf", c["x"], c["y"],
                                             draft=c["draft"]),
                       "sal": cm._grid_interp(cm.CLIM_SO, "so", c["x"], c["y"],
                                              draft=c["draft"])}
        fitted_slope = str(d["melt_slope"]) if "melt_slope" in d else "local"
        if fitted_slope != running:
            PETSc.Sys.Print(f"  WARNING: {os.path.basename(npz_path)} was calibrated "
                            f"under ISMIP7_MELT_SLOPE={fitted_slope} and this check "
                            f"melts under {running}, so no row is the slope K was "
                            f"fitted against")
        if running == "ant":
            fitted_sin = (float(d["sin_alpha_ant"])
                          if fitted_slope == "ant" and "sin_alpha_ant" in d
                          else float("nan"))
            calibration_sin = fitted_sin if np.isfinite(fitted_sin) else sin_alpha_ant()
            if abs(calibration_sin / sin_alpha_ant() - 1.0) > 0.01:
                PETSc.Sys.Print(f"  WARNING: {os.path.basename(npz_path)} was "
                                f"calibrated with sin(alpha) = {calibration_sin:g} "
                                f"and this check melts the forward half with "
                                f"{sin_alpha_ant():g}; the calibration half keeps "
                                f"the file's constant")
            fitted = ", the slope K was fitted against" if fitted_slope == "ant" else ""
            cases.insert(0, (f"calibration half, the constant {calibration_sin:g} "
                             f"on CG1 nodes{fitted}", calibration,
                             np.full_like(calibration["sin_a"], calibration_sin),
                             "nodes", None))
        else:
            cap = float(d["sin_alpha_cap"]) if "sin_alpha_cap" in d else float("nan")
            if not (np.isfinite(cap) and cap > 0):
                cap = cm.default_slope_cap("cg1")
            fitted_on = str(d["geometry_space"]) if "geometry_space" in d else "cg1"
            cases[:0] = [
                (f"calibration half, local, capped at {cap:.0e} on CG1 nodes"
                 + (", the slope K was fitted against"
                    if fitted_slope == "local" and fitted_on == "cg1" else ""),
                 calibration, np.minimum(calibration["sin_a"], cap), "nodes", None),
                ("calibration half, local, uncapped", calibration,
                 calibration["sin_a"], "nodes", None),
            ]
            # A cap inside forcing.compute_sin_alpha would act on this DG0
            # slope, so the capped forward row caps it here.
            cases.append((f"forward half, local, capped at {cap:.0e} on DG0 cells",
                          forward, np.minimum(forward["sin_a"], cap), "cells", None))

    for label, g, sin_a, dofs, melt in cases:
        floating, area = g["floating"], g["area"]
        afl = float(area[floating].sum())
        if melt is None:
            melt = np.where(floating,
                            quadratic_mixed_slope(g["tf"] + g.get("dT", 0.0),
                                                  g["sal"], sin_a, K=g["K"]),
                            0.0)
        over = floating & (melt > bound_m_yr)
        a_over = float(area[over].sum())
        PETSc.Sys.Print(
            f"\n  --- sin_alpha {label} ---\n"
            f"  floating                {afl / 1e6:12.1f} km^2 over "
            f"{int(floating.sum())} {dofs}\n"
            f"  melt max                {melt.max():12.1f} m/yr\n"
            f"  melt p99 (floating)     {np.quantile(melt[floating], 0.99):12.1f} m/yr\n"
            f"  melt area-mean          "
            f"{float((melt * area)[floating].sum()) / afl:12.2f} m/yr\n"
            f"  integrated              "
            f"{float((melt * area)[floating].sum()) * RHO_I / 1e12:12.0f} Gt/yr"
            f"{target}\n"
            f"  {dofs} past the bound    {int(over.sum()):12d}\n"
            f"  area past the bound     {a_over / 1e6:12.1f} km^2 "
            f"({100 * a_over / afl:.3f}% of floating)")

        if not over.any():
            continue
        PETSc.Sys.Print(f"  worst {dofs} (x km, y km, melt m/yr, TF K, "
                        f"draft m, sin_alpha, area km^2):")
        for i in np.argsort(-melt)[:10]:
            PETSc.Sys.Print(
                f"    {g['x'][i] / 1e3:9.1f} {g['y'][i] / 1e3:9.1f} "
                f"{melt[i]:9.1f} {g['tf'][i]:6.2f} {g['draft'][i]:8.1f} "
                f"{sin_a[i]:9.2e} {area[i] / 1e6:8.2f}")
        # The writer's whole-pixel means scale a dof's value by its share of
        # the 8 km pixel, so the median area says how far the grid dilutes the
        # dofs past the bound.
        PETSc.Sys.Print(
            f"  median area of the {dofs} past the bound: "
            f"{np.median(area[over]) / 1e6:.2f} km^2, against "
            f"{ISMIP7_DX ** 2 / 1e6:.0f} km^2 for an 8 km pixel")

    status = 0
    if kind == "deltaT" and "M_obs" in d:
        bad = match_calibration(forward, npz_path, a.match_tol)
        if bad:
            status = 1
            contract = melt_calibration_contract(npz_path) or {}
            fitted_on = contract.get("mesh", "the mesh the offsets were fitted on")
            if contract.get("mesh_build"):
                fitted_on += f" ({contract['mesh_build']})"
            want = contract.get("vertices")
            have = global_size(mesh.coordinates)
            build = (f" This mesh has {have} vertices and the offsets were "
                     f"fitted on {int(want)}, so it is another mesh or another "
                     f"build of it." if want and int(want) != have else "")
            PETSc.Sys.Print(
                f"\n  The forward's melt on this mesh is off the totals its "
                f"offsets were fitted to in basins {bad}.{build} On {fitted_on} "
                f"that means the forward does not apply the melt its "
                f"calibration was fitted to. On another mesh or build it "
                f"measures how far the offsets carry over, and refitting them "
                f"there (calibrate_deltaT.py --K {float(d['K']):.3e}) removes it.")
    if a.ocx is not None:
        status = max(status, compare_with_ocx(a, forward))
    return status


def compare_with_ocx(a, g):
    r"""Melt the forward half with the OCX ocean and set it beside the melt
    from the climatology K was fitted to. Returns the exit status."""
    ocean = ISMIP7Ocean(scenario=OCX, variant=a.ocx)
    years = [int(v) for v in a.ocx_years.split(",") if v.strip()]
    floating, area, sin_a = g["floating"], g["area"], g["sin_a"]
    basin = np.round(cm._grid_interp(cm.IMBIE2_NC, "basinNumber", g["x"], g["y"])).astype(int)
    block, centres = block_ids(g["x"], g["y"])

    def melt_with(tf, sal):
        return np.where(floating, quadratic_mixed_slope(tf + g.get("dT", 0.0), sal,
                                                        sin_a, K=g["K"]), 0.0)

    reference = melt_with(g["tf"], g["sal"])
    PETSc.Sys.Print(
        f"\n=== OCX ocean '{a.ocx}' against the climatology K was fitted to "
        f"(discussion #48) ===\n"
        f"  forward half, uncapped; a region is flagged past "
        f"{100 * a.ocx_tol:.0f}%, regions under {a.ocx_floor:g} Gt/yr ignored")
    # The marker lines a forward logs, so the verdict names the version it was
    # measured on: the reader falls back from its pinned version to the highest
    # one on disk, and upstream regenerates these files (discussion #48).
    for line in describe_forcing_provenance(ocean):
        PETSc.Sys.Print(f"  {line}")
    status = 0
    for year in years:
        other = melt_with(ocean.get_thermal_forcing(year, g["x"], g["y"], draft=g["draft"]),
                          ocean.get_salinity(year, g["x"], g["y"], draft=g["draft"]))
        for name, groups, where in (("IMBIE2 basin", basin, None),
                                    (f"{BLOCK_M / 1e3:.0f} km block", block, centres)):
            ref = melt_by_group(groups, reference, area, floating)
            new = melt_by_group(groups, other, area, floating)
            bad = off_by_more_than(ref, new, a.ocx_tol, a.ocx_floor)
            PETSc.Sys.Print(
                f"\n  --- {year}, by {name}: climatology {sum(ref.values()):.0f} Gt/yr, "
                f"OCX {sum(new.values()):.0f} Gt/yr, {len(bad)} of {len(ref)} flagged ---")
            # every basin, but only the blocks that are off: there are hundreds
            for grp in (sorted(ref) if where is None else
                        sorted(bad, key=lambda b: -abs(new.get(b, 0.0) - ref[b]))):
                at = (f"{grp:6d}" if where is None else
                      f"x {where[grp][0] / 1e3:7.0f} km  y {where[grp][1] / 1e3:7.0f} km")
                ratio = new.get(grp, 0.0) / ref[grp] if ref[grp] > 0 else float("inf")
                PETSc.Sys.Print(
                    f"    {at}  climatology {ref[grp]:8.1f}  OCX {new.get(grp, 0.0):8.1f} Gt/yr"
                    f"  ratio {ratio:5.2f}{'   <-- OFF' if grp in bad else ''}")
            if bad:
                status = 1
    if status:
        PETSc.Sys.Print(
            "\n  The OCX ocean and the climatology disagree by more than the "
            "tolerance somewhere.\n  K is fitted to the climatology, so a "
            "protocol-forced core 11 melts those regions\n  differently from "
            "its own calibration. See discussion #48 before running it.")
    return status


if __name__ == "__main__":
    sys.exit(main())
