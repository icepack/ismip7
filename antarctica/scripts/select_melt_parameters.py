#!/usr/bin/env python3
r"""Our own K05, K50 and K95: the ISMIP7 toolbox objective run on this model's
cells, with a per-basin thermal-forcing offset fitted at every K first
(ocean-forcing protocol, section 2.3).

    --geometry mesh         the forward's DG0 cells (calibrate_melt.forward_geometry)
    --geometry notebook8km  the toolbox notebook's own 8 km grid: the check that
                            this aggregation reproduces the toolbox's own

For every K of the notebook's grid (120 values, 2.5e-6 to 3.0e-4, or up to
--k-max on the same step), variant ``per_k`` fits one offset per IMBIE2 basin
to the observation table (``fit_deltaT``, the calibrate_deltaT.py fit, window
--dt-window, plus or minus 3 K by default), then melts the present-day
climatology, the 12 ocean-model states and the 13 observed states with those
offsets and builds the toolbox's four terms on the cells
(icepack2_tools/melt_selection.py). Variant ``none`` melts without offsets,
the notebook's own order. The vendored ``calculate_objective_function`` then
draws the notebook's weights, 100000 samples per seed; seed 0 is the headline
and seeds 1 to 4 give the spread.

Under ``per_k`` the objective runs on the K that ``melt_selection.TFRule``
admits: every basin fitted inside the window, and present-day thermal forcing
with the offsets at or above --tf-floor on every floating cell, below and
above the --tf-floor-area and --tf-cap-area temperatures on at most their
shares of each basin, and at or below --tf-cap on every cell (off by
default). The objective over every K, the toolbox's handling, is recorded
beside it.

Knobs come from the environment, as for the forward: ISMIP7_DATA_ROOT and
ISMIP7_MELT_OBS_CSV always; under ``mesh`` also ISMIP7_LC, ISMIP7_INV_H5,
ISMIP7_OI_VERSION (default 30_sep) and the melt knobs calibrate_deltaT.py
reads. Serial only. Everything is written under --out:

    ensemble_<variant>.nc      the aggregates, offsets, fit residuals and TF
                               plausibility at every K, with provenance
    selection_<variant>.json   the percentiles per seed and where they fell,
                               on the admitted K and on every K
    deltaT_per_basin_<lc>_K<K>.npz   per_k: the offsets at each selected K,
                               in calibrate_deltaT.py's keys, for
                               ISMIP7_DELTAT_PER_BASIN_NPZ
    tf_present_<lc>.npz        mesh: present-day TF, area and basin of every
                               floating cell, to test another rule offline
    v1_report.json             notebook8km: the checks against the toolbox

    ISMIP7_LC=2000 ISMIP7_INV_H5=<mesh.h5> ISMIP7_MELT_OBS_CSV=<table.csv> \
        python antarctica/scripts/select_melt_parameters.py --geometry mesh \
        --out /abs/dir [--issue30 DIR]
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
from mpi4py import MPI

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, PROJECT)

from icepack2_tools import melt_selection as ms  # noqa: E402
from icepack2_tools.forcing import (  # noqa: E402
    _RHO_I, _oi_climatology_path, load_deltaT_per_basin, melt_slope,
    quadratic_mixed_slope, sin_alpha_ant, _basin_on_mesh,
)

VARIANTS = ("per_k", "none")
# The aggregates the objective's terms are built from.
AGG = ("t1", "t2", "t3_cold", "t3_warm", "t4")
# What the notebook printed (cells 40 and 59), for the notebook8km report.
NOTEBOOK_PERCENTILES = {"p5": 4.75e-5, "p50": 8.5e-5, "p95": 1.375e-4}
NOTEBOOK_TOTALS = {"p5": 878.0, "p50": 1571.0, "p95": 2541.0}
# The ranges a decoded state must lie in wherever it is finite; a fill value
# the reader failed to mask (1e20 in two salinity files) lands outside.
TF_RANGE, SO_RANGE = (-3.0, 12.0), (25.0, 40.0)


def log(msg):
    print(msg, flush=True)


def sha256(path, bufsize=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(bufsize), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit():
    try:
        return subprocess.run(["git", "-C", PROJECT, "rev-parse", "HEAD"],
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return "unknown"


def provenance(inputs, extra):
    r"""The commit, library versions, toolbox pin and the sha256 of every
    input, so a result names what it was computed from."""
    import scipy
    import xarray
    here = os.path.dirname(os.path.abspath(ms.__file__))
    with open(os.path.join(here, "ismip7_parameter_selection_toolbox.source.json")) as f:
        toolbox = json.load(f)
    log(f"  hashing {len(inputs)} inputs")
    out = {
        "commit": git_commit(),
        "numpy": np.__version__, "xarray": xarray.__version__,
        "scipy": scipy.__version__, "toolbox": toolbox,
        "inputs": {name: {"path": path, "sha256": sha256(path)}
                   for name, path in inputs.items()},
    }
    out.update(extra)
    return out


def require(paths):
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError("missing inputs:\n  " + "\n  ".join(missing))


def check_state(path, var, values, bounds):
    v = np.asarray(values, np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        raise ValueError(f"{path}: no finite {var} on the cells")
    lo, hi = float(v.min()), float(v.max())
    if lo < bounds[0] or hi > bounds[1]:
        raise ValueError(f"{path}: {var} spans {lo:g}..{hi:g}, outside "
                         f"{bounds}; a fill value the reader did not mask?")
    return lo, hi


def stack(rows):
    return {name: np.array([r[name] for r in rows]) for name in rows[0]}


def write_ensemble(path, K, agg, per, attrs):
    r"""The aggregates and per-K fit results as one netCDF, p1 first."""
    import xarray as xr
    models = [label for label, _ in ms.NOTEBOOK_MODELS]
    coords = {"p1": K, "basin": np.arange(ms.N_BASINS),
              "bfrn_bin": np.arange(ms.N_BFRN_BINS), "model": models,
              "region": list(ms.REGIONS), "year": list(ms.NOTEBOOK_OBS_YEARS)}
    dv = {
        "t1": (("p1", "basin"), agg["t1"], "Gt/yr"),
        "t2": (("p1", "bfrn_bin"), agg["t2"], "Gt/yr"),
        "t3_cold": (("p1", "model", "basin"), agg["t3_cold"], "kg m-2 yr-1"),
        "t3_warm": (("p1", "model", "basin"), agg["t3_warm"], "kg m-2 yr-1"),
        "t4": (("p1", "region", "year"), agg["t4"], "Gt/yr"),
    }
    units = {"deltaT": "K", "M_dT0": "Gt/yr", "residual": "Gt/yr",
             "sensitivity": "Gt/yr/K", "unrooted": "1", "tf_mean": "degC",
             "tf_min": "degC", "tf_max": "degC", "frac_below_0": "1",
             "frac_above_5": "1", "frac_below_floor": "1",
             "frac_above_cap": "1"}
    for name, values in per.items():
        dv[name] = (("p1", "basin"), np.asarray(values, dtype=np.float64),
                    units[name])
    ds = xr.Dataset({name: xr.DataArray(v, dims=d, attrs={"units": u})
                     for name, (d, v, u) in dv.items()}, coords=coords)
    ds.attrs.update({k: (v if isinstance(v, (str, int, float)) else json.dumps(v))
                     for k, v in attrs.items()})
    ds.to_netcdf(path)


def run_ensemble(K, variant, present, states, cells, bids, M_obs, rho, plaus,
                 window=None, rule=None):
    r"""The four aggregates and the fit at every K for one variant."""
    rows, per = [], []
    ones = np.ones(len(cells.area), bool)
    tf, so = present
    t0 = time.time()
    for i, k in enumerate(K):
        if variant == "per_k":
            dT, M0, resid, sens, flagged = ms.fit_deltaT(
                tf, so, cells.sin_alpha, k, ones, cells.area, cells.basin,
                bids, M_obs, _SERIAL, window=window)
            unrooted = np.isin(bids, [b for b, _ in flagged])
        else:
            dT = np.zeros(len(bids))
            unrooted = np.zeros(len(bids), bool)
        dT_cell = ms.offsets_on_cells(dT, bids, cells.basin)
        a = ms.aggregate(k, dT_cell, present, states, cells, rho)
        if variant == "none":
            M0 = ms.label_totals(ms.melt_rate(tf, so, cells.sin_alpha, k, rho),
                                 cells.area, cells.basin, ms.N_BASINS)
            resid, sens = M0 - M_obs, np.full(len(bids), np.nan)
        row = {"deltaT": dT, "M_dT0": M0, "residual": resid,
               "sensitivity": sens, "unrooted": unrooted.astype(float)}
        row.update(plaus.at(dT, rule) if plaus is not None else {})
        rows.append(a)
        per.append(row)
        if i % 20 == 0 or i == len(K) - 1:
            log(f"    {variant} K {k:.3e}: t1 {a['t1'].sum():7.1f} Gt/yr, "
                f"offsets {np.nanmin(dT):+.2f}..{np.nanmax(dT):+.2f} K, "
                f"{int(unrooted.sum())} basins at the window edge "
                f"({time.time() - t0:.0f} s)")
    return stack(rows), stack(per)


# The driver is serial; fit_deltaT reduces its basin totals over this.
_SERIAL = MPI.COMM_SELF


def select_all(K, agg, targets, seeds, samples, chunk, reso, unrooted):
    r"""Summaries per seed of the objective on these aggregates."""
    terms = ms.toolbox_terms(K, agg, targets)
    terms.update(ms.notebook_weights(terms, targets["t2_weight"], targets["term4"]))
    by_seed = {}
    for seed in seeds:
        t0 = time.time()
        s = ms.select(terms, samples, chunk, seed, reso=reso)
        by_seed[int(seed)] = ms.summarise(s, K, unrooted)
        log(f"    seed {seed}: K05 {by_seed[int(seed)]['p5']:.4e}  "
            f"K50 {by_seed[int(seed)]['p50']:.4e}  "
            f"K95 {by_seed[int(seed)]['p95']:.4e}  ({time.time() - t0:.0f} s)")
    step = float(np.median(np.diff(K))) if len(K) > 1 else float("nan")
    spread = {q: [min(r[q] for r in by_seed.values()),
                  max(r[q] for r in by_seed.values())] for q in ("p5", "p50", "p95")}
    return terms, {
        "headline": by_seed[int(seeds[0])], "by_seed": by_seed,
        "spread": spread,
        "spread_in_steps": {q: (v[1] - v[0]) / step for q, v in spread.items()},
        "sample_size": int(samples), "chunk": int(chunk),
        "seeds": [int(s) for s in seeds],
    }


def k_range(K, ok):
    r"""The span of the K where ``ok`` holds, and whether it has holes."""
    K, ok = np.asarray(K), np.asarray(ok, bool)
    if not ok.any():
        return None
    lo, hi = float(K[ok].min()), float(K[ok].max())
    inside = (K >= lo) & (K <= hi)
    return {"min": lo, "max": hi, "count": int(ok.sum()),
            "contiguous": bool(ok[inside].all())}


def admitted(K, ens, rule):
    r"""The K the rule admits, with each test's failures, logged."""
    ok, tests = rule.admits(ens)
    report = {"rule": rule.as_dict(), "admitted": k_range(K, ok), "tests": {}}
    for name, passed in tests.items():
        failed = np.asarray(K)[~passed]
        report["tests"][name] = {
            "failing": int(failed.size),
            "failing_range": [float(failed.min()), float(failed.max())] if failed.size else None}
        log(f"  {name}: fails at {failed.size} K"
            + (f", {failed.min():.3e}..{failed.max():.3e}" if failed.size else ""))
    span = report["admitted"]
    log("  admitted: " + ("no K" if span is None else
                          f"{span['count']} K, {span['min']:.3e}..{span['max']:.3e}"
                          + ("" if span["contiguous"] else " with holes")))
    return ok, report


# --- geometry: mesh ----------------------------------------------------------

def run_mesh(args):
    import calibrate_melt as cm
    from icepack2_tools.runconfig import geometry_space, raster_sample

    if cm.GEOMETRY != "dg0":
        raise SystemExit("run with ISMIP7_GEOMETRY_SPACE=dg0 (the forward's path)")
    root = cm.DATA_ROOT
    paths = ms.toolbox_paths(root)
    version = os.environ.get("ISMIP7_OI_VERSION", "30_sep")
    clim = {v: _oi_climatology_path(root, v, version) for v in ("tf", "so")}
    bfrn = ms.bfrn_path(root, cm.LC)
    states = ms.state_files(root)
    inputs = {"mesh": cm.INV_H5, "obs_table": cm.OBS_CSV,
              "clim_tf": clim["tf"], "clim_so": clim["so"],
              "imbie2": cm.IMBIE2_NC, "bfrn": bfrn,
              **{k: paths[k] for k in ("term2", "term3_cold", "term3_warm",
                                       "term4", "bfrn8", "shelf")}}
    for kind, label, tf_p, so_p in states:
        inputs[f"{kind}_{label}_tf"], inputs[f"{kind}_{label}_so"] = tf_p, so_p
    require(inputs.values())
    knobs = {"geometry": "mesh", "lc": cm.LC, "oi_version": version,
             "k_max": float(ms.notebook_K_grid(args.k_max)[-1]),
             "geometry_space": geometry_space(), "raster_sample": raster_sample(),
             "melt_slope": melt_slope(),
             "sin_alpha_ant": sin_alpha_ant() if melt_slope() == "ant" else None,
             "sin_alpha_cap": (cm.SIN_ALPHA_CAP if melt_slope() == "local"
                               else float("inf")),
             "rho_i": float(_RHO_I), "dt_window": list(args.window),
             "tf_rule": args.rule.as_dict()}
    log("=== toolbox selection on the forward's DG0 cells ===")
    for k, v in knobs.items():
        log(f"  {k}: {v}")
    prov = provenance(inputs, knobs)

    targets = ms.load_targets(paths, cm.OBS_CSV)
    bids, M_obs = targets["bids"], targets["M_obs"]
    log(f"  observation table: {cm.OBS_CSV}, {M_obs.sum():.1f} Gt/yr over "
        f"{len(bids)} basins")

    mesh = cm._load_mesh()
    if mesh.comm.size != 1:
        raise SystemExit("serial only: run without mpiexec")
    g = cm.forward_geometry(mesh)
    fl = g["floating"]
    x, y, draft, area = g["x"][fl], g["y"][fl], g["draft"][fl], g["area"][fl]
    sin_a = np.asarray(g["sin_a"])[fl]
    if cm.MELT_SLOPE != "ant":
        sin_a = np.minimum(sin_a, cm.SIN_ALPHA_CAP)
    log(f"  {len(x)} floating cells, {area.sum() / 1e6:.0f} km2")

    # Present day at fill 0, as the forward and calibrate_deltaT.py read it.
    tf = cm._grid_interp(clim["tf"], "tf", x, y, draft=draft)
    so = cm._grid_interp(clim["so"], "so", x, y, draft=draft)
    gaps = int(np.isnan(ms.sample_nearest(clim["tf"], "tf", x, y, draft)).sum())
    log(f"  present-day TF {tf.min():.2f}..{tf.max():.2f} degC; floating cells "
        f"with no climatology (read as 0): {gaps}")

    basin = ms.as_labels(ms.sample_nearest(cm.IMBIE2_NC, "basinNumber", x, y))
    stamped = _basin_on_mesh(x, y, cm.IMBIE2_NC)
    differ = int(np.count_nonzero(stamped != basin))
    if differ:
        raise SystemExit(f"{differ} floating cells take another basin from "
                         f"forcing._basin_on_mesh than from the fit's labels")
    bfrn_bin = ms.as_labels(ms.sample_nearest(bfrn, "BFRN_bins", x, y))
    rx, ry, rmap = ms.amundsen_regions(paths["shelf"])
    region = ms.as_labels(ms.sample_grid(rx, ry, rmap, x, y))
    region[region < 0] = 0
    cells = ms.Cells(area, basin, bfrn_bin, region, sin_a)
    area_by_basin = np.bincount(basin[basin >= 0], weights=area[basin >= 0],
                                minlength=ms.N_BASINS) / 1e6
    log("  floating km2 per basin: " + ", ".join(
        f"{b}:{a:.0f}" for b, a in enumerate(area_by_basin)))
    log(f"  cells with no basin {int((basin < 0).sum())}, no buttressing bin "
        f"{int((bfrn_bin < 0).sum())} ({area[bfrn_bin < 0].sum() / 1e6:.0f} km2); "
        f"Pine Island {area[region == 1].sum() / 1e6:.0f} km2, "
        f"Dotson {area[region == 2].sum() / 1e6:.0f} km2")
    # With the offsets in ensemble_per_k.nc, any other thermal forcing rule
    # can be tested from this file alone.
    np.savez_compressed(os.path.join(args.out, f"tf_present_{cm.LC}.npz"),
                        tf=tf, area=area, basin=basin.astype(np.int16),
                        x=x.astype(np.float32), y=y.astype(np.float32))

    sampled = []
    for kind, label, tf_p, so_p in states:
        t = ms.sample_nearest(tf_p, "tf", x, y, draft)
        s = ms.sample_nearest(so_p, "so", x, y, draft)
        tr = check_state(tf_p, "tf", t, TF_RANGE)
        sr = check_state(so_p, "so", s, SO_RANGE)
        log(f"  {kind:4s} {label!s:18s} tf {tr[0]:+.2f}..{tr[1]:+.2f}  so "
            f"{sr[0]:.2f}..{sr[1]:.2f}  cells without data "
            f"{int(np.isnan(t).sum())}")
        sampled.append((kind, label, t, s))

    K = ms.notebook_K_grid(args.k_max)
    plaus = ms.Plausibility(tf, area, basin)
    report = {"basin_area_km2": area_by_basin.tolist(), "floating_cells": int(len(x)),
              "climatology_gaps": gaps}
    status = 0
    for variant in args.variants:
        log(f"\n  --- variant {variant} ---")
        agg, per = run_ensemble(K, variant, (tf, so), sampled, cells, bids,
                                M_obs, float(_RHO_I), plaus, args.window,
                                args.rule)
        ens_path = os.path.join(args.out, f"ensemble_{variant}.nc")
        write_ensemble(ens_path, K, agg, per,
                       {"variant": variant, "provenance": prov})
        ens = _load_ensemble(ens_path)
        Ks = ens["p1"]
        unrooted_any = ens["unrooted"].astype(bool).any(axis=1)
        base = {"variant": variant, "ensemble": os.path.basename(ens_path),
                "rooted_K_range": k_range(Ks, ~unrooted_any),
                "unrooted_K": Ks[unrooted_any].tolist(),
                "provenance": prov, "mesh_report": report}
        if variant != "per_k":
            _, result = select_all(Ks, ens, targets, args.seeds, args.samples,
                                   args.chunk, cm.LC, unrooted_any)
            result.update(base)
            _dump(os.path.join(args.out, f"selection_{variant}.json"), result)
            continue
        log("  thermal forcing rule, " + json.dumps(args.rule.as_dict()))
        ok, rule_report = admitted(Ks, ens, args.rule)
        log("  every K, unfittable K kept with their term 1 penalty (the toolbox's handling):")
        _, every = select_all(Ks, ens, targets, args.seeds[:1], args.samples,
                              args.chunk, cm.LC, unrooted_any)
        base.update({"tf_rule": rule_report,
                     "every_K": {k: every[k] for k in ("headline", "sample_size",
                                                       "chunk", "seeds")}})
        if not ok.any():
            log("  the rule admits no K: nothing selected")
            _dump(os.path.join(args.out, f"selection_{variant}.json"), base)
            status = 3
            continue
        log("  admitted K:")
        _, result = select_all(Ks[ok], {n: ens[n][ok] for n in AGG}, targets,
                               args.seeds, args.samples, args.chunk, cm.LC,
                               unrooted_any[ok])
        result.update(base)
        result["selected"] = write_selected(
            args, cm, result["headline"], Ks, ens, (tf, so), sin_a, x, y,
            area, basin, bids, M_obs, prov, plaus)
        if args.issue30:
            result["issue30"] = compare_issue30(args.issue30, cm.LC, Ks, ens)
        _dump(os.path.join(args.out, f"selection_{variant}.json"), result)
    log("done")
    return status


def _load_ensemble(path):
    import xarray as xr
    with xr.open_dataset(path) as ds:
        return {name: ds[name].values for name in list(ds.data_vars) + ["p1"]}


def write_selected(args, cm, headline, K, ens, present, sin_a, x, y, area,
                   basin, bids, M_obs, prov, plaus):
    r"""The offsets at K05, K50 and K95 (seed 0) in calibrate_deltaT.py's
    keys, each read back through the forward's loader and melted again, with
    the thermal forcing rule's verdict at that K."""
    from icepack2_tools.runconfig import geometry_space
    tf, so = present
    chosen = {}
    for name, q in (("K05", "p5"), ("K50", "p50"), ("K95", "p95")):
        chosen.setdefault(float(headline[q]), []).append(name)
    out = {}
    for k, names in chosen.items():
        hit = np.flatnonzero(np.isclose(K, k, rtol=1e-9, atol=0))
        if hit.size:
            i = int(hit[0])
            dT, M0 = ens["deltaT"][i], ens["M_dT0"][i]
            resid, sens = ens["residual"][i], ens["sensitivity"][i]
            unrooted = ens["unrooted"][i]
            how = "grid"
        else:
            dT, M0, resid, sens, flagged = ms.fit_deltaT(
                tf, so, sin_a, k, np.ones(len(area), bool), area, basin, bids,
                M_obs, _SERIAL, window=args.window)
            unrooted = np.isin(bids, [b for b, _ in flagged]).astype(float)
            how = "refitted off the grid"
        stats = plaus.at(dT, args.rule)
        stats["unrooted"] = unrooted
        admits, tests = args.rule.admits({n: np.asarray(v)[None] for n, v in stats.items()})
        fn = os.path.join(args.out, f"deltaT_per_basin_{cm.LC}_K{k:.3e}.npz")
        ms.write_offsets(
            fn, bids, dT, k, M_obs, M0, resid, sens,
            melt_slope=melt_slope(),
            sin_alpha_ant=(sin_alpha_ant() if melt_slope() == "ant"
                           else float("nan")),
            sin_alpha_cap=(cm.SIN_ALPHA_CAP if melt_slope() == "local"
                           else float("inf")),
            geometry_space=geometry_space(), obs_csv=cm.OBS_CSV,
            imbie2_nc=cm.IMBIE2_NC, inversion=cm.INV_H5,
            selected_as=",".join(names),
            toolbox_commit=prov["toolbox"]["commit"],
            selection_commit=prov["commit"],
            dt_window=np.asarray(args.window, np.float64),
            tf_rule=json.dumps(args.rule.as_dict()))
        field, K_file = load_deltaT_per_basin(fn, x, y, imbie2=cm.IMBIE2_NC)
        melt = quadratic_mixed_slope(tf + field, so, sin_a, K=K_file) * float(_RHO_I)
        totals = ms.label_totals(melt, area, basin, ms.N_BASINS)
        expect = M_obs + np.nan_to_num(resid)
        out[",".join(names)] = {
            "K": k, "how": how, "file": os.path.basename(fn),
            "deltaT": np.asarray(dT).tolist(),
            "reload_max_offset_diff_K": float(np.max(np.abs(
                field - ms.offsets_on_cells(dT, bids, basin)))),
            "reload_max_total_diff_gtyr": float(np.max(np.abs(totals - expect))),
            "total_gtyr": float(totals.sum()),
            "rule_admits": bool(admits[0]),
            "rule_tests": {n: bool(t[0]) for n, t in tests.items()},
            "tf_range": [float(np.nanmin(stats["tf_min"])), float(np.nanmax(stats["tf_max"]))],
            "sensitivity_gt_per_K": float(np.nansum(sens)),
        }
        if not admits[0]:
            log(f"  {','.join(names)} = {k:.4e} fails the rule: "
                + ", ".join(n for n, t in tests.items() if not t[0]))
        log(f"  {','.join(names)} = {k:.4e} ({how}): offsets "
            f"{np.nanmin(dT):+.2f}..{np.nanmax(dT):+.2f} K, reloaded total "
            f"{totals.sum():.1f} Gt/yr, max basin difference "
            f"{out[','.join(names)]['reload_max_total_diff_gtyr']:.2e}")
    return out


def compare_issue30(directory, lc, K, ens):
    r"""The offsets at the grid K equal to issue 30's files, against them."""
    import glob
    out = {}
    for fn in sorted(glob.glob(os.path.join(directory, f"deltaT_per_basin_{lc}_K*.npz"))):
        z = np.load(fn)
        k = float(z["K"])
        hit = np.flatnonzero(np.isclose(K, k, rtol=1e-9, atol=0))
        if not hit.size:
            continue
        i = int(hit[0])
        d = float(np.nanmax(np.abs(ens["deltaT"][i] - z["deltaT_basin"])))
        out[f"{k:.3e}"] = {"max_offset_diff_K": d,
                           "M_dT0_total": float(ens["M_dT0"][i].sum()),
                           "issue30_M_dT0_total": float(np.sum(z["M_dT0"]))}
        log(f"  issue 30 at K {k:.3e}: max offset difference {d:.1e} K, "
            f"uncorrected {ens['M_dT0'][i].sum():.1f} against "
            f"{np.sum(z['M_dT0']):.1f} Gt/yr")
    return out


# --- geometry: notebook8km ---------------------------------------------------

def run_notebook8km(args):
    import xarray as xr

    root = os.environ["ISMIP7_DATA_ROOT"]
    melt_csv = os.environ["ISMIP7_MELT_OBS_CSV"]
    paths = ms.toolbox_paths(root)
    states = ms.state_files(root)
    inputs = {"obs_table": melt_csv,
              **{k: paths[k] for k in ("term2", "term3_cold", "term3_warm",
                                       "term4", "bfrn8", "shelf", "imbie2",
                                       "floating8", "bedmap3", "tf_06nov",
                                       "so_06nov")}}
    for kind, label, tf_p, so_p in states:
        inputs[f"{kind}_{label}_tf"], inputs[f"{kind}_{label}_so"] = tf_p, so_p
    require(inputs.values())
    rho, sin_a = ms.NOTEBOOK_RHO_I, ms.NOTEBOOK_SIN_ALPHA
    knobs = {"geometry": "notebook8km", "rho_i": rho, "sin_alpha": sin_a,
             "k_max": float(ms.notebook_K_grid(args.k_max)[-1])}
    log("=== toolbox selection on the notebook's 8 km grid ===")
    prov = provenance(inputs, knobs)
    targets = ms.load_targets(paths, melt_csv)
    report = {}

    clim_tf = xr.load_dataset(paths["tf_06nov"])
    clim_so = xr.load_dataset(paths["so_06nov"])
    bed = xr.load_dataset(paths["bedmap3"])
    xg, yg = clim_tf["x"].values, clim_tf["y"].values
    if bed["draft"].shape != (len(yg), len(xg)):
        raise SystemExit(f"BedMap3 is {bed['draft'].shape}, the grid "
                         f"{(len(yg), len(xg))}")

    def on_grid(path, var):
        r"""A (y, x) field on the climatology's grid: checked by coordinate
        where the file carries them, paired by position, as the notebook
        pairs it, where it does not."""
        with xr.open_dataset(path) as ds:
            da = ds[var].load()
        if "x" in da.coords and "y" in da.coords:
            if not (np.array_equal(da["x"].values, xg)
                    and np.array_equal(da["y"].values, yg)):
                raise SystemExit(f"{path} is on another grid")
        elif da.shape == (len(yg), len(xg)):
            da = da.assign_coords(y=yg, x=xg)
        else:
            raise SystemExit(f"{path}: {var} is {da.shape}, no coordinates")
        return da.transpose("y", "x")

    fmask = on_grid(paths["floating8"], "mask")
    ff = bed["floating_frac"].values > 0.5
    agree = float(np.mean(ff == (fmask.values > 0.5)))
    flipped = float(np.mean(ff[::-1, :] == (fmask.values > 0.5)))
    report["orientation"] = {"agreement": agree, "flipped_agreement": flipped}
    log(f"  BedMap3 floating_frac > 0.5 against the floating mask: {agree:.4f} "
        f"agreement ({flipped:.4f} flipped in y)")
    if agree < 0.95 or flipped > agree:
        raise SystemExit("BedMap3 does not pair with the grid by position")
    draft = bed["draft"].values
    # The sign is unstated in the file. The conservative remap leaves a few
    # partly floating margin points with a positive mean draft, which the
    # notebook melts at its shallowest level, as the clipped sampler does.
    above = int(np.count_nonzero(~(draft[ff] < 0)))
    report["draft_not_negative_under_shelves"] = above
    log(f"  BedMap3 draft under the shelves: {above} of {int(ff.sum())} points "
        f"not below sea level")
    if above > 0.01 * ff.sum():
        raise SystemExit("BedMap3 draft is not negative under the shelves")

    recipe = ms.notebook_mean_slope(draft, ff, 8000.0, 8000.0)
    report["recipe_sin_alpha"] = recipe
    log(f"  the notebook's slope recipe on this BedMap3: sin(alpha) = "
        f"{recipe:.5e} (its outputs carry {sin_a:.5e})")

    # The notebook's own sampling (cells 10, 12, 15): xarray nearest in z at
    # each point's draft, NaN outside floating_frac > 0.5.
    bed_draft, floating = bed["draft"], bed["floating_frac"] > 0.5

    def notebook_field(path, var):
        with xr.open_dataset(path) as ds:
            da = ds[var].sel(z=bed_draft, method="nearest").load()
        return da.where(floating, np.nan).astype(np.float64)

    fields = {"present": (notebook_field(paths["tf_06nov"], "tf"),
                          notebook_field(paths["so_06nov"], "so"))}
    for kind, label, tf_p, so_p in states:
        fields[(kind, label)] = (notebook_field(tf_p, "tf"), notebook_field(so_p, "so"))
    del clim_tf, clim_so

    # The cells: the grid points the notebook melts, 64 km2 each.
    jy, jx = np.nonzero(ff)
    area = np.full(len(jy), 8000.0 ** 2)
    basins_m = on_grid(paths["imbie2"], "basinNumber")
    bfrn_m = on_grid(paths["bfrn8"], "BFRN_bins").to_dataset(name="BFRN_bins")
    rx, ry, rmap = ms.amundsen_regions(paths["shelf"], xg, yg)
    XX, YY = np.meshgrid(xg, yg)
    rgrid = ms.as_labels(ms.sample_grid(rx, ry, rmap, XX.ravel(), YY.ravel())
                         ).reshape(XX.shape)
    rgrid[rgrid < 0] = 0
    cells = ms.Cells(area, ms.as_labels(basins_m.values[jy, jx]),
                     ms.as_labels(bfrn_m["BFRN_bins"].values[jy, jx]),
                     rgrid[jy, jx], sin_a, agg=fmask.values[jy, jx] > 0.5)

    def at_cells(da):
        return np.asarray(da.values[jy, jx], dtype=np.float64)

    present = tuple(at_cells(f) for f in fields["present"])
    sampled = []
    for kind, label, tf_p, so_p in states:
        t, s = (at_cells(f) for f in fields[(kind, label)])
        check_state(tf_p, "tf", t, TF_RANGE)
        check_state(so_p, "so", s, SO_RANGE)
        sampled.append((kind, label, t, s))

    K = ms.notebook_K_grid(args.k_max)
    log(f"  {len(jy)} grid points melt, {int(cells.basin_agg.__ge__(0).sum())} "
        f"inside the floating mask")
    agg, per = run_ensemble(K, "none", present, sampled, cells,
                            targets["bids"], targets["M_obs"], rho, None)
    ens_path = os.path.join(args.out, "ensemble_none.nc")
    write_ensemble(ens_path, K, agg, per, {"variant": "none", "provenance": prov})
    ens = _load_ensemble(ens_path)
    Ks = ens["p1"]
    terms, result = select_all(Ks, ens, targets, args.seeds, args.samples,
                               args.chunk, 8000, None)

    # V1a: the toolbox's own term functions on the notebook's gridded
    # ensembles (cells 10 to 36), from the same melt.
    log("  V1a: the toolbox's calculate_term1..4 on the gridded ensembles")
    t0 = time.time()
    tb = ms.notebook_terms(
        Ks, rho, sin_a, fields["present"],
        {key: val for key, val in fields.items() if key != "present"},
        fmask, basins_m.rename("basins"), bfrn_m, rgrid, targets,
        inputs["obs_table"])
    log(f"    built in {time.time() - t0:.0f} s")
    v1a = ms.compare_terms(terms, tb)
    for name, v in v1a.items():
        log(f"    {name}: NaN pattern {'equal' if v['same_nan_pattern'] else 'DIFFERS'}, "
            f"max relative difference {v['max_rel_diff']:.2e}")
    ours = ms.select(terms, args.samples, args.chunk, args.seeds[0], reso=8000)
    theirs = ms.select(tb, args.samples, args.chunk, args.seeds[0], reso=8000)
    v1a["min_p1_differing_samples"] = int(np.count_nonzero(ours != theirs))
    v1a["pass"] = bool(all(v["same_nan_pattern"] and v["max_rel_diff"] <= 1e-9
                           for k, v in v1a.items() if k.endswith("_model"))
                       and v1a["min_p1_differing_samples"] <= 10)
    log(f"    min_p1 differs in {v1a['min_p1_differing_samples']} of "
        f"{len(ours)} samples; V1a {'PASS' if v1a['pass'] else 'FAIL'}")
    report["v1a"] = v1a

    # V1c: the mesh path's sampler (scipy nearest, as the forward reads) on
    # the same points against the notebook's xarray selection.
    log("  V1c: the mesh sampler on the grid points")
    with xr.open_dataset(paths["tf_06nov"]) as ds:
        z = ds["z"].values
    mids = np.sort(0.5 * (np.sort(z)[1:] + np.sort(z)[:-1]))
    dcell = draft[jy, jx]
    tie = np.min(np.abs(dcell[:, None] - mids[None, :]), axis=1) < 1e-6
    xc, yc = xg[jx], yg[jy]
    v1c = {"tie_cells": int(tie.sum()), "fields": {}}
    checks = [("present", "tf", paths["tf_06nov"], present[0]),
              ("present", "so", paths["so_06nov"], present[1])]
    for (kind, label, tf_p, so_p), (_, _, t, s) in zip(states, sampled):
        checks += [(f"{kind}_{label}", "tf", tf_p, t), (f"{kind}_{label}", "so", so_p, s)]
    bad = 0
    for name, var, path, ref in checks:
        got = ms.sample_nearest(path, var, xc, yc, dcell).astype(np.float64)
        # the notebook reads float64 files as float64; the sampler loads
        # float32, the forward's precision
        differ = ~((got == ref) | (np.isnan(got) & np.isnan(ref))
                   | np.isclose(got, ref, rtol=1e-6, atol=1e-6))
        off_tie = int(np.count_nonzero(differ & ~tie))
        v1c["fields"][f"{name}_{var}"] = {"differ": int(differ.sum()),
                                          "differ_off_ties": off_tie}
        bad += off_tie
    v1c["pass"] = bad == 0
    log(f"    {v1c['tie_cells']} cells sit on a z midpoint; values differing "
        f"off those: {bad}; V1c {'PASS' if v1c['pass'] else 'FAIL'}")
    report["v1c"] = v1c

    # V1b: against what the notebook printed.
    head = result["headline"]
    totals = {}
    for q, k in (("p5", NOTEBOOK_PERCENTILES["p5"]), ("p50", NOTEBOOK_PERCENTILES["p50"]),
                 ("p95", NOTEBOOK_PERCENTILES["p95"])):
        i = int(np.flatnonzero(np.isclose(Ks, k, rtol=1e-9, atol=0))[0])
        totals[q] = float(np.nansum(ens["t1"][i]))
    report["v1b"] = {
        "totals_at_notebook_K": totals, "notebook_totals": NOTEBOOK_TOTALS,
        "totals_rel_diff": {q: totals[q] / NOTEBOOK_TOTALS[q] - 1 for q in totals},
        "percentiles": {q: head[q] for q in ("p5", "p50", "p95")},
        "notebook_percentiles": NOTEBOOK_PERCENTILES,
        "percentiles_in_steps": {q: (head[q] - NOTEBOOK_PERCENTILES[q]) / 2.5e-6
                                 for q in NOTEBOOK_PERCENTILES},
        "recipe_sin_alpha": recipe, "notebook_sin_alpha": sin_a,
    }
    log("  V1b: totals at the notebook's K " + ", ".join(
        f"{totals[q]:.1f} ({NOTEBOOK_TOTALS[q]:.0f})" for q in totals)
        + "; percentiles " + ", ".join(
        f"{head[q]:.4e} ({NOTEBOOK_PERCENTILES[q]:.4e})" for q in NOTEBOOK_PERCENTILES))
    result.update({"variant": "none", "ensemble": os.path.basename(ens_path),
                   "provenance": prov})
    _dump(os.path.join(args.out, "selection_none.json"), result)
    _dump(os.path.join(args.out, "v1_report.json"), report)
    gate = report["v1a"]["pass"] and report["v1c"]["pass"]
    log(f"done; gate (V1a and V1c) {'PASS' if gate else 'FAIL'}")
    return 0 if gate else 1


def _dump(path, obj):
    def default(o):
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=default)
    log(f"  wrote {path}")


def _bound(text, flag):
    if text.strip().lower() == "none":
        return None
    try:
        return float(text)
    except ValueError:
        raise SystemExit(f"{flag} takes degC or none, not {text!r}")


def _area_bound(words, flag):
    if len(words) == 1 and words[0].strip().lower() == "none":
        return None
    try:
        degc, share = (float(w) for w in words)
    except ValueError:
        raise SystemExit(f"{flag} takes DEGC SHARE or none, not {' '.join(words)!r}")
    if not 0.0 <= share <= 1.0:
        raise SystemExit(f"{flag}: the share {share:g} is not in [0, 1]")
    return degc, share


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--geometry", choices=("mesh", "notebook8km"), required=True)
    ap.add_argument("--out", required=True, help="absolute output directory")
    ap.add_argument("--variants", nargs="+", choices=VARIANTS, default=list(VARIANTS),
                    help="mesh only; notebook8km runs none")
    ap.add_argument("--samples", type=int, default=ms.NOTEBOOK_SAMPLE_SIZE)
    ap.add_argument("--chunk", type=int, default=10000)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--issue30", help="mesh: issue 30's deltaT files for this lc")
    ap.add_argument("--k-max", type=float, default=None,
                    help="extend the notebook's K grid (2.5e-6 to 3.0e-4) upward "
                         "on its step to this K")
    rule = ms.TFRule()

    def shown(v):
        return ["none"] if v is None else [str(x) for x in (v if isinstance(v, tuple) else (v,))]
    ap.add_argument("--dt-window", type=float, default=ms.DT_WINDOW[1],
                    help="search each basin's offset in plus or minus this many K "
                         "(the toolbox uses 2)")
    ap.add_argument("--tf-floor", default=shown(rule.floor)[0], metavar="DEGC|none",
                    help="per_k: every floating cell at or above this")
    ap.add_argument("--tf-floor-area", nargs="+", default=shown(rule.floor_area),
                    metavar="DEGC SHARE|none",
                    help="per_k: at most SHARE of each basin's area below DEGC")
    ap.add_argument("--tf-cap", default=shown(rule.cap)[0], metavar="DEGC|none",
                    help="per_k: every floating cell at or below this")
    ap.add_argument("--tf-cap-area", nargs="+", default=shown(rule.cap_area),
                    metavar="DEGC SHARE|none",
                    help="per_k: at most SHARE of each basin's area above DEGC")
    ap.add_argument("--keep-unfitted", action="store_true",
                    help="per_k: admit a K at which some basin cannot reach its "
                         "observed total inside the window")
    args = ap.parse_args()
    if not os.path.isabs(args.out):
        raise SystemExit("--out must be absolute")
    if not args.dt_window > 0:
        raise SystemExit("--dt-window must be positive")
    args.window = (-args.dt_window, args.dt_window)
    args.rule = ms.TFRule(floor=_bound(args.tf_floor, "--tf-floor"),
                          floor_area=_area_bound(args.tf_floor_area, "--tf-floor-area"),
                          cap=_bound(args.tf_cap, "--tf-cap"),
                          cap_area=_area_bound(args.tf_cap_area, "--tf-cap-area"),
                          rooted=not args.keep_unfitted)
    os.makedirs(args.out, exist_ok=True)
    if args.geometry == "mesh":
        return run_mesh(args) or 0
    return run_notebook8km(args)


if __name__ == "__main__":
    sys.exit(main())
