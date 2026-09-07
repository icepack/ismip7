#!/usr/bin/env python3
r"""Overlay our (projection - CTRL) sea-level contribution on the actual
ISMIP6-Antarctica ensemble (Seroussi et al. 2020 computed scalars).

check_ismip6_track.py audits a single run against observed budget envelopes;
this script answers the sharper question "does our forced response sit inside
the ISMIP6 ensemble?" using the per-model
``computed_ivaf_minus_ctrl_proj_AIS_<group>_<model>_<exp>.nc`` archive
(local mirror: ``~/data/ismip6``; ivaf in m^3, already differenced against
each model's ctrl_proj, time in calendar years).

Our signal is built the same way: VAF(t) - VAF(t0) from the projection CSV,
minus the same difference from the CTRL CSV, sign-flipped to a sea-level
contribution in mm SLE (positive = sea-level rise). Both CSVs come from
run_simulation's timeseries output (vaf_mm_sle is already mm SLE on the
standard 3.625e14 m^2 ocean).

The experiment pool is auto-selected from the projection filename's scenario
(ssp126 -> the RCP2.6 members; ssp370/ssp585 -> the SSP5-8.5 + RCP8.5
members, per Seroussi et al. 2020). The mapping is deliberately best-effort:
the point is magnitude and envelope membership, not exact forcing
equivalence - pass --exps to override the pool explicitly.

Usage:
    python compare_ismip6.py <proj.csv> <ctrl.csv> [--exps exp05,exp07]
                             [--ismip6-dir ~/data/ismip6]

Exit 0 if our trajectory lies inside the ensemble min-max envelope at the
last common year, 1 if outside, 2 on usage/data errors.
"""

import argparse
import csv
import glob
import os
import sys

import numpy as np

OCEAN_AREA_M2 = 3.625e14  # ISMIP6 convention, matches simulation.py's 362.5e12 km^2

# ISMIP6-AIS experiment -> emission-scenario grouping (Seroussi et al. 2020,
# The Cryosphere, "ISMIP6 Antarctica projections"). These are the standard
# core projections; the split by emission scenario gives a broad-consistency
# envelope for our ISMIP7 scenarios. NB: verify against the protocol before
# using for anything quantitative; override with --exps to be explicit.
ISMIP6_POOLS = {
    "rcp26":  ["exp03", "exp07"],                              # NorESM1-M RCP2.6
    "rcp85":  ["exp01", "exp02", "exp04", "exp05", "exp06",    # RCP8.5 (CMIP5),
               "exp08", "exp09", "exp10"],                     #   various ESMs/melt
    "ssp585": ["exp11", "exp12", "exp13"],                     # SSP5-8.5 (CMIP6)
}
# ISMIP7 scenario -> the closest-emission ISMIP6 pool. ssp585 uses the CMIP6
# SSP5-8.5 members (direct analog) plus the RCP8.5 tier for a fuller envelope;
# ssp370 (no direct analog) uses the same high-emission envelope as an upper
# reference; ssp126 uses RCP2.6.
SCENARIO_TO_POOL = {
    "ssp126": ISMIP6_POOLS["rcp26"],
    "ssp370": ISMIP6_POOLS["ssp585"] + ISMIP6_POOLS["rcp85"],
    "ssp585": ISMIP6_POOLS["ssp585"] + ISMIP6_POOLS["rcp85"],
}


def pool_for(proj_fn):
    r"""(scenario, exp-list) auto-selected from the projection filename, or
    (None, None) if no ssp* token is recognized."""
    base = os.path.basename(proj_fn).lower()
    for scen, pool in SCENARIO_TO_POOL.items():
        if scen in base:
            return scen, pool
    return None, None


def load_csv_vaf(fn):
    with open(fn) as f:
        rows = list(csv.DictReader(f))
    if not rows or "vaf_mm_sle" not in rows[0]:
        print(f"{fn}: not a run_simulation timeseries (no vaf_mm_sle)")
        sys.exit(2)
    yr = np.array([float(r["year"]) for r in rows])
    vaf = np.array([float(r["vaf_mm_sle"]) for r in rows])
    return yr, vaf


def our_slr(proj_fn, ctrl_fn):
    r"""(years, mm SLE contribution): -(dVAF_proj - dVAF_ctrl), both from
    their own first output row, CTRL interpolated onto projection years."""
    yp, vp = load_csv_vaf(proj_fn)
    yc, vc = load_csv_vaf(ctrl_fn)
    y0, y1 = max(yp[0], yc[0]), min(yp[-1], yc[-1])
    m = (yp >= y0) & (yp <= y1)
    if not m.any():
        print(f"no overlap between {proj_fn} ({yp[0]}-{yp[-1]}) "
              f"and {ctrl_fn} ({yc[0]}-{yc[-1]})")
        sys.exit(2)
    yrs = yp[m]
    dproj = vp[m] - vp[m][0]
    dctrl = np.interp(yrs, yc, vc) - np.interp(yrs[0], yc, vc)
    return yrs, -(dproj - dctrl)


def load_ensemble(ismip6_dir, exps):
    r"""{(group, model, exp): (years, mm SLE)} from the computed scalars."""
    try:
        import xarray as xr
    except ImportError:
        print("xarray required for the ISMIP6 archive")
        sys.exit(2)
    out = {}
    for exp in exps:
        pat = os.path.join(ismip6_dir, "*", "*", exp,
                           f"computed_ivaf_minus_ctrl_proj_AIS_*_{exp}.nc")
        for fn in sorted(glob.glob(pat)):
            parts = fn.split(os.sep)
            group, model = parts[-4], parts[-3]
            try:
                ds = xr.open_dataset(fn, decode_times=False)
                rhoi = float(ds["rhoi"][0]) if "rhoi" in ds else 917.0
                rhow = float(ds["rhow"][0]) if "rhow" in ds else 1000.0
                yrs = np.asarray(ds["time"].values, dtype=float)
                slr = -np.asarray(ds["ivaf"].values, dtype=float) \
                    * (rhoi / rhow) / OCEAN_AREA_M2 * 1e3
                ds.close()
            except Exception as e:
                print(f"  (skipping {os.path.basename(fn)}: {e})")
                continue
            good = np.isfinite(slr)
            if good.sum() >= 2:
                out[(group, model, exp)] = (yrs[good], slr[good])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("proj_csv")
    ap.add_argument("ctrl_csv")
    ap.add_argument("--exps", default=None,
                    help="ISMIP6 experiments to pool (comma list); default: "
                         "auto-select the emission-scenario analog from the "
                         "projection filename")
    ap.add_argument("--ismip6-dir",
                    default=os.path.expanduser("~/data/ismip6"))
    args = ap.parse_args()

    yrs, ours = our_slr(args.proj_csv, args.ctrl_csv)
    if args.exps:
        exps = [e.strip() for e in args.exps.split(",") if e.strip()]
        scen = "(explicit --exps)"
    else:
        scen, exps = pool_for(args.proj_csv)
        if exps is None:
            print(f"no ssp scenario recognized in {os.path.basename(args.proj_csv)}; "
                  f"pass --exps to choose the ISMIP6 pool")
            sys.exit(2)
    ens = load_ensemble(args.ismip6_dir, exps)
    if not ens:
        print(f"no ISMIP6 members found under {args.ismip6_dir} for {exps}")
        sys.exit(2)

    # ensemble members interpolated onto our years (their time base is
    # 2016..2101; clip our window to it)
    e0 = min(y[0] for y, _ in ens.values())
    e1 = max(y[-1] for y, _ in ens.values())
    m = (yrs >= e0) & (yrs <= e1)
    if not m.any():
        print(f"our window {yrs[0]:.1f}-{yrs[-1]:.1f} has no overlap with "
              f"the ensemble ({e0:.0f}-{e1:.0f})")
        sys.exit(2)
    yrs, ours = yrs[m], ours[m]
    grid = np.array([
        np.interp(yrs, y, s, left=np.nan, right=np.nan)
        for y, s in ens.values()
    ])

    print(f"ISMIP6 ensemble comparison: {os.path.basename(args.proj_csv)} "
          f"- {os.path.basename(args.ctrl_csv)}")
    print(f"  scenario pool {scen}: {len(ens)} members from "
          f"{sorted(set(k[2] for k in ens))}, overlap {yrs[0]:.1f}-{yrs[-1]:.1f}\n")
    print(f"  {'year':>7} {'ours':>8} {'ens p5':>8} {'median':>8} "
          f"{'p95':>8} {'min':>8} {'max':>8}   [mm SLE vs ctrl]")
    marks = [yrs[0]] + [y for y in (2020, 2025, 2050, 2075, 2100)
                        if yrs[0] < y <= yrs[-1]] + [yrs[-1]]
    for y in sorted(set(round(v, 1) for v in marks)):
        i = int(np.argmin(np.abs(yrs - y)))
        col = grid[:, i]
        col = col[np.isfinite(col)]
        if not col.size:
            continue
        print(f"  {yrs[i]:>7.1f} {ours[i]:>8.2f} "
              f"{np.percentile(col, 5):>8.2f} {np.percentile(col, 50):>8.2f} "
              f"{np.percentile(col, 95):>8.2f} "
              f"{col.min():>8.2f} {col.max():>8.2f}")

    col = grid[:, -1]
    col = col[np.isfinite(col)]
    if not col.size:
        print(f"\n  no ensemble member covers year {yrs[-1]:.1f}; cannot "
              f"evaluate envelope membership")
        sys.exit(2)
    lo, hi = float(col.min()), float(col.max())
    med = float(np.median(col))
    spread = max(hi - lo, 1e-9)
    ov = float(ours[-1])
    inside = lo <= ov <= hi
    print(f"\n  at {yrs[-1]:.1f}: ours {ov:+.2f} mm vs ensemble "
          f"[{lo:+.2f}, {hi:+.2f}] (median {med:+.2f}, {col.size} members) -> "
          f"{'INSIDE ensemble envelope' if inside else 'OUTSIDE ensemble envelope'}")
    # SLR here is the sea-level contribution (+ = ice loss). Give a directional
    # parametrization hint when outside, so 'check how we parametrize' is actionable.
    if not inside:
        if ov < lo:
            gap = (lo - ov) / spread
            print(f"  BELOW the ensemble by {lo - ov:.2f} mm ({gap:.1f}x the "
                  f"ensemble spread): too little sea-level contribution -> the "
                  f"sheet is gaining/holding too much mass. Check SMB (the ~9% "
                  f"aSMB unit inflation over-adds snowfall), or melt/discharge "
                  f"being too weak (K_melt, friction).")
        else:
            gap = (ov - hi) / spread
            print(f"  ABOVE the ensemble by {ov - hi:.2f} mm ({gap:.1f}x the "
                  f"ensemble spread): too much sea-level contribution -> too "
                  f"much loss. Check melt (K_melt too high), or numerical "
                  f"discharge/instability at the front.")
    else:
        pctl = 100.0 * float((col < ov).sum()) / col.size
        print(f"  broadly consistent (~{pctl:.0f}th percentile of the ensemble).")
    sys.exit(0 if inside else 1)


if __name__ == "__main__":
    main()
