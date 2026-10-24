#!/usr/bin/env python3
r"""One table for a calving campaign: each member's front flux and front score.

Reads the members' time series (``<run>_timeseries.csv``) and, where present,
their ``score_front.py`` JSONs, and prints a Markdown table, one row per
member. The front flux has two parts the budget books separately: the ice
leaving through the mesh boundary (``outflux_gtyr``; on a mesh without a
buffer the boundary is the front, so this is nearly all of it) and the ice
a calving law removed inside the mesh (``calv_gt`` per step, converted to a
rate with the rows' spacing). Their sum is compared with the observed front
flux, BedMachine's 500 m front under MEaSUREs, 1,138 to 1,297 Gt/yr
(icepack/ismip7#167). Mass and VAF rates are fitted over the run after its
first year, to leave out the start-up transient.

    python antarctica/scripts/campaign_table.py --tag calv75 \
        --results antarctica/results --scores scores/
"""
import argparse
import csv
import glob
import json
import os
import re

import numpy as np

OBS_FRONT_FLUX = (1138.0, 1297.0)  # Gt/yr, BedMachine front under MEaSUREs


def read_series(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    cols = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]
            if k not in ("collapse_flagged_cells", "collapse_removed_cells", "collapse_held_cells")}
    return cols


def rate_after_first_year(t, y):
    m = t >= t[0] + 1.0
    if m.sum() < 2:
        return float("nan")
    return float(np.polyfit(t[m], y[m], 1)[0])


def summarize(cols):
    t = cols["year"]
    dt = float(np.median(np.diff(t))) if len(t) > 1 else float("nan")
    calv_rate = cols["calv_gt"] / dt
    out = {
        "t0": float(t[0]), "t1": float(t[-1]), "years": float(t[-1] - t[0]), "n": int(len(t)),
        "outflux": float(cols["outflux_gtyr"].mean()),
        "calv": float(calv_rate.mean()),
        "melt": float(cols["melt_gtyr"].mean()),
        "smb": float(cols["smb_gtyr"].mean()),
        "dmdt": rate_after_first_year(t, cols["mass_gt"]),
        "dvaf": rate_after_first_year(t, cols["vaf_mm_sle"]),
        "resid_max": float(np.abs(cols["resid_gt"]).max()),
        "clamp": float(cols["clamp_gt"].sum() / max(t[-1] - t[0], dt)),
    }
    out["front_total"] = out["outflux"] + out["calv"]
    return out


def read_score(path):
    with open(path) as f:
        d = json.load(f)
    s = d["summary"]
    last = [r for r in d["states"] if r.get("epoch") is not None]
    inland = [r["inland_ice_free_km2"] for r in last if "inland_ice_free_km2" in r]
    return {"rms": s["mean_rms_m"], "bias": s["mean_bias_m"], "obs_rms": s.get("mean_obs_rms_m"),
            "epochs": len(s["epochs"]),
            "inland_change": (inland[-1] - inland[0]) if len(inland) > 1 else None}


def fmt(x, nd=0):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "-"
    return f"{x:,.{nd}f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True, help="campaign tag, e.g. calv75")
    ap.add_argument("--results", default="antarctica/results")
    ap.add_argument("--scores", default=None, help="directory of <tag>_<member>.json from score_front.py")
    ap.add_argument("--experiment", default="hist_cesm2_waccm")
    args = ap.parse_args()
    pat = os.path.join(args.results, f"{args.experiment}_{args.tag}_*_timeseries.csv")
    files = sorted(glob.glob(pat))
    if not files:
        raise SystemExit(f"no time series match {pat}")
    rows = []
    for path in files:
        m = re.search(rf"{re.escape(args.tag)}_(.+?)_\d+_timeseries\.csv$", os.path.basename(path))
        member = m.group(1) if m else os.path.basename(path)
        s = summarize(read_series(path))
        sc = None
        if args.scores:
            sp = os.path.join(args.scores, f"{args.tag}_{member}.json")
            if os.path.exists(sp):
                sc = read_score(sp)
        rows.append((member, s, sc))
    lo, hi = OBS_FRONT_FLUX
    print(f"Campaign {args.tag}, {args.experiment}: front flux against the observed {lo:,.0f} to {hi:,.0f} Gt/yr "
          f"(BedMachine front under MEaSUREs); rates after the first year.\n")
    print("| member | years | boundary outflux (Gt/yr) | law calving (Gt/yr) | front total (Gt/yr) | vs observed | "
          "dM/dt (Gt/yr) | dVAF/dt (mm SLE/yr) | max resid (Gt) | front rms (km) | bias (km) | obs->model rms (km) | "
          "epochs | inland ice-free change (km2) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for member, s, sc in rows:
        ft = s["front_total"]
        vs = "inside" if lo <= ft <= hi else ("below" if ft < lo else "above")
        score = (fmt(sc["rms"] / 1e3, 1), fmt(sc["bias"] / 1e3, 1), fmt(sc["obs_rms"] / 1e3, 1) if sc["obs_rms"] else "-",
                 str(sc["epochs"]), fmt(sc["inland_change"])) if sc else ("-",) * 5
        print(f"| {member} | {s['t0']:.2f} to {s['t1']:.2f} ({s['years']:.2f}) | {fmt(s['outflux'])} | {fmt(s['calv'], 1)} | "
              f"{fmt(ft)} | {vs} | {fmt(s['dmdt'])} | {fmt(s['dvaf'], 2)} | {s['resid_max']:.1e} | "
              + " | ".join(score) + " |")


if __name__ == "__main__":
    main()
