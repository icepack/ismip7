#!/usr/bin/env python3
r"""Is a run "on the ISMIP6 track"? Audit a timeseries CSV against observed
Antarctic budget envelopes and ISMIP6-class control-run drift.

The ISMIP7 initialization is only trustworthy if the CTRL it launches behaves
like an ISMIP6-era control: component fluxes at observed magnitudes, a small
total drift, and no numerical runaway. This script encodes those envelopes so
"the simulation is on the right track" is a checkable claim rather than a
vibe. Reference numbers:

- Total mass balance dM/dt: IMBIE (Shepherd et al. 2018, Nature) AIS
  1992-2017 = -109 +/- 56 Gt/yr; 2012-2017 = -219 +/- 43 Gt/yr. A control
  may also legitimately drift positive after init, so the envelope is wide.
- Ice-shelf basal melt: Rignot et al. 2013 = 1325 +/- 235 Gt/yr;
  Adusumilli et al. 2020 steady-state = ~1100 Gt/yr; the Paolo/Adusumilli
  IMBIE2 target used for our K calibration integrates 865 Gt/yr.
- Calving flux: Rignot et al. 2013 = 1265 +/- 140 Gt/yr (front discharge
  on a buffered mesh = outflux + fixed-front calving tally).
- SMB (grounded+shelves): RACMO2 ~2300-2740 Gt/yr (van Wessem et al. 2018).
- ISMIP6 ctrl_proj drift (Seroussi et al. 2020, TC): multi-model control
  VAF drift is small vs the forced signal; we require |dVAF/dt| <= 2 mm
  SLE/yr after the first-year init transient (the 2016+ window; the
  post-continuation transient is real but must average out fast).
- Stability: the Jul 2026 split-step blow-up doubled outflux per step.
  Any sustained year-over-year outflux growth factor >= 1.5, or outflux
  beyond 6000 Gt/yr, is an automatic FAIL regardless of the other rows.

Usage:
    python check_ismip6_track.py <timeseries.csv> [--dt DT]
    python check_ismip6_track.py            # newest *_timeseries.csv in results/

Exit code 0 iff no FAIL rows (WARNs allowed), so gates can chain on it.
"""

import csv
import glob
import os
import sys

import numpy as np

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "results")

# (name, lo, hi, fail_lo, fail_hi, units) - WARN outside [lo, hi],
# FAIL outside [fail_lo, fail_hi].
ENVELOPES = {
    "smb":       ("SMB",                 2000.0, 2900.0, 1500.0, 3500.0, "Gt/yr"),
    "melt":      ("shelf basal melt",     600.0, 1800.0,  200.0, 3000.0, "Gt/yr"),
    "discharge": ("front discharge",      700.0, 2400.0,  200.0, 6000.0, "Gt/yr"),
    "dmdt":      ("dM/dt (post-2016)",   -400.0,  200.0, -1000.0, 500.0, "Gt/yr"),
    "dvafdt":    ("dVAF/dt (post-2016)",   -2.0,    2.0,   -5.0,    5.0, "mm SLE/yr"),
    "resid":     ("budget residual",       -0.5,    0.5,   -5.0,    5.0, "Gt/yr"),
}


def load(csv_fn):
    with open(csv_fn) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"{csv_fn}: empty timeseries")
    cols = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return cols


def main():
    args = []
    dt_arg = None
    skip_next = False
    for i, a in enumerate(sys.argv[1:]):
        if skip_next:
            skip_next = False
            continue
        if a == "--dt":
            dt_arg = float(sys.argv[i + 2])
            skip_next = True
            continue
        if not a.startswith("--"):
            args.append(a)
    if args:
        csv_fn = args[0]
    else:
        cands = sorted(glob.glob(os.path.join(RESULTS_DIR, "*_timeseries.csv")),
                       key=os.path.getmtime)
        if not cands:
            raise SystemExit("no *_timeseries.csv in results/")
        csv_fn = cands[-1]

    c = load(csv_fn)
    required = ("year", "mass_gt", "vaf_mm_sle", "smb_gtyr", "melt_gtyr",
                "outflux_gtyr", "calv_gt", "resid_gt")
    missing = [k for k in required if k not in c]
    if missing:
        print(f"{csv_fn}: missing budget column(s): {', '.join(missing)} "
              f"(legacy timeseries format?)")
        sys.exit(2)
    yr = c["year"]
    dt = dt_arg if dt_arg else (float(np.median(np.diff(yr))) if len(yr) > 1 else 1.0)
    n_yr1 = max(1, int(round(1.0 / dt)))          # steps in the init-transient year

    # calv/clamp/resid columns are per-STEP Gt; convert to rates.
    calv_rate = c["calv_gt"] / dt
    resid_rate = c["resid_gt"] / dt
    discharge = c["outflux_gtyr"] + calv_rate     # total front discharge

    post = slice(n_yr1, None) if len(yr) > n_yr1 else slice(None)
    dm = np.diff(c["mass_gt"], prepend=c["mass_gt"][0]) / dt
    dvaf = np.diff(c["vaf_mm_sle"], prepend=c["vaf_mm_sle"][0]) / dt

    vals = {
        "smb":       float(np.mean(c["smb_gtyr"])),
        "melt":      float(np.mean(c["melt_gtyr"])),
        "discharge": float(np.mean(discharge[post])),
        "dmdt":      float(np.mean(dm[post])),
        "dvafdt":    float(np.mean(dvaf[post])),
        "resid":     float(np.abs(resid_rate).max()),
    }

    # Runaway detector: sustained growth of the front discharge. Compare
    # year-block means so dt does not matter; the Jul 2026 blow-up grows
    # ~2x per 0.1-yr step (~1000x/yr) and trips this in year one.
    nblk = max(1, int(round(1.0 / dt)))
    blocks = [np.mean(discharge[i:i + nblk]) for i in range(0, len(discharge), nblk)]
    growth = [blocks[i + 1] / max(abs(blocks[i]), 1e-9) for i in range(len(blocks) - 1)]
    runaway = (max(np.abs(discharge)) > 6000.0
               or any(g >= 1.5 and blocks[i + 1] > 1000.0 for i, g in enumerate(growth)))

    print(f"ISMIP6-track audit: {os.path.basename(csv_fn)}")
    print(f"  {len(yr)} steps, {yr[0]:.1f}->{yr[-1]:.1f}, dt={dt:.3g} yr\n")
    print(f"  {'quantity':<22} {'run':>10}   {'obs/ISMIP6 envelope':<22} verdict")
    n_fail = 0
    for key, (name, lo, hi, flo, fhi, unit) in ENVELOPES.items():
        v = vals[key]
        if flo <= v <= fhi:
            verdict = "PASS" if lo <= v <= hi else "WARN"
        else:
            verdict, n_fail = "FAIL", n_fail + 1
        print(f"  {name:<22} {v:>10.1f}   [{lo:>7.1f}, {hi:>7.1f}] {unit:<9} {verdict}")
    v = "FAIL" if runaway else "PASS"
    n_fail += int(runaway)
    print(f"  {'no discharge runaway':<22} {max(np.abs(discharge)):>10.1f}   "
          f"[peak < 6000, growth<1.5x/yr] {v}")

    print(f"\n  {'ON TRACK' if n_fail == 0 else 'OFF TRACK'} "
          f"({n_fail} FAIL row{'s' if n_fail != 1 else ''})")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
