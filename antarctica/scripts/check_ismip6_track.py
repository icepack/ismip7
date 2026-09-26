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
  Adusumilli et al. 2020 steady-state = ~1100 Gt/yr; the IMBIE2 target our
  K calibration uses integrates 1067.4 Gt/yr (Paolo/Davison/Adusumilli; the
  older Paolo/Adusumilli table integrates 865 Gt/yr).
- Calving flux: Rignot et al. 2013 = 1265 +/- 140 Gt/yr (front discharge
  on a buffered mesh = outflux + fixed-front calving tally).
- SMB (grounded+shelves): RACMO2 ~2300-2740 Gt/yr (van Wessem et al. 2018).
- ISMIP6 ctrl_proj drift (Seroussi et al. 2020, TC): multi-model control
  VAF drift is small vs the forced signal; we require |dVAF/dt| <= 2 mm
  SLE/yr after the first-year init transient (the 2016+ window; the
  post-continuation transient is real but must average out fast).
- Stability: the Jul 2026 split-step blow-up doubled outflux per step.
  A numerical runaway is an automatic FAIL regardless of the other rows: front
  discharge growing by >= 1.5x in each of two consecutive years, or a year
  whose median discharge exceeds 6000 Gt/yr. Single-step spikes do not count:
  an emptying event books a step of discharge and the budget closes, and the
  Jul 2026 blow-up that this exists to catch grew ~2x per 0.1-yr step.

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

_ANT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(_ANT))

from icepack2_tools.timeseries import row_steps  # noqa: E402

RESULTS_DIR = os.path.join(_ANT, "results")

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


def runaway_detected(discharge, dt):
    r"""True when the front discharge is running away, not merely spiking.

    ``dt`` is the step, one value or one per row. Blocks of one model year,
    cut by elapsed time, make the test independent of dt. Two clauses, both
    sustained: the block MEDIAN above 6000 Gt/yr (a mean is moved by one step, a median is
    not), or a growth factor >= 1.5 in each of two consecutive years to above
    1000 Gt/yr. Growth uses block means on purpose: a year whose mean is lifted
    1.5x by several large steps counts as growth. The median clause is what
    makes a single step harmless, and the growth clause is the stricter of the
    two (it keeps the 285-year n=4 CTRL, whose growth years carry 2 and 4 steps
    above 6000 Gt/yr out of 10, as FAIL). A trailing block shorter than half a
    year is dropped, since one step would set its median. The Jul 2026 blow-up
    (~2x per 0.1-yr step) fails the median clause in year one; the growth
    clause needs three year blocks and cannot fire before year three. A single
    emptying step, which the peak-of-any-step clause used to flag on cores 2,
    3 and 7 of the July matrix, fails neither.
    """
    discharge = np.abs(np.asarray(discharge, float))
    steps = np.broadcast_to(np.asarray(dt, float), discharge.shape)
    year = np.floor(np.cumsum(steps) - 1e-6).astype(int)
    edges = np.flatnonzero(np.diff(year)) + 1
    blocks = np.split(discharge, edges)
    spans = np.split(steps, edges)
    if len(blocks) > 1 and spans[-1].sum() < 0.5 - 1e-6:
        blocks.pop()
    means = [np.mean(b) for b in blocks]
    medians = [np.median(b) for b in blocks]
    growth = [means[i + 1] / max(means[i], 1e-9) for i in range(len(means) - 1)]
    sustained = any(growth[i] >= 1.5 and growth[i + 1] >= 1.5 and means[i + 2] > 1000.0
                    for i in range(len(growth) - 1))
    return bool(max(medians) > 6000.0 or sustained)


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
    if dt_arg:
        dt = np.full(len(yr), dt_arg)
    else:
        try:
            dt = np.asarray(row_steps(yr))
        except ValueError as e:
            print(f"{csv_fn}: {e}; pass --dt")
            sys.exit(2)

    # calv/clamp/resid columns are per-STEP Gt; convert to rates.
    calv_rate = c["calv_gt"] / dt
    resid_rate = c["resid_gt"] / dt
    discharge = c["outflux_gtyr"] + calv_rate     # total front discharge

    post = np.cumsum(dt) > 1.0 + 1e-6              # past the init-transient year
    if not post.any():
        post = slice(None)
    dm = np.diff(c["mass_gt"], prepend=c["mass_gt"][0]) / dt
    dvaf = np.diff(c["vaf_mm_sle"], prepend=c["vaf_mm_sle"][0]) / dt

    vals = {
        "smb":       float(np.average(c["smb_gtyr"], weights=dt)),
        "melt":      float(np.average(c["melt_gtyr"], weights=dt)),
        "discharge": float(np.average(discharge[post], weights=dt[post])),
        "dmdt":      float(np.average(dm[post], weights=dt[post])),
        "dvafdt":    float(np.average(dvaf[post], weights=dt[post])),
        "resid":     float(np.abs(resid_rate).max()),
    }

    runaway = runaway_detected(discharge, dt)

    print(f"ISMIP6-track audit: {os.path.basename(csv_fn)}")
    dt_lo, dt_hi = dt.min(), dt.max()
    dt_txt = f"{dt_lo:.3g}" if dt_hi - dt_lo <= 1e-6 * dt_hi else f"{dt_lo:.3g}..{dt_hi:.3g}"
    print(f"  {len(yr)} steps, {yr[0]:.1f}->{yr[-1]:.1f}, dt={dt_txt} yr\n")
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
          f"[yr-median < 6000, growth<1.5x for 2 yr] {v}")

    print(f"\n  {'ON TRACK' if n_fail == 0 else 'OFF TRACK'} "
          f"({n_fail} FAIL row{'s' if n_fail != 1 else ''})")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
