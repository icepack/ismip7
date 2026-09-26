#!/usr/bin/env python3
r"""Overlay the budget timeseries of several forward runs on one figure.

Reads each run's ``*_timeseries.csv`` (the per-step mass-budget audit that
``run_simulation`` writes: mass, VAF, SMB, melt, outflux, calving, apparent
MB) and plots mass change, VAF change and the flux components against time,
one line per run. This is the quickest way to see whether two runs that
differ in ONE knob (a MAP, ``ISMIP7_APPARENT_MB``, a front treatment) part
ways, and when.

Usage:
    python antarctica/scripts/compare_runs.py LABEL=PATH [LABEL=PATH ...]
        [--out PNG] [--smooth YEARS]

PATH is a ``*_timeseries.csv`` or the run prefix it belongs to. Fluxes are
box-smoothed over ``--smooth`` years (default 1) because the per-step values
jitter from step to step. The step of each row is read from the year column
(``icepack2_tools.timeseries.row_steps``). Read-only.

Units in the CSV are not uniform: the ``*_gtyr`` columns are rates [Gt/yr],
while ``calv_gt``, ``clamp_gt`` and ``resid_gt`` are the mass moved by ONE
step [Gt] and must be divided by dt before they can be added to a rate.
"""

import argparse
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
from icepack2_tools.timeseries import row_steps  # noqa: E402

COLS = ["year", "vaf_mm_sle", "mass_gt", "smb_gtyr", "melt_gtyr",
        "outflux_gtyr", "calv_gt", "clamp_gt", "resid_gt", "amb_gtyr"]


def load(path):
    if not path.endswith(".csv"):
        path = f"{path}_timeseries.csv"
    data = {c: [] for c in COLS}
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                vals = [float(r.get(c, "nan")) for c in COLS]
            except ValueError:
                continue
            for c, v in zip(COLS, vals):
                data[c].append(v)
    return {c: np.asarray(v) for c, v in data.items()}


def smooth(y, n):
    """Centred running mean over n samples; the window shrinks at the ends so
    the series is not dragged toward zero there."""
    if n <= 1 or len(y) < n:
        return y
    c = np.concatenate([[0.0], np.cumsum(y)])
    i = np.arange(len(y))
    lo = np.maximum(i - n // 2, 0)
    hi = np.minimum(i + (n - n // 2), len(y))
    return (c[hi] - c[lo]) / (hi - lo)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs", nargs="+", help="LABEL=PATH")
    ap.add_argument("--out", default=None)
    ap.add_argument("--smooth", type=float, default=1.0)
    args = ap.parse_args()

    runs = []
    for spec in args.runs:
        label, sep, path = spec.rpartition("=")
        if not sep:
            path, label = spec, os.path.basename(spec)
        runs.append((label, load(path)))

    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    ax = axes.ravel()
    for label, d in runs:
        t = d["year"]
        try:
            dt = np.asarray(row_steps(t))
        except ValueError as e:
            print(f"{label}: {e}; not plotted")
            continue
        n = max(1, int(round(args.smooth / dt.mean())))
        ax[0].plot(t, d["mass_gt"] - d["mass_gt"][0], label=label)
        ax[1].plot(t, d["vaf_mm_sle"] - d["vaf_mm_sle"][0], label=label)
        dmdt = np.gradient(d["mass_gt"], np.cumsum(dt))
        ax[2].plot(t, smooth(dmdt, n), label=label)
        ax[3].plot(t, smooth(d["melt_gtyr"], n), label=label)
        ax[4].plot(t, smooth(d["outflux_gtyr"] + d["calv_gt"] / dt, n),
                   label=label)
        ax[5].plot(t, smooth(d["smb_gtyr"] + d["amb_gtyr"], n), label=label)
    titles = ["mass change [Gt]", "VAF change [mm SLE]", "dM/dt [Gt/yr]",
              "shelf melt [Gt/yr]", "front discharge + calving [Gt/yr]",
              "SMB + apparent MB [Gt/yr]"]
    for a, ttl in zip(ax, titles):
        a.set_title(ttl)
        a.grid(alpha=0.3)
        a.set_xlabel("year")
    ax[2].axhline(0, color="k", lw=0.8)
    ax[2].axhspan(-400, 200, color="0.85", zorder=0,
                  label="ISMIP6-track envelope")
    ax[0].legend()
    ax[2].legend(fontsize=8)
    fig.suptitle("forward-run budget comparison", fontsize=14)
    fig.tight_layout()
    out = args.out or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "figs", "compare_runs.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == "__main__":
    main()
