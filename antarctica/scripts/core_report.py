#!/usr/bin/env python3
r"""Generate the tracked per-core run report (antarctica/reports/).

Each ISMIP7 core experiment lands in the repo as one commit whose payload is
a small markdown report: run configuration (the ISMIP7_* environment as seen
by this process - run it in the same shell/env as the run itself - with the
run-shaping knobs resolved to their effective values, see effective_env), budget
rows at marker years from the timeseries CSV, the observational audit
(check_ismip6_track.py) and, for projections, the ISMIP6-ensemble overlay
(compare_ismip6.py), plus provenance (git SHA, log path, checkpoint files).
Results h5/CSVs stay gitignored; the report is the reviewable record.

Validity banners are written by --superseded REASON, which emits a blockquote
directly under the title. Use it whenever regenerating a run that is known
invalid: this report is the ONLY committed record (timeseries and checkpoints
are gitignored), so regenerating without the banner silently reinstates an
invalid result as the provenance.
antarctica/reports/MATRIX_STATUS.md owns the full invalidation detail - the
per-core banner is only a short pointer to it.

Usage (from the run shell, so the env is the run env):
    python core_report.py --core 1 --name hist_cesm2_waccm \
        --csv results/hist_cesm2_waccm_32000_timeseries.csv \
        --log results/logs/core1_....log \
        [--ctrl-csv results/ctrl..._timeseries.csv --exps exp01,...]
"""

import argparse
import csv
import os
import subprocess
import sys
from datetime import date

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_ANT = os.path.dirname(_SCRIPTS)
_PROJECT = os.path.dirname(_ANT)
sys.path.insert(0, _PROJECT)

from icepack2_tools.climatology import (
    CLIM_POOL_MARKER, clim_scenario, clim_start, clim_end,
)
from icepack2_tools.runconfig import (
    N_FLOW_DEFAULT, friction, geometry_space, lc, lc_coarse,
)


def effective_env():
    r"""The ISMIP7_* environment as the run sees it, with the run-shaping
    knobs resolved to their EFFECTIVE values rather than only the ones that
    happen to be exported.

    A knob left at its default is absent from ``os.environ``, so recording the
    environment as-set makes two reports byte-identical even when a default
    has since been flipped underneath them - which is exactly how the
    ISMIP7_N_FLOW ambiguity documented in reports/MATRIX_STATUS.md arose, and
    what flipping ISMIP7_CLIM_SCENARIO to ssp126 would otherwise repeat. The
    report is the only committed record of a run, so it has to state the value
    the run used. Defaulted entries are marked so the distinction between
    "exported" and "resolved" is not lost either.

    ISMIP7_GEOMETRY_SPACE is the realized case: MATRIX_STATUS.md records that
    every existing core ran under CG1 while the default is now DG0, and the
    two are NOT interchangeable. So this resolves the whole set of run-shaping
    knobs, not a sample of it.
    """
    env = {k: v for k, v in os.environ.items()
           if k.startswith("ISMIP7_") or k == "OMP_NUM_THREADS"}
    resolved = {
        "ISMIP7_LC": str(lc()),
        "ISMIP7_LC_COARSE": str(lc_coarse()),
        "ISMIP7_FRICTION": friction(),
        "ISMIP7_GEOMETRY_SPACE": geometry_space(),
        "ISMIP7_N_FLOW": N_FLOW_DEFAULT,
        "ISMIP7_CLIM_SCENARIO": clim_scenario(),
        "ISMIP7_CLIM_START": str(clim_start()),
        "ISMIP7_CLIM_END": str(clim_end()),
    }
    for k, v in resolved.items():
        if k not in env:
            env[k] = f"{v}    # default (not exported)"
    return dict(sorted(env.items()))


def climatology_pool(log_path):
    r"""The reference-climate pool the run actually built, lifted out of its
    log.

    The env block records which WINDOW was requested; only the run knows how
    much of it existed on disk. A core re-referenced over 15 of the 30 years
    is a different baseline from one re-referenced over all 30, and the run
    warns rather than refusing, so that fact has to reach the only committed
    record of the run.

    Always returns at least one line. Emitting nothing when the log is absent
    would leave "full window", "half window" and "nobody passed --log"
    indistinguishable in the report, which defeats the point of recording it.
    """
    if not log_path or not os.path.isfile(log_path):
        return [f"{CLIM_POOL_MARKER} NOT RECORDED (log `{log_path}` is not "
                f"readable from here; re-run with --log pointing at the run "
                f"log to capture it)"]
    lines = []
    try:
        with open(log_path, errors="replace") as f:
            for line in f:
                if CLIM_POOL_MARKER in line:
                    stripped = line.strip()
                    if stripped not in lines:
                        lines.append(stripped)
    except OSError as e:
        return [f"{CLIM_POOL_MARKER} NOT RECORDED (log `{log_path}` could not "
                f"be read: {e})"]
    if not lines:
        return [f"{CLIM_POOL_MARKER} none reported in `{log_path}` (expected "
                f"for a CTRL running on the RACMO climatology, which builds "
                f"no ESM pool)"]
    return lines


def sh(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    return (r.stdout + r.stderr).strip(), r.returncode


def csv_marker_rows(fn, n=8):
    with open(fn) as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        return "(empty timeseries)"
    head, body = rows[0], rows[1:]
    idx = sorted({0, len(body) - 1}
                 | {int(round(i * (len(body) - 1) / (n - 1))) for i in range(n)})
    out = [" | ".join(head), " | ".join("---" for _ in head)]
    out += [" | ".join(body[i]) for i in idx]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", type=int, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--log", default="(see results/logs)")
    ap.add_argument("--ctrl-csv", default=None,
                    help="CTRL timeseries for the ISMIP6-ensemble overlay")
    ap.add_argument("--exps", default="exp01,exp02,exp03,exp04,exp05")
    ap.add_argument("--notes", default="")
    ap.add_argument("--superseded", default=None, metavar="REASON",
                    help="Write a SUPERSEDED banner directly under the title. "
                         "Use whenever the run is known invalid: this report "
                         "is the only committed record of the run (the "
                         "timeseries and checkpoints are gitignored), so a "
                         "regenerated report without the banner silently "
                         "reinstates an invalid result as the provenance. "
                         "Keep REASON short and point at "
                         "reports/MATRIX_STATUS.md for the detail.")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    sha, _ = sh(["git", "-C", _ANT, "rev-parse", "--short", "HEAD"])
    env = effective_env()
    audit, audit_rc = sh([sys.executable,
                          os.path.join(_SCRIPTS, "check_ismip6_track.py"),
                          args.csv])
    ens, ens_rc = ("", None)
    if args.ctrl_csv:
        ens, ens_rc = sh([sys.executable,
                          os.path.join(_SCRIPTS, "compare_ismip6.py"),
                          args.csv, args.ctrl_csv, "--exps", args.exps])

    out = args.out or os.path.join(
        _ANT, "reports", f"core{args.core:02d}_{args.name}_32km.md")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write(f"# Core {args.core}: {args.name} (32 km)\n\n")
        # Banner goes directly under the title so a reader cannot miss it.
        # Deliberately a FLAG rather than a hard-coded string: runs made after
        # the ice-front fixes are valid, and a baked-in banner would mislabel
        # them (which is why the docs pass declined to emit one from code).
        if args.superseded:
            # Every line of REASON is quoted: a multi-line reason (the shape
            # every hand-written banner in reports/ already has) would
            # otherwise break out of the blockquote after the first line and
            # render as ordinary body text, splitting the banner in two.
            banner = args.superseded.strip().splitlines() or [""]
            banner[0] = f"**SUPERSEDED.** {banner[0]}"
            banner += ["", "See `reports/MATRIX_STATUS.md` for which results "
                           "are currently valid."]
            for line in banner:
                f.write(("> " + line).rstrip() + "\n")
            f.write("\n")
        f.write(f"- date: {date.today().isoformat()}\n")
        f.write(f"- git: {sha}\n")
        f.write(f"- log: `{args.log}`\n")
        f.write(f"- timeseries: `{args.csv}` (gitignored; this report is "
                f"the tracked record)\n")
        f.write(f"- observational audit: "
                f"{'ON TRACK' if audit_rc == 0 else 'OFF TRACK'}\n")
        for line in climatology_pool(args.log):
            f.write(f"- {line}\n")
        if ens_rc is not None:
            f.write(f"- ISMIP6 ensemble: "
                    f"{'inside envelope' if ens_rc == 0 else 'outside envelope'}"
                    f" (pool: {args.exps})\n")
        if args.notes:
            f.write(f"\n{args.notes}\n")
        f.write("\n## Run environment\n\n```\n")
        for k, v in env.items():
            f.write(f"{k}={v}\n")
        f.write("```\n\n## Budget at marker years\n\n")
        f.write(csv_marker_rows(args.csv))
        f.write("\n\n## Observational audit\n\n```\n" + audit + "\n```\n")
        if ens:
            f.write("\n## ISMIP6 ensemble overlay\n\n```\n" + ens + "\n```\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
