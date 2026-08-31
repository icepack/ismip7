#!/usr/bin/env python3
r"""Generate the tracked per-core run report (antarctica/reports/).

Each ISMIP7 core experiment lands in the repo as one commit whose payload is
a small markdown report: run configuration (the ISMIP7_* environment as seen
by this process - run it in the same shell/env as the run itself), budget
rows at marker years from the timeseries CSV, the observational audit
(check_ismip6_track.py) and, for projections, the ISMIP6-ensemble overlay
(compare_ismip6.py), plus provenance (git SHA, log path, checkpoint files).
Results h5/CSVs stay gitignored; the report is the reviewable record.

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
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    sha, _ = sh(["git", "-C", _ANT, "rev-parse", "--short", "HEAD"])
    env = {k: v for k, v in sorted(os.environ.items())
           if k.startswith("ISMIP7_") or k == "OMP_NUM_THREADS"}
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
        f.write(f"- date: {date.today().isoformat()}\n")
        f.write(f"- git: {sha}\n")
        f.write(f"- log: `{args.log}`\n")
        f.write(f"- timeseries: `{args.csv}` (gitignored; this report is "
                f"the tracked record)\n")
        f.write(f"- observational audit: "
                f"{'ON TRACK' if audit_rc == 0 else 'OFF TRACK'}\n")
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
