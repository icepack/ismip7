#!/usr/bin/env python3
"""
Aggregate per-run timing JSON records into TIMING_MATRIX.md.

Reads results/timing/timing_*.json and writes one table per solver/timing tag.
Qualification probes and incomplete/non-five-year records are excluded.

Usage:
    python scripts/build_timing_matrix.py
"""

import glob
import json
import os
from datetime import datetime, timezone

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIMING_DIR = os.path.join(_ROOT, "results", "timing")
OUT_FN = os.path.join(_ROOT, "TIMING_MATRIX.md")


def load_records():
    records = []
    for path in sorted(glob.glob(os.path.join(TIMING_DIR, "timing_*.json"))):
        with open(path) as f:
            records.append(json.load(f))
    return records


def matrix_records(records):
    return [
        record for record in records
        if record.get("timing_kind", "matrix") == "matrix"
        and record.get("t_start") == 2015.0
        and record.get("t_end") == 2020.0
        and record.get("dt") == 1.0
        and record.get("nsteps") == 5
        and record.get("completed_steps") == 5
        and record.get("t_final") == 2020.0
    ]


def render_table(records, per_step=False):
    core_counts = sorted({record["ncores"] for record in records})
    rows = sorted(
        {(record["lc"], record["lc_coarse"]) for record in records},
        key=lambda value: (value[0], value[1]),
    )
    index = {
        (record["lc"], record["lc_coarse"], record["ncores"]): record
        for record in records
    }
    unit = "s/step" if per_step else "s"
    key = "seconds_per_step" if per_step else "run_seconds"
    digits = 2 if per_step else 1
    header = (
        "| LC (m) | LC_coarse (m) | Vertices | Cells | "
        + " | ".join(f"{count} cores ({unit})" for count in core_counts)
        + " |"
    )
    lines = [header, "|" + "|".join(["---"] * (4 + len(core_counts))) + "|"]
    for lc, lc_coarse in rows:
        candidates = [
            index.get((lc, lc_coarse, count)) for count in core_counts
        ]
        sample = next((record for record in candidates if record), None)
        vertices = sample["vertices"] if sample else "—"
        cells = sample["cells"] if sample else "—"
        timings = [
            f"{record[key]:.{digits}f}" if record else "—"
            for record in candidates
        ]
        lines.append(
            f"| {lc} | {lc_coarse} | {vertices} | {cells} | "
            + " | ".join(timings)
            + " |"
        )
    return lines


def main():
    records = matrix_records(load_records())
    if not records:
        raise SystemExit(
            f"No timing records found in {TIMING_DIR}/. "
            "Run `make timing` (or at least the transient target) first."
        )

    lines = [
        "# Antarctica timing matrix",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "Only completed 5-year transient runs (`2015`–`2020`, `dt=1.0`) "
        "with `timing_kind=matrix` are included. Qualification probes are "
        "reported by their gate and excluded here.",
    ]

    groups = {}
    for record in records:
        tag = record.get("timing_tag", "legacy-untagged")
        groups.setdefault(tag, []).append(record)

    for tag, group in sorted(groups.items()):
        sample = group[0]
        mode = sample.get("diagnostic_solver_mode", sample["linear_solver"])
        lines.extend(["", f"## `{tag}`", "", f"Diagnostic solver: `{mode}`", ""])
        lines.extend(render_table(group))
        lines.extend(["", "### Per-step timing", ""])
        lines.extend(render_table(group, per_step=True))

    lines.append("")
    with open(OUT_FN, "w") as f:
        f.write("\n".join(lines))
    print(f"Wrote {OUT_FN} ({len(records)} completed timing records)")


if __name__ == "__main__":
    main()
