#!/usr/bin/env python3
"""
Aggregate one campaign's per-run timing records into TIMING_MATRIX.md.

Reads results/timing/timing_<tag>_*.json and writes partial timing/status tables.
Qualification probes and incomplete/non-five-step records are excluded.

Usage:
    python scripts/build_timing_matrix.py --tag <timing-tag>
"""

import argparse
import glob
import json
import math
import os
from datetime import datetime, timezone

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIMING_DIR = os.path.join(_ROOT, "results", "timing")
OUT_FN = os.path.join(_ROOT, "TIMING_MATRIX.md")
MATRIX_T_START = 2015.0
MATRIX_STEPS = 5
MATRIX_DT_2500 = 0.25
MATRIX_REFERENCE_LC = 2500.0
STATUS_LABELS = {
    "submitted": "QUEUED",
    "running": "RUNNING",
    "finished": "NO RECORD",
    "failed": "FAILED",
    "submission_failed": "SUBMIT FAIL",
    "not_runnable": "NOT RUNNABLE",
}


def load_records(tag=None):
    records = []
    errors = []
    name = "timing_*.json" if tag is None else f"timing_{glob.escape(tag)}_*.json"
    for path in sorted(glob.glob(os.path.join(TIMING_DIR, name))):
        try:
            with open(path) as f:
                records.append(json.load(f))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"{os.path.basename(path)}: {exc}")
    return records, errors


def _parse_int_list(value):
    try:
        values = tuple(int(item) for item in value.split(":") if item)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if not values:
        raise argparse.ArgumentTypeError("expected a colon-separated integer list")
    return values


def _reason_diverged(reason):
    text = str(reason)
    if text.startswith("DIVERGED"):
        return True
    try:
        return int(text) < 0
    except ValueError:
        return False


def _diverged_reasons(record, summary_name):
    reasons = record.get(summary_name, {}).get("reason_counts", {})
    return {reason: count for reason, count in reasons.items()
            if _reason_diverged(reason)}


def _is_matrix_record(record):
    try:
        dt = MATRIX_DT_2500 * float(record["lc"]) / MATRIX_REFERENCE_LC
        t_end = MATRIX_T_START + MATRIX_STEPS * dt
        return (
            record.get("timing_kind", "matrix") == "matrix"
            and math.isclose(
                float(record["t_start"]), MATRIX_T_START, abs_tol=1e-12
            )
            and math.isclose(float(record["dt"]), dt, abs_tol=1e-12)
            and record.get("nsteps") == MATRIX_STEPS
            and record.get("completed_steps") == MATRIX_STEPS
            and math.isclose(float(record["t_end"]), t_end, abs_tol=1e-12)
            and math.isclose(float(record["t_final"]), t_end, abs_tol=1e-12)
            and not _diverged_reasons(record, "diagnostic_solve_summary")
            and not _diverged_reasons(record, "transport_solve_summary")
        )
    except (KeyError, TypeError, ValueError):
        return False


def matrix_records(records):
    return [record for record in records if _is_matrix_record(record)]


def render_table(records, rows, core_counts, per_step=False):
    index = {
        (record["lc"], record["lc_coarse"], record["ncores"]): record
        for record in records
    }
    unit = "s/step" if per_step else "s"
    key = "seconds_per_step" if per_step else "run_seconds"
    digits = 2 if per_step else 1
    header = (
        "| LC (m) | LC_coarse (m) | dt (yr) | Vertices | Cells | "
        + " | ".join(f"{count} cores ({unit})" for count in core_counts)
        + " |"
    )
    lines = [header, "|" + "|".join(["---"] * (5 + len(core_counts))) + "|"]
    for lc, lc_coarse in rows:
        candidates = [
            index.get((lc, lc_coarse, count)) for count in core_counts
        ]
        sample = next((record for record in candidates if record), None)
        dt = f"{sample['dt']:.3g}" if sample else "—"
        vertices = sample["vertices"] if sample else "—"
        cells = sample["cells"] if sample else "—"
        timings = [
            f"{record[key]:.{digits}f}" if record else "—"
            for record in candidates
        ]
        lines.append(
            f"| {lc} | {lc_coarse} | {dt} | {vertices} | {cells} | "
            + " | ".join(timings)
            + " |"
        )
    return lines


def _read_status(tag, lc, lc_coarse, ncores):
    path = os.path.join(
        TIMING_DIR, f"status_{tag}_{lc}_{lc_coarse}_{ncores}.txt"
    )
    try:
        with open(path) as stream:
            tokens = stream.read().strip().split()
    except OSError:
        return None
    if not tokens:
        return {"state": "invalid", "path": path}
    fields = {"state": tokens[0], "path": path}
    for token in tokens[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return fields


def _record_issue(record):
    requested = record.get("nsteps", "?")
    completed = record.get("completed_steps", "?")
    if completed != MATRIX_STEPS or requested != MATRIX_STEPS:
        return f"incomplete record ({completed}/{requested} steps)"
    diagnostic = _diverged_reasons(record, "diagnostic_solve_summary")
    transport = _diverged_reasons(record, "transport_solve_summary")
    if diagnostic or transport:
        return (
            "completed only after a failed solve/rescue "
            f"(diagnostic={diagnostic or '{}'}, transport={transport or '{}'})"
        )
    return (
        "record does not match the current five-step resolution-scaled "
        "timestep policy"
    )


def _classify(record, status):
    if record is not None and _is_matrix_record(record):
        return "OK", "completed and accepted"
    if record is not None:
        issue = _record_issue(record)
        label = "INCOMPLETE" if "incomplete" in issue else "EXCLUDED"
        return label, issue
    if status is None:
        return "NOT RUN", "no status stamp or timing record"
    state = status["state"]
    label = STATUS_LABELS.get(state, "UNKNOWN")
    details = [state.replace("_", " ")]
    for key in ("job_id", "exit_code", "timestamp"):
        if key in status:
            details.append(f"{key}={status[key]}")
    if state == "running":
        details.append("or terminated before writing a completion stamp")
    return label, ", ".join(details)


def render_status_table(rows, core_counts, classifications):
    header = (
        "| LC (m) | LC_coarse (m) | dt (yr) | "
        + " | ".join(f"{count} cores" for count in core_counts)
        + " |"
    )
    lines = [header, "|" + "|".join(["---"] * (3 + len(core_counts))) + "|"]
    for lc, lc_coarse in rows:
        dt = MATRIX_DT_2500 * lc / MATRIX_REFERENCE_LC
        labels = [
            classifications[(lc, lc_coarse, ncores)][0]
            for ncores in core_counts
        ]
        lines.append(
            f"| {lc} | {lc_coarse} | {dt:.3g} | "
            + " | ".join(labels)
            + " |"
        )
    return lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", required=True)
    parser.add_argument("--lcs", type=_parse_int_list,
                        default=_parse_int_list("500:1000:2000:2500:5000"))
    parser.add_argument("--ratios", type=_parse_int_list,
                        default=_parse_int_list("20:10"))
    parser.add_argument("--cores", type=_parse_int_list,
                        default=_parse_int_list("16:32:64"))
    args = parser.parse_args()

    all_records, load_errors = load_records(args.tag)
    tag_records = [
        record for record in all_records
        if record.get("timing_tag") == args.tag
    ]
    record_index = {
        (record.get("lc"), record.get("lc_coarse"), record.get("ncores")): record
        for record in tag_records
    }
    rows = sorted(
        {(lc, ratio * lc) for lc in args.lcs for ratio in args.ratios},
        key=lambda value: (value[0], value[1]),
    )
    core_counts = sorted(args.cores)
    classifications = {}
    for lc, lc_coarse in rows:
        for ncores in core_counts:
            key = (lc, lc_coarse, ncores)
            status = _read_status(args.tag, *key)
            classifications[key] = _classify(record_index.get(key), status)

    expected_keys = set(classifications)
    records = [
        record for record in matrix_records(tag_records)
        if (record["lc"], record["lc_coarse"], record["ncores"])
        in expected_keys
    ]
    counts = {}
    for label, _ in classifications.values():
        counts[label] = counts.get(label, 0) + 1

    lines = [
        "# Antarctica timing matrix",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        f"Campaign tag: `{args.tag}`. Expected lanes: {len(classifications)}; "
        f"accepted: {counts.get('OK', 0)}.",
        "",
        "Only completed five-step transient runs with `timing_kind=matrix` "
        "are included. The timestep scales with fine resolution as "
        "`dt = 0.25 × LC / 2500` years. Qualification probes are reported "
        "by their gate and excluded here.",
        "",
        "## Run status",
        "",
    ]
    lines.extend(render_status_table(rows, core_counts, classifications))
    lines.extend([
        "",
        "`OK` records are included below. `EXCLUDED` finished only through a "
        "failed solve/rescue or used the wrong run policy. `QUEUED` means the "
        "job was submitted but its wrapper did not start; it can also mean it "
        "was canceled before starting. `RUNNING` can also mean Slurm terminated "
        "the job before it wrote its final stamp.",
    ])

    non_ok = [
        (key, value) for key, value in classifications.items()
        if value[0] != "OK"
    ]
    if non_ok or load_errors:
        lines.extend(["", "### Status notes", ""])
        for (lc, lc_coarse, ncores), (label, detail) in non_ok:
            lines.append(
                f"- LC={lc}, LC_coarse={lc_coarse}, {ncores} cores: "
                f"**{label}** — {detail}."
            )
        for error in load_errors:
            lines.append(f"- Unreadable timing record: `{error}`")

    lines.extend(["", "## Successful timings", ""])
    if records:
        sample = records[0]
        mode = (
            sample.get("diagnostic_solver_mode")
            or sample.get("linear_solver", "unknown")
        )
        lines.extend([f"Diagnostic solver: `{mode}`", ""])
        lines.extend(render_table(records, rows, core_counts))
        lines.extend(["", "### Per-step timing", ""])
        lines.extend(render_table(records, rows, core_counts, per_step=True))
    else:
        lines.append("No accepted timing records are available yet.")

    lines.append("")
    with open(OUT_FN, "w") as f:
        f.write("\n".join(lines))
    print(
        f"Wrote {OUT_FN} ({len(records)} accepted of "
        f"{len(classifications)} expected timing records)"
    )


if __name__ == "__main__":
    main()
