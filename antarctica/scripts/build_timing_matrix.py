#!/usr/bin/env python3
"""Render accepted timings and explicit failure states for one campaign."""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from timing_campaign import (
    CAMPAIGN_TAG,
    DISPLAY_CORES,
    MASS_RESIDUAL_TOL_GT,
    MATRIX_DT_2500,
    MATRIX_REFERENCE_LC,
    MATRIX_STEPS,
    MATRIX_T_START,
    diverged_reasons,
    mesh_rows,
    planned_lanes,
    timing_status_basename,
    validate_timing_record,
)

_ROOT = Path(__file__).resolve().parents[1]


def _load_records(timing_dir, tag):
    records = []
    errors = []
    pattern = timing_dir / f"timing_{glob.escape(tag)}_*.json"
    for path_text in sorted(glob.glob(os.fspath(pattern))):
        path = Path(path_text)
        try:
            with open(path) as stream:
                records.append(json.load(stream))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"{path.name}: {exc}")
    return records, errors


def _read_status(timing_dir, tag, lane):
    path = timing_dir / timing_status_basename(tag, *lane)
    try:
        tokens = path.read_text().strip().split()
    except OSError:
        return None
    if not tokens:
        return {"state": "invalid", "path": os.fspath(path)}
    status = {"state": tokens[0], "path": os.fspath(path)}
    for token in tokens[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            status[key] = value
    return status


def _legacy_valid(record):
    try:
        lc = int(record["lc"])
        dt = MATRIX_DT_2500 * lc / MATRIX_REFERENCE_LC
        t_end = MATRIX_T_START + MATRIX_STEPS * dt
        return (
            record.get("timing_kind", "matrix") == "matrix"
            and math.isclose(float(record["t_start"]), MATRIX_T_START)
            and math.isclose(float(record["dt"]), dt)
            and int(record["nsteps"]) == MATRIX_STEPS
            and int(record["completed_steps"]) == MATRIX_STEPS
            and math.isclose(float(record["t_end"]), t_end)
            and math.isclose(float(record["t_final"]), t_end)
            and not diverged_reasons(record, "diagnostic_solve_summary")
            and not diverged_reasons(record, "transport_solve_summary")
        )
    except (KeyError, TypeError, ValueError):
        return False


def _failure_class(category):
    if category in {
        "diagnostic_convergence",
        "transport_convergence",
        "transport_mass_budget",
        "step_mass_budget",
    }:
        return "NUMERICAL FAILURE"
    if category in {"external_termination", "oom", "slurm_oom"}:
        return "OOM / EXTERNAL"
    if category in {"incomplete_output", "stale_running"}:
        return "INCOMPLETE OUTPUT"
    return "FAILED"


def _stale_running(status, now=None):
    if status.get("state") != "running" or "timestamp" not in status:
        return False
    try:
        stamp = datetime.fromisoformat(
            status["timestamp"].replace("Z", "+00:00")
        )
    except ValueError:
        return False
    now = now or datetime.now(timezone.utc)
    return (now - stamp).total_seconds() > 24 * 3600


def _legacy_log_evidence(status):
    job_id = status.get("job_id")
    if not job_id:
        return None
    paths = [_ROOT / f"ant_timing_{job_id}{suffix}"
             for suffix in (".txt", ".err")]
    text = ""
    for path in paths:
        try:
            text += path.read_text(errors="replace")
        except OSError:
            pass
    if "Transport mass residual" in text or "Step mass residual" in text:
        return "NUMERICAL FAILURE", "mass-budget exception in archived Slurm log"
    if "DIVERGED_" in text or "rescue-continuation" in text:
        return "NUMERICAL FAILURE", "solver divergence/rescue in archived Slurm log"
    if any(path.is_file() for path in paths):
        return "INCOMPLETE OUTPUT", "archived Slurm log ended without a record"
    return None


def _classify(record, status, lane, strict):
    if record is not None:
        if strict:
            valid, detail = validate_timing_record(
                record, lc=lane[0], lc_coarse=lane[1], ncores=lane[2]
            )
        else:
            valid = _legacy_valid(record)
            detail = "completed legacy five-step timing run"
        if valid:
            return "OK", detail
        failure = record.get("failure") or {}
        category = failure.get("category")
        if category:
            phase = failure.get("phase", "unknown phase")
            return _failure_class(category), f"{category} at {phase}"
        diagnostic = diverged_reasons(record, "diagnostic_solve_summary")
        transport = diverged_reasons(record, "transport_solve_summary")
        if diagnostic or transport:
            return (
                "NUMERICAL FAILURE",
                f"diagnostic={diagnostic or '{}'}, transport={transport or '{}'}",
            )
        completed = record.get("completed_steps", "?")
        requested = record.get("nsteps", "?")
        return "INCOMPLETE OUTPUT", f"{completed}/{requested} steps; {detail}"

    if status is None:
        return "NOT RUN", "no record or status stamp"
    state = status["state"]
    if state == "blocked_by_scout":
        return "BLOCKED BY SCOUT", "scout did not pass"
    if state == "failed":
        if not strict:
            evidence = _legacy_log_evidence(status)
            if evidence is not None and evidence[0] == "NUMERICAL FAILURE":
                return evidence
        category = status.get("category")
        exit_code = status.get("exit_code")
        if category == "external_termination" or exit_code in {"137", "9"}:
            return "OOM / EXTERNAL", (
                f"external termination, exit_code={exit_code or 'unknown'}"
            )
        if category:
            return _failure_class(category), category
        return "FAILED", f"failed, exit_code={exit_code or 'unknown'}"
    if state == "finished":
        return "INCOMPLETE OUTPUT", "finished stamp but no readable record"
    if state in {"running", "submitted", "submitting"}:
        if not strict:
            job_id = status.get("job_id")
            if job_id and any(
                (_ROOT / f"ant_timing_{job_id}{suffix}").is_file()
                for suffix in (".txt", ".err")
            ):
                return (
                    "INCOMPLETE OUTPUT",
                    "archived job output exists but no completion record",
                )
        if _stale_running(status):
            return (
                "INCOMPLETE OUTPUT",
                "stale running stamp older than 24 hours and no record",
            )
        return state.upper(), f"{state}; job_id={status.get('job_id', 'unknown')}"
    if state == "not_runnable":
        return "NOT RUNNABLE", status.get("reason", "missing prerequisite")
    if state == "submission_failed":
        return "SUBMIT FAILED", "sbatch submission failed"
    return "UNKNOWN", state


def _status_table(rows, classifications):
    header = (
        "| LC (m) | LC_coarse (m) | dt (yr) | "
        + " | ".join(f"{cores} cores" for cores in DISPLAY_CORES)
        + " |"
    )
    lines = [header, "|" + "|".join(["---"] * 6) + "|"]
    for lc, lc_coarse in rows:
        labels = [
            classifications[(lc, lc_coarse, cores)][0]
            for cores in DISPLAY_CORES
        ]
        dt = MATRIX_DT_2500 * lc / MATRIX_REFERENCE_LC
        lines.append(
            f"| {lc} | {lc_coarse} | {dt:.3g} | "
            + " | ".join(labels)
            + " |"
        )
    return lines


def _timing_table(rows, record_index, per_step=False):
    key = "seconds_per_step" if per_step else "run_seconds"
    unit = "s/step" if per_step else "s"
    digits = 2 if per_step else 1
    header = (
        "| LC (m) | LC_coarse (m) | Vertices | Cells | "
        + " | ".join(f"{cores} cores ({unit})" for cores in DISPLAY_CORES)
        + " |"
    )
    lines = [header, "|" + "|".join(["---"] * 7) + "|"]
    for lc, lc_coarse in rows:
        candidates = [
            record_index.get((lc, lc_coarse, cores))
            for cores in DISPLAY_CORES
        ]
        sample = next((record for record in candidates if record), None)
        vertices = sample.get("vertices", "—") if sample else "—"
        cells = sample.get("cells", "—") if sample else "—"
        values = [
            f"{record[key]:.{digits}f}" if record is not None else "—"
            for record in candidates
        ]
        lines.append(
            f"| {lc} | {lc_coarse} | {vertices} | {cells} | "
            + " | ".join(values)
            + " |"
        )
    return lines


def render(tag, timing_dir, output, legacy_full_matrix=False):
    rows = mesh_rows()
    configured = (
        {(lc, lc_coarse, cores) for lc, lc_coarse in rows
         for cores in DISPLAY_CORES}
        if legacy_full_matrix
        else set(planned_lanes())
    )
    displayed = {
        (lc, lc_coarse, cores)
        for lc, lc_coarse in rows
        for cores in DISPLAY_CORES
    }
    records, load_errors = _load_records(timing_dir, tag)
    all_index = {
        (record.get("lc"), record.get("lc_coarse"), record.get("ncores")): record
        for record in records
        if record.get("timing_tag") == tag
    }
    classifications = {}
    accepted = {}
    strict = not legacy_full_matrix and tag == CAMPAIGN_TAG
    for lane in sorted(displayed):
        if lane not in configured:
            classifications[lane] = (
                "NOT PLANNED",
                "excluded by the selected 20-lane policy",
            )
            continue
        status = _read_status(timing_dir, tag, lane)
        classification = _classify(all_index.get(lane), status, lane, strict)
        classifications[lane] = classification
        if classification[0] == "OK":
            accepted[lane] = all_index[lane]

    counts = {}
    for lane in configured:
        label = classifications[lane][0]
        counts[label] = counts.get(label, 0) + 1
    policy_text = (
        "This is the archived cold-start 30-lane campaign. Its generic exit "
        "stamps are supplemented with the copied Slurm logs so numerical "
        "failures remain distinguishable from incomplete output."
        if legacy_full_matrix
        else
        "The primary timing covers only the five-step transient loop; setup "
        "and checkpoint loading are recorded separately. The timestep is "
        "`0.25 × LC / 2500` years. Strict lanes use `scpc_mumps`, an "
        "exact-mesh prepared cache, and no rescue or subcycle recovery."
    )
    lines = [
        "# Antarctica timing matrix",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        f"Campaign tag: `{tag}`. Configured lanes: {len(configured)}; "
        f"accepted: {counts.get('OK', 0)}.",
        "",
        policy_text,
        "",
        "## Run status",
        "",
    ]
    lines.extend(_status_table(rows, classifications))
    non_ok = [
        (lane, classifications[lane])
        for lane in sorted(configured)
        if classifications[lane][0] != "OK"
    ]
    if non_ok or load_errors:
        lines.extend(["", "### Status notes", ""])
        for (lc, lc_coarse, cores), (label, detail) in non_ok:
            lines.append(
                f"- LC={lc}, LC_coarse={lc_coarse}, {cores} cores: "
                f"**{label}** — {detail}."
            )
        for error in load_errors:
            lines.append(f"- Unreadable timing record: `{error}`")
    lines.extend(["", "## Successful timings", ""])
    if accepted:
        lines.extend(_timing_table(rows, accepted))
        lines.extend(["", "### Per-step timing", ""])
        lines.extend(_timing_table(rows, accepted, per_step=True))
    else:
        lines.append("No accepted timing records are available yet.")
    lines.extend([
        "",
        f"Scout mass-residual acceptance limit: {MASS_RESIDUAL_TOL_GT:g} Gt.",
        "",
    ])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines))
    print(
        f"Wrote {output} ({len(accepted)} accepted of "
        f"{len(configured)} configured timing records)"
    )
    return classifications


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--timing-dir", type=Path, default=_ROOT / "results/timing")
    parser.add_argument("--output", type=Path, default=_ROOT / "TIMING_MATRIX.md")
    parser.add_argument("--legacy-full-matrix", action="store_true")
    args = parser.parse_args()
    nested = _ROOT / "antarctica/results"
    if nested.is_dir() and any(path.is_file() for path in nested.rglob("*")):
        raise SystemExit(
            f"Refusing to build matrix while nested results tree exists: {nested}. "
            "Run `make sync-results` first."
        )
    render(
        args.tag,
        args.timing_dir.resolve(),
        args.output.resolve(),
        legacy_full_matrix=args.legacy_full_matrix,
    )


if __name__ == "__main__":
    main()
