#!/usr/bin/env python3
r"""Fail-loud gate for the transient solver qualification probes.

The timing driver records requested and completed steps separately because an
exhausted rescue ladder saves and returns cleanly.  This checker requires the
whole interval, the intended number of CSV rows, no diverged diagnostic solve,
and a mass-budget residual that rounds to zero in the persisted CSV.
"""

import argparse
import csv
import json
import math


def _is_diverged(reason):
    if reason.startswith("DIVERGED"):
        return True
    try:
        return int(reason) < 0
    except ValueError:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--min-steps", type=int, required=True)
    # CSV stores four decimal places, so 5e-5 accepts only rows persisted as
    # 0.0000 (including signed zero), not a small but nonzero rounded residual.
    parser.add_argument("--residual-tol-gt", type=float, default=5e-5)
    args = parser.parse_args()

    with open(args.record) as stream:
        record = json.load(stream)
    if record.get("timing_kind") != "qualification":
        raise SystemExit(
            f"FAIL: {args.record} is not a qualification record "
            f"(timing_kind={record.get('timing_kind')!r})"
        )
    if record.get("nsteps", 0) < args.min_steps:
        raise SystemExit(
            f"FAIL: record requested {record.get('nsteps', 0)} steps; "
            f"need at least {args.min_steps}"
        )
    if record.get("completed_steps", 0) < args.min_steps:
        raise SystemExit(
            f"FAIL: record completed {record.get('completed_steps', 0)} steps; "
            f"need at least {args.min_steps}"
        )

    summary = record.get("diagnostic_solve_summary", {})
    reasons = summary.get("reason_counts", {})
    diverged = {name: count for name, count in reasons.items()
                if _is_diverged(name)}
    if diverged:
        raise SystemExit(f"FAIL: diagnostic solve divergence recorded: {diverged}")

    transport_summary = record.get("transport_solve_summary", {})
    transport_reasons = transport_summary.get("reason_counts", {})
    transport_diverged = {
        name: count for name, count in transport_reasons.items()
        if _is_diverged(name)
    }
    if transport_diverged:
        raise SystemExit(
            f"FAIL: transport KSP divergence recorded: {transport_diverged}"
        )
    if transport_summary.get("count", 0) < args.min_steps:
        raise SystemExit(
            "FAIL: record contains only "
            f"{transport_summary.get('count', 0)} transport solves; "
            f"need at least {args.min_steps}"
        )
    transport_mass_resid = transport_summary.get("mass_residual_gt_max")
    if (
        transport_mass_resid is None
        or not math.isfinite(transport_mass_resid)
        or transport_mass_resid > args.residual_tol_gt
    ):
        raise SystemExit(
            "FAIL: max |transport mass residual| is "
            f"{transport_mass_resid!r} Gt; limit is "
            f"{args.residual_tol_gt:.6g} Gt"
        )

    with open(args.csv, newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < args.min_steps:
        raise SystemExit(
            f"FAIL: timeseries contains {len(rows)} completed steps; "
            f"need at least {args.min_steps}"
        )

    last_year = float(rows[-1]["year"])
    expected_year = float(record["t_end"])
    if not math.isclose(last_year, expected_year, abs_tol=0.51 * record["dt"]):
        raise SystemExit(
            f"FAIL: timeseries stops at {last_year}, expected {expected_year}"
        )

    residuals = [abs(float(row["resid_gt"])) for row in rows]
    max_residual = max(residuals, default=math.inf)
    if not math.isfinite(max_residual) or max_residual > args.residual_tol_gt:
        raise SystemExit(
            f"FAIL: max |mass residual| is {max_residual:.6g} Gt; "
            f"limit is {args.residual_tol_gt:.6g} Gt"
        )

    mode = record.get("diagnostic_solver_mode", "unknown")
    print(
        f"PASS: {mode}: {len(rows)} transient steps, "
        f"{summary.get('count', 0)} diagnostic solves, "
        f"max |transport/step mass residual|="
        f"{transport_mass_resid:.3g}/{max_residual:.3g} Gt"
    )


if __name__ == "__main__":
    main()
