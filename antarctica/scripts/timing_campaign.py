#!/usr/bin/env python3
"""Shared policy and validation for the cached transient timing campaign.

This module deliberately has no Firedrake dependency.  The Makefile helpers,
Slurm wrappers, report builder, and unit tests all import the same lane and
record policy so the submission matrix cannot drift away from the report.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.fspath(Path(__file__).resolve().parents[2]))

from icepack2_tools.runconfig import TARGET_MESH_GEOMETRY_METHOD


RECORD_SCHEMA_VERSION = 3
CACHE_SCHEMA_VERSION = 3
CACHE_ROLE = "timing-initial-state"
CACHE_REQUIRED_FIELDS = (
    "log_friction",
    "log_fluidity",
    "velocity_obs",
    "thickness",
    "bed",
    "surface",
    "fluidity_prior",
    "velocity",
    "membrane_stress",
    "basal_stress",
    "H_init",
    "phi_eff",
    "C_w0",
    "N_ref",
)

LCS = (500, 1000, 2000, 2500, 5000)
RATIOS = (10, 20)
DISPLAY_CORES = (16, 32, 64)
CORES_BY_LC = {
    500: (32, 64),
    1000: (16, 32, 64),
    2000: (16, 32),
    2500: (16, 32),
    5000: (16,),
}

SOLVER_MODE = "scpc_mumps"
SOURCE_TAG = "dg0_logvelnet"
SOURCE_INVERSION_BASENAME = (
    "inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5"
)
CAMPAIGN_TAG = (
    "scpc_mumps_5step_dt0p25at2500_dg0_logvelnet_cached_strict_v3"
)
AMB_PROBE_TAG = f"{CAMPAIGN_TAG}_ambdiv_probe"
CACHE_TAG = "scpc_mumps_dg0_logvelnet_v3"

MATRIX_T_START = 2015.0
MATRIX_STEPS = 5
MATRIX_DT_2500 = 0.25
MATRIX_REFERENCE_LC = 2500.0
BUFFER_M = 20000
MASS_RESIDUAL_TOL_GT = 5.0e-5

MEMORY_BY_LC = {
    500: "240G",
    1000: "96G",
    2000: "80G",
    2500: "64G",
    5000: "32G",
}


def mesh_rows(lcs=LCS, ratios=RATIOS):
    """Return unique ``(lc, lc_coarse)`` rows in display order."""
    return tuple(sorted(
        {(int(lc), int(lc) * int(ratio))
         for lc in lcs for ratio in ratios},
        key=lambda item: (item[0], item[1]),
    ))


def planned_lanes():
    """Return the selected 20 ``(lc, lc_coarse, ncores)`` lanes."""
    return tuple(
        (lc, lc_coarse, ncores)
        for lc, lc_coarse in mesh_rows()
        for ncores in CORES_BY_LC[lc]
    )


def scout_lanes():
    """Return one lowest-retained-core scout for every target mesh."""
    return tuple(
        (lc, lc_coarse, CORES_BY_LC[lc][0])
        for lc, lc_coarse in mesh_rows()
    )


def scaling_lanes():
    scouts = set(scout_lanes())
    return tuple(lane for lane in planned_lanes() if lane not in scouts)


def expected_dt(lc):
    return MATRIX_DT_2500 * float(lc) / MATRIX_REFERENCE_LC


def expected_t_end(lc):
    return MATRIX_T_START + MATRIX_STEPS * expected_dt(lc)


def mesh_basename(lc, lc_coarse, buffer_m=BUFFER_M):
    return f"antarctica_{int(lc_coarse)}_{int(lc)}_buffered{int(buffer_m)}.msh"


def cache_stem(lc, lc_coarse, buffer_m=BUFFER_M):
    return (
        f"initial_state_{CACHE_TAG}_{int(lc)}_{int(lc_coarse)}"
        f"_buffered{int(buffer_m)}"
    )


def cache_paths(cache_dir, lc, lc_coarse, buffer_m=BUFFER_M):
    stem = cache_stem(lc, lc_coarse, buffer_m)
    root = Path(cache_dir)
    return root / f"{stem}.h5", root / f"{stem}.json"


def timing_record_basename(tag, lc, lc_coarse, ncores):
    return f"timing_{tag}_{int(lc)}_{int(lc_coarse)}_{int(ncores)}.json"


def timing_status_basename(tag, lc, lc_coarse, ncores):
    return f"status_{tag}_{int(lc)}_{int(lc_coarse)}_{int(ncores)}.txt"


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def solver_configuration_fingerprint(configuration):
    """Fingerprint setup-relevant solver settings, excluding recovery policy.

    Cache construction deliberately enables adaptive continuation while strict
    transient lanes disable rescue/subcycling.  Those recovery settings must
    differ, so the cache identity covers the nonlinear/linear operators and
    tolerances that define the prepared state, not the transient recovery
    policy.
    """
    selected = {
        "diagnostic_mode": configuration.get("diagnostic_mode"),
        "diagnostic_petsc_options": configuration.get(
            "diagnostic_petsc_options"
        ),
        "transport_petsc_options": configuration.get(
            "transport_petsc_options"
        ),
        "snes_atol_policy": configuration.get("snes_atol_policy"),
    }
    encoded = json.dumps(
        selected, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def atomic_write_json(path, payload):
    """Write JSON by rename so interrupted jobs never publish half a record."""
    path = os.fspath(path)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + ".tmp.",
                               dir=directory)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise


def atomic_write_status(path, state, **fields):
    tokens = [state]
    for key, value in fields.items():
        if value is None or value == "":
            continue
        text = str(value).replace(" ", "_").replace("\n", "_")
        tokens.append(f"{key}={text}")
    path = os.fspath(path)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + ".tmp.",
                               dir=directory)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(" ".join(tokens) + "\n")
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise


def reason_diverged(reason):
    text = str(reason)
    if text.startswith("DIVERGED"):
        return True
    try:
        return int(text) < 0
    except ValueError:
        return False


def diverged_reasons(record, summary_name):
    reasons = record.get(summary_name, {}).get("reason_counts", {})
    return {
        str(reason): count
        for reason, count in reasons.items()
        if reason_diverged(reason)
    }


def validate_cache_manifest(
    manifest,
    *,
    lc,
    lc_coarse,
    cache_path=None,
    source_sha256=None,
    mesh_sha256=None,
    solver_fingerprint=None,
):
    """Return ``(valid, detail)`` for a target-mesh initial-state cache."""
    expected = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_role": CACHE_ROLE,
        "lc": int(lc),
        "lc_coarse": int(lc_coarse),
        "buffer_m": BUFFER_M,
        "diagnostic_solver_mode": SOLVER_MODE,
        "friction": "budd",
        "geometry_space": "dg0",
        "n_flow": 3.0,
        "a4_factor": 1.0,
        "t_yr": MATRIX_T_START,
        "mesh_basename": mesh_basename(lc, lc_coarse),
        "source_inversion_basename": SOURCE_INVERSION_BASENAME,
        "geometry_source_method": TARGET_MESH_GEOMETRY_METHOD,
    }
    for key, value in expected.items():
        actual = manifest.get(key)
        if isinstance(value, float):
            try:
                matches = math.isclose(float(actual), value, abs_tol=1e-12)
            except (TypeError, ValueError):
                matches = False
        else:
            matches = actual == value
        if not matches:
            return False, f"cache {key}={actual!r}; expected {value!r}"
    for key in ("source_inversion_sha256", "source_mesh_sha256"):
        if not manifest.get(key):
            return False, f"cache {key} is missing"
    for key in ("geometry_source", "geometry_source_basename"):
        if not manifest.get(key):
            return False, f"cache {key} is missing"
    if os.path.basename(os.fspath(manifest["geometry_source"])) != manifest[
        "geometry_source_basename"
    ]:
        return False, "cache geometry-source basename is inconsistent"
    checkpoint_fields = manifest.get("checkpoint_fields")
    if not isinstance(checkpoint_fields, list):
        return False, "cache checkpoint_fields is missing"
    missing_fields = sorted(
        set(CACHE_REQUIRED_FIELDS) - set(checkpoint_fields)
    )
    if missing_fields:
        return False, (
            "cache checkpoint is incomplete; missing fields: "
            + ", ".join(missing_fields)
        )
    if cache_path is not None:
        actual = os.path.realpath(os.fspath(manifest.get("cache_path", "")))
        expected_path = os.path.realpath(os.fspath(cache_path))
        if actual != expected_path:
            return False, f"cache path={actual!r}; expected {expected_path!r}"
        if not os.path.isfile(expected_path):
            return False, f"cache file is missing: {expected_path}"
    if source_sha256 is not None \
            and manifest.get("source_inversion_sha256") != source_sha256:
        return False, "cache source-inversion checksum is stale"
    if mesh_sha256 is not None \
            and manifest.get("source_mesh_sha256") != mesh_sha256:
        return False, "cache source-mesh checksum is stale"
    if solver_fingerprint is not None and manifest.get(
        "solver_configuration_fingerprint"
    ) != solver_fingerprint:
        return False, "cache solver configuration is stale"
    return True, "cache provenance matches"


def validate_timing_record(
    record,
    *,
    lc=None,
    lc_coarse=None,
    ncores=None,
    require_success=True,
    timing_kind="matrix",
    timing_tag=CAMPAIGN_TAG,
    apparent_mb_mode=None,
):
    """Validate the strict five-step contract and return ``(valid, detail)``."""
    if record.get("record_schema_version") != RECORD_SCHEMA_VERSION:
        return False, "record schema is not the cached-strict schema"
    if require_success and record.get("run_status") != "success":
        failure = record.get("failure") or {}
        detail = failure.get("category", record.get("run_status", "failed"))
        phase = failure.get("phase")
        return False, detail + (f" at {phase}" if phase else "")
    if record.get("timing_kind") != timing_kind:
        return False, f"timing_kind={record.get('timing_kind')!r}"
    if record.get("timing_tag") != timing_tag:
        return False, f"timing_tag={record.get('timing_tag')!r}"
    if apparent_mb_mode is not None:
        if record.get("apparent_mb_mode") != apparent_mb_mode:
            return False, (
                "apparent_mb_mode="
                f"{record.get('apparent_mb_mode')!r}"
            )
        try:
            uncapped = math.isclose(
                float(record["apparent_mb_cap_m_per_yr"]),
                0.0,
                abs_tol=1e-12,
            )
        except (KeyError, TypeError, ValueError):
            uncapped = False
        if not uncapped:
            return False, "apparent-MB probe was capped"
    elif timing_kind == "matrix" and record.get("apparent_mb_mode") is not None:
        return False, "matrix timing unexpectedly used apparent MB"
    if record.get("diagnostic_solver_mode") != SOLVER_MODE:
        return False, (
            "diagnostic_solver_mode="
            f"{record.get('diagnostic_solver_mode')!r}"
        )

    checks = (("lc", lc), ("lc_coarse", lc_coarse), ("ncores", ncores))
    for key, expected in checks:
        if expected is not None and record.get(key) != int(expected):
            return False, f"record {key}={record.get(key)!r}; expected {expected}"

    try:
        record_lc = int(record["lc"])
        interval_ok = (
            math.isclose(float(record["t_start"]), MATRIX_T_START,
                         abs_tol=1e-12)
            and math.isclose(float(record["dt"]), expected_dt(record_lc),
                             abs_tol=1e-12)
            and int(record["nsteps"]) == MATRIX_STEPS
            and int(record["completed_steps"]) == MATRIX_STEPS
            and math.isclose(float(record["t_end"]), expected_t_end(record_lc),
                             abs_tol=1e-12)
            and math.isclose(float(record["t_final"]), expected_t_end(record_lc),
                             abs_tol=1e-12)
        )
    except (KeyError, TypeError, ValueError):
        interval_ok = False
    if not interval_ok:
        return False, "record did not complete the required five-step interval"

    diagnostic = diverged_reasons(record, "diagnostic_solve_summary")
    if diagnostic:
        return False, f"diagnostic divergence: {diagnostic}"
    labels = [stat.get("label", "")
              for stat in record.get("diagnostic_solves", [])]
    if len(labels) < MATRIX_STEPS:
        return False, "record has fewer than five diagnostic solves"
    if any("rescue" in label or "trust-region" in label for label in labels):
        return False, "record used the rescue ladder"
    if any(label != f"step-{index}-direct"
           for index, label in enumerate(labels, 1)):
        return False, f"unexpected diagnostic solve sequence: {labels}"

    transport = diverged_reasons(record, "transport_solve_summary")
    if transport:
        return False, f"transport divergence: {transport}"
    transport_summary = record.get("transport_solve_summary", {})
    if int(transport_summary.get("count", 0)) != MATRIX_STEPS:
        return False, "record does not contain exactly five transport solves"
    transport_labels = [
        stat.get("label", "") for stat in record.get("transport_solves", [])
    ]
    expected_transport_labels = [
        f"step-{index}-substep-1/1"
        for index in range(1, MATRIX_STEPS + 1)
    ]
    if transport_labels != expected_transport_labels:
        return False, f"unexpected transport solve sequence: {transport_labels}"
    residual = transport_summary.get("mass_residual_gt_max")
    try:
        residual_ok = (
            residual is not None
            and math.isfinite(float(residual))
            and float(residual) <= MASS_RESIDUAL_TOL_GT
        )
    except (TypeError, ValueError):
        residual_ok = False
    if not residual_ok:
        return False, f"transport mass residual {residual!r} exceeds tolerance"
    step_residual = record.get("step_mass_residual_gt_max")
    try:
        step_residual_ok = (
            step_residual is not None
            and math.isfinite(float(step_residual))
            and float(step_residual) <= MASS_RESIDUAL_TOL_GT
        )
    except (TypeError, ValueError):
        step_residual_ok = False
    if not step_residual_ok:
        return False, f"step mass residual {step_residual!r} exceeds tolerance"

    cache = record.get("cache_validation", {})
    if cache.get("status") != "valid":
        return False, f"cache validation={cache.get('status')!r}"
    if record.get("rescue_enabled") is not False:
        return False, "strict timing record did not disable rescue"
    if record.get("solver_configuration", {}).get("subcycles") != [1]:
        return False, "strict timing record did not restrict subcycles to [1]"
    if record.get("timing_scope") != "transient_loop_only":
        return False, f"timing_scope={record.get('timing_scope')!r}"
    if timing_kind == "matrix":
        return True, "completed strict cached five-step timing run"
    return True, "completed strict cached five-step run"


assert len(mesh_rows()) == 10
assert len(scout_lanes()) == 10
assert len(scaling_lanes()) == 10
assert len(planned_lanes()) == 20


def _print_lanes(lanes):
    for lane in lanes:
        print(" ".join(str(value) for value in lane))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("rows", "planned", "scout", "scale"):
        subparsers.add_parser(command)
    validate = subparsers.add_parser("validate-cache")
    validate.add_argument("--manifest", required=True)
    validate.add_argument("--cache", required=True)
    validate.add_argument("--lc", required=True, type=int)
    validate.add_argument("--lc-coarse", required=True, type=int)
    validate.add_argument("--source-sha256")
    validate.add_argument("--mesh-sha256")
    validate.add_argument("--solver-fingerprint")
    args = parser.parse_args()

    if args.command == "rows":
        _print_lanes(mesh_rows())
    elif args.command == "planned":
        _print_lanes(planned_lanes())
    elif args.command == "scout":
        _print_lanes(scout_lanes())
    elif args.command == "scale":
        _print_lanes(scaling_lanes())
    elif args.command == "validate-cache":
        try:
            with open(args.manifest) as stream:
                manifest = json.load(stream)
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"INVALID: {exc}") from exc
        valid, detail = validate_cache_manifest(
            manifest,
            lc=args.lc,
            lc_coarse=args.lc_coarse,
            cache_path=args.cache,
            source_sha256=args.source_sha256,
            mesh_sha256=args.mesh_sha256,
            solver_fingerprint=args.solver_fingerprint,
        )
        if not valid:
            raise SystemExit(f"INVALID: {detail}")
        print(f"VALID: {detail}")


if __name__ == "__main__":
    main()
