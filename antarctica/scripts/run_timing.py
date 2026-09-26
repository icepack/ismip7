#!/usr/bin/env python3
"""Run one transient timing lane and always publish a structured record.

Matrix lanes must restart from an exact-mesh prepared cache. Setup and
operator construction are measured separately; the primary timing begins at
the marker immediately before the transient loop in ``simulation``.
"""

from __future__ import annotations

import json
import os
import re
import sys
from collections import Counter
from pathlib import Path
from time import perf_counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.fspath(Path(__file__).resolve().parents[2]))

from firedrake import COMM_WORLD
from firedrake.petsc import PETSc
from mpi4py import MPI

from icepack2_tools.runconfig import buffer_m as _buffer_m, mesh_override
from icepack2_tools.solverconfig import (
    diagnostic_solver_label,
    diagnostic_solver_mode,
    fixed_front_enabled,
    rescue_enabled,
    solver_provenance,
)
from simulation import lc, run_simulation, setup_model
from timing_campaign import (
    CACHE_SOLVER_MODE,
    CONTRACTS,
    LANE_SOLVER_MODES,
    MAP_CHECK_KIND,
    MAP_CHECK_STEPS,
    RECORD_SCHEMA_VERSION,
    host_provenance,
    parse_campaign_tag,
    atomic_write_json,
    atomic_write_status,
    solver_configuration_fingerprint,
    validate_cache_manifest,
    validate_map_check_manifest,
)
from icepack2_tools.runconfig import friction as _friction

_ROOT = Path(__file__).resolve().parents[1]
# A MAP check keeps its records beside its own MAP (ISMIP7_TIMING_DIR); the
# campaign's stay in results/timing.
TIMING_DIR = Path(os.environ.get("ISMIP7_TIMING_DIR") or _ROOT / "results" / "timing")

T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2020"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "5"))
DIAGNOSTIC_LINEAR_SOLVER = diagnostic_solver_mode()
LINEAR_SOLVER_LABEL = diagnostic_solver_label(DIAGNOSTIC_LINEAR_SOLVER)
# Taken at import, before anything compiles: whether this lane starts from an
# empty kernel cache is part of how its step time should be read.
HOST_PROVENANCE = host_provenance()
TIMING_KIND = os.environ.get("ISMIP7_TIMING_KIND", "matrix")
APPARENT_MB_MODE = os.environ.get("ISMIP7_APPARENT_MB")
APPARENT_MB_CAP = float(os.environ.get("ISMIP7_AMB_CAP", "0"))
FIXED_FRONT = fixed_front_enabled()
# Filled by _validate_lane_contract() for matrix/probe kinds.
LANE_CONTRACT = {}


def _safe_tag(env_name, default):
    raw = os.environ.get(env_name, default)
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", raw).strip("-.")
    if not value:
        raise ValueError(f"{env_name} must contain a filename-safe character")
    return value


TIMING_TAG = _safe_tag("ISMIP7_TIMING_TAG", DIAGNOSTIC_LINEAR_SOLVER)
EXPERIMENT_NAME = _safe_tag("ISMIP7_TIMING_EXPERIMENT", "timing")
RESTART_FROM = os.environ.get("ISMIP7_RESTART") or None
CACHE_MANIFEST = os.environ.get("ISMIP7_TIMING_CACHE_MANIFEST") or None
STATUS_PATH = os.environ.get("ISMIP7_TIMING_STATUS") or None


def _validate_map_check_contract():
    """A map-check lane runs the strict contract under a campaign solver, on
    the MAP's own mesh (cold start) or on a prepared map-check cache."""
    if not TIMING_TAG.startswith("mapcheck_"):
        raise RuntimeError(
            f"timing tag {TIMING_TAG!r} does not name a map check"
        )
    if DIAGNOSTIC_LINEAR_SOLVER not in LANE_SOLVER_MODES:
        raise RuntimeError(
            f"a map check times one of {LANE_SOLVER_MODES}, not "
            f"{DIAGNOSTIC_LINEAR_SOLVER!r}"
        )
    contract = CONTRACTS["strict"]
    if APPARENT_MB_MODE != contract["apparent_mb_mode"]:
        raise RuntimeError(
            f"a map check runs strict: ISMIP7_APPARENT_MB must be unset, got "
            f"{APPARENT_MB_MODE!r}"
        )
    if abs(APPARENT_MB_CAP - contract["apparent_mb_cap_m_per_yr"]) > 1e-12:
        raise RuntimeError(
            f"a map check runs strict: ISMIP7_AMB_CAP must be 0, got {APPARENT_MB_CAP:g}"
        )
    if FIXED_FRONT != contract["fixed_front"]:
        raise RuntimeError(
            "a map check runs strict: ISMIP7_FIXED_FRONT must be 0/unset"
        )
    steps = int(os.environ.get("ISMIP7_MAP_CHECK_STEPS", str(MAP_CHECK_STEPS)))
    LANE_CONTRACT.update(contract, name="strict", steps=steps, dt_2500=None)
    return contract


def _validate_lane_contract():
    """The tag names the contract; the environment must match it exactly."""
    if TIMING_KIND == MAP_CHECK_KIND:
        return _validate_map_check_contract()
    if TIMING_KIND not in {"matrix", "cache_probe"}:
        return None
    try:
        spec = parse_campaign_tag(TIMING_TAG)
    except ValueError as exc:
        raise RuntimeError(
            f"timing tag {TIMING_TAG!r} does not name a lane contract: {exc}"
        ) from exc
    is_probe = spec["lane"] == "probe"
    if is_probe != (TIMING_KIND == "cache_probe"):
        raise RuntimeError(
            f"timing tag {TIMING_TAG!r} is a "
            f"{'probe' if is_probe else 'matrix'} tag but ISMIP7_TIMING_KIND="
            f"{TIMING_KIND!r}"
        )
    if spec["solver"] != DIAGNOSTIC_LINEAR_SOLVER:
        raise RuntimeError(
            f"timing tag {TIMING_TAG!r} names solver {spec['solver']!r} but "
            f"ISMIP7_DIAGNOSTIC_LINEAR_SOLVER={DIAGNOSTIC_LINEAR_SOLVER!r}"
        )
    contract = CONTRACTS[spec["contract"]]
    if APPARENT_MB_MODE != contract["apparent_mb_mode"]:
        raise RuntimeError(
            f"contract {spec['contract']!r} requires ISMIP7_APPARENT_MB="
            f"{contract['apparent_mb_mode']!r}, got {APPARENT_MB_MODE!r}"
        )
    if abs(APPARENT_MB_CAP - contract["apparent_mb_cap_m_per_yr"]) > 1e-12:
        raise RuntimeError(
            f"contract {spec['contract']!r} requires ISMIP7_AMB_CAP="
            f"{contract['apparent_mb_cap_m_per_yr']:g}, got {APPARENT_MB_CAP:g}"
        )
    if FIXED_FRONT != contract["fixed_front"]:
        raise RuntimeError(
            f"contract {spec['contract']!r} requires ISMIP7_FIXED_FRONT="
            f"{'1' if contract['fixed_front'] else '0/unset'}"
        )
    LANE_CONTRACT.update(contract, name=spec["contract"], steps=spec["steps"],
                         dt_2500=spec["dt_2500"])
    return contract


def _budget_summary(steps):
    """Totals over the per-step mass budget published by run_simulation."""
    if not steps:
        return None
    def total(key):
        return float(sum(float(step.get(key, 0.0) or 0.0) for step in steps))
    return {
        "steps": steps,
        "calving_gt_total": total("calv_gt"),
        "clamp_gt_total": total("clamp_gt"),
        "limit_gt_total": total("limit_gt"),
        "outflux_gt_total": total("out_gt"),
        "apparent_mb_gt_per_yr": float(steps[0].get("amb_gt_per_yr", 0.0) or 0.0),
        "resid_gt_max": max(abs(float(step.get("resid_gt", 0.0) or 0.0)) for step in steps),
    }


def _summary(stats, reason_key, iteration_key, residual_key=None):
    summary = {
        "count": len(stats),
        "reason_counts": dict(Counter(str(stat[reason_key]) for stat in stats)),
        f"{iteration_key}_total": sum(int(stat[iteration_key]) for stat in stats),
        f"{iteration_key}_max": max(
            (int(stat[iteration_key]) for stat in stats), default=0
        ),
        "seconds_total": sum(float(stat.get("seconds", 0.0)) for stat in stats),
    }
    if residual_key is not None:
        residuals = [
            abs(float(stat[residual_key]))
            for stat in stats
            if stat.get(residual_key) is not None
        ]
        summary[f"{residual_key}_max"] = max(residuals, default=None)
    return summary


def _diagnostic_summary(stats):
    summary = _summary(stats, "snes_reason", "snes_iterations")
    summary.update({
        "linear_iterations_total": sum(
            int(stat["linear_iterations"]) for stat in stats
        ),
        "linear_iterations_max": max(
            (int(stat["linear_iterations"]) for stat in stats), default=0
        ),
    })
    # SCPC's own count of the work on the condensed system. It includes the
    # line search's solves, which ``linear_iterations`` (SNES's) does not.
    for key in ("condensed_solves", "condensed_iterations"):
        if stats and all(key in stat for stat in stats):
            summary[f"{key}_total"] = sum(int(stat[key]) for stat in stats)
    return summary


def _transport_summary(stats):
    summary = _summary(
        stats, "ksp_reason", "ksp_iterations", "mass_residual_gt"
    )
    summary["residual_norm_max"] = max(
        (float(stat["residual_norm"]) for stat in stats), default=0.0
    )
    return summary


def _normal(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "item"):
        return value.item()
    return value


def _load_and_validate_map_check_cache():
    """A map-check lane restarts from its prepared cache (the transferred
    state) or cold-starts from the MAP on the MAP's own mesh."""
    if mesh_override():
        raise RuntimeError(
            "A map-check lane solves on the mesh its MAP or cache carries; "
            "ISMIP7_MESH must be unset or `checkpoint`"
        )
    source_basename = os.path.basename(os.environ.get("ISMIP7_INVERSION", ""))
    if not RESTART_FROM:
        return {
            "status": "not_required",
            "detail": "cold start from the MAP on its own mesh",
            "source_inversion_basename": source_basename,
        }
    if not CACHE_MANIFEST:
        raise RuntimeError(
            "A map-check lane restarting from a cache needs "
            "ISMIP7_TIMING_CACHE_MANIFEST"
        )
    with open(CACHE_MANIFEST) as stream:
        manifest = json.load(stream)
    fingerprint = solver_configuration_fingerprint(
        solver_provenance(CACHE_SOLVER_MODE)
    )
    valid, detail = validate_map_check_manifest(
        manifest,
        lc=lc,
        lc_coarse=int(os.environ["ISMIP7_LC_COARSE"]),
        buffer_m=int(round(_buffer_m())),
        friction=_friction(),
        source_basename=source_basename,
        mesh_name=manifest.get("mesh_basename"),
        cache_path=RESTART_FROM,
        solver_fingerprint=fingerprint,
    )
    if not valid:
        raise RuntimeError(detail)
    return {
        "status": "valid",
        "detail": detail,
        "cache_path": str(Path(RESTART_FROM).resolve()),
        "manifest_path": str(Path(CACHE_MANIFEST).resolve()),
        "cache_schema_version": manifest["cache_schema_version"],
        "source_inversion": manifest["source_inversion"],
        "source_inversion_basename": manifest["source_inversion_basename"],
        "source_inversion_sha256": manifest["source_inversion_sha256"],
        "source_mesh_sha256": manifest["source_mesh_sha256"],
        "geometry_source": manifest["geometry_source"],
        "geometry_source_basename": manifest["geometry_source_basename"],
        "geometry_source_method": manifest["geometry_source_method"],
        "friction_gate": manifest.get("friction_gate"),
        "transfer_fill": manifest.get("transfer_fill"),
        "solver_configuration_fingerprint": fingerprint,
        "manifest": manifest,
    }


def _load_and_validate_cache():
    if TIMING_KIND == MAP_CHECK_KIND:
        return _load_and_validate_map_check_cache()
    if TIMING_KIND not in {"matrix", "cache_probe"}:
        return {"status": "not_required"}
    if mesh_override():
        raise RuntimeError(
            "Cached timing must load the mesh embedded in its cache; "
            "ISMIP7_MESH must be unset or `checkpoint`"
        )
    if not RESTART_FROM or not CACHE_MANIFEST:
        raise RuntimeError(
            "Cached timing requires ISMIP7_RESTART and "
            "ISMIP7_TIMING_CACHE_MANIFEST"
        )
    with open(CACHE_MANIFEST) as stream:
        manifest = json.load(stream)
    lc_coarse = int(os.environ["ISMIP7_LC_COARSE"])
    # The cache is judged against the solver it was PREPARED under, in this
    # lane's environment; the solver the lane times is its tag's, and every
    # campaign solver starts from the same prepared state.
    fingerprint = solver_configuration_fingerprint(
        solver_provenance(CACHE_SOLVER_MODE)
    )
    valid, detail = validate_cache_manifest(
        manifest,
        lc=lc,
        lc_coarse=lc_coarse,
        cache_path=RESTART_FROM,
        solver_fingerprint=fingerprint,
    )
    if not valid:
        raise RuntimeError(detail)
    return {
        "status": "valid",
        "detail": detail,
        "cache_path": str(Path(RESTART_FROM).resolve()),
        "manifest_path": str(Path(CACHE_MANIFEST).resolve()),
        "cache_schema_version": manifest["cache_schema_version"],
        "source_inversion": manifest["source_inversion"],
        "source_inversion_basename": manifest["source_inversion_basename"],
        "source_inversion_sha256": manifest["source_inversion_sha256"],
        "source_mesh_sha256": manifest["source_mesh_sha256"],
        "geometry_source": manifest["geometry_source"],
        "geometry_source_basename": manifest["geometry_source_basename"],
        "geometry_source_method": manifest["geometry_source_method"],
        "friction_gate": manifest["friction_gate"],
        "solver_configuration_fingerprint": fingerprint,
        "manifest": manifest,
    }


def _validate_loaded_cache(ctx, cache_validation):
    if TIMING_KIND not in {"matrix", "cache_probe", MAP_CHECK_KIND}:
        return
    if "manifest" not in cache_validation:
        return  # a map-check cold start has no cache to hold to
    manifest = cache_validation["manifest"]
    attrs = {
        key: _normal(value)
        for key, value in ctx.get("checkpoint_metadata", {}).items()
    }
    expected = {
        "timing_cache_schema_version": manifest["cache_schema_version"],
        "timing_cache_role": manifest["cache_role"],
        "source_inversion": manifest["source_inversion"],
        "source_inversion_sha256": manifest["source_inversion_sha256"],
        "source_mesh_sha256": manifest["source_mesh_sha256"],
        "geometry_source": manifest["geometry_source"],
        "geometry_source_method": manifest["geometry_source_method"],
        "diagnostic_solver_mode": manifest["diagnostic_solver_mode"],
        "solver_configuration_fingerprint": manifest[
            "solver_configuration_fingerprint"
        ],
        "geometry_space": manifest["geometry_space"],
        "n_flow": manifest["n_flow"],
        "a4_factor": manifest["a4_factor"],
        "friction_gate": manifest.get("friction_gate"),
    }
    for key, expected_value in expected.items():
        actual = attrs.get(key)
        if isinstance(expected_value, float):
            try:
                matches = abs(float(actual) - expected_value) <= 1.0e-12
            except (TypeError, ValueError):
                matches = False
        else:
            matches = actual == expected_value
        if not matches:
            raise RuntimeError(
                f"loaded cache {key}={actual!r}; expected {expected_value!r}"
            )
    context_expected = {
        "lc": manifest["lc"],
        "lc_coarse": manifest["lc_coarse"],
        "buffer_m": manifest["buffer_m"],
        "t_yr": manifest["t_yr"],
        "mesh_basename": manifest["mesh_basename"],
        "friction": manifest["friction"],
        "geometry_space": manifest["geometry_space"],
    }
    context_actual = {
        "lc": ctx.get("lc"),
        "lc_coarse": ctx.get("lc_coarse"),
        "buffer_m": int(round(float(ctx.get("buffer_m", -1)))),
        "t_yr": ctx.get("t_restart"),
        "mesh_basename": ctx.get("mesh_basename"),
        "friction": ctx.get("friction"),
        "geometry_space": "dg0" if ctx.get("geom_dg") else "cg1",
    }
    for key, expected_value in context_expected.items():
        if context_actual[key] != expected_value:
            raise RuntimeError(
                f"loaded cache {key}={context_actual[key]!r}; "
                f"expected {expected_value!r}"
            )


def _write_record(record, target_lc_coarse, ncores):
    if COMM_WORLD.rank == 0:
        out_fn = TIMING_DIR / (
            f"timing_{TIMING_TAG}_{lc}_{target_lc_coarse}_{ncores}.json"
        )
        atomic_write_json(out_fn, record)
        failure = record.get("failure") or {}
        if STATUS_PATH:
            if record["run_status"] == "success":
                atomic_write_status(
                    STATUS_PATH,
                    "finished",
                    exit_code=0,
                    completed_steps=record["completed_steps"],
                )
            else:
                atomic_write_status(
                    STATUS_PATH,
                    "failed",
                    category=failure.get("category", "unknown_exception"),
                    phase=failure.get("phase", "unknown"),
                    completed_steps=record["completed_steps"],
                )
        PETSc.Sys.Print(f"Timing record -> {out_fn}")
    # Do not let a non-root rank re-raise and trigger launcher cleanup while
    # rank 0 is still atomically publishing the catchable-failure record.
    COMM_WORLD.barrier()


def main():
    _validate_lane_contract()
    TIMING_DIR.mkdir(parents=True, exist_ok=True)
    ncores = COMM_WORLD.size
    target_lc_coarse = int(os.environ.get("ISMIP7_LC_COARSE", "0"))
    target_buffer_m = _buffer_m()
    overall_t0 = perf_counter()
    configuration = solver_provenance()
    ctx = None
    results = []
    cache_validation = {"status": "not_checked"}
    caught = None
    effective_t_start = T_START
    global_vertices = None
    global_cells = None
    activity = "cache_validation"

    PETSc.Sys.Print(
        f"Timing run: lc={lc} lc_coarse={target_lc_coarse} "
        f"ncores={ncores} solver={LINEAR_SOLVER_LABEL}/transport-gmres "
        f"t={T_START}->{T_END} dt={DT} restart={RESTART_FROM or 'none'}"
    )

    try:
        cache_validation = _load_and_validate_cache()
        activity = "setup"
        ctx = setup_model(
            restart_from=RESTART_FROM,
            # A validated pristine cache may build its apparent-MB reference
            # under any contract that asks for one (production closure).
            allow_timing_cache_a_ref=(
                LANE_CONTRACT.get("apparent_mb_mode") == "div"
            ),
        )
        activity = "loaded_cache_validation"
        try:
            _validate_loaded_cache(ctx, cache_validation)
        except Exception as exc:
            cache_validation["status"] = "invalid"
            cache_validation["detail"] = str(exc)
            raise
        activity = "setup"
        mesh = ctx["mesh"]
        target_lc_coarse = int(ctx["lc_coarse"])
        target_buffer_m = float(ctx["buffer_m"])
        restart_time = ctx.get("t_restart")
        if RESTART_FROM is not None and restart_time is None:
            raise RuntimeError(
                f"Timing restart {RESTART_FROM} has no t_yr attribute"
            )
        effective_t_start = T_START if restart_time is None else restart_time
        global_vertices = int(ctx["Q"].dim())
        global_cells = int(mesh.comm.allreduce(mesh.cell_set.size, op=MPI.SUM))
        PETSc.Sys.Print(
            f"Timing compute mesh: {global_vertices} global vertices, "
            f"{global_cells} global cells"
        )
        activity = "transient_setup_or_loop"
        results = run_simulation(
            ctx,
            experiment_name=EXPERIMENT_NAME,
            t_start=T_START,
            t_end=T_END,
            dt=DT,
            output_interval=OUTPUT_INTERVAL,
            checkpoint_interval=10**9,
            forcing_callback=None,
        )
        requested_steps = int(round((T_END - effective_t_start) / DT))
        if len(results) != requested_steps:
            ctx["failure"] = {
                "category": "incomplete_output",
                "phase": "transient_loop",
                "exception_type": "RuntimeError",
                "message": f"completed {len(results)} of {requested_steps} steps",
            }
            raise RuntimeError(ctx["failure"]["message"])
    except Exception as exc:
        caught = exc
        if ctx is None:
            failure = {
                "category": (
                    "cache_invalid"
                    if activity == "cache_validation"
                    else "setup_failure"
                ),
                "phase": activity,
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
        else:
            failure = ctx.get("failure") or {
                "category": (
                    "cache_invalid"
                    if activity == "loaded_cache_validation"
                    else (
                        "unknown_exception"
                        if "transient_t0" in ctx
                        else "setup_failure"
                    )
                ),
                "phase": (
                    "transient_loop" if "transient_t0" in ctx else activity
                ),
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
    else:
        failure = None

    # run_simulation publishes its live list into ctx before entering the
    # timestep loop.  If a later diagnostic/transport solve raises, its return
    # assignment above never executes, so recover that list here rather than
    # falsely recording zero completed steps.
    if ctx is not None:
        results = ctx.get("results", results)

    now = perf_counter()
    transient_t0 = ctx.get("transient_t0") if ctx is not None else None
    setup_seconds = (
        transient_t0 - overall_t0
        if transient_t0 is not None
        else now - overall_t0
    )
    transient_seconds = (
        ctx.get("transient_seconds", now - transient_t0)
        if transient_t0 is not None
        else 0.0
    )
    requested_steps = max(0, int(round((T_END - effective_t_start) / DT)))
    solver_stats = ctx.get("solver_stats", []) if ctx is not None else []
    transport_stats = ctx.get("transport_stats", []) if ctx is not None else []
    final_year = results[-1][0] if results else effective_t_start
    step_mass_residuals = [abs(float(row[8])) for row in results]

    record = {
        "record_schema_version": RECORD_SCHEMA_VERSION,
        "run_status": "failed" if failure else "success",
        "failure": failure,
        "lc": int(lc),
        "lc_coarse": target_lc_coarse,
        "buffer_m": target_buffer_m,
        "mesh_basename": ctx.get("mesh_basename", "") if ctx else "",
        "mesh_input": os.environ.get("ISMIP7_MESH", ""),
        "boundary_ids_input": os.environ.get("ISMIP7_BNDIDS", ""),
        "inversion_input": os.environ.get("ISMIP7_INVERSION", ""),
        "restart_input": RESTART_FROM or "",
        # Which MAP the cache's mixed state descends from: the per-mesh
        # short invert, or the transferred 2.5 km source for meshes that are
        # not re-inverted. The matrix report keeps the two apart.
        "initial_state_source": cache_validation.get(
            "source_inversion_basename", ""
        ),
        "cache_validation": {
            key: value
            for key, value in cache_validation.items()
            if key != "manifest"
        },
        "ncores": ncores,
        "host": HOST_PROVENANCE,
        "vertices": global_vertices,
        "cells": global_cells,
        "t_start": effective_t_start,
        "t_end": T_END,
        "dt": DT,
        "nsteps": requested_steps,
        "completed_steps": len(results),
        "t_final": final_year,
        "timing_kind": TIMING_KIND,
        "timing_tag": TIMING_TAG,
        "experiment_name": EXPERIMENT_NAME,
        "apparent_mb_mode": APPARENT_MB_MODE,
        "apparent_mb_cap_m_per_yr": APPARENT_MB_CAP,
        "fixed_front": FIXED_FRONT,
        "front_hmin_m": ctx.get("front_hmin") if ctx else None,
        "fixed_front_cells": ctx.get("fixed_front_cells") if ctx else None,
        "contract": LANE_CONTRACT.get("name"),
        "matrix_steps": LANE_CONTRACT.get("steps"),
        "matrix_dt_2500": LANE_CONTRACT.get("dt_2500"),
        "tripwire": ctx.get("tripwire") if ctx else None,
        "budget": _budget_summary(ctx.get("step_budget", []) if ctx else []),
        "diagnostic_solver_mode": DIAGNOSTIC_LINEAR_SOLVER,
        "linear_solver": LINEAR_SOLVER_LABEL,
        "transport_solver": "gmres-bjacobi-ilu",
        "solver_configuration": configuration,
        "rescue_enabled": rescue_enabled(),
        "timing_scope": "transient_loop_only",
        "setup_seconds": setup_seconds,
        "transient_seconds": transient_seconds,
        "run_seconds": transient_seconds,
        "seconds_per_step": transient_seconds / max(len(results), 1),
        "diagnostic_solve_summary": _diagnostic_summary(solver_stats),
        "diagnostic_solves": solver_stats,
        "transport_solve_summary": _transport_summary(transport_stats),
        "transport_solves": transport_stats,
        "step_mass_residual_gt_max": max(step_mass_residuals, default=None),
        "field_extrema": ctx.get("field_stats", []) if ctx else [],
    }
    if TIMING_KIND == MAP_CHECK_KIND:
        record["map_check"] = {
            "friction": ctx.get("friction") if ctx else _friction(),
            "mesh_role": "transferred" if RESTART_FROM else "native",
            "transfer_fill": (
                cache_validation.get("transfer_fill")
                if RESTART_FROM else (ctx.get("transfer_fill") if ctx else None)
            ),
        }
    _write_record(record, target_lc_coarse, ncores)

    if caught is not None:
        raise caught
    PETSc.Sys.Print(
        f"Timing: setup={setup_seconds:.1f}s, "
        f"transient={transient_seconds:.1f}s "
        f"({record['seconds_per_step']:.2f}s/step)"
    )


if __name__ == "__main__":
    main()
