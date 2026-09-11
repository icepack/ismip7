#!/usr/bin/env python3
"""Advance the cached timing campaign by one idempotent Slurm stage."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.fspath(Path(__file__).resolve().parents[2]))

from icepack2_tools.solverconfig import (
    SNES_DIVERGENCE_TOL_DEFAULT,
    solver_provenance,
)
from timing_campaign import (
    BUFFER_M,
    CAMPAIGN_TAG,
    CACHE_TAG,
    MEMORY_BY_LC,
    SOLVER_MODE,
    atomic_write_status,
    cache_paths,
    expected_dt,
    expected_t_end,
    mesh_basename,
    mesh_rows,
    scaling_lanes,
    scout_lanes,
    sha256_file,
    solver_configuration_fingerprint,
    timing_record_basename,
    timing_status_basename,
    validate_cache_manifest,
    validate_timing_record,
)

_ROOT = Path(__file__).resolve().parents[1]
_LIVE_SLURM_STATES = {
    "CONFIGURING",
    "COMPLETING",
    "PENDING",
    "REQUEUED",
    "REQUEUE_FED",
    "REQUEUE_HOLD",
    "RESIZING",
    "RUNNING",
    "STAGE_OUT",
    "SUSPENDED",
}
_TERMINAL_SLURM_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "SPECIAL_EXIT",
    "STOPPED",
    "TIMEOUT",
}


def read_status(path):
    try:
        tokens = Path(path).read_text().strip().split()
    except OSError:
        return None
    if not tokens:
        return {"state": "invalid"}
    parsed = {"state": tokens[0]}
    for token in tokens[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            parsed[key] = value
    return parsed


def read_record(path):
    try:
        with open(path) as stream:
            return json.load(stream), None
    except FileNotFoundError:
        return None, None
    except (OSError, json.JSONDecodeError) as exc:
        return None, str(exc)


def _slurm_outcome_from_outputs(queue_output, accounting_output):
    """Classify one Slurm allocation from machine-readable command output."""
    if queue_output.strip():
        return "active", None, None
    rows = [line for line in accounting_output.splitlines() if line.strip()]
    if not rows:
        return "unknown", None, None
    fields = rows[0].split("|")
    raw_state = fields[0].strip()
    state = raw_state.split()[0].rstrip("+").upper() if raw_state else ""
    exit_code = fields[1].strip() if len(fields) > 1 else ""
    if state in _LIVE_SLURM_STATES:
        return "active", state, exit_code
    if state in _TERMINAL_SLURM_STATES:
        return "terminal", state, exit_code
    return "unknown", state or None, exit_code or None


def _slurm_outcome(status):
    """Return ``(kind, Slurm state, exit code)`` for an active-looking stamp.

    Unknown scheduler state is treated conservatively by callers: never submit
    a duplicate merely because accounting is delayed or unavailable.
    """
    if not status or status.get("state") not in {
        "submitting", "submitted", "running"
    }:
        return "not_active_stamp", None, None
    job_id = status.get("job_id")
    if not job_id or shutil.which("squeue") is None:
        return "unknown", None, None
    queue = subprocess.run(
        ["squeue", "-h", "-j", job_id, "-o", "%i"],
        check=False,
        capture_output=True,
        text=True,
    )
    if queue.returncode == 0 and queue.stdout.strip():
        return "active", None, None
    # A completed job may already be purged from squeue, which can be reported
    # as either an empty successful query or an "invalid job id" error.  In
    # both cases accounting, not the squeue return code, is authoritative.
    if shutil.which("sacct") is None:
        return "unknown", None, None
    accounting = subprocess.run(
        [
            "sacct", "-n", "-X", "-j", job_id,
            "-o", "State,ExitCode", "-P",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if accounting.returncode != 0:
        return "unknown", None, None
    return _slurm_outcome_from_outputs(queue.stdout, accounting.stdout)


def _timestamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def scaling_decision(scout_state):
    if scout_state == "passed":
        return "submit"
    if scout_state == "failed":
        return "block"
    return "wait"


def _parse_mesh(value):
    parts = value.replace(":", "/").split("/")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            "mesh must be written as LC/LC_coarse, for example 2500/25000"
        )
    try:
        pair = tuple(int(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("mesh values must be integers") from exc
    if pair not in mesh_rows():
        choices = ", ".join(f"{lc}/{lcc}" for lc, lcc in mesh_rows())
        raise argparse.ArgumentTypeError(
            f"mesh {value!r} is not in the timing policy; choose {choices}"
        )
    return pair


class CampaignManager:
    def __init__(self, args):
        self.root = Path(args.root).resolve()
        self.mesh_dir = self.root / "mesh"
        self.cache_dir = Path(args.cache_dir).resolve()
        self.timing_dir = Path(args.timing_dir).resolve()
        self.logs_dir = self.root / "results" / "logs"
        self.inversion = Path(args.inversion).resolve()
        self.partition = args.partition
        self.walltime = args.walltime
        self.dry_run = args.dry_run
        self.force = args.force
        self.assume_valid_caches = args.assume_valid_caches
        self.only_mesh = args.only_mesh
        self.monitor = args.monitor
        self.submit_failures = 0
        self.source_sha256 = None
        self.mesh_checksums = {}
        self.solver_configuration = solver_provenance()
        self.solver_fingerprint = solver_configuration_fingerprint(
            self.solver_configuration
        )
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timing_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def selected_rows(self):
        rows = mesh_rows()
        if self.only_mesh is None:
            return rows
        return tuple(row for row in rows if row == self.only_mesh)

    def selected_lanes(self, lanes):
        if self.only_mesh is None:
            return lanes
        return tuple(lane for lane in lanes if lane[:2] == self.only_mesh)

    def reconcile_status(self, status_path):
        """Turn a dead Slurm allocation's live-looking stamp into a failure."""
        status = read_status(status_path)
        kind, slurm_state, exit_code = _slurm_outcome(status)
        if kind in {"active", "unknown"}:
            return status, True
        if kind != "terminal":
            return status, False
        category = (
            "incomplete_output"
            if slurm_state == "COMPLETED"
            else "external_termination"
        )
        reconciled = {
            "state": "failed",
            "category": category,
            "phase": "scheduler",
            "slurm_state": slurm_state,
            "exit_code": exit_code,
            "job_id": status.get("job_id"),
            "timestamp": _timestamp(),
        }
        if not self.dry_run:
            atomic_write_status(
                status_path,
                reconciled.pop("state"),
                **reconciled,
            )
            reconciled["state"] = "failed"
        return reconciled, False

    def source_checksum(self):
        if self.source_sha256 is None:
            if not self.inversion.is_file():
                raise FileNotFoundError(self.inversion)
            self.source_sha256 = sha256_file(self.inversion)
        return self.source_sha256

    def mesh_path(self, lc, lc_coarse):
        return self.mesh_dir / mesh_basename(lc, lc_coarse)

    def mesh_checksum(self, lc, lc_coarse):
        key = (lc, lc_coarse)
        if key not in self.mesh_checksums:
            path = self.mesh_path(lc, lc_coarse)
            if not path.is_file():
                raise FileNotFoundError(path)
            self.mesh_checksums[key] = sha256_file(path)
        return self.mesh_checksums[key]

    def boundary_path(self, lc, lc_coarse):
        return self.mesh_dir / (
            f"boundary_ids_antarctica_{lc_coarse}_{lc}_buffered{BUFFER_M}.json"
        )

    def cache_status_path(self, lc, lc_coarse):
        return self.timing_dir / (
            f"status_cache_{CACHE_TAG}_{lc}_{lc_coarse}.txt"
        )

    def cache_validation(self, lc, lc_coarse):
        if self.dry_run and self.assume_valid_caches:
            return True, "assumed valid for dry-run command inspection"
        cache, manifest_path = cache_paths(
            self.cache_dir, lc, lc_coarse
        )
        try:
            with open(manifest_path) as stream:
                manifest = json.load(stream)
        except (OSError, json.JSONDecodeError) as exc:
            return False, f"cache manifest unavailable: {exc}"
        try:
            source_sha256 = self.source_checksum()
            mesh_sha256 = self.mesh_checksum(lc, lc_coarse)
        except FileNotFoundError as exc:
            return False, f"cache source unavailable: {exc}"
        return validate_cache_manifest(
            manifest,
            lc=lc,
            lc_coarse=lc_coarse,
            cache_path=cache,
            source_sha256=source_sha256,
            mesh_sha256=mesh_sha256,
            solver_fingerprint=self.solver_fingerprint,
        )

    def _submit(self, job_name, ncores, memory, exports, script, status_path):
        command = [
            "sbatch",
            "--parsable",
            f"--job-name={job_name}",
            f"--ntasks-per-node={ncores}",
            f"--partition={self.partition}",
            f"--time={self.walltime}",
            f"--mem={memory}",
            "--export=ALL," + ",".join(
                f"{key}={value}" for key, value in exports.items()
            ),
            os.fspath(script),
        ]
        print("DRY RUN:" if self.dry_run else "SUBMIT:", shlex.join(command))
        if self.dry_run:
            return
        if shutil.which("sbatch") is None:
            raise RuntimeError("sbatch is unavailable; use --dry-run off-cluster")
        atomic_write_status(status_path, "submitting", timestamp=_timestamp())
        result = subprocess.run(
            command, check=False, capture_output=True, text=True
        )
        if result.returncode != 0:
            atomic_write_status(
                status_path,
                "submission_failed",
                exit_code=result.returncode,
                timestamp=_timestamp(),
            )
            self.submit_failures += 1
            print(result.stderr.strip(), file=sys.stderr)
            return
        job_id = result.stdout.strip().split(";", 1)[0]
        atomic_write_status(
            status_path, "submitted", job_id=job_id, timestamp=_timestamp()
        )

    def prepare(self):
        try:
            source_sha256 = self.source_checksum()
        except FileNotFoundError as exc:
            raise SystemExit(f"Cannot prepare caches: {exc}") from exc
        print(f"Source inversion sha256: {source_sha256}")
        for lc, lc_coarse in self.selected_rows():
            valid, detail = self.cache_validation(lc, lc_coarse)
            if valid and not self.force:
                print(f"CACHE OK {lc}/{lc_coarse}: {detail}")
                continue
            status_path = self.cache_status_path(lc, lc_coarse)
            status, active = self.reconcile_status(status_path)
            if not self.force and active:
                print(f"CACHE ACTIVE {lc}/{lc_coarse}: {status['state']}")
                continue
            if not self.force and status and status.get("state") == "failed":
                print(
                    f"CACHE FAILED {lc}/{lc_coarse}: retry with "
                    "FORCE_TIMING=1 after inspection"
                )
                continue
            mesh = self.mesh_path(lc, lc_coarse)
            boundary = self.boundary_path(lc, lc_coarse)
            missing = [path for path in (mesh, boundary) if not path.is_file()]
            if missing:
                message = ",".join(os.fspath(path) for path in missing)
                print(f"CACHE NOT RUNNABLE {lc}/{lc_coarse}: missing {message}")
                if not self.dry_run:
                    atomic_write_status(
                        status_path, "not_runnable", reason="missing_input"
                    )
                continue
            cache, manifest = cache_paths(self.cache_dir, lc, lc_coarse)
            raw = cache.with_suffix(".parallel.h5")
            ncores = 32 if lc == 500 else 16
            exports = {
                "ISMIP7_LC": lc,
                "ISMIP7_LC_COARSE": lc_coarse,
                "ISMIP7_BUFFER_M": BUFFER_M,
                "ISMIP7_MESH": mesh,
                "ISMIP7_BNDIDS": boundary,
                "ISMIP7_INVERSION": self.inversion,
                "ISMIP7_FRICTION": "budd",
                "ISMIP7_GEOMETRY_SPACE": "dg0",
                "ISMIP7_N_FLOW": "3.0",
                "ISMIP7_A4_FACTOR": "1.0",
                "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": SOLVER_MODE,
                "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
                "ISMIP7_RESCUE_ENABLED": "1",
                "ISMIP7_TIMING_CACHE_RAW": raw,
                "ISMIP7_TIMING_CACHE": cache,
                "ISMIP7_TIMING_CACHE_MANIFEST": manifest,
                "ISMIP7_TIMING_CACHE_STATUS": status_path,
            }
            self._submit(
                f"timing_cache_{lc}_{lc_coarse}",
                ncores,
                MEMORY_BY_LC[lc],
                exports,
                self.root / "scripts/batch_runners/timing_prepare.script",
                status_path,
            )

    def _lane_paths(self, lc, lc_coarse, ncores):
        return (
            self.timing_dir
            / timing_record_basename(CAMPAIGN_TAG, lc, lc_coarse, ncores),
            self.timing_dir
            / timing_status_basename(CAMPAIGN_TAG, lc, lc_coarse, ncores),
        )

    def lane_result(self, lc, lc_coarse, ncores):
        record_path, status_path = self._lane_paths(lc, lc_coarse, ncores)
        record, error = read_record(record_path)
        if error:
            return "invalid", error
        if record is not None:
            valid, detail = validate_timing_record(
                record, lc=lc, lc_coarse=lc_coarse, ncores=ncores
            )
            return ("passed" if valid else "failed"), detail
        status, active = self.reconcile_status(status_path)
        if active:
            return "active", status["state"]
        if status and status.get("state") in {
            "failed", "finished", "submission_failed", "not_runnable"
        }:
            return "failed", status.get("category", status["state"])
        return "missing", "no record or active status"

    def _submit_lane(self, lane):
        lc, lc_coarse, ncores = lane
        state, detail = self.lane_result(*lane)
        if state in {"passed", "active"} and not self.force:
            print(f"LANE {state.upper()} {lc}/{lc_coarse}/{ncores}: {detail}")
            return
        if state == "failed" and not self.force:
            print(f"LANE FAILED {lc}/{lc_coarse}/{ncores}: {detail}")
            return
        valid, cache_detail = self.cache_validation(lc, lc_coarse)
        if not valid:
            print(f"LANE WAITING CACHE {lc}/{lc_coarse}/{ncores}: {cache_detail}")
            return
        cache, manifest = cache_paths(self.cache_dir, lc, lc_coarse)
        boundary = self.boundary_path(lc, lc_coarse)
        _, status_path = self._lane_paths(lc, lc_coarse, ncores)
        if not boundary.is_file() and not (
            self.dry_run and self.assume_valid_caches
        ):
            print(
                f"LANE NOT RUNNABLE {lc}/{lc_coarse}/{ncores}: "
                f"missing {boundary}"
            )
            if not self.dry_run:
                atomic_write_status(
                    status_path,
                    "not_runnable",
                    reason="missing_boundary_ids",
                )
            return
        dt = expected_dt(lc)
        exports = {
            "ISMIP7_LC": lc,
            "ISMIP7_LC_COARSE": lc_coarse,
            "ISMIP7_BUFFER_M": BUFFER_M,
            "ISMIP7_BNDIDS": boundary,
            "ISMIP7_INVERSION": self.inversion,
            "ISMIP7_RESTART": cache,
            "ISMIP7_TIMING_CACHE_MANIFEST": manifest,
            "ISMIP7_FRICTION": "budd",
            "ISMIP7_GEOMETRY_SPACE": "dg0",
            "ISMIP7_N_FLOW": "3.0",
            "ISMIP7_A4_FACTOR": "1.0",
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": SOLVER_MODE,
            "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
            "ISMIP7_RESCUE_ENABLED": "0",
            "ISMIP7_SUBCYCLES": "1",
            "ISMIP7_T_END": f"{expected_t_end(lc):.12g}",
            "ISMIP7_DT": f"{dt:.12g}",
            "ISMIP7_TIMING_KIND": "matrix",
            "ISMIP7_TIMING_TAG": CAMPAIGN_TAG,
            "ISMIP7_TIMING_EXPERIMENT": (
                f"timing_{CAMPAIGN_TAG}_lcc{lc_coarse}_n{ncores}"
            ),
            "ISMIP7_TIMING_STATUS": status_path,
        }
        if self.monitor:
            exports.update({
                "ISMIP7_SNES_MONITOR": "1",
                "ISMIP7_SNES_LOG": self.logs_dir / (
                    f"timing_snes_{CAMPAIGN_TAG}_{lc}_{lc_coarse}"
                    f"_{ncores}.log"
                ),
            })
        self._submit(
            f"timing_{lc}_{lc_coarse}_{ncores}",
            ncores,
            MEMORY_BY_LC[lc],
            exports,
            self.root / "scripts/batch_runners/timing_transient.script",
            status_path,
        )

    def scout(self):
        for lane in self.selected_lanes(scout_lanes()):
            self._submit_lane(lane)

    def scale(self):
        scouts = {(lc, lc_coarse): ncores
                  for lc, lc_coarse, ncores in scout_lanes()}
        for lane in self.selected_lanes(scaling_lanes()):
            lc, lc_coarse, ncores = lane
            scout_ncores = scouts[(lc, lc_coarse)]
            scout_state, detail = self.lane_result(
                lc, lc_coarse, scout_ncores
            )
            _, status_path = self._lane_paths(*lane)
            decision = scaling_decision(scout_state)
            if decision == "block":
                print(
                    f"BLOCKED BY SCOUT {lc}/{lc_coarse}/{ncores}: {detail}"
                )
                if not self.dry_run:
                    atomic_write_status(
                        status_path,
                        "blocked_by_scout",
                        scout_cores=scout_ncores,
                        reason=detail,
                    )
                continue
            if decision == "wait":
                print(
                    f"WAITING FOR SCOUT {lc}/{lc_coarse}/{ncores}: "
                    f"{scout_state} ({detail})"
                )
                continue
            self._submit_lane(lane)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "scout", "scale"))
    parser.add_argument("--root", default=_ROOT)
    parser.add_argument("--cache-dir", default=_ROOT / "results/timing/cache")
    parser.add_argument("--timing-dir", default=_ROOT / "results/timing")
    parser.add_argument(
        "--inversion",
        default=(
            _ROOT / "mesh/"
            "inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5"
        ),
    )
    parser.add_argument("--partition", default="general")
    parser.add_argument("--walltime", default="12:00:00")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--assume-valid-caches",
        action="store_true",
        help="dry-run only: print downstream commands before caches exist",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--only-mesh",
        type=_parse_mesh,
        help="operate only on one LC/LC_coarse pair, for example 2500/25000",
    )
    parser.add_argument(
        "--monitor",
        action="store_true",
        help="enable per-lane SNES/KSP logs under results/logs",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.assume_valid_caches and not args.dry_run:
        raise SystemExit("--assume-valid-caches is allowed only with --dry-run")
    os.environ["ISMIP7_DIAGNOSTIC_LINEAR_SOLVER"] = SOLVER_MODE
    # A timing campaign is a fixed solver configuration, not an ambient-shell
    # experiment.  Pin PETSc's actual PETSC_UNLIMITED sentinel so an exported
    # legacy -1 cannot silently restore the 1e4 DIVERGED_DTOL cutoff.
    os.environ["ISMIP7_SNES_DIVERGENCE_TOL"] = SNES_DIVERGENCE_TOL_DEFAULT
    os.environ["ISMIP7_FRICTION"] = "budd"
    os.environ["ISMIP7_GEOMETRY_SPACE"] = "dg0"
    os.environ["ISMIP7_N_FLOW"] = "3.0"
    os.environ["ISMIP7_A4_FACTOR"] = "1.0"
    manager = CampaignManager(args)
    getattr(manager, args.stage)()
    if manager.submit_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
