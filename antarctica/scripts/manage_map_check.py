#!/usr/bin/env python3
r"""Drive one released MAP through its checks, stage by stage.

    make -C antarctica map-check MAP_CHECK_FRICTION=regularized_coulomb|budd
    make -C antarctica map-check-dry-run ...

One MAP under one law, on its own mesh (native) and transferred onto the
production mesh: fetch and verify, score the t = 0 state on each mesh, prepare
the transferred cache and audit its fill, run the strict 10-step lane on each
mesh, run a 10-year control on each, audit the controls, summarise. Every
stage stamps a status file under ``results/map_check/<map stem>/``, and each
call prints the stage table and submits the first stage that is runnable
(``--stages all`` submits every runnable one, ``--stages NAME`` one). Nothing
here touches the timing campaign's caches, tags or records: the cache role,
the cache stem and the lane tag are all named apart (timing_campaign).

The reasoning, the acceptance rules and the results live in
``antarctica/MAP_CHECK.md``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, os.fspath(_SCRIPTS))
sys.path.insert(0, os.fspath(_SCRIPTS.parents[1]))

from check_solver_qualification import qualification_verdict  # noqa: E402
from icepack2_tools.runconfig import DT_DEFAULT, MESH_FROM_CHECKPOINT  # noqa: E402
from icepack2_tools.solverconfig import (  # noqa: E402
    SNES_DIVERGENCE_TOL_DEFAULT,
    solver_provenance,
)
from manage_timing_campaign import (  # noqa: E402
    SlurmStageRunner,
    _timestamp,
    read_record,
    read_status,
)
from timing_campaign import (  # noqa: E402
    CACHE_SOLVER_MODE,
    MAP_CHECK_CACHE_ROLE,
    MAP_CHECK_KIND,
    MAP_CHECK_LAW_TAGS,
    MAP_CHECK_STEPS,
    MASS_RESIDUAL_TOL_GT,
    MEMORY_BY_LC,
    TRIPWIRE_DEFAULTS,
    atomic_write_json,
    atomic_write_status,
    expected_dt,
    lane_solver,
    map_check_cache_stem,
    map_check_run_tag,
    map_check_stem,
    map_check_tag,
    mesh_basename,
    sha256_file,
    solver_configuration_fingerprint,
    validate_map_check_manifest,
)

# The stages in the order the table prints them and `--stages next` tries them.
STAGES = (
    "fetch",
    "repack",
    "score_native",
    "prepare_transfer",
    "audit_cache",
    "score_transfer",
    "lane_transfer",
    "lane_native",
    "control_transfer",
    "control_native",
    "audit_controls",
    "summary",
)
DEPENDENCIES = {
    "fetch": (),
    "repack": ("fetch",),
    "score_native": ("repack",),
    "prepare_transfer": ("repack",),
    "audit_cache": ("prepare_transfer",),
    "score_transfer": ("prepare_transfer",),
    "lane_transfer": ("prepare_transfer",),
    "lane_native": ("repack",),
    "control_transfer": ("lane_transfer",),
    "control_native": ("lane_native",),
    "audit_controls": ("control_transfer", "control_native"),
    "summary": (),
}
MESH_OF_STAGE = {
    "score_native": "native", "lane_native": "native", "control_native": "native",
    "prepare_transfer": "transferred", "audit_cache": "transferred",
    "score_transfer": "transferred", "lane_transfer": "transferred",
    "control_transfer": "transferred",
}
# Stamp states that mean the allocation ended without what the stage needs.
_TERMINAL_STATES = {"failed", "submission_failed", "not_runnable", "invalid"}
_RELEASE_API = "https://api.github.com/repos/{owner}/{repo}/releases/tags/{tag}"
_RELEASE_URL_RE = re.compile(
    r"https://github\.com/(?P<owner>[^/]+)/(?P<repo>[^/]+)/releases/download/"
    r"(?P<tag>[^/]+)/(?P<asset>[^/]+)$"
)


def md5_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.md5()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def release_md5(url, basename, fetch=None):
    """The md5 the release notes list beside ``basename``, or None."""
    match = _RELEASE_URL_RE.match(url)
    if not match:
        return None
    api = _RELEASE_API.format(**match.groupdict())
    if fetch is None:
        def fetch(target):
            import urllib.request
            with urllib.request.urlopen(target, timeout=60) as response:
                return response.read().decode("utf-8")
    try:
        body = json.loads(fetch(api)).get("body", "")
    except Exception:  # network, JSON, anything: the caller then wants --md5
        return None
    for line in body.splitlines():
        if basename in line:
            found = re.findall(r"\b[0-9a-f]{32}\b", line)
            if found:
                return found[-1]
    return None


def _normal(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "item"):
        return value.item()
    return value


def read_map_attrs(path):
    """Root attributes of a Firedrake checkpoint via h5py, or None when h5py
    is not importable here."""
    try:
        import h5py
    except ImportError:
        return None
    with h5py.File(path, "r") as handle:
        attrs = {key: _normal(value) for key, value in handle["/"].attrs.items()}
        fields = sorted(
            name for name in handle.get("topologies", {}) or []
        )
    attrs["_topologies"] = fields
    return attrs


def _finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


class MapCheckManager(SlurmStageRunner):
    def __init__(self, args):
        # absolute(), never resolve(): the Makefile test runs a sandbox whose
        # scripts/ is a symlink into the checkout, and the sandbox is the root.
        self.root = Path(args.root).absolute()
        self.mesh_dir = self.root / "mesh"
        self.results_dir = self.root / "results"
        self.root_dir = Path(args.root_dir).absolute()
        self.friction = args.friction
        if self.friction not in MAP_CHECK_LAW_TAGS:
            raise SystemExit(
                f"--friction must be one of {tuple(MAP_CHECK_LAW_TAGS)}, "
                f"not {self.friction!r}"
            )
        self.url = args.url
        self.md5 = args.md5
        self.check_map_attrs = bool(args.check_map_attrs)
        self.repack = bool(args.repack)
        self.map_path = Path(args.map).absolute()
        self.basename = self.map_path.name
        self.stem = map_check_stem(self.map_path)
        self.download_path = self.root_dir / "maps" / "download" / self.basename
        self.fetch_target = self.download_path if self.repack else self.map_path
        self.stem_dir = self.root_dir / self.stem
        # Source mesh triple (the MAP's own) and the target's.
        self.lc, self.lc_coarse, self.buffer_m = (
            int(args.lc), int(args.lc_coarse), int(round(float(args.buffer_m)))
        )
        self.target_lc, self.target_lc_coarse, self.target_buffer_m = (
            int(args.target_lc), int(args.target_lc_coarse),
            int(round(float(args.target_buffer_m))),
        )
        self.target_mesh = Path(args.target_mesh).absolute()
        self.target_bndids = (
            Path(args.target_bndids).absolute() if args.target_bndids
            else self.mesh_dir / f"boundary_ids_{self.target_mesh.stem}.json"
        )
        self.native_bndids = self.mesh_dir / (
            "boundary_ids_"
            + mesh_basename(self.lc, self.lc_coarse, self.buffer_m)[:-4]
            + ".json"
        )
        self.solver = lane_solver(args.solver)
        self.cores = int(args.cores)
        self.diag_cores = int(args.diag_cores)
        self.diag_mem = args.diag_mem
        self.control_mem = args.control_mem
        self.control_time = args.control_time
        self.t_end = float(args.t_end)
        self.k_npz = args.k_npz or None
        self.steps = MAP_CHECK_STEPS
        self.queue = args.queue
        self.partition = args.partition
        self.constraint = args.constraint
        self.walltime = args.walltime
        self.dry_run = args.dry_run
        self.assume_passed = bool(args.assume_passed) and self.dry_run
        self.force = args.force
        self.controls_after_failed_lane = bool(
            getattr(args, "controls_after_failed_lane", False)
        )
        self.stages = args.stages
        if self.assume_passed and self.stages == "next":
            # A dry run that assumes every dependency passed is asking to see
            # the whole ladder, one composed request per stage.
            self.stages = "all"
        self.submit_failures = 0
        self.tag = map_check_tag(self.stem, self.friction, self.solver)
        cache_stem = map_check_cache_stem(
            self.stem, self.friction, self.target_lc, self.target_lc_coarse,
            self.target_buffer_m,
        )
        self.cache_dir = self.stem_dir / "cache"
        self.cache_path = self.cache_dir / f"{cache_stem}.h5"
        self.cache_raw = self.cache_dir / f"{cache_stem}.parallel.h5"
        self.cache_manifest = self.cache_dir / f"{cache_stem}.json"
        self.solver_fingerprint = solver_configuration_fingerprint(
            solver_provenance(CACHE_SOLVER_MODE)
        )
        self._sha = {}
        if not self.dry_run:
            self.stem_dir.mkdir(parents=True, exist_ok=True)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self.fetch_target.parent.mkdir(parents=True, exist_ok=True)

    # --- names -------------------------------------------------------------

    def status_path(self, stage):
        return self.stem_dir / f"status_{stage}.txt"

    def score_json(self, role):
        return self.stem_dir / f"score_{'native' if role == 'native' else 'transfer'}.json"

    def census_path(self):
        return self.stem_dir / "census_native.txt"

    def lane_lc(self, role):
        return self.lc if role == "native" else self.target_lc

    def lane_lc_coarse(self, role):
        return self.lc_coarse if role == "native" else self.target_lc_coarse

    def lane_dt(self, role):
        return expected_dt(self.lane_lc(role))

    def control_dt(self, role):
        r"""The step of the 10-year control. The control on the production
        mesh runs the production step (issue 20); the one on the MAP's own
        mesh, and both strict lanes, keep the timing matrix's rule, so the
        lanes stay comparable with the matrix rows."""
        if role == "native":
            return self.lane_dt(role)
        return float(DT_DEFAULT)

    def lane_experiment(self, role):
        return f"{self.tag}_lcc{self.lane_lc_coarse(role)}_n{self.cores}"

    def lane_record(self, role):
        return self.stem_dir / (
            f"timing_{self.tag}_{self.lane_lc(role)}_{self.lane_lc_coarse(role)}"
            f"_{self.cores}.json"
        )

    def lane_csv(self, role):
        return self.results_dir / (
            f"{self.lane_experiment(role)}_{self.lane_lc(role)}_timeseries.csv"
        )

    def run_tag(self, role):
        return map_check_run_tag(self.stem, self.friction, str(self.lane_lc(role)))

    def control_experiment(self, role):
        return f"ctrl2015_cesm2_waccm_{self.run_tag(role)}"

    def control_csv(self, role):
        return self.results_dir / (
            f"{self.control_experiment(role)}_{self.lane_lc(role)}_timeseries.csv"
        )

    def control_final(self, role):
        return self.results_dir / (
            f"{self.control_experiment(role)}_{self.lane_lc(role)}_final.h5"
        )

    def control_job_name(self, role):
        return f"mc_{MAP_CHECK_LAW_TAGS[self.friction]}_ctrl{self.lane_lc(role)}_{self.stem[-13:]}"

    def sha(self, path):
        key = os.fspath(path)
        if key not in self._sha:
            self._sha[key] = sha256_file(path) if Path(path).is_file() else None
        return self._sha[key]

    # --- environment -------------------------------------------------------

    def common_exports(self):
        return {
            "ISMIP7_INVERSION": self.map_path,
            "ISMIP7_FRICTION": self.friction,
            "ISMIP7_GEOMETRY_SPACE": "dg0",
            "ISMIP7_N_FLOW": "3.0",
            "ISMIP7_A4_FACTOR": "1.0",
            "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
        }

    def mesh_exports(self, role, mesh=MESH_FROM_CHECKPOINT):
        """The mesh triple and sidecar of one role. ``mesh`` is the compute
        mesh: the sentinel (solve on the MAP's or the cache's own mesh) or the
        target .msh (a transfer)."""
        if role == "native":
            return {
                "ISMIP7_LC": self.lc,
                "ISMIP7_LC_COARSE": self.lc_coarse,
                "ISMIP7_BUFFER_M": self.buffer_m,
                "ISMIP7_BNDIDS": self.native_bndids,
                "ISMIP7_MESH": MESH_FROM_CHECKPOINT,
            }
        return {
            "ISMIP7_LC": self.target_lc,
            "ISMIP7_LC_COARSE": self.target_lc_coarse,
            "ISMIP7_BUFFER_M": self.target_buffer_m,
            "ISMIP7_BNDIDS": self.target_bndids,
            "ISMIP7_MESH": mesh,
        }

    # --- verdicts from artifacts ------------------------------------------

    def verdict_fetch(self):
        record, _ = read_record(self.stem_dir / "fetch.json")
        if record is None:
            return None, ""
        if record.get("ok"):
            return "passed", record.get("detail", "")
        return "failed", record.get("detail", "")

    def verdict_repack(self):
        if not self.repack:
            return ("passed", "no repack asked for") if self.fetch_target.is_file() else (None, "")
        if self.map_path.is_file() and (
            not self.download_path.is_file()
            or self.map_path.stat().st_mtime >= self.download_path.stat().st_mtime
        ):
            return "passed", os.fspath(self.map_path)
        return None, ""

    def verdict_score(self, role):
        record, error = read_record(self.score_json(role))
        if error:
            return "failed", error
        if record is None:
            return None, ""
        if isinstance(record, list):
            record = record[0]
        if not _finite(record.get("ratio")):
            return "failed", f"ratio={record.get('ratio')!r}"
        a_prior = record.get("a_prior_range") or [None, None]
        if not (_finite(a_prior[0]) and float(a_prior[0]) > 0.0):
            return "failed", f"fluidity prior minimum {a_prior[0]!r} is not positive"
        fill = record.get("transfer_fill") or {}
        filled = {k: v["missing"] for k, v in fill.items() if v.get("missing")}
        if role == "native" and filled:
            return "failed", f"a native score filled dofs: {filled}"
        bands = " ".join(f"{b['ratio']:.2f}" for b in record.get("bands", []))
        detail = f"ratio {record['ratio']:.2f} (bands {bands})"
        if role != "native":
            detail += f", filled {sum(filled.values())} dofs"
        misfit = record.get("initial_misfit")
        if _finite(misfit):
            detail += f", misfit0 {float(misfit):.3e}"
        return "passed", detail

    def verdict_prepare(self):
        manifest, error = read_record(self.cache_manifest)
        if error:
            return "failed", error
        if manifest is None:
            return None, ""
        valid, detail = validate_map_check_manifest(
            manifest,
            lc=self.target_lc,
            lc_coarse=self.target_lc_coarse,
            buffer_m=self.target_buffer_m,
            friction=self.friction,
            source_basename=self.basename,
            mesh_name=self.target_mesh.name,
            cache_path=self.cache_path,
            source_sha256=self.sha(self.map_path),
            mesh_sha256=self.sha(self.target_mesh),
            solver_fingerprint=self.solver_fingerprint,
        )
        if not valid:
            return "failed", detail
        fill = manifest.get("transfer_fill") or {}
        filled = {k: v.get("missing", 0) for k, v in fill.items() if v.get("missing")}
        return "passed", (
            f"filled {filled}" if filled else "no dof outside the source mesh"
        )

    def verdict_json_present(self, path):
        record, error = read_record(path)
        if error:
            return "failed", error
        if record is None:
            return None, ""
        return "passed", os.fspath(path)

    def verdict_lane(self, role):
        record, error = read_record(self.lane_record(role))
        if error:
            return "failed", error
        if record is None:
            return None, ""
        rows = []
        if self.lane_csv(role).is_file():
            with open(self.lane_csv(role), newline="") as stream:
                rows = list(csv.DictReader(stream))
        ok, detail = qualification_verdict(
            record, rows, self.steps, MASS_RESIDUAL_TOL_GT, kinds=(MAP_CHECK_KIND,)
        )
        if not ok:
            return "failed", detail
        if record.get("timing_tag") != self.tag:
            return "failed", f"timing_tag={record.get('timing_tag')!r}"
        if record.get("initial_state_source") != self.basename:
            return "failed", (
                f"initial_state_source={record.get('initial_state_source')!r}"
            )
        if record.get("tripwire"):
            return "failed", f"tripwire {record['tripwire']}"
        if record.get("rescue_enabled") is not False:
            return "failed", "rescue was enabled"
        summary = record.get("diagnostic_solve_summary", {})
        newton = summary.get("snes_iterations_total", 0) / max(record.get("completed_steps", 1), 1)
        return "passed", (
            f"{record.get('seconds_per_step', float('nan')):.1f} s/step on "
            f"{record.get('ncores')} ranks, {newton:.1f} Newton/step, "
            f"resid max {record.get('step_mass_residual_gt_max')}"
        )

    def verdict_control(self, role):
        path = self.control_csv(role)
        if not path.is_file():
            return None, ""
        with open(path, newline="") as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            return None, "empty timeseries"
        dt = self.control_dt(role)
        last = float(rows[-1]["year"])
        resid = max(abs(float(row["resid_gt"])) for row in rows)
        if resid > MASS_RESIDUAL_TOL_GT:
            return "failed", f"mass residual {resid:.3g} Gt at year {last}"
        if last < self.t_end - 0.5 * dt:
            return None, f"at {last:.2f} of {self.t_end:.0f}"
        if not self.control_final(role).is_file():
            return None, f"reached {last:.2f}, final state not written yet"
        vaf = [float(row["vaf_mm_sle"]) for row in rows]
        years = [float(row["year"]) for row in rows]
        span = years[-1] - years[0]
        trend = (vaf[-1] - vaf[0]) / span if span > 0 else float("nan")
        return "passed", f"reached {last:.2f}, dVAF/dt {trend:+.2f} mm SLE/yr over the run"

    def verdict(self, stage):
        if stage == "fetch":
            return self.verdict_fetch()
        if stage == "repack":
            return self.verdict_repack()
        if stage in ("score_native", "score_transfer"):
            return self.verdict_score(MESH_OF_STAGE[stage])
        if stage == "prepare_transfer":
            return self.verdict_prepare()
        if stage == "audit_cache":
            return self.verdict_json_present(self.stem_dir / "cache_audit.json")
        if stage in ("lane_transfer", "lane_native"):
            return self.verdict_lane(MESH_OF_STAGE[stage])
        if stage in ("control_transfer", "control_native"):
            return self.verdict_control(MESH_OF_STAGE[stage])
        if stage == "audit_controls":
            return self.verdict_json_present(self.stem_dir / "audit_controls.json")
        if stage == "summary":
            return (("passed", "written") if (self.stem_dir / "summary.md").is_file()
                    else (None, ""))
        raise KeyError(stage)

    # --- states ------------------------------------------------------------

    def _job_named_active(self, name):
        """Whether Slurm holds a job of this name (a chained forward changes
        id, never name). None when squeue is not here to ask."""
        if shutil.which("squeue") is None:
            return None
        result = subprocess.run(
            ["squeue", "-h", "-n", name, "-o", "%i"],
            check=False, capture_output=True, text=True,
        )
        if result.returncode != 0:
            return None
        return bool(result.stdout.strip())

    def inputs_missing(self, stage):
        role = MESH_OF_STAGE.get(stage)
        missing = []
        if stage in ("repack", "score_native", "lane_native", "control_native",
                     "prepare_transfer") and not self.map_path.is_file() \
                and not (stage == "repack" and self.download_path.is_file()):
            missing.append(os.fspath(self.map_path))
        if role == "native" and not self.native_bndids.is_file():
            missing.append(os.fspath(self.native_bndids))
        if role == "transferred":
            if not self.target_bndids.is_file():
                missing.append(os.fspath(self.target_bndids))
            if stage in ("prepare_transfer", "control_transfer") \
                    and not self.target_mesh.is_file():
                missing.append(os.fspath(self.target_mesh))
        return missing

    def state(self, stage, states):
        """``(state, detail)`` with ``states`` the already computed states of
        earlier stages (the table is filled in order)."""
        verdict, detail = self.verdict(stage)
        # --force reaches the named stages alone: a failed lane stays failed
        # when its control is the one being forced, so that the control's
        # dependency reads as failed (and the switch below may lift it) rather
        # than as pending.
        forced = self.force and stage in self.forced()
        if verdict == "passed" and not forced:
            return "passed", detail
        if verdict == "failed" and not forced:
            return "failed", detail
        status = read_status(self.status_path(stage))
        if status is not None and not forced:
            state = status.get("state")
            if state in {"submitting", "submitted", "running"}:
                if stage.startswith("control_"):
                    live = self._job_named_active(
                        self.control_job_name(MESH_OF_STAGE[stage])
                    )
                    if live:
                        return "active", f"{state} job={status.get('job_id', '?')} (chain live)"
                    if live is None:
                        return "active", f"{state} job={status.get('job_id', '?')}"
                    status, active = self.reconcile_status(self.status_path(stage))
                    if active:
                        return "active", f"{status.get('state')} job={status.get('job_id', '?')}"
                    return "failed", (
                        f"{status.get('category', status.get('state'))}"
                        + (f" {detail}" if detail else "")
                    )
                status, active = self.reconcile_status(self.status_path(stage))
                if active:
                    return "active", f"{status.get('state')} job={status.get('job_id', '?')}"
                return "failed", (
                    f"{status.get('category', status.get('state'))}"
                    + (f" {detail}" if detail else "")
                )
            if state in _TERMINAL_STATES:
                return "failed", status.get("category") or status.get("reason") or state
            if state == "finished" and verdict is None:
                if stage.startswith("control_"):
                    return "active", detail or "chain running"
                return "failed", "finished without its artifacts"
        dep_states = [states.get(dep, ("pending", ""))[0] for dep in DEPENDENCIES[stage]]
        if self.assume_passed:
            dep_states = ["passed" for _ in dep_states]
        # A control cold-starts from the MAP under the production configuration
        # (apparent mass balance, fixed front), which the strict lane does not
        # run. When the lane failed on the state itself, the switch lets the
        # control run anyway, for the cost per simulated year and the drift
        # under production settings; the summary says which lane failed.
        lanes_failed = []
        if self.controls_after_failed_lane and stage.startswith("control_"):
            for index, dep in enumerate(DEPENDENCIES[stage]):
                if dep.startswith("lane_") and dep_states[index] == "failed":
                    dep_states[index] = "passed"
                    lanes_failed.append(dep)
        if any(s in {"failed", "blocked", "not_runnable"} for s in dep_states):
            return "blocked", "a dependency failed"
        if any(s != "passed" for s in dep_states):
            return "waiting", "for " + ", ".join(
                dep for dep, s in zip(DEPENDENCIES[stage], dep_states) if s != "passed"
            )
        missing = self.inputs_missing(stage)
        if missing and not self.assume_passed:
            return "not_runnable", "missing " + ", ".join(missing)
        if lanes_failed:
            detail = (
                ", ".join(lanes_failed) + " failed; the control runs the "
                "production configuration regardless (--controls-after-failed-lane)"
            )
        return "pending", detail

    def forced(self):
        return set(STAGES) if self.stages == "all" else {self.stages}

    def table(self):
        states = {}
        for stage in STAGES:
            states[stage] = self.state(stage, states)
        return states

    def print_table(self, states):
        print(f"MAP  {self.basename}  law {self.friction}  solver {self.solver}")
        print(f"root {self.stem_dir}")
        print(f" {'#':>2}  {'stage':<17} {'mesh':<19} {'state':<12} detail")
        for number, stage in enumerate(STAGES):
            role = MESH_OF_STAGE.get(stage)
            mesh = "-"
            if role == "native":
                mesh = f"{self.lc}/{self.lc_coarse}/b{self.buffer_m}"
            elif role == "transferred":
                mesh = f"{self.target_lc}/{self.target_lc_coarse}/b{self.target_buffer_m}"
            state, detail = states[stage]
            print(f" {number:>2}  {stage:<17} {mesh:<19} {state:<12} {detail}")

    # --- stage runners -------------------------------------------------------

    def run_fetch(self):
        status_path = self.status_path("fetch")
        target = self.fetch_target
        expected = self.md5
        if expected == "auto":
            expected = release_md5(self.url, self.basename)
        if self.dry_run:
            print(
                f"DRY RUN: curl -L --fail -o {target} {self.url}; md5 must be "
                f"{expected or '(from the release notes)'}"
            )
            return
        if not target.is_file():
            atomic_write_status(status_path, "running", phase="fetch", timestamp=_timestamp())
            tmp = target.with_suffix(target.suffix + ".part")
            result = subprocess.run(
                ["curl", "-L", "--fail", "-o", os.fspath(tmp), self.url],
                check=False,
            )
            if result.returncode != 0:
                atomic_write_status(status_path, "failed", category="download",
                                    exit_code=result.returncode, timestamp=_timestamp())
                self.submit_failures += 1
                return
            os.replace(tmp, target)
        digest = md5_file(target)
        record = {"path": os.fspath(target), "url": self.url, "md5": digest,
                  "md5_expected": expected, "size": target.stat().st_size,
                  "timestamp": _timestamp()}
        problems = []
        if expected is None:
            problems.append("no md5 to check against (pass --md5)")
        elif digest != expected:
            problems.append(f"md5 {digest} differs from the release's {expected}")
        if self.check_map_attrs and not problems:
            attrs = read_map_attrs(target)
            record["attrs"] = attrs
            if attrs is None:
                record["attrs_note"] = "h5py not importable here; attributes unchecked"
            else:
                want = {
                    "friction": self.friction,
                    "geometry_space": "dg0",
                    "mesh_basename": mesh_basename(self.lc, self.lc_coarse, self.buffer_m),
                }
                for key, value in want.items():
                    got = attrs.get(key)
                    if got is not None and str(got) != str(value):
                        problems.append(f"{key}={got!r}, expected {value!r}")
                if attrs.get("n_flow") is not None and abs(float(attrs["n_flow"]) - 3.0) > 1e-9:
                    problems.append(f"n_flow={attrs['n_flow']!r}, expected 3")
                for key, value in (("lc", self.lc), ("lc_coarse", self.lc_coarse),
                                   ("buffer_m", self.buffer_m)):
                    got = attrs.get(key)
                    if got is not None and int(round(float(got))) != value:
                        problems.append(f"{key}={got!r}, expected {value}")
        record["ok"] = not problems
        record["detail"] = (
            "; ".join(problems) if problems
            else f"md5 ok, {record['size'] / 1e6:.0f} MB"
            + ("" if not record.get("attrs") else
               f", attrs friction={record['attrs'].get('friction')} "
               f"mesh={record['attrs'].get('mesh_basename')}")
        )
        atomic_write_json(self.stem_dir / "fetch.json", record)
        atomic_write_status(
            status_path, "finished" if record["ok"] else "failed",
            category="verification" if problems else "ok", timestamp=_timestamp(),
        )
        print(f"FETCH {'PASSED' if record['ok'] else 'FAILED'}: {record['detail']}")

    def run_repack(self):
        self._submit(
            f"mc_repack_{self.stem[-13:]}", 1, "32G",
            {"TIMING_REDISTRIBUTE_INPUT": self.download_path,
             "TIMING_REDISTRIBUTE_OUTPUT": self.map_path},
            self.root / "scripts/batch_runners/timing_redistribute.script",
            self.status_path("repack"),
        )

    def run_score(self, role):
        exports = self.common_exports()
        exports.update(self.mesh_exports(role))
        exports.update({
            "ISMIP7_MAP": self.map_path,
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": CACHE_SOLVER_MODE,
            "ISMIP7_RESCUE_ENABLED": "1",
            "ISMIP7_MAP_CHECK_SCORE_JSON": self.score_json(role),
            "ISMIP7_MAP_CHECK_STATUS": self.status_path(f"score_{'native' if role == 'native' else 'transfer'}"),
        })
        if role == "native":
            if self.friction == "budd":
                exports["ISMIP7_MAP_CHECK_CENSUS"] = self.census_path()
            walltime, mem = self.walltime, self.diag_mem
        else:
            exports["ISMIP7_MAP_CHECK_RESTART"] = self.cache_path
            exports["ISMIP7_TIMING_CACHE_MANIFEST"] = self.cache_manifest
            walltime, mem = "02:00:00", MEMORY_BY_LC.get(self.target_lc, self.diag_mem)
        saved = self.walltime
        self.walltime = walltime
        try:
            self._submit(
                f"mc_score_{'nat' if role == 'native' else 'tra'}_{self.stem[-13:]}",
                self.diag_cores, mem, exports,
                self.root / "scripts/batch_runners/map_check_score.script",
                exports["ISMIP7_MAP_CHECK_STATUS"],
            )
        finally:
            self.walltime = saved

    def run_prepare(self):
        status_path = self.status_path("prepare_transfer")
        exports = self.common_exports()
        exports.update(self.mesh_exports("transferred", mesh=self.target_mesh))
        exports.update({
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": CACHE_SOLVER_MODE,
            "ISMIP7_RESCUE_ENABLED": "1",
            "ISMIP7_TIMING_CACHE_ROLE": MAP_CHECK_CACHE_ROLE,
            "ISMIP7_TIMING_CACHE_RAW": self.cache_raw,
            "ISMIP7_TIMING_CACHE": self.cache_path,
            "ISMIP7_TIMING_CACHE_MANIFEST": self.cache_manifest,
            "ISMIP7_TIMING_CACHE_STATUS": status_path,
        })
        self._submit(
            f"mc_prepare_{self.stem[-13:]}", self.diag_cores,
            MEMORY_BY_LC.get(self.target_lc, self.diag_mem), exports,
            self.root / "scripts/batch_runners/timing_prepare.script",
            status_path,
        )

    def run_audit_cache(self):
        status_path = self.status_path("audit_cache")
        exports = self.common_exports()
        exports.update(self.mesh_exports("transferred"))
        exports.update({
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": CACHE_SOLVER_MODE,
            "ISMIP7_TIMING_CACHE": self.cache_path,
            "ISMIP7_TIMING_CACHE_MANIFEST": self.cache_manifest,
            "ISMIP7_TIMING_CACHE_AUDIT_OUTPUT": self.stem_dir / "cache_audit.json",
            "ISMIP7_TIMING_CACHE_AUDIT_STATUS": status_path,
        })
        saved = self.walltime
        self.walltime = "02:00:00"
        try:
            self._submit(
                f"mc_audit_{self.stem[-13:]}", self.diag_cores,
                MEMORY_BY_LC.get(self.target_lc, self.diag_mem), exports,
                self.root / "scripts/batch_runners/timing_cache_audit.script",
                status_path,
            )
        finally:
            self.walltime = saved

    def run_lane(self, role):
        stage = f"lane_{'native' if role == 'native' else 'transfer'}"
        status_path = self.status_path(stage)
        dt = self.lane_dt(role)
        exports = self.common_exports()
        exports.update(self.mesh_exports(role))
        exports.update({
            "ISMIP7_RESTART": self.cache_path if role != "native" else "",
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": self.solver,
            "ISMIP7_RESCUE_ENABLED": "0",
            "ISMIP7_SUBCYCLES": "1",
            "ISMIP7_T_END": f"{2015.0 + self.steps * dt:.12g}",
            "ISMIP7_DT": f"{dt:.12g}",
            "ISMIP7_TIMING_KIND": MAP_CHECK_KIND,
            "ISMIP7_TIMING_TAG": self.tag,
            "ISMIP7_TIMING_EXPERIMENT": self.lane_experiment(role),
            "ISMIP7_TIMING_DIR": self.stem_dir,
            "ISMIP7_TIMING_STATUS": status_path,
            "ISMIP7_MAP_CHECK_STEPS": self.steps,
            "ISMIP7_AMB_CAP": "0",
            "ISMIP7_OUTPUT_INTERVAL": "1",
        })
        if role != "native":
            exports["ISMIP7_TIMING_CACHE_MANIFEST"] = self.cache_manifest
        exports.update(TRIPWIRE_DEFAULTS)
        self._submit(
            f"mc_lane_{'nat' if role == 'native' else 'tra'}_{self.stem[-13:]}",
            self.cores, MEMORY_BY_LC.get(self.lane_lc(role), self.diag_mem),
            exports,
            self.root / "scripts/batch_runners/timing_transient.script",
            status_path,
        )

    def run_control(self, role):
        stage = f"control_{'native' if role == 'native' else 'transfer'}"
        status_path = self.status_path(stage)
        exports = self.common_exports()
        exports.update(self.mesh_exports(
            role, mesh=MESH_FROM_CHECKPOINT if role == "native" else self.target_mesh
        ))
        exports.update({
            "ISMIP7_EXPERIMENT": "control",
            "ISMIP7_T_END": f"{self.t_end:g}",
            "ISMIP7_DT": f"{self.control_dt(role):.12g}",
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": self.solver,
            "ISMIP7_RUN_TAG": self.run_tag(role),
            "ISMIP7_OUTPUT": "1",
            "ISMIP7_OUTPUT_INTERVAL": "2",
        })
        if self.k_npz:
            exports["ISMIP7_K_PER_BASIN_NPZ"] = self.k_npz
        exports.update(TRIPWIRE_DEFAULTS)
        saved = self.walltime
        self.walltime = self.control_time
        try:
            self._submit(
                self.control_job_name(role), self.cores, self.control_mem,
                exports, None, status_path, kind="projection",
            )
        finally:
            self.walltime = saved

    def run_audit_controls(self):
        status_path = self.status_path("audit_controls")
        exports = self.common_exports()
        exports.update({
            "ISMIP7_MAP_CHECK_CSV_NATIVE": self.control_csv("native"),
            "ISMIP7_MAP_CHECK_CSV_TRANSFER": self.control_csv("transferred"),
            "ISMIP7_MAP_CHECK_FINAL_NATIVE": self.control_final("native"),
            "ISMIP7_MAP_CHECK_FINAL_TRANSFER": self.control_final("transferred"),
            "ISMIP7_MAP_CHECK_DT_NATIVE": f"{self.control_dt('native'):.12g}",
            "ISMIP7_MAP_CHECK_DT_TRANSFER": f"{self.control_dt('transferred'):.12g}",
            "ISMIP7_MAP_CHECK_AUDIT_JSON": self.stem_dir / "audit_controls.json",
            "ISMIP7_MAP_CHECK_PNG": self.stem_dir / "compare_controls.png",
            "ISMIP7_MAP_CHECK_STATUS": status_path,
            "ISMIP7_MAP_CHECK_NATIVE_LC": self.lc,
            "ISMIP7_MAP_CHECK_NATIVE_LC_COARSE": self.lc_coarse,
            "ISMIP7_MAP_CHECK_NATIVE_BUFFER_M": self.buffer_m,
            "ISMIP7_MAP_CHECK_NATIVE_BNDIDS": self.native_bndids,
            "ISMIP7_MAP_CHECK_TRANSFER_LC": self.target_lc,
            "ISMIP7_MAP_CHECK_TRANSFER_LC_COARSE": self.target_lc_coarse,
            "ISMIP7_MAP_CHECK_TRANSFER_BUFFER_M": self.target_buffer_m,
            "ISMIP7_MAP_CHECK_TRANSFER_BNDIDS": self.target_bndids,
        })
        if self.k_npz:
            exports["ISMIP7_K_PER_BASIN_NPZ"] = self.k_npz
        saved = self.walltime
        self.walltime = "02:00:00"
        try:
            self._submit(
                f"mc_auditc_{self.stem[-13:]}", self.diag_cores, "64G", exports,
                self.root / "scripts/batch_runners/map_check_audit.script",
                status_path,
            )
        finally:
            self.walltime = saved

    def run_summary(self, states):
        lines = [
            f"# MAP check: {self.basename}", "",
            f"Law `{self.friction}`, lane solver `{self.solver}`, "
            f"native mesh {self.lc}/{self.lc_coarse} buffer {self.buffer_m} m, "
            f"transferred onto `{self.target_mesh.name}`. Generated {_timestamp()}.", "",
            "| # | stage | mesh | state | detail |", "|---|---|---|---|---|",
        ]
        for number, stage in enumerate(STAGES):
            if stage == "summary":
                continue
            role = MESH_OF_STAGE.get(stage)
            mesh = ("-" if role is None else
                    f"{self.lc}/{self.lc_coarse}" if role == "native"
                    else f"{self.target_lc}/{self.target_lc_coarse}")
            state, detail = states[stage]
            lines.append(f"| {number} | {stage} | {mesh} | {state} | {detail} |")
        for role in ("native", "transferred"):
            record, _ = read_record(self.score_json(role))
            if record:
                if isinstance(record, list):
                    record = record[0]
                lines += ["", f"## t = 0 on the {role} mesh", "",
                          f"Q(u_model)/Q(u_obs) {record.get('ratio', float('nan')):.3f}; "
                          f"initial misfit {record.get('initial_misfit')}; "
                          f"fluidity prior range {record.get('a_prior_range')}.", "",
                          "| band (m/yr) | model (Gt/yr) | u_obs (Gt/yr) | ratio |",
                          "|---|---|---|---|"]
                for band in record.get("bands", []):
                    lines.append(
                        f"| {band['band']} | {band['q_model']:.0f} | "
                        f"{band['q_obs']:.0f} | {band['ratio']:.2f} |"
                    )
        for role in ("native", "transferred"):
            record, _ = read_record(self.lane_record(role))
            if record:
                summary = record.get("diagnostic_solve_summary", {})
                lines += ["", f"## Strict {self.steps}-step lane, {role} mesh", "",
                          f"{record.get('seconds_per_step', float('nan')):.2f} s per step on "
                          f"{record.get('ncores')} ranks under `{record.get('diagnostic_solver_mode')}`, "
                          f"dt {record.get('dt')} yr, {summary.get('snes_iterations_total', 0)} "
                          f"Newton iterations over {record.get('completed_steps')} steps, "
                          f"mass residual max {record.get('step_mass_residual_gt_max')} Gt, "
                          f"tripwire {record.get('tripwire')}."]
        ran_after_failure = [
            role for role, suffix in (("native", "native"), ("transferred", "transfer"))
            if states[f"lane_{suffix}"][0] == "failed"
            and states[f"control_{suffix}"][0]
            not in {"waiting", "blocked", "pending", "not_runnable"}
        ]
        if ran_after_failure:
            lines += ["", "## Controls after a failed lane", "",
                      "The strict lane failed on the " + " and ".join(ran_after_failure)
                      + " mesh and the control on that mesh ran anyway "
                      "(`--controls-after-failed-lane`). A control cold-starts from "
                      "the MAP under the production configuration, apparent mass "
                      "balance and a fixed front, which cancels the t = 0 thickness "
                      "tendency the strict contract exposes; its cost and drift are "
                      "production numbers and no stability verdict on the state."]
        audit, _ = read_record(self.stem_dir / "audit_controls.json")
        if audit:
            lines += ["", "## 10-year controls", "", "```",
                      json.dumps(audit, indent=2, sort_keys=True)[:6000], "```"]
        text = "\n".join(lines) + "\n"
        if self.dry_run:
            print(f"DRY RUN: would write {self.stem_dir / 'summary.md'}")
            return
        self.stem_dir.mkdir(parents=True, exist_ok=True)
        tmp = self.stem_dir / "summary.md.tmp"
        tmp.write_text(text)
        os.replace(tmp, self.stem_dir / "summary.md")
        print(f"SUMMARY -> {self.stem_dir / 'summary.md'}")

    def run_stage(self, stage, states):
        if stage == "fetch":
            self.run_fetch()
        elif stage == "repack":
            self.run_repack()
        elif stage in ("score_native", "score_transfer"):
            self.run_score(MESH_OF_STAGE[stage])
        elif stage == "prepare_transfer":
            self.run_prepare()
        elif stage == "audit_cache":
            self.run_audit_cache()
        elif stage in ("lane_transfer", "lane_native"):
            self.run_lane(MESH_OF_STAGE[stage])
        elif stage in ("control_transfer", "control_native"):
            self.run_control(MESH_OF_STAGE[stage])
        elif stage == "audit_controls":
            self.run_audit_controls()
        elif stage == "summary":
            self.run_summary(states)
        else:
            raise KeyError(stage)

    def run(self):
        states = self.table()
        self.print_table(states)
        runnable = [s for s in STAGES if states[s][0] == "pending" or
                    (self.force and s in self.forced() and states[s][0] in
                     {"passed", "failed", "pending"})]
        if self.stages == "next":
            chosen = [s for s in runnable if s != "summary"][:1] or (
                ["summary"] if "summary" in runnable else []
            )
        elif self.stages == "all":
            chosen = runnable
        else:
            if self.stages not in STAGES:
                raise SystemExit(f"--stages must be next, all or one of {STAGES}")
            chosen = [self.stages] if self.stages in runnable else []
            if not chosen:
                print(f"{self.stages}: {states[self.stages][0]} ({states[self.stages][1]}); nothing to submit")
        for stage in chosen:
            print(f"--- {stage}")
            self.run_stage(stage, states)
        if not chosen and self.stages in ("next", "all"):
            active = [s for s in STAGES if states[s][0] == "active"]
            print("nothing to submit" + (f"; active: {', '.join(active)}" if active else ""))
        # The summary always reflects the table as it stands.
        if "summary" not in chosen and not self.dry_run:
            self.run_summary(states)
        return 1 if self.submit_failures else 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="print the stage table and submit what is runnable")
    run.add_argument("--map", required=True, help="the MAP the checks read (the repack)")
    run.add_argument("--url", required=True, help="release asset URL of the MAP")
    run.add_argument("--md5", default="auto",
                     help="md5 the release lists, or `auto` to read the release notes")
    run.add_argument("--friction", required=True, choices=sorted(MAP_CHECK_LAW_TAGS))
    run.add_argument("--lc", default=2000)
    run.add_argument("--lc-coarse", default=5000)
    run.add_argument("--buffer-m", default=0)
    run.add_argument("--target-mesh", required=True)
    run.add_argument("--target-lc", default=1000)
    run.add_argument("--target-lc-coarse", default=10000)
    run.add_argument("--target-buffer-m", default=20000)
    run.add_argument("--target-bndids", default=None)
    run.add_argument("--root", default=_SCRIPTS.parent,
                     help="the antarctica/ directory to run from")
    run.add_argument("--root-dir", default=None,
                     help="where the checks keep their files (default results/map_check)")
    run.add_argument("--stages", default="next", help="next | all | <stage>")
    run.add_argument("--solver", default="scpc_gamg")
    run.add_argument("--cores", default=64, help="ranks of the lanes and controls")
    run.add_argument("--diag-cores", default=16, help="ranks of the scores, prepare, audits")
    run.add_argument("--diag-mem", default="96G")
    run.add_argument("--control-mem", default="240G")
    run.add_argument("--control-time", default="12:00:00")
    run.add_argument("--t-end", default=2025)
    run.add_argument("--k-npz", default=None, help="ISMIP7_K_PER_BASIN_NPZ of the controls")
    run.add_argument("--queue", default="short")
    run.add_argument("--partition", default=None)
    run.add_argument("--constraint", default=None)
    run.add_argument("--walltime", default="12:00:00")
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--assume-passed", action="store_true",
                     help="dry run only: show every stage's request as if its "
                          "dependencies had passed")
    run.add_argument("--force", action="store_true",
                     help="resubmit the named stage(s) whatever their state")
    run.add_argument("--controls-after-failed-lane", action="store_true",
                     help="run a control although the strict lane on its mesh "
                          "failed (production configuration, cost and drift)")
    run.add_argument("--check-map-attrs", type=int, default=1)
    run.add_argument("--repack", type=int, default=1)
    args = parser.parse_args(argv)
    if args.root_dir is None:
        args.root_dir = Path(args.root) / "results" / "map_check"
    return args


def main(argv=None):
    args = parse_args(argv)
    manager = MapCheckManager(args)
    return manager.run()


if __name__ == "__main__":
    raise SystemExit(main())
