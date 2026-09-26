r"""The map-check stage machine (antarctica/scripts/manage_map_check.py).

One released MAP through its checks, with each stage's state read off its
artifacts and stamps, and its requests composed through submit.sh like the
timing campaign's. Pure Python: no Firedrake, a fake sbatch where a
submission is exercised.
"""
import csv
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ANT = REPO / "antarctica"
sys.path.insert(0, str(ANT / "scripts"))

import manage_map_check as mmc  # noqa: E402
from timing_campaign import parse_campaign_tag  # noqa: E402

BASENAME = "inversion_icepack2_rc_n3_dg0_logvelnet_2000_int5000_bilap_snap20260922_0948.h5"
URL = f"https://github.com/icepack/ismip7/releases/download/maps-2km-snap-2026-09-22/{BASENAME}"

FAKE_SBATCH = '''#!/bin/bash
printf 'CWD: %s\\n' "$PWD" >> "$SBATCH_CALLS"
env | grep '^ISMIP7_' | sort | sed 's/^/ENV: /' >> "$SBATCH_CALLS"
printf 'ARG: %s\\n' "$@" >> "$SBATCH_CALLS"
echo "4343;cluster"
'''


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    root = tmp_path / "repo" / "antarctica"
    (root / "mesh").mkdir(parents=True)
    (root / "results").mkdir()
    (root / "scripts").symlink_to(ANT / "scripts")
    for name in ("boundary_ids_antarctica_10000_1000_buffered20000.json",
                 "boundary_ids_antarctica_5000_2000_buffered0.json"):
        shutil.copy(ANT / "mesh" / name, root / "mesh" / name)
    (root / "mesh" / "antarctica_10000_1000_buffered20000.msh").write_text("mesh\n")
    download = root / "results/map_check/maps/download" / BASENAME
    download.parent.mkdir(parents=True)
    download.write_text("a released MAP\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "sbatch").write_text(FAKE_SBATCH)
    (bin_dir / "sbatch").chmod(0o755)
    (bin_dir / "hostname").write_text("#!/bin/bash\necho nowhere.example\n")
    (bin_dir / "hostname").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setenv("ISMIP7_LOCAL_ENV", os.devnull)
    monkeypatch.setenv("SBATCH_CALLS", str(tmp_path / "sbatch_calls.txt"))
    monkeypatch.setenv("ISMIP7_FIREDRAKE", str(download))
    for name in ("ISMIP7_SITE", "ISMIP7_REPO", "ISMIP7_CORES_PER_NODE",
                 "ISMIP7_MEM_PER_NODE", "ISMIP7_MESH", "ISMIP7_FRICTION"):
        monkeypatch.delenv(name, raising=False)
    return tmp_path


def manager(sandbox, *extra, dry_run=False):
    root = sandbox / "repo" / "antarctica"
    download = root / "results/map_check/maps/download" / BASENAME
    argv = [
        "run", "--root", str(root),
        "--map", str(root / "results/map_check/maps" / BASENAME),
        "--url", URL, "--md5", md5(download), "--friction", "regularized_coulomb",
        "--target-mesh", str(root / "mesh/antarctica_10000_1000_buffered20000.msh"),
        "--check-map-attrs", "0", "--k-npz", "/k/K_issue11_mesh2500.npz",
        *extra,
    ]
    if dry_run:
        argv.append("--dry-run")
    return mmc.MapCheckManager(mmc.parse_args(argv))


def calls(sandbox):
    path = sandbox / "sbatch_calls.txt"
    return path.read_text().splitlines() if path.exists() else []


def states(m):
    return {stage: state for stage, (state, _detail) in m.table().items()}


def test_the_tag_and_the_cache_are_named_apart_from_the_campaign(sandbox):
    m = manager(sandbox)
    assert m.tag.startswith("mapcheck_rc_scpc_gamg_10step_")
    with pytest.raises(ValueError):
        parse_campaign_tag(m.tag)
    assert not m.cache_path.name.startswith("initial_state_")
    assert m.run_tag("native") == "mapcheck_rc_snap20260922_0948_2000"
    assert m.run_tag("transferred") == "mapcheck_rc_snap20260922_0948_1000"


def test_nothing_done_means_fetch_is_the_one_runnable_stage(sandbox):
    got = states(manager(sandbox))
    assert got["fetch"] == "pending"
    assert got["summary"] == "pending"
    assert all(got[s] == "waiting" for s in mmc.STAGES if s not in ("fetch", "summary"))


def test_fetch_verifies_the_md5_and_then_repack_is_next(sandbox, capsys):
    m = manager(sandbox)
    m.run_fetch()
    record = json.loads((m.stem_dir / "fetch.json").read_text())
    assert record["ok"] and record["md5"] == m.md5
    got = states(m)
    assert got["fetch"] == "passed" and got["repack"] == "pending"
    assert got["score_native"] == "waiting" and got["lane_native"] == "waiting"


def test_a_wrong_md5_fails_the_fetch_and_blocks_everything_after(sandbox):
    m = manager(sandbox, "--md5", "0" * 32)
    m.run_fetch()
    got = states(m)
    assert got["fetch"] == "failed"
    assert got["repack"] == "blocked" and got["lane_native"] == "blocked"


def test_a_repacked_map_unlocks_the_three_independent_stages(sandbox):
    m = manager(sandbox)
    m.run_fetch()
    m.map_path.write_text("repacked\n")
    got = states(m)
    assert got["repack"] == "passed"
    assert {got[s] for s in ("score_native", "prepare_transfer", "lane_native")} == {"pending"}
    assert got["audit_cache"] == "waiting" and got["control_native"] == "waiting"


def test_the_first_pending_stage_is_submitted_through_the_site_file(sandbox, monkeypatch, capsys):
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    m = manager(sandbox)
    m.run_fetch()
    assert m.run() == 0
    seen = calls(sandbox)
    assert seen and seen[-1] == "ARG: scripts/batch_runners/timing_redistribute.script"
    assert "ARG: --ntasks-per-node=1" in seen and "ARG: -A" in seen
    status = mmc.read_status(m.status_path("repack"))
    assert status["state"] == "submitted" and status["job_id"] == "4343"
    # A second call finds the repack active (no squeue here to say otherwise)
    # and submits nothing more.
    (sandbox / "sbatch_calls.txt").unlink()
    m.run()
    assert calls(sandbox) == []
    assert "nothing to submit; active: repack" in capsys.readouterr().out


def test_a_control_goes_through_submit_sh_projection(sandbox, monkeypatch):
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    m = manager(sandbox)
    m.run_control("native")
    seen = calls(sandbox)
    assert seen[-1] == "ARG: antarctica/scripts/batch_runners/projection.sbatch"
    assert "ARG: --ntasks-per-node=64" in seen and "ARG: --mem=240G" in seen
    assert f"ARG: {m.control_job_name('native')}" in seen
    # submit.sh sets the KEY=VALUE pairs in sbatch's environment and passes a
    # bare --export=ALL, so the job's settings are read off the fake's record.
    env = [line[len("ENV: "):] for line in seen if line.startswith("ENV: ")]
    for pair in ("ISMIP7_EXPERIMENT=control", "ISMIP7_MESH=checkpoint",
                 "ISMIP7_FRICTION=regularized_coulomb", "ISMIP7_LC=2000",
                 "ISMIP7_T_END=2025", "ISMIP7_OUTPUT=1",
                 "ISMIP7_RUN_TAG=mapcheck_rc_snap20260922_0948_2000",
                 "ISMIP7_K_PER_BASIN_NPZ=/k/K_issue11_mesh2500.npz"):
        assert pair in env, pair
    assert any(line.startswith("ISMIP7_TRIPWIRE_U_MAX=") for line in env)
    assert "ARG: --export=ALL" in seen
    assert not any(line.startswith("ARG: --queue") for line in seen)
    status = mmc.read_status(m.status_path("control_native"))
    # The projection form prints its composed line on standard output ahead
    # of the id; the stamp carries the id alone.
    assert status["state"] == "submitted" and status["job_id"] == "4343"
    with pytest.raises(ValueError):
        m._submit("x", 1, "1G", {}, None, m.status_path("control_native"),
                  dependency="afterok:1", kind="projection")


def test_the_production_mesh_control_runs_the_production_step(sandbox, monkeypatch):
    r"""Issue 20: the control on the production mesh is the production
    configuration, step included. The strict lanes and the control on the
    MAP's own mesh keep the matrix's rule, and the audit divides by the step
    each control ran."""
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    m = manager(sandbox)
    assert m.lane_dt("transferred") == pytest.approx(0.05)
    assert m.control_dt("transferred") == pytest.approx(0.025)
    assert m.control_dt("native") == m.lane_dt("native") == pytest.approx(0.1)

    def env():
        seen = calls(sandbox)
        (sandbox / "sbatch_calls.txt").unlink()
        return [line[len("ENV: "):] for line in seen if line.startswith("ENV: ")]

    m.run_control("transferred")
    assert "ISMIP7_DT=0.025" in env()
    m.run_control("native")
    assert "ISMIP7_DT=0.1" in env()
    m.run_audit_controls()
    got = env()
    assert "ISMIP7_MAP_CHECK_DT_TRANSFER=0.025" in got
    assert "ISMIP7_MAP_CHECK_DT_NATIVE=0.1" in got


def _lane_record(m, role, **overrides):
    record = {
        "timing_kind": "map_check", "timing_tag": m.tag, "nsteps": 10,
        "completed_steps": 10, "t_end": 2015.0 + 10 * m.lane_dt(role),
        "dt": m.lane_dt(role), "ncores": 64, "seconds_per_step": 30.0,
        "initial_state_source": m.basename, "tripwire": None,
        "rescue_enabled": False, "step_mass_residual_gt_max": 0.0,
        "diagnostic_solve_summary": {"count": 10, "reason_counts": {"2": 10},
                                     "snes_iterations_total": 120},
        "transport_solve_summary": {"count": 10, "reason_counts": {"2": 10},
                                    "mass_residual_gt_max": 0.0},
    }
    record.update(overrides)
    m.stem_dir.mkdir(parents=True, exist_ok=True)
    m.lane_record(role).write_text(json.dumps(record))
    m.lane_csv(role).parent.mkdir(parents=True, exist_ok=True)
    with open(m.lane_csv(role), "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["year", "resid_gt"])
        for k in range(1, 11):
            writer.writerow([f"{2015.0 + k * m.lane_dt(role):.4f}", "0.0000"])


def test_a_lane_passes_the_qualification_rule_and_names_its_source(sandbox):
    m = manager(sandbox)
    _lane_record(m, "native")
    state, detail = m.verdict("lane_native")
    assert state == "passed" and "30.0 s/step" in detail and "12.0 Newton/step" in detail
    _lane_record(m, "native", initial_state_source="somebody_elses_map.h5")
    assert m.verdict("lane_native")[0] == "failed"
    _lane_record(m, "native", tripwire={"step": 3})
    assert m.verdict("lane_native")[0] == "failed"
    _lane_record(m, "native", completed_steps=9)
    assert m.verdict("lane_native")[0] == "failed"


def test_a_failed_lane_blocks_its_control_unless_the_switch_says_otherwise(sandbox):
    r"""Quartz, 22 September: both Budd lanes tripped the runaway tripwire at
    step 1 on the snapshot's own hotspots, and the controls, whose production
    configuration cancels that t = 0 tendency, were what the mesh question
    still needed. The switch runs the control on a failed lane's mesh and
    says so; the other control keeps waiting for a lane that has not run."""
    m = manager(sandbox)
    m.run_fetch()
    m.map_path.write_text("repacked\n")
    _lane_record(m, "native", tripwire={"step": 1}, completed_steps=1)
    assert states(m)["lane_native"] == "failed"
    assert states(m)["control_native"] == "blocked"
    told = manager(sandbox, "--controls-after-failed-lane")
    got = told.table()
    assert got["control_native"][0] == "pending"
    assert "lane_native failed" in got["control_native"][1]
    assert got["control_transfer"][0] == "waiting"
    assert got["audit_controls"][0] == "waiting"
    # Forcing the control alone leaves the lane failed (a resubmission with a
    # smaller request found the lane read as pending and the control waiting).
    forced = manager(sandbox, "--controls-after-failed-lane", "--force",
                     "--stages", "control_native")
    got = forced.table()
    assert got["lane_native"][0] == "failed"
    assert got["control_native"][0] == "pending"


def _control_csv(m, role, last_year, resid="0.0000"):
    m.control_csv(role).parent.mkdir(parents=True, exist_ok=True)
    dt = m.control_dt(role)
    with open(m.control_csv(role), "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["year", "vaf_mm_sle", "resid_gt"])
        year = 2015.0
        vaf = 0.0
        while year < last_year - 0.5 * dt:
            year += dt
            vaf += 0.01 * dt
            writer.writerow([f"{year:.4f}", f"{vaf:.4f}", resid])


def test_a_control_is_done_when_its_series_reaches_the_end_and_its_state_exists(sandbox):
    m = manager(sandbox)
    _control_csv(m, "transferred", 2020.0)
    state, detail = m.verdict("control_transfer")
    assert state is None and "of 2025" in detail
    _control_csv(m, "transferred", 2025.0)
    assert m.verdict("control_transfer")[0] is None
    m.control_final("transferred").write_text("final\n")
    state, detail = m.verdict("control_transfer")
    assert state == "passed" and "dVAF/dt +0.01" in detail
    _control_csv(m, "transferred", 2025.0, resid="0.5000")
    assert m.verdict("control_transfer")[0] == "failed"


def test_a_dry_run_that_assumes_its_dependencies_shows_every_stage(sandbox, monkeypatch, capsys):
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    m = manager(sandbox, "--assume-passed", dry_run=True)
    assert m.run() == 0
    out = capsys.readouterr().out
    for stage in mmc.STAGES:
        assert f"--- {stage}" in out
    assert out.count("DRY RUN: site iu_quartz: ") == 10
    assert "ISMIP7_MESH=checkpoint" in out
    assert f"ISMIP7_MESH={m.target_mesh}" in out
    assert "--ntasks-per-node=64" in out and "--ntasks-per-node=16" in out
    assert "ISMIP7_TIMING_KIND=map_check" in out
    assert calls(sandbox) == []
    assert not (m.stem_dir / "summary.md").exists()


def test_release_notes_give_the_md5_beside_the_basename():
    body = ("| `x.h5` | RC | 60 | 07:30 | `30d64e08ec9dd651d8a67eb00da3e30d` |\n"
            "| `y.h5` | Budd | 120 | 09:05 | `1c5d1031651f873e386c66aea34a69fe` |\n")
    fetch = lambda url: json.dumps({"body": body})  # noqa: E731
    assert mmc.release_md5(URL.replace(BASENAME, "y.h5"), "y.h5", fetch) == (
        "1c5d1031651f873e386c66aea34a69fe")
    assert mmc.release_md5(URL, "z.h5", fetch) is None
    assert mmc.release_md5("https://example.org/not-a-release", "y.h5", fetch) is None
