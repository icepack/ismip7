r"""The staged initialization routine's state machine.

``antarctica/scripts/initialization_routine.sh`` turns the document's stages
into submissions of the batch runners, each under the environment it recorded
at ``init`` plus the stage's own knobs, and hands itself on from the link that
finishes a stage. What it must get right is the order, the inputs (which MAP
or state each stage starts from, under which name), the knobs that may not
leak between stages, and idempotence: a call that finds nothing new to start
starts nothing. Those decisions run here against a stub ``sbatch`` that
records what it was given, with the stage products faked on disk.
"""

import os
import subprocess
import sys
from pathlib import Path

import h5py
import pytest

REPO = Path(__file__).resolve().parent.parent

SBATCH_STUB = (
    "#!/bin/bash\n"
    'printf "ARGV: %s\\n" "$*" >> "$SBATCH_CALLS"\n'
    'env | grep -E "^(ISMIP7_|SLURM_)" | sort | sed "s/^/ENV: /" >> "$SBATCH_CALLS"\n'
    'echo "END" >> "$SBATCH_CALLS"\n'
    'echo "$((RANDOM + 1000))"\n'
)


@pytest.fixture
def sandbox(tmp_path):
    r"""A checkout whose scripts are the repository's and whose results
    directory is the sandbox's own."""
    ant = tmp_path / "antarctica"
    ant.mkdir()
    (ant / "scripts").symlink_to(REPO / "antarctica" / "scripts")
    (tmp_path / "icepack2_tools").symlink_to(REPO / "icepack2_tools")
    (ant / "results").mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "sbatch").write_text(SBATCH_STUB)
    (bin_dir / "sbatch").chmod(0o755)
    (bin_dir / "python").write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (bin_dir / "python").chmod(0o755)
    (tmp_path / "out").mkdir()
    return tmp_path


def base_env(sandbox):
    calls = sandbox / "sbatch_calls.txt"
    return {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "SBATCH_CALLS": str(calls),
        "USER": "tester",
        # the objective the routine records
        "ISMIP7_LC": "32000",
        "ISMIP7_MESH": "/meshes/antarctica_320000_32000.msh",
        "ISMIP7_MISFIT_SCALE": "14000",
        "ISMIP7_DHDT_WEIGHT": "1",
        "ISMIP7_GRAD_PRECOND": "mass_consistent",
        "PYTHONPATH": str(sandbox),
        # a stage knob left in the caller's shell must not be recorded
        "ISMIP7_INVERT": "both",
        "ISMIP7_APPARENT_MB": "0",
        "ISMIP7_ROUTINE_LAUNCH": "slurm",
        "ISMIP7_ROUTINE_SBATCH_INV": "-p scavenge --time=01:00:00 --mem=240G",
        "ISMIP7_ROUTINE_SBATCH_FWD": "-p scavenge --time=01:00:00 --mem=180G",
        "ISMIP7_ROUTINE_ITERS": "40",
        "ISMIP7_ROUTINE_RELAX_YR": "20",
        "ISMIP7_ROUTINE_WARM_START": str(sandbox / "warm.h5"),
    }


def routine(sandbox, *args, env=None):
    merged = base_env(sandbox)
    merged.update(env or {})
    script = sandbox / "antarctica" / "scripts" / "initialization_routine.sh"
    proc = subprocess.run(
        ["bash", str(script), str(sandbox / "out"), *args], env=merged,
        cwd=str(sandbox), capture_output=True, text=True, timeout=120)
    return proc.returncode, proc.stdout + proc.stderr


def submits(sandbox):
    r"""One dict per recorded sbatch call: argv and the ISMIP7 environment,
    then the record is cleared."""
    calls = sandbox / "sbatch_calls.txt"
    if not calls.exists():
        return []
    out = []
    for block in calls.read_text().split("END\n"):
        if not block.strip():
            continue
        lines = block.strip().splitlines()
        env = dict(line[5:].split("=", 1) for line in lines[1:] if line.startswith("ENV: "))
        out.append({"argv": lines[0][6:], "env": env})
    calls.unlink()
    return out


def finish_inversion(path, its=40, **attrs):
    with h5py.File(path, "w") as h:
        h["/"].attrs["objective_iteration"] = its
        h["/"].attrs["objective_total"] = 1.0e6
        h["/"].attrs["misfit_scale"] = 14000.0
        h["/"].attrs["fluidity_control"] = "floating"
        for key, value in attrs.items():
            h["/"].attrs[key] = value
    Path(str(path) + ".done").touch()


def finish_forward(path, t_yr):
    with h5py.File(path, "w") as h:
        h["/"].attrs["t_yr"] = float(t_yr)
        h["/"].attrs["dt_yr"] = 0.05
        h["/"].attrs["stalled"] = 0
        h["/"].attrs["fssa_reference"] = "none"


def by_name(subs, name):
    found = [s for s in subs if f"-J {name} " in s["argv"] or s["argv"].endswith(f"-J {name}")]
    assert len(found) == 1, f"{name}: {[s['argv'] for s in subs]}"
    return found[0]


def test_the_routine_runs_itself_stage_by_stage(sandbox):
    out = sandbox / "out"
    results = sandbox / "antarctica" / "results"
    (sandbox / "warm.h5").write_text("warm\n")

    rc, log = routine(sandbox, "init")
    assert rc == 0, log
    env_text = (out / "routine.env").read_text()
    assert 'export ISMIP7_MISFIT_SCALE="14000"' in env_text
    assert "ISMIP7_INVERT" not in env_text, "stage knobs are never recorded"
    assert "ISMIP7_APPARENT_MB" not in env_text
    assert 'ISMIP7_ROUTINE_TAG=out' in env_text

    # init starts the first stage and nothing else: everything after it waits
    subs = submits(sandbox)
    assert [s["argv"].split(" -J ")[1].split()[0] for s in subs] == ["rt_out_shelf1"]
    s1 = subs[0]["env"]
    assert s1["ISMIP7_INVERT"] == "phi"
    assert s1["ISMIP7_PHI_BOUNDS"] == "-2.3,2.3"
    assert s1["ISMIP7_FLUIDITY_CONTROL"] == "floating"
    assert s1["ISMIP7_MAXITER"] == "40" and s1["ISMIP7_CHAIN_ITER_BUDGET"] == "40"
    assert s1["ISMIP7_MAP_OUT"] == str(out / "map_shelf1.h5")
    assert s1["ISMIP7_WARM_START"] == str(sandbox / "warm.h5")
    assert s1["ISMIP7_WARM_START_STRICT"] == "0"
    assert s1["ISMIP7_CHAIN_THEN_CMD"].endswith(f"initialization_routine.sh {out} next")
    assert "ISMIP7_DHDT_REWEIGHT" not in s1
    assert "SLURM_JOB_ID" not in s1
    assert "inversion.sbatch" in subs[0]["argv"]
    assert "--time=01:00:00 --mem=240G" in subs[0]["argv"]

    # a second call finds the stage started and starts nothing
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    assert submits(sandbox) == []

    # the friction stage warm-starts from the shelf MAP and re-weights the fit
    finish_inversion(out / "map_shelf1.h5")
    rc, log = routine(sandbox, "next", env={"ISMIP7_INVERT": "phi", "SLURM_JOB_ID": "77"})
    assert rc == 0, log
    (s2,) = submits(sandbox)
    assert "-J rt_out_friction1" in s2["argv"]
    assert s2["env"]["ISMIP7_INVERT"] == "theta"
    assert s2["env"]["ISMIP7_DHDT_REWEIGHT"] == "24"
    assert "ISMIP7_PHI_BOUNDS" not in s2["env"]
    assert s2["env"]["ISMIP7_WARM_START"] == str(out / "map_shelf1.h5")
    assert s2["env"]["ISMIP7_DHDT_WEIGHT"] == "1", "the recorded objective reaches every stage"
    assert "SLURM_JOB_ID" not in s2["env"], "the calling job's scheduler variables stay behind"

    finish_inversion(out / "map_friction1.h5")
    routine(sandbox, "next")
    (s3,) = submits(sandbox)
    assert "-J rt_out_shelf2" in s3["argv"] and s3["env"]["ISMIP7_MAXITER"] == "20"
    finish_inversion(out / "map_shelf2.h5")
    routine(sandbox, "next")
    (s4,) = submits(sandbox)
    assert "-J rt_out_friction2" in s4["argv"] and s4["env"]["ISMIP7_CHAIN_ITER_BUDGET"] == "20"

    # the relaxation: the control driver, apparent MB off, front pinned
    finish_inversion(out / "map_friction2.h5", its=20, misfit_norm="sigma")
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    (relax,) = submits(sandbox)
    assert "-J rt_out_relax" in relax["argv"] and "projection.sbatch" in relax["argv"]
    e = relax["env"]
    assert e["ISMIP7_EXPERIMENT"] == "control" and e["ISMIP7_ESM"] == "CESM2-WACCM"
    assert e["ISMIP7_INVERSION"] == str(out / "map_friction2.h5")
    assert e["ISMIP7_RUN_TAG"] == "out_relax" and e["ISMIP7_T_END"] == "2035"
    assert e["ISMIP7_APPARENT_MB"] == "0" and e["ISMIP7_OUTPUT"] == "0"
    assert e["ISMIP7_FIXED_FRONT"] == "1" and e["ISMIP7_CALVING"] == "none"
    assert "ISMIP7_INVERT" not in e and "ISMIP7_MAP_OUT" not in e
    assert "--mem=180G" in relax["argv"]

    # short of its end year the relaxation is still in progress
    final = results / "ctrl2015_cesm2_waccm_out_relax_32000_final.h5"
    finish_forward(final, 2030)
    routine(sandbox, "next")
    assert submits(sandbox) == []
    rc, log = routine(sandbox, "status")
    assert rc == 0 and "relax                in progress (t=2030.0 of 2035)" in log, log

    # at its end the state is stamped and the re-inversion starts from it,
    # under the name the relaxation module requires
    finish_forward(final, 2035)
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    state = out / "relaxed_state.h5"
    assert state.exists()
    with h5py.File(state, "r") as h:
        a = h["/"].attrs
        assert int(a["relaxation_end_state"]) == 1
        assert float(a["relax_t_start"]) == 2015.0 and float(a["relax_t_end"]) == 2035.0
        assert os.path.basename(a["relax_source_map"]) == "map_friction2.h5"
        assert float(a["misfit_scale"]) == 14000.0, "the MAP's objective travels with the state"
        assert float(a["relax_dt"]) == 0.05
    (r1,) = submits(sandbox)
    assert "-J rt_out_r_shelf1" in r1["argv"]
    assert r1["env"]["ISMIP7_WARM_START"] == str(state)
    assert r1["env"]["ISMIP7_MAP_OUT"] == str(out / "map_friction2_relax2015.h5")
    assert r1["env"]["ISMIP7_INVERT"] == "phi"

    finish_inversion(out / "map_friction2_relax2015.h5")
    routine(sandbox, "next")
    (r2,) = submits(sandbox)
    assert "-J rt_out_r_friction1" in r2["argv"]
    assert r2["env"]["ISMIP7_WARM_START"] == str(out / "map_friction2_relax2015.h5")
    finish_inversion(out / "map_r_friction1.h5")
    routine(sandbox, "next")
    (r3,) = submits(sandbox)
    assert "-J rt_out_r_shelf2" in r3["argv"]
    finish_inversion(out / "map_r_shelf2.h5")
    routine(sandbox, "next")
    (r4,) = submits(sandbox)
    assert "-J rt_out_r_friction2" in r4["argv"]

    # the last MAP is the forward MAP; the historicals (one per ESM) and OCX
    # start from it with the apparent MB off and the submission's front
    finish_inversion(out / "map_r_friction2.h5", its=20)
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    assert os.readlink(out / "map_forward.h5") == "map_r_friction2.h5"
    subs = submits(sandbox)
    names = sorted(s["argv"].split(" -J ")[1].split()[0] for s in subs)
    assert names == ["rt_out_hist_cesm_waccm", "rt_out_hist_mri_esm2", "rt_out_ocx"]
    h = by_name(subs, "rt_out_hist_mri_esm2")["env"]
    assert h["ISMIP7_EXPERIMENT"] == "hist_mri_esm2" and h["ISMIP7_ESM"] == "MRI-ESM2-0"
    assert h["ISMIP7_INVERSION"] == str(out / "map_forward.h5")
    assert h["ISMIP7_RUN_TAG"] == "out" and h["ISMIP7_T_START"] == "2002"
    assert h["ISMIP7_APPARENT_MB"] == "0" and h["ISMIP7_OUTPUT"] == "1"
    assert h["ISMIP7_CALVING"] == "hfb" and h["ISMIP7_FRONT_ADVANCE"] == "none"
    assert h["ISMIP7_FRACTURE"] == "none"
    assert h["ISMIP7_GEOMETRY_BACKDATE"] == "12", "the forced 2002 year starts from the 2003 geometry"
    assert "ISMIP7_FIXED_FRONT" not in h and "ISMIP7_T_END" not in h
    assert "ISMIP7_INVERT" not in h and "ISMIP7_FLUIDITY_CONTROL" not in h

    # one historical at its end releases its own control and projections only
    finish_forward(results / "hist_cesm2_waccm_out_32000_final.h5", 2015)
    finish_forward(results / "hist_mri_esm2_0_out_32000_final.h5", 2010)
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    subs = submits(sandbox)
    names = sorted(s["argv"].split(" -J ")[1].split()[0] for s in subs)
    assert names == ["rt_out_control_cesm_waccm", "rt_out_ssp126_cesm_waccm",
                     "rt_out_ssp370_cesm_waccm", "rt_out_ssp585_cesm_waccm"]
    c = by_name(subs, "rt_out_control_cesm_waccm")["env"]
    assert c["ISMIP7_EXPERIMENT"] == "control" and c["ISMIP7_ESM"] == "CESM2-WACCM"
    # no reference: the historical carries none and one cannot be built
    # from its evolved endpoint (simulation.setup_model refuses)
    assert c["ISMIP7_APPARENT_MB"] == "0" and c["ISMIP7_RUN_TAG"] == "out"
    assert c["ISMIP7_FRACTURE"] == "none" and "ISMIP7_T_START" not in c
    assert "ISMIP7_GEOMETRY_BACKDATE" not in c
    p = by_name(subs, "rt_out_ssp585_cesm_waccm")["env"]
    assert p["ISMIP7_EXPERIMENT"] == "ssp585_cesm_waccm"
    assert p["ISMIP7_APPARENT_MB"] == "0" and p["ISMIP7_FRACTURE"] == "mask_front"
    # the law's own grounded rule everywhere unless the routine says otherwise
    assert "ISMIP7_CALVING_FLOATING_ONLY" not in p and "ISMIP7_CALVING_FLOATING_ONLY" not in c

    finish_forward(results / "hist_mri_esm2_0_out_32000_final.h5", 2015)
    routine(sandbox, "next")
    names = sorted(s["argv"].split(" -J ")[1].split()[0] for s in submits(sandbox))
    assert names == ["rt_out_control_mri_esm2", "rt_out_ssp126_mri_esm2",
                     "rt_out_ssp370_mri_esm2", "rt_out_ssp585_mri_esm2"]

    rc, log = routine(sandbox, "next")
    assert rc == 0 and submits(sandbox) == []
    rc, log = routine(sandbox, "status")
    assert rc == 0, log
    assert "r_friction2          done (20 its" in log
    assert "hist_cesm_waccm      done (t=2015.0)" in log
    assert "ssp585_mri_esm2      started, no checkpoint yet" in log


def test_init_refuses_to_overwrite_a_recorded_routine(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, _ = routine(sandbox, "init")
    assert rc == 0
    submits(sandbox)
    rc, log = routine(sandbox, "init")
    assert rc != 0 and "exists" in log
    assert submits(sandbox) == []


def test_slurm_mode_needs_its_resources(sandbox):
    rc, log = routine(sandbox, "init", env={"ISMIP7_ROUTINE_SBATCH_FWD": ""})
    assert rc != 0 and "ISMIP7_ROUTINE_SBATCH_FWD" in log


def test_forwards_alone_run_from_the_map_given(sandbox):
    r"""Without an inversion stage the forwards start from
    ISMIP7_ROUTINE_WARM_START: the routine as a plain forward suite."""
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={
        "ISMIP7_ROUTINE_STAGES": "hist control", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM"})
    assert rc == 0, log
    assert os.readlink(sandbox / "out" / "map_forward.h5") == str(sandbox / "warm.h5")
    (h,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in h["argv"]
    assert h["env"]["ISMIP7_INVERSION"] == str(sandbox / "out" / "map_forward.h5")


def test_resume_submits_a_started_stage_again(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    routine(sandbox, "init")
    submits(sandbox)
    rc, log = routine(sandbox, "resume", "shelf1")
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_shelf1" in again["argv"]
    assert again["env"]["ISMIP7_CHAIN_ITER_BUDGET"] == "40"
    rc, log = routine(sandbox, "resume", "friction1")
    assert rc != 0 and "not ready" in log


def test_a_2003_historical_keeps_the_default_backdate(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={
        "ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
        "ISMIP7_ROUTINE_HIST_START": "2003"})
    assert rc == 0, log
    (h,) = submits(sandbox)
    assert h["env"]["ISMIP7_T_START"] == "2003"
    assert "ISMIP7_GEOMETRY_BACKDATE" not in h["env"]


@pytest.mark.parametrize("precond", ["none", "mass", ""])
def test_init_refuses_inversion_stages_off_the_tao_path(sandbox, precond):
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_GRAD_PRECOND": precond}
    rc, log = routine(sandbox, "init", env=env)
    assert rc != 0 and "ISMIP7_GRAD_PRECOND" in log, log
    assert submits(sandbox) == []
    assert not (sandbox / "out" / "routine.env").exists(), "a refused init records nothing"


def test_forwards_alone_need_no_tao_metric(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={
        "ISMIP7_GRAD_PRECOND": "none", "ISMIP7_ROUTINE_STAGES": "hist",
        "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM"})
    assert rc == 0, log
    assert len(submits(sandbox)) == 1


def test_a_lock_whose_holder_is_gone_is_taken(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    routine(sandbox, "init")
    submits(sandbox)
    dead = subprocess.Popen(["true"])
    dead.wait()
    lock = sandbox / "out" / ".routine.lock"
    lock.mkdir()
    (lock / "holder").write_text(f"{dead.pid} {os.uname().nodename}\n")
    finish_inversion(sandbox / "out" / "map_shelf1.h5")
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    (s,) = submits(sandbox)
    assert "-J rt_out_friction1" in s["argv"]
    assert not lock.exists()


def test_a_live_lock_is_waited_for(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    routine(sandbox, "init")
    submits(sandbox)
    lock = sandbox / "out" / ".routine.lock"
    lock.mkdir()
    (lock / "holder").write_text(f"{os.getpid()} {os.uname().nodename}\n")
    finish_inversion(sandbox / "out" / "map_shelf1.h5")
    with pytest.raises(subprocess.TimeoutExpired):
        merged = base_env(sandbox)
        subprocess.run(
            ["bash", str(sandbox / "antarctica" / "scripts" / "initialization_routine.sh"),
             str(sandbox / "out"), "next"],
            env=merged, cwd=str(sandbox), capture_output=True, text=True, timeout=8)
    assert submits(sandbox) == []
    assert (lock / "holder").read_text().split()[0] == str(os.getpid())


def test_watch_resubmits_a_dead_chain_and_cancels_a_silent_link(sandbox):
    r"""A started stage with no job left in the queue is resubmitted (a link
    killed by the scheduler queues no successor); a running forward link whose
    log has been silent longer than the stall limit is cancelled and its node
    excluded from the stage's next submission."""
    import time
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    queue = sandbox / "queue.txt"
    (bin_dir / "squeue").write_text(
        "#!/bin/bash\n"
        'name=""; while [ $# -gt 0 ]; do [ "$1" = -n ] && name=$2; shift; done\n'
        '[ -f "$QUEUE_FILE" ] && grep "^$name " "$QUEUE_FILE" | cut -d" " -f2-\n'
        "exit 0\n")
    (bin_dir / "scancel").write_text('#!/bin/bash\necho "scancel $*" >> "$SBATCH_CALLS"\n')
    for name in ("squeue", "scancel"):
        (bin_dir / name).chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist control", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "QUEUE_FILE": str(queue), "ISMIP7_ROUTINE_WATCH_ONCE": "1",
           "ISMIP7_ROUTINE_STALL_MIN": "15"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in h["argv"]

    # the chain died: nothing queued under the stage's name
    queue.write_text("")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert "--exclude" not in again["argv"]
    assert "no job in the queue; resubmitting" in log

    # a running link whose log went silent 20 minutes ago
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    silent = logs / "ismip7_fwd_4242.out"
    silent.write_text("t=2003.0\n")
    old = time.time() - 20 * 60
    os.utime(silent, (old, old))
    queue.write_text("rt_out_hist_cesm_waccm 4242 RUNNING bb5u20c1 0:30:00\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert "silent for 20 min" in log
    calls = (sandbox / "sbatch_calls.txt").read_text()
    assert "scancel 4242" in calls
    assert (sandbox / "out" / "exclude_nodes").read_text().strip() == "bb5u20c1"
    # the cancelled job is gone on the next pass: the stage is resubmitted
    # without that node
    queue.write_text("")
    (sandbox / "sbatch_calls.txt").unlink()
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "--exclude=bb5u20c1" in again["argv"]

    # a running link with a fresh log is left alone
    silent.write_text("t=2004.0\n")
    queue.write_text("rt_out_hist_cesm_waccm 4243 RUNNING bb5u21c1 0:10:00\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == []
    assert "scancel 4243" not in (sandbox / "sbatch_calls.txt").read_text() if (sandbox / "sbatch_calls.txt").exists() else True


def test_watch_stamps_the_node_exclusions_onto_a_pending_link_that_has_none(sandbox):
    r"""A chain submitted before the watchdog excluded a node hands every
    successor an empty exclusion list (site_core.sh copies the predecessor's),
    so the watch stamps the stage's current list onto a pending link that
    carries none and leaves one that already has a list alone."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    queue = sandbox / "queue.txt"
    (bin_dir / "squeue").write_text(
        "#!/bin/bash\n"
        'name=""; while [ $# -gt 0 ]; do [ "$1" = -n ] && name=$2; shift; done\n'
        '[ -f "$QUEUE_FILE" ] && grep "^$name " "$QUEUE_FILE" | cut -d" " -f2-\n'
        "exit 0\n")
    (bin_dir / "scontrol").write_text(
        "#!/bin/bash\n"
        'if [ "$1" = show ]; then\n'
        '  if grep -qx "$3" "$NULL_FILE"; then echo "JobId=$3 ExcNodeList=(null) Features=x"; '
        'else echo "JobId=$3 ExcNodeList=bb1u1c1 Features=x"; fi\n'
        'else echo "scontrol $*" >> "$SBATCH_CALLS"; fi\n')
    for name in ("squeue", "scontrol"):
        (bin_dir / name).chmod(0o755)
    nulls = sandbox / "nulls.txt"
    nulls.write_text("4300\n")
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "QUEUE_FILE": str(queue), "NULL_FILE": str(nulls), "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (sandbox / "sbatch_calls.txt").unlink(missing_ok=True)
    queue.write_text("rt_out_hist_cesm_waccm 4300 PENDING (null) 0:00 Priority\n"
                     "rt_out_hist_cesm_waccm 4301 PENDING (null) 0:00 Priority\n")

    # nothing excluded yet: nothing to stamp
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert not (sandbox / "sbatch_calls.txt").exists()

    (sandbox / "out" / "exclude_nodes").write_text("bb6u15g1\nbb6u13g1\nbb6u15g1\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    calls = (sandbox / "sbatch_calls.txt").read_text()
    assert calls == "scontrol update JobId=4300 ExcNodeList=bb6u13g1,bb6u15g1\n"
    assert "job 4300 (hist_cesm_waccm) pending without node exclusions" in log
    assert "4301" not in log


def test_the_routine_can_confine_the_law_to_floating_ice_on_the_projections(sandbox):
    r"""``ISMIP7_ROUTINE_CALVING_FLOATING_ONLY=ssp`` hands the projections
    ``ISMIP7_CALVING_FLOATING_ONLY=1`` and nothing else (the historical and
    control keep the law's own rule, as their checkpoints record it);
    ``all`` gates every forward; anything else is refused."""
    (sandbox / "warm.h5").write_text("warm\n")
    results = sandbox / "antarctica" / "results"
    env = {"ISMIP7_ROUTINE_STAGES": "hist control proj",
           "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM", "ISMIP7_ROUTINE_SCENARIOS": "ssp126",
           "ISMIP7_ROUTINE_CALVING_FLOATING_ONLY": "ssp"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    assert "ISMIP7_CALVING_FLOATING_ONLY" not in h["env"]
    assert "ISMIP7_ROUTINE_CALVING_FLOATING_ONLY=ssp" in (sandbox / "out" / "routine.env").read_text()
    finish_forward(results / "hist_cesm2_waccm_out_32000_final.h5", 2015)
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    subs = submits(sandbox)
    assert "ISMIP7_CALVING_FLOATING_ONLY" not in by_name(subs, "rt_out_control_cesm_waccm")["env"]
    assert by_name(subs, "rt_out_ssp126_cesm_waccm")["env"]["ISMIP7_CALVING_FLOATING_ONLY"] == "1"

    env["ISMIP7_ROUTINE_CALVING_FLOATING_ONLY"] = "all"
    (sandbox / "out" / "routine.env").unlink()
    (sandbox / "out" / "control_cesm_waccm.jobs").unlink()
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    assert by_name(submits(sandbox), "rt_out_control_cesm_waccm")["env"]["ISMIP7_CALVING_FLOATING_ONLY"] == "1"

    env["ISMIP7_ROUTINE_CALVING_FLOATING_ONLY"] = "projections"
    (sandbox / "out" / "routine.env").unlink()
    rc, log = routine(sandbox, "init", env=env)
    assert rc != 0 and "ISMIP7_ROUTINE_CALVING_FLOATING_ONLY must be none, ssp or all" in log


def test_watch_does_nothing_when_the_scheduler_cannot_be_reached(sandbox):
    r"""An unreachable controller answers nothing; that is not an empty queue.
    The 14 October 2026 outage resubmitted every stage on top of its live chain."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text('#!/bin/bash\necho "squeue: error: Unable to contact slurm controller" >&2\nexit 1\n')
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert "squeue failed; skipping" in log
    assert submits(sandbox) == [], "nothing may be resubmitted on a scheduler failure"



def test_watch_parks_a_stage_whose_driver_failed(sandbox):
    r"""A link that exits with a driver error ends its chain like a scheduler
    kill does, but a resubmission only repeats the error: on 15 October 2026
    each control and projection of two routines was resubmitted about 60
    times, two minutes apart, on an a_ref refusal. The watchdog parks such a
    stage with the error, status shows it, and resume clears it."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    queue = sandbox / "queue.txt"
    (bin_dir / "squeue").write_text(
        "#!/bin/bash\n"
        'name=""; while [ $# -gt 0 ]; do [ "$1" = -n ] && name=$2; shift; done\n'
        '[ -f "$QUEUE_FILE" ] && grep "^$name " "$QUEUE_FILE" | cut -d" " -f2-\n'
        "exit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "QUEUE_FILE": str(queue), "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    jid = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    (logs / f"ismip7_fwd_{jid}.out").write_text(
        "=== forward exit 1 ===\nchain: stopping, driver exited 1\n")
    (logs / f"ismip7_fwd_{jid}.err").write_text(
        "Traceback\nRuntimeError: ISMIP7_APPARENT_MB is set but restart checkpoint x has no a_ref_mb\n")
    queue.write_text("")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == [], "a driver failure is not resubmitted"
    assert "stopped on a driver failure" in log
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert "no a_ref_mb" in blocked
    rc, log = routine(sandbox, "watch", env=env)
    assert submits(sandbox) == [], "a parked stage stays parked"
    rc, log = routine(sandbox, "status", env=env)
    assert "BLOCKED" in log and "no a_ref_mb" in log
    rc, log = routine(sandbox, "resume", "hist_cesm_waccm", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists()


def _stalled_link(sandbox, resumed, stalled):
    r"""A stage whose newest forward link resumed at ``resumed`` and stalled
    at ``stalled``, as projection.sbatch and the driver log it."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text("#!/bin/bash\nexit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    jid = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    (logs / f"ismip7_fwd_{jid}.out").write_text(
        f"  Restart: evolved geometry + frozen anchors loaded (t_yr={resumed}, friction=budd)\n"
        "  Step 1: rescue ladder + subcycles exhausted, saving and stopping\n"
        "=== 2026-10-16T12:11:54-05:00  forward exit 0 ===\n"
        f"chain: the solver stalled (rescue ladder + subcycles exhausted) at {stalled}.\n"
        "chain: NOT resubmitting. A stall is not a run that needs more time;\n")
    return env, jid


def test_watch_parks_a_stage_that_stalled_where_it_resumed(sandbox):
    r"""A forward whose solver stalled saves its state, exits 0 and does not
    chain. The watchdog used to read that as a link the scheduler took and
    resubmit it; a link that stalls where it resumed makes no progress, and on
    16 October 2026 the seq0 ssp126 projection was relaunched link after link,
    each stalling in its first step at 2019.95. Parked, with the year."""
    env, jid = _stalled_link(sandbox, "2019.95", "2019.95")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == [], "a stall without progress is not resubmitted"
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert blocked.strip() == f"job {jid}: the solver stalled at 2019.95, where the link resumed"
    rc, log = routine(sandbox, "resume", "hist_cesm_waccm", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]


def test_watch_resubmits_a_stage_that_advanced_before_it_stalled(sandbox):
    r"""A link that advanced and then stalled is resubmitted: the next link
    starts the substep controller afresh at the stall, and the seq0 ssp585
    projection's restart passed the step its link (2015.0 to 2021.25) had
    stalled in."""
    env, _ = _stalled_link(sandbox, "2015.0", "2021.25")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists()


def test_a_scheduler_kill_is_still_resubmitted(sandbox):
    r"""A link the scheduler took (no driver failure in its log) is resubmitted."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text("#!/bin/bash\nexit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    routine(sandbox, "init", env=env)
    submits(sandbox)
    jid = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    (logs / f"ismip7_fwd_{jid}.out").write_text("  t=2005.3 ...\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]


SACCT_LINKS_STUB = (
    "#!/bin/bash\n"
    'for a in "$@"; do case "$a" in --name=*) cat "$LINKS"; exit 0 ;; esac; done\n'
    'echo "FAILED"\n'
)
PBLAS = "PBLAS ERROR 'Illegal IY = 0, IY must be at least 1'\n"
ABORTED = "=== forward exit 137 ===\nchain: stopping, driver exited 137\n"
CLEAN = "  t=2004.2 ...\n=== forward exit 0 ===\nchain: at 2004.25 of 2015.0, resubmitting\n"


def _solver_abort_stage(sandbox, strict=True, sacct=SACCT_LINKS_STUB):
    r"""A started historical stage whose accounting lists the job ids in
    sandbox/links (what sacct answers by job name inside its window), and
    the stage's first link."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text("#!/bin/bash\nexit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    (bin_dir / "sacct").write_text(sacct)
    (bin_dir / "sacct").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1", "ISMIP7_ROUTINE_STRICT": "1" if strict else "0",
           "LINKS": str(sandbox / "links")}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    (sandbox / "logs").mkdir(exist_ok=True)
    return env, _first_link(sandbox)


def _first_link(sandbox):
    return (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]


def _link(sandbox, jid, out, err=None):
    logs = sandbox / "logs"
    (logs / f"ismip7_fwd_{jid}.out").write_text(out)
    if err is not None:
        (logs / f"ismip7_fwd_{jid}.err").write_text(err)


def _links(sandbox, *jids):
    (sandbox / "links").write_text("".join(f"{j}\n" for j in jids))


def test_a_strict_park_names_the_link_that_actually_ended(sandbox):
    r"""A chain's own successors are submitted by projection.sbatch, not by
    the routine, so OUT/<stage>.jobs holds only the first link. The watchdog
    read that line as the stage's newest link and on 21 October 2026 parked
    seq0f's MRI historical "COMPLETED without a successor", naming a link
    that had ended cleanly two links before the one that aborted in MUMPS.
    The accounting knows every link by job name. The newest link's solver
    abort is retried once; the retry's own abort, the second in a row, parks
    the stage naming both links, though the accounting's window (opened at
    the rotated OUT/<stage>.jobs) no longer lists the first."""
    env, first = _solver_abort_stage(sandbox)
    _links(sandbox, first, "9999777")
    _link(sandbox, first, CLEAN)
    _link(sandbox, "9999777", "Restart state must be re-solved: ||F||=1.51e+08\n" + ABORTED,
          PBLAS + "MPI_ABORT was invoked on rank 5 in communicator MPI_COMM_WORLD\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert (sandbox / "out" / "hist_cesm_waccm.solver_retry").read_text().strip() == "9999777"

    retry = _first_link(sandbox)
    _links(sandbox, retry)
    _link(sandbox, retry, ABORTED,
          "PBLAS ERROR 'Array subscript out of bounds: IY = 1793, DESCY[M_] = 908'\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == [], "a second solver abort in a row is not resubmitted"
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert blocked.startswith(f"job {retry}:"), blocked
    assert "twice in a row" in blocked and "after 9999777" in blocked
    assert f"ismip7_fwd_{retry}.err" in blocked

    rc, log = routine(sandbox, "resume", "hist_cesm_waccm", env=env)
    assert rc == 0, log
    assert not (sandbox / "out" / "hist_cesm_waccm.solver_retry").exists()


def test_a_link_that_aborted_in_the_linear_solver_is_resubmitted_once(sandbox):
    r"""A ScaLAPACK or MUMPS abort in a link's first factorizations is a node
    fault (every October 2026 production abort was on a shared node), so the
    stage is resubmitted once, under strict mode too, and parks only on a
    second such abort in a row (the test above). A driver exit with another
    error still parks at once."""
    env, first = _solver_abort_stage(sandbox)
    _links(sandbox, first, "9999777")
    _link(sandbox, first, CLEAN)
    _link(sandbox, "9999777", ABORTED,
          "[ 6] MUMPS error in solve: INFOG(1)=-3, INFO(2)=3\n"
          "petsc4py.PETSc.Error: error code 76\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists(), log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert "aborted in the linear solver" in log and "resubmitting once" in log

    # the same link ending on an ordinary driver error parks at once
    (sandbox / "out" / "hist_cesm_waccm.jobs").write_text(f"{first} 2026-10-24T10:00:00-05:00\n")
    (sandbox / "out" / "hist_cesm_waccm.solver_retry").unlink()
    _link(sandbox, "9999777", ABORTED, "RuntimeError: ISMIP7_FRACTURE=mask_front needs a mask\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert blocked.startswith("job 9999777: RuntimeError"), blocked


def test_an_abort_after_a_clean_link_gets_its_own_retry(sandbox):
    r"""The retry is spent only by an abort in a row: once a later link of
    the stage ended cleanly (the retry chained on to a successor), the next
    isolated solver abort is resubmitted once again."""
    env, first = _solver_abort_stage(sandbox, strict=False)
    _links(sandbox, first)
    _link(sandbox, first, ABORTED, PBLAS)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert len(submits(sandbox)) == 1
    marker = sandbox / "out" / "hist_cesm_waccm.solver_retry"
    assert marker.read_text().strip() == first

    retry = _first_link(sandbox)
    _links(sandbox, retry, "9999777")
    _link(sandbox, retry, CLEAN)
    _link(sandbox, "9999777", ABORTED, "[ 3] MUMPS error in solve: INFOG(1)=-3\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert len(submits(sandbox)) == 1, log
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists(), log
    assert marker.read_text().strip() == "9999777"


def test_a_solver_abort_retry_holds_without_the_accounting(sandbox):
    r"""Where sacct answers nothing, the newest link is the stage's
    OUT/<stage>.jobs entry, and the retry is still granted once: the
    resubmitted link's abort parks the stage naming both links."""
    env, first = _solver_abort_stage(sandbox, strict=False, sacct="#!/bin/bash\nexit 0\n")
    _link(sandbox, first, ABORTED, PBLAS)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert len(submits(sandbox)) == 1
    assert "resubmitting once" in log
    assert "in the linear solver on an unknown node;" in log, log

    retry = _first_link(sandbox)
    _link(sandbox, retry, ABORTED, PBLAS)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == []
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert blocked.startswith(f"job {retry}:") and f"after {first}" in blocked, blocked


def test_a_strict_routine_resubmits_a_link_the_watchdog_cancelled(sandbox):
    r"""The accounting prints a cancelled link as ``CANCELLED+`` (the "by
    <uid>" is cut off). A strict routine must read that as the scheduler's
    doing and resubmit, as it does for ``TIMEOUT``: on 22 October 2026 the
    watchdog cancelled dg1f's hung ssp126 link and then parked the stage on
    its own cancellation."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text("#!/bin/bash\nexit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    (bin_dir / "sacct").write_text(
        "#!/bin/bash\n"
        'for a in "$@"; do case "$a" in --name=*) echo "$FIRST"; exit 0 ;; esac; done\n'
        'echo "CANCELLED+"\n')
    (bin_dir / "sacct").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1", "ISMIP7_ROUTINE_STRICT": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    first = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    env["FIRST"] = first
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    (logs / f"ismip7_fwd_{first}.out").write_text("Restart state must be re-solved: ||F||=1.5e+08\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists(), log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert "no job in the queue; resubmitting" in log



def _watched_stage(sandbox, sacct_state=None, **extra):
    r"""A started historical under the watchdog's fakes: ``squeue`` answers
    from QUEUE_FILE (``<name> <row>`` lines), ``scancel`` records its
    arguments, and ``sacct``, when a state is given, names the first link as
    the newest and reports that state for it."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    queue = sandbox / "queue.txt"
    cancels = sandbox / "scancel_calls.txt"
    (bin_dir / "squeue").write_text(
        "#!/bin/bash\n"
        'name=""; while [ $# -gt 0 ]; do [ "$1" = -n ] && name=$2; shift; done\n'
        '[ -f "$QUEUE_FILE" ] && grep "^$name " "$QUEUE_FILE" | cut -d" " -f2-\n'
        "exit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    (bin_dir / "scancel").write_text('#!/bin/bash\necho "$*" >> "$SCANCEL_CALLS"\n')
    (bin_dir / "scancel").chmod(0o755)
    if sacct_state is not None:
        (bin_dir / "sacct").write_text(
            "#!/bin/bash\n"
            'for a in "$@"; do case "$a" in --name=*) echo "$FIRST"; exit 0 ;; esac; done\n'
            f'echo "{sacct_state}"\n')
        (bin_dir / "sacct").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "QUEUE_FILE": str(queue), "SCANCEL_CALLS": str(cancels),
           "ISMIP7_ROUTINE_WATCH_ONCE": "1", **extra}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    first = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    env["FIRST"] = first
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    return env, first, queue, cancels, logs


def _frozen_log(logs, jid, zeros):
    (logs / f"ismip7_fwd_{jid}.out").write_text(
        "Restart state re-solved\n"
        + "".join(f"  momentum: snes_its=0 fnorm=1.2e+03 step {k}\n" for k in range(zeros))
        + "  momentum: snes_its=4 fnorm=3.1e-02 step x\n")


def test_watch_parks_a_stage_whose_velocity_is_frozen(sandbox):
    r"""A link whose momentum solves "converge" in zero iterations carries the
    restart's velocity forward unchanged (the frozen controls of 19-21 October
    2026). At ISMIP7_ROUTINE_FROZEN_MAX such solves the watchdog cancels the
    running link and parks the stage, strict or not; below it, nothing."""
    env, first, queue, cancels, logs = _watched_stage(sandbox)
    queue.write_text(f"rt_out_hist_cesm_waccm {first} RUNNING bb6u1g1 10:00 None\n")
    _frozen_log(logs, first, 4)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert not cancels.exists(), "four zero-iteration solves are below the default of five"
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists()
    _frozen_log(logs, first, 5)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert cancels.read_text().split() == [first]
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert blocked.startswith(f"job {first}: 5 zero-iteration momentum solves")
    queue.write_text("")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == [], "a frozen stage stays parked"


def test_a_frozen_count_of_zero_disables_the_check(sandbox):
    env, first, queue, cancels, logs = _watched_stage(sandbox, ISMIP7_ROUTINE_FROZEN_MAX="0")
    queue.write_text(f"rt_out_hist_cesm_waccm {first} RUNNING bb6u1g1 10:00 None\n")
    _frozen_log(logs, first, 0)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    _frozen_log(logs, first, 50)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert not cancels.exists()
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists()


@pytest.mark.parametrize("value", ["-1", "five", "2.5"])
def test_a_bad_frozen_count_is_refused(sandbox, value):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={"ISMIP7_ROUTINE_FROZEN_MAX": value})
    assert rc != 0 and "ISMIP7_ROUTINE_FROZEN_MAX" in log
    assert submits(sandbox) == []


@pytest.mark.parametrize("state", ["COMPLETED", "FAILED", "OUT_OF_MEMORY"])
def test_a_strict_routine_parks_a_link_the_scheduler_did_not_end(sandbox, state):
    r"""Strict mode resubmits only a link the scheduler ended; a link that
    ended on its own without queueing a successor is parked with its state,
    even when its log shows no driver error or stall."""
    env, first, queue, cancels, logs = _watched_stage(
        sandbox, state, ISMIP7_ROUTINE_STRICT="1")
    (logs / f"ismip7_fwd_{first}.out").write_text("  t=2005.3 ...\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == []
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert blocked.startswith(f"job {first}: ended {state} without a successor"), blocked


@pytest.mark.parametrize("state", ["RUNNING", "COMPLETING", "PENDING", "REQUEUED"])
def test_a_strict_routine_waits_for_a_link_the_accounting_holds_live(sandbox, state):
    r"""squeue can drop a link before slurmdbd records its end, and a requeue
    passes through states that are not an end. Strict mode neither resubmits
    (a duplicate chain) nor parks such a stage; a later pass decides."""
    env, first, queue, cancels, logs = _watched_stage(
        sandbox, state, ISMIP7_ROUTINE_STRICT="1")
    (logs / f"ismip7_fwd_{first}.out").write_text("  t=2005.3 ...\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == []
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists()
    assert "has not ended in the accounting" in log


def test_watch_cancels_a_successor_that_can_never_start_and_resubmits(sandbox):
    r"""A successor behind a predecessor that ended in a way its dependency
    excludes sits PENDING on DependencyNeverSatisfied for ever and made the
    stage look alive for a day (seq1 control_mri_esm2, 19-20 October 2026).
    The watchdog cancels it and resubmits the stage."""
    env, first, queue, cancels, logs = _watched_stage(sandbox)
    queue.write_text("rt_out_hist_cesm_waccm 7777001 PENDING (null) 0:00 DependencyNeverSatisfied\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert cancels.read_text().split() == ["7777001"]
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]


def test_the_forwards_run_the_routine_s_calving_law(sandbox):
    r"""The forwards take the routine's calving law, its parameters and its
    advance rule, recorded at init; the defaults are the submission front
    (hfb at 150 kPa, retreat-only), the collapse mask stays with the
    projections, and the legacy pinned front stays off under a law."""
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_CALVING": "vonmises",
           "ISMIP7_ROUTINE_CALVING_PARAMS": "sigma_max_fl=0.275,sigma_max_gr=1.0",
           "ISMIP7_ROUTINE_FRONT_ADVANCE": "free"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    e = h["env"]
    assert e["ISMIP7_CALVING"] == "vonmises"
    assert e["ISMIP7_CALVING_PARAMS"] == "sigma_max_fl=0.275,sigma_max_gr=1.0"
    assert e["ISMIP7_FRONT_ADVANCE"] == "free"
    assert e["ISMIP7_FRACTURE"] == "none"
    assert "ISMIP7_FIXED_FRONT" not in e
    env_text = (sandbox / "out" / "routine.env").read_text()
    assert "ISMIP7_ROUTINE_CALVING=vonmises" in env_text


def test_a_law_without_parameters_runs_on_its_defaults(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_CALVING": "nocalving", "ISMIP7_ROUTINE_CALVING_PARAMS": "",
           "ISMIP7_ROUTINE_FRONT_ADVANCE": "free"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    assert h["env"]["ISMIP7_CALVING"] == "nocalving"
    assert "ISMIP7_CALVING_PARAMS" not in h["env"]


def test_a_routine_recorded_before_the_calving_knobs_runs_the_submission_front(sandbox):
    r"""The routines started on 13-14 October recorded no calving knobs; their
    forwards ran the submission front (ismip7_front_routine.sh) and must keep
    it, not be handed an empty law."""
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    envfile = sandbox / "out" / "routine.env"
    envfile.write_text("".join(
        line for line in envfile.read_text().splitlines(keepends=True)
        if not line.startswith(("export ISMIP7_ROUTINE_CALVING",
                                "export ISMIP7_ROUTINE_FRONT_ADVANCE"))))
    rc, log = routine(sandbox, "resume", "hist_cesm_waccm", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    e = h["env"]
    assert e["ISMIP7_CALVING"] == "hfb"
    assert e["ISMIP7_CALVING_PARAMS"] == "sigma_max=0.15"
    assert e["ISMIP7_FRONT_ADVANCE"] == "none"
    assert "ISMIP7_FIXED_FRONT" not in e


def test_an_unknown_advance_rule_is_refused_at_init(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={"ISMIP7_ROUTINE_FRONT_ADVANCE": "sometimes"})
    assert rc != 0 and "ISMIP7_ROUTINE_FRONT_ADVANCE" in log
    assert submits(sandbox) == []
