r"""The self-chaining forward runner's verdict and resubmit decisions.

``antarctica/scripts/batch_runners/projection.sbatch`` is the only thing
standing between a five-day projection and a chain that either loops forever
or reports a stall as a success. Its logic cannot be exercised on NOTS without
burning a queue slot, so it runs here against a Slurm shim: stub ``srun``,
``sbatch``, ``scontrol`` and ``module`` on PATH, and a stub driver that writes
a checkpoint with the ``t_yr`` and ``stalled`` attributes the real drivers
write and prints the same ``Saved:`` / ``Restart:`` lines. Everything else,
including every decision under test, is the shipped script.

The stub driver reads its knobs through ``icepack2_tools.runconfig``, exactly
as the real drivers do, so a value the runner mangles on the way through shows
up here as the wrong behaviour rather than as a string comparison.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SBATCH = REPO / "antarctica" / "scripts" / "batch_runners" / "projection.sbatch"
JOB_ID = "424242"

# The stub driver: a real driver's externally visible behaviour, and nothing
# else. It resolves its knobs through runconfig (so a bad value aborts the run
# the way experiment.py does), writes the final checkpoint, and reports it.
DRIVER = r'''
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from icepack2_tools.runconfig import (
    apparent_mb_mode, auto_resume, ismip7_output, smb_elevation_feedback,
)
from icepack2_tools.solverconfig import diagnostic_solver_mode

try:
    resume = auto_resume()
    amb = apparent_mb_mode()
    out = ismip7_output()
    feedback = smb_elevation_feedback()
except ValueError as exc:
    print(f"driver: {exc}")
    sys.exit(2)
print(f"driver: auto_resume={resume} apparent_mb={amb} ismip7_output={out}")
print(f"driver: smb_elevation_feedback={feedback}")
print(f"driver: solver={diagnostic_solver_mode()} dt={os.environ.get('ISMIP7_DT')}")
# srun takes SLURM_EXPORT_ENV as its own --export, so a list there would
# reach the real driver's environment through the launcher.
print(f"driver: srun export list={os.environ.get('SLURM_EXPORT_ENV', 'none')}")

start = os.environ.get("FAKE_START_YEAR", "")
if start:
    print(f"Restart: evolved geometry from checkpoint, t_yr={start}")

# run_simulation prints this for every driver, and the chain reads the run's
# end year back out of it rather than re-declaring the per-experiment default
t_end = float(os.environ.get("ISMIP7_T_END", "2301"))
print(f"Time-stepping: {start or 2015.0}->{t_end}, dt=0.1yr, 10 steps")

rc = int(os.environ.get("FAKE_RC", "0"))
if rc:
    sys.exit(rc)

import h5py

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(out, exist_ok=True)
path = os.path.join(out, "fake_experiment_2500_final.h5")
with h5py.File(path, "w") as h:
    h["/"].attrs["t_yr"] = float(os.environ["FAKE_T_YR"])
    h["/"].attrs["stalled"] = int(os.environ.get("FAKE_STALLED", "0"))
print(f"Saved: {path}")
'''

STUBS = {
    # Slurm's launcher: run the stub driver instead of the real one, with the
    # environment the runner handed it.
    "srun": '#!/bin/bash\nexec "$FAKE_PYTHON" "$FAKE_DRIVER"\n',
    # Record the resubmit so the test can see whether one happened and what it
    # carried. --export=ALL means the successor inherits this environment, so
    # dump the ISMIP7 part of it too.
    "sbatch": (
        "#!/bin/bash\n"
        'printf "ARGV: %s\\n" "$*" >> "$SBATCH_CALLS"\n'
        'env | grep -E "^(ISMIP7_|SLURM_GET_USER_ENV=|SLURM_EXPORT_ENV=)" | sort'
        ' | sed "s/^/ENV: /" >> "$SBATCH_CALLS"\n'
        'echo "Submitted batch job 999999"\n'
    ),
    "scontrol": (
        "#!/bin/bash\n"
        'case "$1 $2" in\n'
        '  "show job") echo "JobId=%s TimeLimit=1-00:00:00 Partition=commons'
        ' Features=cascadelake MinMemoryNode=240G" ;;\n'
        '  "show node") echo "CPUTot=40 RealMemory=187135 ActiveFeatures=cascadelake" ;;\n'
        "esac\n"
    ) % JOB_ID,
    "module": "#!/bin/bash\nexit 0\n",
}


@pytest.fixture
def sandbox(tmp_path):
    r"""A submit directory that looks like the repository to the runner."""
    for name in ("antarctica", "icepack2_tools"):
        (tmp_path / name).symlink_to(REPO / name)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in STUBS.items():
        stub = bin_dir / name
        stub.write_text(body)
        stub.chmod(0o755)
    # The runner calls `python`; give it this interpreter.
    (bin_dir / "python").write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (bin_dir / "python").chmod(0o755)
    (tmp_path / "driver.py").write_text(DRIVER)
    # ismip7_activate sources this and refuses to run if it is unreadable.
    (tmp_path / "activate").write_text("# stub venv\n")
    return tmp_path


def run_job(sandbox, **env):
    r"""Submit the runner once. Returns (exit code, job log, resubmits)."""
    log = sandbox / "logs" / f"ismip7_fwd_{JOB_ID}.out"
    log.parent.mkdir(exist_ok=True)
    calls = sandbox / "sbatch_calls.txt"
    base = {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "SCRATCH": str(sandbox),
        "SLURM_JOB_ID": JOB_ID,
        "SLURM_NTASKS": "12",
        "SLURM_JOB_NUM_NODES": "1",
        "SLURM_JOB_PARTITION": "commons",
        "SLURM_JOB_NAME": "ismip7_fwd",
        "SLURM_SUBMIT_DIR": str(sandbox),
        # sites/local.sh takes every setting from this environment, so the
        # runner logic is exercised without a scheduler or a site file.
        "ISMIP7_SITE": "local",
        # The sandbox reaches the real batch_runners through a symlink, so keep
        # this checkout's own sites/local.env out of the test.
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_FIREDRAKE": str(sandbox / "activate"),
        "ISMIP7_REPO": str(sandbox),
        "FAKE_PYTHON": sys.executable,
        "FAKE_DRIVER": str(sandbox / "driver.py"),
        "SBATCH_CALLS": str(calls),
        "FAKE_T_YR": "2050",
    }
    base.update({k: v for k, v in env.items() if v is not None})
    with open(log, "w") as fh:
        proc = subprocess.run(
            ["bash", str(SBATCH)], env=base, cwd=str(sandbox),
            stdout=fh, stderr=subprocess.STDOUT,
        )
    return proc.returncode, log.read_text(), (calls.read_text() if calls.exists() else "")


def test_short_run_chains_once(sandbox):
    r"""The case the chain exists for: stopped by the wall clock, short of
    t_end, having advanced. One successor, and it carries neither the wall
    budget nor a restart."""
    rc, log, calls = run_job(sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000")
    assert rc == 0, log
    assert "resubmitting" in log
    assert calls.count("ARGV:") == 1
    assert "--dependency=afterok:424242" in calls
    assert "ENV: ISMIP7_WALL_STOP_MIN" not in calls
    assert "ENV: ISMIP7_RESTART" not in calls


def test_the_successor_is_given_this_job_s_allocation(sandbox):
    r"""The job script carries no resource directives, so a successor
    submitted from inside a job has to be told the running job's own
    allocation. Without it the link lands on the cluster's defaults: wrong
    partition, wrong wall limit, and a default of one task, which the runner
    refuses outright - a 285-year projection would stop after its first link.
    --hint=nomultithread and -J are not readable back from scontrol, so they
    have to be restated. Without the first, Slurm packs the ranks onto
    hyperthreads and halves per-rank memory bandwidth for the rest of the
    chain. Without the second, sbatch falls back to the script filename and
    two concurrent chains become indistinguishable in squeue."""
    rc, log, calls = run_job(sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000")
    assert rc == 0, log
    argv = [line for line in calls.splitlines() if line.startswith("ARGV:")]
    assert len(argv) == 1, calls
    for flag in ("-p commons", "-C cascadelake", "--time=1-00:00:00",
                 "--mem=240G", "-N 1", "-n 12", "--cpus-per-task=1",
                 "--hint=nomultithread", "-J ismip7_fwd"):
        assert flag in argv[0], f"successor lost {flag}: {argv[0]}"


def test_a_list_form_first_link_leaves_nothing_behind(sandbox):
    r"""sbatch gives a job submitted with an --export list two variables that
    a bare --export=ALL would pass to every successor. SLURM_GET_USER_ENV=1
    has slurmd rebuild the login environment when the successor starts and
    hold the job when that fails. SLURM_EXPORT_ENV holds the list, which srun
    takes as its own --export, so the ISMIP7_RESTART the runner unsets before
    resubmitting would come back in every later link's driver and send it to
    the first link's starting checkpoint."""
    listed = "ALL,ISMIP7_SITE=local,ISMIP7_RESTART=/results/hist_final.h5"
    rc, log, calls = run_job(sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000",
                             ISMIP7_RESTART="/results/hist_final.h5",
                             SLURM_GET_USER_ENV="1", SLURM_EXPORT_ENV=listed)
    assert rc == 0, log
    assert "driver: srun export list=none" in log
    argv = [line for line in calls.splitlines() if line.startswith("ARGV:")]
    assert len(argv) == 1, calls
    assert [word for word in argv[0].split() if word.startswith("--export")] == ["--export=ALL"]
    assert "SLURM_GET_USER_ENV" not in calls and "SLURM_EXPORT_ENV" not in calls
    assert "ENV: ISMIP7_RESTART" not in calls


def test_the_successor_is_given_the_site_s_extra_flags(sandbox):
    r"""ISMIP7_SBATCH_EXTRA carries what a site insists on for every
    submission, a QOS for instance. scontrol does not hand it back in a form
    the chain reads, so a successor without it would be refused by exactly the
    partition the first link was accepted on."""
    rc, log, calls = run_job(sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000",
                             ISMIP7_SBATCH_EXTRA="--qos=long --mail-type=FAIL")
    assert rc == 0, log
    argv = [line for line in calls.splitlines() if line.startswith("ARGV:")]
    assert len(argv) == 1, calls
    assert "--qos=long --mail-type=FAIL" in argv[0], argv[0]


def test_reaching_t_end_finishes(sandbox):
    r"""At t_end there is nothing left to chain."""
    rc, log, calls = run_job(
        sandbox, FAKE_T_YR="2300", FAKE_START_YEAR="2250", ISMIP7_T_END="2300")
    assert rc == 0, log
    assert "finished" in log
    assert calls == ""


def test_stall_is_a_failure_even_with_the_chain_off(sandbox):
    r"""The verdict owns the exit code. A stalled run reports failure whether
    or not a successor was ever going to be submitted, so Slurm does not
    record it COMPLETED and mail it as a success."""
    for off in ({"ISMIP7_CHAIN": "0"}, {"ISMIP7_AUTO_RESUME": "0"}):
        rc, log, calls = run_job(
            sandbox, FAKE_T_YR="2037", FAKE_STALLED="1", FAKE_START_YEAR="2000",
            **off)
        assert rc == 1, log
        assert "the solver stalled" in log
        assert "would cold-start" not in log
        assert calls == "", f"a stalled run must not chain: {calls}"
        (sandbox / "sbatch_calls.txt").unlink(missing_ok=True)


def test_driver_exit_stops_the_chain(sandbox):
    r"""A non-zero driver exit is a real failure, never retried."""
    rc, log, calls = run_job(sandbox, FAKE_RC="7")
    assert rc == 7, log
    assert "driver exited 7" in log
    assert calls == ""


@pytest.mark.parametrize("value", ["0", ""])
def test_auto_resume_off_stops_the_chain(sandbox, value):
    r"""With resume off a successor would cold-start and repeat these years
    forever. The empty string is one of the off spellings, and
    ``submit.sh projection ISMIP7_AUTO_RESUME=`` delivers it (set empty in
    sbatch's environment), so it has to survive the runner's defaulting."""
    rc, log, calls = run_job(
        sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000",
        ISMIP7_AUTO_RESUME=value)
    assert rc == 0, log
    assert "driver: auto_resume=False" in log
    assert "auto-resume is off" in log
    assert calls == "", f"resume is off, so nothing may be resubmitted: {calls}"


def test_unset_auto_resume_defaults_on(sandbox):
    r"""Unset is the documented default: resume, and chain."""
    rc, log, calls = run_job(sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000")
    assert rc == 0, log
    assert "driver: auto_resume=True" in log
    assert calls.count("ARGV:") == 1


def test_unset_output_writes_the_submission(sandbox):
    r"""Every experiment this runner offers is a core experiment, and the
    yearly fields and scalars cannot be recovered afterwards: a projection
    that reaches 2300 without them has to be run again. So unset means on."""
    rc, log, _ = run_job(sandbox, FAKE_T_YR="2301")
    assert rc == 0, log
    assert "driver: auto_resume=True apparent_mb=balance ismip7_output=True" in log


def test_unset_smb_feedback_is_on_in_the_driver(sandbox):
    r"""The runner exports nothing for the feedback; the driver's own default
    (runconfig) turns it on."""
    rc, log, _ = run_job(sandbox, FAKE_T_YR="2301")
    assert rc == 0, log
    assert "driver: smb_elevation_feedback=True" in log


@pytest.mark.parametrize("value", ["0", ""])
def test_smb_feedback_off_reaches_the_driver(sandbox, value):
    r"""``submit.sh projection ISMIP7_SMB_ELEVATION_FEEDBACK=0`` (or empty) is
    how an A/B run turns the feedback off; the runner passes it through."""
    rc, log, _ = run_job(sandbox, FAKE_T_YR="2301",
                         ISMIP7_SMB_ELEVATION_FEEDBACK=value)
    assert rc == 0, log
    assert "driver: smb_elevation_feedback=False" in log


def test_a_bad_smb_feedback_value_aborts_the_run(sandbox):
    rc, log, calls = run_job(sandbox, ISMIP7_SMB_ELEVATION_FEEDBACK="yes")
    assert rc == 2, log
    assert "ISMIP7_SMB_ELEVATION_FEEDBACK must be 1 to enable" in log
    assert calls == ""


@pytest.mark.parametrize("value", ["0", ""])
def test_output_off_reaches_the_driver(sandbox, value):
    r"""Both documented off spellings survive the runner's defaulting, so a
    pipeline exercise can still say so; ``submit.sh projection ISMIP7_OUTPUT=``
    delivers the empty one."""
    rc, log, _ = run_job(sandbox, FAKE_T_YR="2301", ISMIP7_OUTPUT=value)
    assert rc == 0, log
    assert "ismip7_output=False" in log


def test_a_bad_output_value_aborts_the_run(sandbox):
    r"""runconfig owns the value set; the runner must not normalise a typo
    into a decision about whether a submission gets written."""
    rc, log, calls = run_job(sandbox, ISMIP7_OUTPUT="yes")
    assert rc == 2, log
    assert "ISMIP7_OUTPUT must be 1 to enable" in log
    assert calls == ""


def test_a_bad_auto_resume_value_aborts_the_run(sandbox):
    r"""The runner hands the value to the driver as given; runconfig owns the
    value set. Normalising in the shell would turn an ``=off`` meant to
    disable resume into a 1 and silently continue an old trajectory."""
    rc, log, calls = run_job(sandbox, ISMIP7_AUTO_RESUME="off")
    assert rc == 2, log
    assert "ISMIP7_AUTO_RESUME must be an integer flag" in log
    assert calls == ""


@pytest.mark.parametrize(
    "value,expected",
    [("1", "balance"), ("balance", "balance"), ("div", "div"),
     ("0", "None"), ("off", "None"), ("none", "None"), ("", "None")],
)
def test_apparent_mb_reaches_the_driver_as_given(sandbox, value, expected):
    r"""Including the empty string, which the runner must not overwrite with
    its default: a run asked to drop the correction would otherwise get the
    full balanced one, with the t=0 tendency forced to zero."""
    _, log, _ = run_job(
        sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000",
        ISMIP7_APPARENT_MB=value)
    assert f"apparent_mb={expected}" in log


def test_no_progress_stops_the_chain(sandbox):
    r"""A clean exit that advanced no years is not the wall-clock case: the
    setup alone spent the budget, and a successor under the same budget would
    do the same. Report it and end the chain, with no retry state carried
    between jobs."""
    rc, log, calls = run_job(
        sandbox, FAKE_T_YR="2000", FAKE_START_YEAR="2000")
    assert rc == 1, log
    assert "no progress" in log
    assert calls == ""


def test_no_reported_checkpoint_does_not_chain(sandbox):
    r"""No verdict, no resubmit: a job whose driver never reported a final
    checkpoint has no year to compare against t_end, and chaining blind would
    repeat the years it cannot see."""
    (sandbox / "driver.py").write_text(
        DRIVER.replace('print(f"Saved: {path}")', 'print("driver: no checkpoint")'))
    rc, log, calls = run_job(sandbox, FAKE_T_YR="2050")
    assert rc == 0, log
    assert "could not read this job's final checkpoint" in log
    assert calls == ""


def test_the_runner_names_the_production_solver_and_step(sandbox):
    r"""The forward runner is where production leaves the full-Jacobian
    reference: solverconfig's own default has to stay full_mumps, because the
    inversion reads it too. The successor inherits the choice with the rest of
    the environment, so a chain cannot change solver between links."""
    rc, log, calls = run_job(sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000")
    assert rc == 0, log
    assert "driver: solver=scpc_gamg dt=0.05" in log
    assert "    solver=scpc_gamg dt=0.05 " in log
    assert "ENV: ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_gamg" in calls
    assert "ENV: ISMIP7_DT=0.05" in calls


def test_a_named_solver_and_step_win(sandbox):
    r"""The reference stays one variable away."""
    rc, log, _ = run_job(
        sandbox, FAKE_T_YR="2050", FAKE_START_YEAR="2000",
        ISMIP7_DIAGNOSTIC_LINEAR_SOLVER="full_mumps", ISMIP7_DT="0.1")
    assert rc == 0, log
    assert "driver: solver=full_mumps dt=0.1" in log
