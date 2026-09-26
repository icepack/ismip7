r"""The self-chaining inversion runner's done-marker and warm-start decisions,
and the objective a resumed link minimises.

``antarctica/scripts/batch_runners/inversion.sbatch`` decides, after every
link of a multi-day chain, whether the MAP on disk is finished or whether the
successor should warm-start and spend another budget on it. Getting that wrong
costs a queue slot either way: job 1339342 re-inverted the MAP job 1339328 had
already finished, and a marker written on a MAP whose write never completed
would leave a truncated checkpoint with no job to rebuild it.

Like ``test_projection_chain.py``, this runs the shipped script against a Slurm
shim: stub ``srun``, ``sbatch``, ``scontrol`` and ``module`` on PATH, and a stub
driver with the real driver's externally visible behaviour - it writes the
checkpoint at ``ISMIP7_MAP_OUT``, prints the same ``Saved MAP:`` line only after
that write returns, and touches ``<map>.done`` itself. Everything else,
including every decision under test, is the shipped script.

A resumed link also has to minimise the objective the link before it did
(issue 68). The weight driver below resolves its log-velocity weight with the
inversion's own ``resolve_log_vel_weight`` from the attributes of an HDF5
checkpoint, and the last tests hold that function to its cases.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from icepack2_tools.optimization import recorded_objective, resolve_log_vel_weight

REPO = Path(__file__).resolve().parent.parent
SBATCH = REPO / "antarctica" / "scripts" / "batch_runners" / "inversion.sbatch"
JOB_ID = "424243"

# The stub driver, in the real driver's order: the checkpoint write first, then
# the "Saved MAP:" line, then the marker. FAKE_DIE_AFTER names the step a wall
# clock kill lands on, which is the whole point of the ordering.
DRIVER = r'''
import os
import sys

map_out = os.environ["ISMIP7_MAP_OUT"]
print(f"driver: warm_start={os.environ.get('ISMIP7_WARM_START', 'none')}")
print(f"driver: maxiter={os.environ['ISMIP7_MAXITER']}")

die = os.environ.get("FAKE_DIE_AFTER", "")
print("Optimization finished: CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH")
if die == "optimizer":
    sys.exit(137)

with open(map_out, "w") as fh:                   # save_map: mode "w", truncating
    fh.write("checkpoint\n")
if die == "map_write":
    sys.exit(137)

print(f"Saved MAP: {map_out} (misfit_norm=sigma)")
if die == "marker":
    sys.exit(137)

open(map_out + ".done", "w").close()
'''

STUBS = {
    "srun": '#!/bin/bash\nexec "$FAKE_PYTHON" "$FAKE_DRIVER"\n',
    "sbatch": (
        "#!/bin/bash\n"
        'printf "ARGV: %s\\n" "$*" >> "$SBATCH_CALLS"\n'
        'env | grep -E "^(ISMIP7_|SLURM_GET_USER_ENV=|SLURM_EXPORT_ENV=)" | sort'
        ' | sed "s/^/ENV: /" >> "$SBATCH_CALLS"\n'
        'echo "999999"\n'
    ),
    "scontrol": (
        "#!/bin/bash\n"
        'case "$1 $2" in\n'
        '  "show job") echo "JobId=%s TimeLimit=3-00:00:00 Features=(null)'
        ' MinMemoryNode=240G" ;;\n'
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
    (bin_dir / "python").write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (bin_dir / "python").chmod(0o755)
    (tmp_path / "driver.py").write_text(DRIVER)
    (tmp_path / "activate").write_text("# stub venv\n")
    return tmp_path


def map_out(sandbox):
    return sandbox / "inversion_map.h5"


def run_job(sandbox, job_id=JOB_ID, **env):
    r"""Submit the runner once. Returns (exit code, job log, resubmits)."""
    log = sandbox / "logs" / f"ismip7_inv_{job_id}.out"
    log.parent.mkdir(exist_ok=True)
    calls = sandbox / "sbatch_calls.txt"
    base = {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "SCRATCH": str(sandbox),
        "SLURM_JOB_ID": job_id,
        "SLURM_NTASKS": "32",
        "SLURM_JOB_NUM_NODES": "1",
        "SLURM_JOB_PARTITION": "long",
        "SLURM_SUBMIT_DIR": str(sandbox),
        # sites/local.sh takes every setting from this environment, so the
        # runner logic is exercised without a scheduler or a site file.
        "ISMIP7_SITE": "local",
        # The sandbox reaches the real batch_runners through a symlink, so keep
        # this checkout's own sites/local.env out of the test.
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_FIREDRAKE": str(sandbox / "activate"),
        "ISMIP7_REPO": str(sandbox),
        "ISMIP7_MAP_OUT": str(map_out(sandbox)),
        "FAKE_PYTHON": sys.executable,
        "FAKE_DRIVER": str(sandbox / "driver.py"),
        "SBATCH_CALLS": str(calls),
    }
    base.update({k: v for k, v in env.items() if v is not None})
    with open(log, "w") as fh:
        proc = subprocess.run(
            ["bash", str(SBATCH)], env=base, cwd=str(sandbox),
            stdout=fh, stderr=subprocess.STDOUT,
        )
    return proc.returncode, log.read_text(), (calls.read_text() if calls.exists() else "")


def test_a_finished_inversion_is_marked_done_and_its_successor_exits(sandbox):
    r"""The 1339328/1339342 case. A run that saved its MAP leaves the marker,
    and the successor Slurm already queued must find nothing to do: no driver
    launched, no further link queued."""
    rc, log, calls = run_job(sandbox)
    assert rc == 0, log
    assert map_out(sandbox).exists()
    assert (sandbox / "inversion_map.h5.done").exists()
    assert calls.count("ARGV:") == 1
    assert "--dependency=afterany:424243" in calls

    rc, log, calls = run_job(sandbox, job_id="424244")
    assert rc == 0, log
    assert "already finished; nothing to do" in log
    assert "driver:" not in log, "the successor must not launch the inversion"
    assert calls.count("ARGV:") == 1, f"no third link may be queued: {calls}"


def test_a_kill_during_the_map_write_leaves_no_marker(sandbox):
    r"""The regression the marker rule was re-keyed for. The optimizer returned
    and said so, then the wall clock killed the job inside save_map, which
    truncates the periodic checkpoint. Nothing may be marked done on that: the
    log's "Optimization finished" is not evidence that a MAP reached disk."""
    rc, log, _ = run_job(sandbox, FAKE_DIE_AFTER="map_write")
    assert rc == 137, log
    assert "Optimization finished" in log
    assert not (sandbox / "inversion_map.h5.done").exists()


@pytest.mark.parametrize("die", ["optimizer", "map_write"])
def test_the_successor_warm_starts_from_an_unfinished_map(sandbox, die):
    r"""No marker means the budget is unspent, so the successor hands the
    checkpoint to the driver rather than cold-starting it."""
    map_out(sandbox).write_text("checkpoint\n")
    rc, _, _ = run_job(sandbox, FAKE_DIE_AFTER=die)
    assert rc == 137
    assert not (sandbox / "inversion_map.h5.done").exists()

    rc, log, calls = run_job(sandbox, job_id="424244")
    assert rc == 0, log
    assert f"driver: warm_start={map_out(sandbox)}" in log
    assert calls.count("ARGV:") == 2, "each unfinished link queues a successor"


def test_the_fallback_marks_a_saved_map_the_driver_could_not(sandbox):
    r"""The runner's rule is the fallback for a driver that saved the MAP but
    never got to write the marker itself."""
    rc, log, _ = run_job(sandbox, FAKE_DIE_AFTER="marker")
    assert rc == 137, log
    assert map_out(sandbox).exists()
    assert (sandbox / "inversion_map.h5.done").exists()

    rc, log, _ = run_job(sandbox, job_id="424244")
    assert rc == 0, log
    assert "already finished; nothing to do" in log


def test_an_existing_marker_queues_nothing(sandbox):
    r"""The marker is checked before the successor is queued, so a chain that
    is already finished cannot extend itself."""
    (sandbox / "inversion_map.h5.done").touch()
    rc, log, calls = run_job(sandbox)
    assert rc == 0, log
    assert "already finished; nothing to do" in log
    assert calls == "", f"nothing may be queued once the marker exists: {calls}"


def test_the_successor_carries_its_depth_in_sbatch_s_environment(sandbox):
    r"""The list form --export=ALL,ISMIP7_CHAIN_DEPTH=N has Slurm rebuild the
    login environment when the successor starts and hold the job when that
    fails, as IU Quartz held 10555984. The depth goes in sbatch's environment
    under a bare ALL, and what a list-form first link left in the job
    (SLURM_GET_USER_ENV=1 and the list itself in SLURM_EXPORT_ENV) stays out
    of the successor."""
    rc, log, calls = run_job(
        sandbox, FAKE_DIE_AFTER="optimizer", ISMIP7_CHAIN_DEPTH="1",
        SLURM_GET_USER_ENV="1",
        SLURM_EXPORT_ENV="ALL,ISMIP7_SITE=local,ISMIP7_MAXITER=200")
    assert rc == 137, log
    argv = [line for line in calls.splitlines() if line.startswith("ARGV:")]
    assert len(argv) == 1, calls
    assert [word for word in argv[0].split() if word.startswith("--export")] == ["--export=ALL"]
    assert "ENV: ISMIP7_CHAIN_DEPTH=2" in calls.splitlines()
    assert "SLURM_GET_USER_ENV" not in calls and "SLURM_EXPORT_ENV" not in calls
    assert "chain: successor 999999 queued (depth 2 of 4)" in log


def _successor_argv(calls):
    argv = [line for line in calls.splitlines() if line.startswith("ARGV:")]
    assert len(argv) == 1, calls
    return argv[0].split()


def test_a_multinode_link_hands_its_layout_to_the_successor(sandbox):
    r"""A link that asked for 2 nodes x 32 tasks queues its successor with the
    same layout. Without it the successor's srun, reading the link's
    SLURM_NTASKS_PER_NODE from the exported environment, asked an allocation
    laid out by Slurm for more than it held (NOTS 1614035)."""
    rc, log, calls = run_job(
        sandbox, FAKE_DIE_AFTER="optimizer",
        SLURM_JOB_NUM_NODES="2", SLURM_NTASKS="64", SLURM_NTASKS_PER_NODE="32")
    assert rc == 137, log
    words = _successor_argv(calls)
    assert words[words.index("-N") + 1] == "2"
    assert words[words.index("-n") + 1] == "64"
    assert "--ntasks-per-node=32" in words


def test_a_link_without_a_layout_request_leaves_it_to_slurm(sandbox):
    r"""A link that named no per-node layout does not invent one for its
    successor."""
    rc, log, calls = run_job(sandbox, FAKE_DIE_AFTER="optimizer")
    assert rc == 137, log
    words = _successor_argv(calls)
    assert not [word for word in words if word.startswith("--ntasks-per-node")]


def test_the_chain_depth_cap_stops_it(sandbox):
    r"""The cap is the only bound on a chain whose links keep failing early."""
    rc, log, calls = run_job(
        sandbox, FAKE_DIE_AFTER="optimizer",
        ISMIP7_CHAIN_DEPTH="4", ISMIP7_CHAIN_MAX="4")
    assert rc == 137, log
    assert "reached ISMIP7_CHAIN_MAX=4, no successor" in log
    assert calls == ""


# The real driver's weight, in the real driver's order: read the warm start's
# attributes, resolve the weight, write it into the checkpoint as save_map
# does. FAKE_DERIVED stands in for the chi^2 / log ratio at the state a link
# starts from, which is all the solve contributes to the decision.
WEIGHT_DRIVER = r'''
import os
import sys

import h5py

from icepack2_tools.optimization import recorded_objective, resolve_log_vel_weight


class Checkpoint:
    """CheckpointFile's has_attr/get_attr, which read h5py's root attributes."""

    def __init__(self, handle):
        self.handle = handle

    def has_attr(self, path, key):
        return key in self.handle[path].attrs

    def get_attr(self, path, key):
        return self.handle[path].attrs[key]


map_out = os.environ["ISMIP7_MAP_OUT"]
warm = os.environ.get("ISMIP7_WARM_START", "")
recorded = None
if warm:
    with h5py.File(warm, "r") as handle:
        recorded = recorded_objective(Checkpoint(handle))
misfit_norm = os.environ.get("ISMIP7_MISFIT_NORM", "sigma").lower()
eps = float(os.environ.get("ISMIP7_LOG_VEL_EPS", "1.0"))
weight, source, note = resolve_log_vel_weight(
    os.environ.get("ISMIP7_LOG_VEL_WEIGHT", "0"), float(os.environ["FAKE_DERIVED"]),
    recorded, misfit_norm=misfit_norm, eps=eps)
print(f"driver: log_vel_weight={weight!r} source={source}")

with h5py.File(map_out, "w") as handle:
    handle["/"].attrs["misfit_norm"] = misfit_norm
    handle["/"].attrs["log_vel_weight"] = float(weight)
    handle["/"].attrs["log_vel_weight_source"] = source
    handle["/"].attrs["log_vel_eps"] = eps
if os.environ.get("FAKE_DIE_AFTER") == "map_write":
    sys.exit(137)

print(f"Saved MAP: {map_out} (misfit_norm={misfit_norm} log_vel_weight={weight:g})")
open(map_out + ".done", "w").close()
'''


def test_two_links_minimise_one_log_velocity_weight(sandbox):
    r"""Issue 68, in the numbers IU's 32 km chain measured. Under the runner's
    default ISMIP7_LOG_VEL_WEIGHT=auto the first link derived 17452 at its
    cold start, and the successor, warm-starting from the checkpoint of
    evaluation 140, re-derived 2495 there: two links of one MAP minimising
    different objectives. The successor holds the weight the checkpoint
    records, and the MAP it finishes carries that weight."""
    h5py = pytest.importorskip("h5py")
    driver = sandbox / "weight_driver.py"
    driver.write_text(WEIGHT_DRIVER)

    rc, log, _ = run_job(sandbox, FAKE_DRIVER=str(driver),
                         FAKE_DERIVED="17451.9897644", FAKE_DIE_AFTER="map_write")
    assert rc == 137, log
    assert "driver: log_vel_weight=17451.9897644 source=derived" in log
    assert not (sandbox / "inversion_map.h5.done").exists()

    rc, log, _ = run_job(sandbox, job_id="424244", FAKE_DRIVER=str(driver),
                         FAKE_DERIVED="2495.0")
    assert rc == 0, log
    assert f"warm-starting theta/phi from the periodic checkpoint {map_out(sandbox)}" in log
    assert "driver: log_vel_weight=17451.9897644 source=warm_start" in log
    with h5py.File(map_out(sandbox), "r") as handle:
        assert float(handle["/"].attrs["log_vel_weight"]) == 17451.9897644
        assert handle["/"].attrs["log_vel_weight_source"] == "warm_start"


SIGMA = {"misfit_norm": "sigma", "log_vel_eps": 1.0}


def resolve(requested, derived, recorded):
    return resolve_log_vel_weight(requested, derived, recorded,
                                  misfit_norm="sigma", eps=1.0)


@pytest.mark.parametrize("requested", ["auto", "AUTO", " auto "])
def test_auto_holds_the_weight_the_warm_start_records(requested):
    recorded = {**SIGMA, "log_vel_weight": 17451.9897644}
    assert resolve(requested, 2495.0, recorded) == (17451.9897644, "warm_start", "")


def test_auto_derives_the_weight_at_a_cold_start():
    assert resolve("auto", 17451.9897644, None) == (17451.9897644, "derived", "")


@pytest.mark.parametrize("recorded, why", [
    # A MAP written before the log term existed, or an adapted-mesh checkpoint.
    ({"misfit_norm": "sigma"}, "it records no log-velocity weight"),
    ({}, "it records no log-velocity weight"),
    # Inverted without the log term: holding its 0 would drop the term auto asks for.
    ({**SIGMA, "log_vel_weight": 0.0}, "it records weight 0"),
    # The weight scales a chi^2 in other units, or a log term with another floor.
    ({**SIGMA, "log_vel_weight": 1.0e3, "misfit_norm": "none"},
     "misfit_norm none there, sigma here"),
    ({**SIGMA, "log_vel_weight": 1.0e3, "log_vel_eps": 10.0},
     "log_vel_eps 10.0 there, 1 here"),
])
def test_auto_derives_when_the_warm_start_holds_no_comparable_weight(recorded, why):
    weight, source, note = resolve("auto", 2495.0, recorded)
    assert (weight, source) == (2495.0, "derived")
    assert note == f"the warm start's weight is not used: {why}"


def test_a_number_is_used_as_given():
    recorded = {**SIGMA, "log_vel_weight": 17451.9897644}
    # IU resubmitted its second link with the recorded weight by hand.
    assert resolve("17451.9897644", None, recorded) == (17451.9897644, "requested", "")
    assert resolve("0", None, None) == (0.0, "requested", "")
    weight, source, note = resolve("2495", None, recorded)
    assert (weight, source) == (2495.0, "requested")
    assert note == "it replaces the weight 17452 the warm start was minimised under"


def test_a_checkpoint_file_returns_what_save_map_wrote(tmp_path):
    r"""The attributes come back through firedrake.CheckpointFile itself, as a
    numpy float and a str, and hold an auto weight."""
    fd = pytest.importorskip("firedrake")
    path = str(tmp_path / "map.h5")
    with fd.CheckpointFile(path, "w") as chk:
        chk.save_mesh(fd.UnitSquareMesh(2, 2))
        chk.set_attr("/", "misfit_norm", "sigma")
        chk.set_attr("/", "log_vel_weight", 17451.9897644)
        chk.set_attr("/", "log_vel_eps", 1.0)
    with fd.CheckpointFile(path, "r") as chk:
        recorded = recorded_objective(chk)
    assert recorded == {**SIGMA, "log_vel_weight": 17451.9897644}
    assert resolve("auto", 2495.0, recorded) == (17451.9897644, "warm_start", "")
