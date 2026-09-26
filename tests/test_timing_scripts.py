r"""The timing campaign's job scripts, run against a Slurm shim.

They used to carry one cluster's account, modules and venv. They now source
`site_core.sh` and nothing else, so what is pinned here is what that must and
must not do to a lane:

- the site's Firedrake is activated and ranks start through `ismip7_mpirun`;
- no model default leaks in. `site_env.sh` exports ISMIP7_MESH, which
  `run_timing.py` refuses on a matrix lane, and defaults ISMIP7_FRICTION to a
  law the campaign does not run. A lane's model settings are its exports alone;
- the kernel cache persists between jobs, because a lane's seconds_per_step has
  no warm-up excluded;
- the status file still walks running -> finished/failed as the manager expects.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BR = "scripts/batch_runners"

DRIVER = r'''
import os, sys
with open(os.environ["FAKE_ENV_DUMP"], "a") as fh:
    fh.write("ARGV " + " ".join(sys.argv[1:]) + "\n")
    for key in sorted(os.environ):
        if key.startswith(("ISMIP7_", "PYOP2_", "FIREDRAKE_", "OMP_", "OPENBLAS_",
                           "XDG_", "MPLCONFIGDIR", "FAKE_VENV", "LOOPY_")):
            fh.write(f"{key}={os.environ[key]}\n")
status = os.environ.get("FAKE_WRITE_STATUS")
if status:
    open(os.environ["ISMIP7_TIMING_STATUS"], "w").write(status + "\n")
sys.exit(int(os.environ.get("FAKE_RC", "0")))
'''

STUBS = {
    # `srun -n N python -u script args...`: record N, hand the rest to the driver.
    "srun": ('#!/bin/bash\n'
             'echo "SRUN $*" >> "$FAKE_ENV_DUMP"\n'
             'shift 2; shift 2\n'
             'exec "$FAKE_PYTHON" "$FAKE_DRIVER" "$@"\n'),
    "scontrol": "#!/bin/bash\nexit 0\n",
    "module": "#!/bin/bash\nexit 0\n",
}


@pytest.fixture
def sandbox(tmp_path):
    r"""A submit directory that looks like antarctica/ to the job script."""
    root = tmp_path / "repo" / "antarctica"
    root.mkdir(parents=True)
    (root / "scripts").symlink_to(REPO / "antarctica" / "scripts")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in STUBS.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    (tmp_path / "driver.py").write_text(DRIVER)
    (tmp_path / "activate").write_text("export FAKE_VENV_ACTIVE=1\n")
    return tmp_path


def run_script(sandbox, name, **env):
    root = sandbox / "repo" / "antarctica"
    base = {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "SLURM_JOB_ID": "777",
        "SLURM_NTASKS": "16",
        "SLURM_SUBMIT_DIR": str(root),
        "ISMIP7_SITE": "local",
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_FIREDRAKE": str(sandbox / "activate"),
        "ISMIP7_REPO": str(sandbox / "repo"),
        "FAKE_PYTHON": sys.executable,
        "FAKE_DRIVER": str(sandbox / "driver.py"),
        "FAKE_ENV_DUMP": str(sandbox / "env_dump.txt"),
    }
    base.update(env)
    proc = subprocess.run(["bash", str(root / BR / name)], env=base, cwd=str(sandbox),
                          capture_output=True, text=True)
    dump = sandbox / "env_dump.txt"
    return proc, (dump.read_text().splitlines() if dump.exists() else [])


def test_a_lane_runs_in_the_site_s_environment_with_only_its_own_model_settings(sandbox):
    status = sandbox / "status.txt"
    proc, seen = run_script(sandbox, "timing_transient.script",
                            ISMIP7_TIMING_STATUS=str(status), ISMIP7_LC="500")
    assert proc.returncode == 0, proc.stderr
    assert "SRUN -n 16 python -u scripts/run_timing.py" in seen
    assert "FAKE_VENV_ACTIVE=1" in seen
    assert "OMP_NUM_THREADS=1" in seen and "OPENBLAS_NUM_THREADS=1" in seen
    assert "ISMIP7_LC=500" in seen
    leaked = [line for line in seen if line.startswith((
        "ISMIP7_MESH=", "ISMIP7_FRICTION=", "ISMIP7_LC_COARSE=", "ISMIP7_MAP_DEFAULT=",
        "ISMIP7_GEOMETRY_SPACE=", "ISMIP7_N_FLOW=", "ISMIP7_DATA_ROOT="))]
    assert leaked == []
    assert "site    local" in proc.stdout and "friction=" not in proc.stdout
    assert status.read_text().startswith("finished exit_code=0 job_id=777 ")


def test_the_kernel_cache_persists_between_lanes(sandbox):
    r"""By default Firedrake's own location (no PYOP2_CACHE_DIR, and no empty
    per-job directory left behind); a named one when the site gives it."""
    proc, seen = run_script(sandbox, "timing_transient.script")
    assert proc.returncode == 0, proc.stderr
    assert not any(line.startswith("PYOP2_CACHE_DIR=") for line in seen)
    assert not (sandbox / ".pyop2_cache" / "777").exists()
    jit = sandbox / "jit"
    proc, seen = run_script(sandbox, "timing_transient.script",
                            ISMIP7_TIMING_JIT_CACHE=str(jit))
    assert proc.returncode == 0, proc.stderr
    assert f"PYOP2_CACHE_DIR={jit}/pyop2" in seen
    assert f"FIREDRAKE_TSFC_KERNEL_CACHE_DIR={jit}/tsfc" in seen
    assert (jit / "pyop2").is_dir()


def test_loopys_dict_stays_per_lane_while_the_kernel_cache_is_shared(sandbox):
    r"""Rice 1559476 died 40 s in with a pytools NoSuchEntryError, racing a job
    submitted a second earlier. `make timing-scout` submits one lane per mesh at
    once, so a shared loopy dict puts the timing lanes in exactly that state:
    the warm tree the lanes were measured against comes back for PyOP2 and TSFC
    alone, and loopy keeps the per-job directory ismip7_activate gave it, with
    its dict switched off (Quartz 10569250 and 10569252 hung on its lock)."""
    per_job = sandbox / ".pyop2_cache" / "xdg" / "777"
    jit = sandbox / "jit"
    for env in ({}, {"ISMIP7_TIMING_JIT_CACHE": str(jit)}):
        (sandbox / "env_dump.txt").unlink(missing_ok=True)
        proc, seen = run_script(sandbox, "timing_transient.script", **env)
        assert proc.returncode == 0, proc.stderr
        assert f"XDG_CACHE_HOME={per_job}" in seen
        assert "LOOPY_NO_CACHE=1" in seen
        assert per_job.is_dir()
    assert not (jit / "xdg").exists()
    assert not (sandbox / ".pyop2_cache" / "777").exists()


def test_a_cache_the_site_s_own_environment_names_is_the_one_kept(sandbox):
    r"""Found on Quartz: IU's firedrake modulefile sets PYOP2_CACHE_DIR and the
    TSFC cache to a scratch directory, and that is the warm cache every lane
    there was measured with. ismip7_activate replaces it with a per-job one, so
    simply unsetting that afterwards sent the first ported scout to an empty
    default cache: 17 s of compilation inside step 1 and 37.2 s/step against a
    32.8 reference. What the modules or venv set has to come back."""
    site_cache = sandbox / "scratch" / "firedrake.cache"
    (sandbox / "activate").write_text(
        "export FAKE_VENV_ACTIVE=1\n"
        f"export PYOP2_CACHE_DIR={site_cache}\n"
        f"export FIREDRAKE_TSFC_KERNEL_CACHE_DIR={site_cache}\n")
    proc, seen = run_script(sandbox, "timing_transient.script")
    assert proc.returncode == 0, proc.stderr
    assert f"PYOP2_CACHE_DIR={site_cache}" in seen
    assert f"FIREDRAKE_TSFC_KERNEL_CACHE_DIR={site_cache}" in seen
    # The site's XDG_CACHE_HOME is not restored with them: see
    # test_loopys_dict_stays_per_lane_while_the_kernel_cache_is_shared.
    assert f"XDG_CACHE_HOME={sandbox / '.pyop2_cache' / 'xdg' / '777'}" in seen
    assert not (sandbox / ".pyop2_cache" / "777").exists()
    # A cache named for the timing lanes still wins over the site's.
    jit = sandbox / "jit"
    proc, seen = run_script(sandbox, "timing_transient.script",
                            ISMIP7_TIMING_JIT_CACHE=str(jit))
    assert f"PYOP2_CACHE_DIR={jit}/pyop2" in seen[len(seen) // 2:]


def test_a_killed_lane_is_stamped_failed_and_a_specific_stamp_is_kept(sandbox):
    status = sandbox / "status.txt"
    proc, _ = run_script(sandbox, "timing_transient.script",
                         ISMIP7_TIMING_STATUS=str(status), FAKE_RC="137")
    assert proc.returncode == 137
    assert status.read_text().startswith(
        "failed category=external_termination exit_code=137 job_id=777 ")
    proc, _ = run_script(sandbox, "timing_transient.script",
                         ISMIP7_TIMING_STATUS=str(status), FAKE_RC="1",
                         FAKE_WRITE_STATUS="failed category=tripwire phase=step_3")
    assert proc.returncode == 1
    assert status.read_text() == "failed category=tripwire phase=step_3\n"


def test_prepare_repacks_on_one_rank_through_the_same_launcher(sandbox):
    cache = sandbox / "cache"
    cache.mkdir()
    status = cache / "status.txt"
    final, manifest = cache / "state.h5", cache / "state.json"
    final.write_text("x"); manifest.write_text("{}")
    proc, seen = run_script(
        sandbox, "timing_prepare.script",
        ISMIP7_TIMING_CACHE_RAW=str(cache / "raw.h5"), ISMIP7_TIMING_CACHE=str(final),
        ISMIP7_TIMING_CACHE_MANIFEST=str(manifest), ISMIP7_TIMING_CACHE_STATUS=str(status),
        ISMIP7_TIMING_CACHE_PRISTINE=str(cache / "state.prepare.h5"))
    assert proc.returncode == 0, proc.stderr
    launches = [line for line in seen if line.startswith("SRUN ")]
    assert launches[0] == "SRUN -n 16 python -u scripts/prepare_timing_cache.py"
    assert launches[1].startswith("SRUN -n 1 python -u scripts/redistribute_checkpoint.py --input ")
    assert status.read_text().startswith("finished phase=published exit_code=0 ")
    assert (cache / "state.prepare.h5").exists()


@pytest.mark.parametrize("name", sorted(
    p.name for p in (REPO / "antarctica" / BR).glob("*.script")))
def test_no_job_script_names_a_cluster_a_person_or_a_resource(name):
    text = (REPO / "antarctica" / BR / name).read_text()
    for needle in ("/N/u/", "r00905", "@iu.edu", "module load", "srun "):
        assert needle not in text, f"{name} still carries {needle!r}"
    directives = [line for line in text.splitlines() if line.startswith("#SBATCH")]
    assert len(directives) == 2
    assert all(line.split()[1] in ("-o", "-e") for line in directives), directives
    assert ". scripts/batch_runners/site_core.sh" in text
    assert ". scripts/batch_runners/site_env.sh" not in text


def test_the_score_script_scores_a_map_and_then_its_transferred_state(sandbox):
    r"""map_check_score.script: score_map.py --json on the MAP's own mesh,
    --restart for a prepared cache, the census only when asked for."""
    root = sandbox / "repo" / "antarctica"
    map_path = sandbox / "map.h5"
    map_path.write_text("map\n")
    status = sandbox / "status.txt"
    out = sandbox / "score.json"
    proc, seen = run_script(sandbox, "map_check_score.script",
                            ISMIP7_MAP=str(map_path), ISMIP7_MAP_CHECK_SCORE_JSON=str(out),
                            ISMIP7_MAP_CHECK_STATUS=str(status), ISMIP7_MESH="checkpoint",
                            ISMIP7_FRICTION="regularized_coulomb")
    assert proc.returncode == 0, proc.stderr
    assert f"SRUN -n 16 python -u scripts/score_map.py --json {out} {map_path}" in seen
    assert not any("check_budd_map.py" in line for line in seen)
    assert "ISMIP7_MESH=checkpoint" in seen and "ISMIP7_FRICTION=regularized_coulomb" in seen
    assert status.read_text().startswith("finished phase=score exit_code=0 job_id=777 ")
    (sandbox / "env_dump.txt").unlink()
    census = sandbox / "census.txt"
    proc, seen = run_script(sandbox, "map_check_score.script",
                            ISMIP7_MAP=str(map_path), ISMIP7_MAP_CHECK_SCORE_JSON=str(out),
                            ISMIP7_MAP_CHECK_STATUS=str(status),
                            ISMIP7_MAP_CHECK_RESTART=str(sandbox / "cache.h5"),
                            ISMIP7_MAP_CHECK_CENSUS=str(census), ISMIP7_FRICTION="budd")
    assert proc.returncode == 0, proc.stderr
    assert (f"SRUN -n 16 python -u scripts/score_map.py --json {out} "
            f"--restart {sandbox / 'cache.h5'} {map_path}") in seen
    assert f"SRUN -n 16 python -u scripts/check_budd_map.py {map_path}" in seen
    assert "ISMIP7_CHECK_FRICTION=budd" in seen
    assert census.is_file()
    assert status.read_text().startswith("finished phase=score exit_code=0 ")
    proc, _ = run_script(sandbox, "map_check_score.script",
                         ISMIP7_MAP=str(map_path), ISMIP7_MAP_CHECK_SCORE_JSON=str(out),
                         ISMIP7_MAP_CHECK_STATUS=str(status), FAKE_RC="3")
    assert proc.returncode == 3
    assert status.read_text().startswith("failed phase=score exit_code=3 ")


def test_the_audit_script_collects_every_part_whatever_each_returns(sandbox):
    r"""map_check_audit.script keeps each tool's exit code and output in one
    JSON; a tool that fails on its input is recorded, never fatal."""
    for name in ("native.csv", "transfer.csv", "native_final.h5", "transfer_final.h5"):
        (sandbox / name).write_text("year,mass_gt\n2015.0,1\n")
    out = sandbox / "audit_controls.json"
    status = sandbox / "status.txt"
    env = {
        "ISMIP7_MAP_CHECK_CSV_NATIVE": str(sandbox / "native.csv"),
        "ISMIP7_MAP_CHECK_CSV_TRANSFER": str(sandbox / "transfer.csv"),
        "ISMIP7_MAP_CHECK_FINAL_NATIVE": str(sandbox / "native_final.h5"),
        "ISMIP7_MAP_CHECK_FINAL_TRANSFER": str(sandbox / "transfer_final.h5"),
        "ISMIP7_MAP_CHECK_AUDIT_JSON": str(out),
        "ISMIP7_MAP_CHECK_PNG": str(sandbox / "overlay.png"),
        "ISMIP7_MAP_CHECK_STATUS": str(status),
        "ISMIP7_MAP_CHECK_NATIVE_LC": "2000", "ISMIP7_MAP_CHECK_NATIVE_LC_COARSE": "5000",
        "ISMIP7_MAP_CHECK_NATIVE_BUFFER_M": "0", "ISMIP7_MAP_CHECK_NATIVE_BNDIDS": "/b0.json",
        "ISMIP7_MAP_CHECK_TRANSFER_LC": "1000", "ISMIP7_MAP_CHECK_TRANSFER_LC_COARSE": "10000",
        "ISMIP7_MAP_CHECK_TRANSFER_BUFFER_M": "20000", "ISMIP7_MAP_CHECK_TRANSFER_BNDIDS": "/b20.json",
    }
    proc, seen = run_script(sandbox, "map_check_audit.script", **env)
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(out.read_text())
    assert set(payload) == {"track", "overlay", "region_budget"}
    assert payload["track"]["native"]["exit_code"] is not None
    assert f"SRUN -n 16 python -u scripts/region_budget.py {sandbox / 'native_final.h5'} --csv {sandbox / 'native.csv'}" in seen
    assert "ISMIP7_LC=2000" in seen and "ISMIP7_LC=1000" in seen and "ISMIP7_MESH=checkpoint" in seen
    assert status.read_text().startswith("finished phase=collect exit_code=0 ")


def test_the_melt_bound_job_ends_with_the_check_s_own_status(sandbox):
    # The check is serial and starts under plain `python`, so stub that the way
    # srun is stubbed: drop `-u` and the script path, hand the rest to the driver.
    stub = sandbox / "bin" / "python"
    stub.write_text('#!/bin/bash\n[ "$1" = -u ] && shift\nshift\n'
                    'exec "$FAKE_PYTHON" "$FAKE_DRIVER" "$@"\n')
    stub.chmod(0o755)
    state = sandbox / "state_1000_final.h5"
    state.write_text("")
    k_npz = sandbox / "K.npz"
    k_npz.write_text("")

    # exit 1 is the verdict "flagged"; the job has to fail with it so that
    # sacct and an afterok dependency both see it
    proc, seen = run_script(sandbox, "check_melt_bound.script", FAKE_RC="1",
                            ISMIP7_LC="1000", ISMIP7_INV_H5=str(state),
                            ISMIP7_K_PER_BASIN_NPZ=str(k_npz))
    assert proc.returncode == 1, proc.stderr
    assert f"ARGV --ocx main --npz {k_npz}" in seen
    assert f"ISMIP7_INV_H5={state}" in seen and "OMP_NUM_THREADS=1" in seen
    assert "check_melt_bound exit status: 1" in proc.stdout

    # with no override the check reads the tracked calibration itself, and an
    # offsets file is passed on the same way a legacy K is
    proc, seen = run_script(sandbox, "check_melt_bound.script",
                            ISMIP7_INV_H5=str(state), ISMIP7_OCX_OCEAN="warm")
    assert proc.returncode == 0, proc.stderr
    assert [ln for ln in seen if ln.startswith("ARGV")][-1] == "ARGV --ocx warm"
    proc, seen = run_script(sandbox, "check_melt_bound.script",
                            ISMIP7_INV_H5=str(state), ISMIP7_OCX_OCEAN="none",
                            ISMIP7_DELTAT_PER_BASIN_NPZ=str(k_npz))
    assert proc.returncode == 0, proc.stderr
    assert [ln for ln in seen if ln.startswith("ARGV")][-1] == f"ARGV --npz {k_npz}"
    # with neither, the check runs on its defaults; site_core.sh's set -u must
    # not trip over the empty argument list (bash before 4.4 did)
    proc, seen = run_script(sandbox, "check_melt_bound.script",
                            ISMIP7_INV_H5=str(state), ISMIP7_OCX_OCEAN="none")
    assert proc.returncode == 0, proc.stderr
    assert [ln for ln in seen if ln.startswith("ARGV")][-1].rstrip() == "ARGV"

    # a missing mesh source or calibration, or two calibrations, stops the
    # job before the check starts
    n_before = len(seen)
    proc, seen = run_script(sandbox, "check_melt_bound.script", ISMIP7_LC="1000")
    assert proc.returncode != 0 and "ISMIP7_INV_H5 is required" in proc.stderr
    proc, seen = run_script(sandbox, "check_melt_bound.script",
                            ISMIP7_INV_H5=str(state),
                            ISMIP7_K_PER_BASIN_NPZ=str(sandbox / "absent.npz"))
    assert proc.returncode == 2 and "melt calibration not found" in proc.stderr
    proc, seen = run_script(sandbox, "check_melt_bound.script",
                            ISMIP7_INV_H5=str(state),
                            ISMIP7_K_PER_BASIN_NPZ=str(k_npz),
                            ISMIP7_DELTAT_PER_BASIN_NPZ=str(k_npz))
    assert proc.returncode == 2 and "both set" in proc.stderr
    assert len(seen) == n_before


def test_the_scalar_processing_job_keeps_the_tool_clear_and_ends_with_the_verdict(sandbox):
    r"""scalar_processing.script: the organisers' tool runs from its own venv
    with the Firedrake environment's PYTHONPATH scrubbed; it reads params.nc
    from the tree, as the organisers' run of it will, and the job stops
    without one; a run that is its own reference has to name its stamped year
    and a paired run must not; the job ends with the tool's status when that
    fails and the comparison's otherwise."""
    venv = sandbox / "tools" / "bin"
    venv.mkdir(parents=True)
    fake = ('#!/bin/bash\n'
            'echo "$(basename "$0") $* PYTHONPATH=${PYTHONPATH-unset}" >> "$FAKE_ENV_DUMP"\n'
            'case "$(basename "$0")" in\n'
            '    ismip7-scalars) [ "$1" = --version ] && { echo "ismip7-scalars 0.1.0"; exit 0; }\n'
            '                    exit "${FAKE_TOOL_RC:-0}" ;;\n'
            '    python) echo 0.5.1 ;;\n'
            'esac\n')
    for name in ("ismip7-scalars", "ismip7-scalars-set-params", "python"):
        (venv / name).write_text(fake)
        (venv / name).chmod(0o755)
    # the comparison runs under the job's own python: drop `-u` and the script
    # path, hand the rest to the driver
    stub = sandbox / "bin" / "python"
    stub.write_text('#!/bin/bash\n[ "$1" = -u ] && shift\nshift\n'
                    'exec "$FAKE_PYTHON" "$FAKE_DRIVER" "$@"\n')
    stub.chmod(0o755)
    tree = sandbox / "tree"
    (tree / "AIS" / "RICE" / "icepack2" / "CORE" / "C007").mkdir(parents=True)
    params = tree / "AIS" / "RICE" / "icepack2" / "params.nc"
    params.write_text("")
    grids = sandbox / "grids"
    grids.mkdir()
    (grids / "af2_AIS_08000m_v1.nc").write_text("")
    native = sandbox / "run_ismip7_scalars.csv"
    native.write_text("")
    env = dict(ISMIP7_TOOLS_VENV=str(sandbox / "tools"), ISMIP7_SCALAR_TREE=str(tree),
               ISMIP7_SCALAR_EXPERIMENT="ssp585", ISMIP7_SCALAR_CONFIGID="C007",
               ISMIP7_SCALAR_OUT=str(sandbox / "out"), ISMIP7_SCALAR_DATAPATH=str(grids),
               PYTHONPATH="/firedrake/site-packages")

    proc, seen = run_script(sandbox, "scalar_processing.script", **env)
    assert proc.returncode != 0 and "ISMIP7_SCALAR_REFYEAR is required" in proc.stderr
    proc, seen = run_script(sandbox, "scalar_processing.script", ISMIP7_SCALAR_HIST="historical",
                            ISMIP7_SCALAR_HIST_CONFIGID="C001", ISMIP7_SCALAR_REFYEAR="2016", **env)
    assert proc.returncode == 2 and "refused" in proc.stderr
    assert not [ln for ln in seen if ln.startswith("ismip7-scalars --region")]

    # exit 1 is the comparison's verdict; the job has to end with it
    proc, seen = run_script(sandbox, "scalar_processing.script", ISMIP7_SCALAR_REFYEAR="2016",
                            ISMIP7_SCALAR_NATIVE_CSV=str(native), FAKE_RC="1", **env)
    assert proc.returncode == 1, proc.stderr
    tool = next(ln for ln in seen if ln.startswith("ismip7-scalars --region"))
    assert "--hist ssp585 --hist-configid C007 --refyear 2016" in tool
    assert f"--modelpath {tree}/AIS --outpath {sandbox}/out/tool" in tool
    assert "--params-path" not in tool and tool.endswith("PYTHONPATH=unset")
    assert not [ln for ln in seen if ln.startswith("ismip7-scalars-set-params")]
    argv = next(ln for ln in seen if ln.startswith("ARGV"))
    assert f"--submission {tree}/AIS/RICE/icepack2/CORE/C007" in argv
    assert f"--params {params}" in argv
    assert "--refyear 2016" in argv and f"--native-csv {native}" in argv
    assert "af2_AIS_08000m_v1.nc sha256" in argv and "--native-af2" not in argv
    assert "compare_scalars exit status: 1" in proc.stdout

    # a paired run hands the comparison the historical's folder, and a run
    # whose native scalars carry the area factor says so
    proc, seen = run_script(sandbox, "scalar_processing.script", ISMIP7_SCALAR_HIST="historical",
                            ISMIP7_SCALAR_HIST_CONFIGID="C001", ISMIP7_SCALAR_NATIVE_AF2="1", **env)
    assert proc.returncode == 0, proc.stderr
    argv = [ln for ln in seen if ln.startswith("ARGV")][-1]
    assert f"--hist-submission {tree}/AIS/RICE/icepack2/CORE/C001" in argv
    assert "--native-af2" in argv

    # a failed tool ends the job with its status and no comparison
    n = len([ln for ln in seen if ln.startswith("ARGV")])
    proc, seen = run_script(sandbox, "scalar_processing.script", ISMIP7_SCALAR_REFYEAR="2016",
                            FAKE_TOOL_RC="2", **env)
    assert proc.returncode == 2 and "ismip7-scalars exit status: 2" in proc.stdout
    assert len([ln for ln in seen if ln.startswith("ARGV")]) == n

    # a tree without params.nc stops the job before the tool runs, and says
    # how to write it
    params.unlink()
    n = len([ln for ln in seen if ln.startswith("ismip7-scalars --region")])
    proc, seen = run_script(sandbox, "scalar_processing.script", ISMIP7_SCALAR_REFYEAR="2016", **env)
    assert proc.returncode == 2 and f"not found: {params}" in proc.stderr
    assert "ismip7-scalars-set-params --region AIS --group RICE" in proc.stderr
    assert len([ln for ln in seen if ln.startswith("ismip7-scalars --region")]) == n
