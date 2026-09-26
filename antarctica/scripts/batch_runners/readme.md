# Slurm batch runners

Scheduler job scripts for the ISMIP7 icepack2 runs, on any cluster. Every wall
time and memory figure here was measured.

## How it works

| | |
|---|---|
| `sites/<cluster>.sh` | everything that differs between clusters: Firedrake venv (or container image), module loads, partitions by job length, account, node constraint, per-node limits, and paths that hold a checkout and 300 GB of forcing. Tracked. |
| `sites/local.env` | everything that differs between two people on one cluster: the account to charge, a private build or image, the work filesystem, mail flags. Ignored by git; `sites/local.env.example` is the blank. Read before the site file, so it wins. |
| `site_core.sh` | picks the site file (`ISMIP7_SITE`, else hostname), refuses to run while a required field is empty, and provides `ismip7_activate`, `ismip7_mpirun`, `ismip7_chain_resources` and `ismip7_banner`. Chooses no model configuration. |
| `site_env.sh` | `site_core.sh` plus the model defaults the inversion and projection runners share (`ISMIP7_LC`, `ISMIP7_MESH`, `ISMIP7_FRICTION`, ...). The timing campaign's scripts source `site_core.sh` alone, because a lane's model settings are its exports and nothing else. |
| `submit.sh` | composes the scheduler command from the site file, so `submit.sh inversion ISMIP7_LC=2000` and `make -C antarctica timing` are the same line everywhere. |

The job scripts carry log paths and nothing else. A `#SBATCH` line is parsed
before any shell runs, so it cannot read a site file, and a header that
disagrees with the command line is rejected outright ("Requested node
configuration is not available"). A job submitted by hand without resources
stops immediately and says so.

```bash
submit.sh inversion  ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_5000_2000_buffered0.msh
submit.sh projection ISMIP7_EXPERIMENT=ssp585_cesm_waccm
submit.sh smoke                                  # minutes, debug partition
submit.sh inversion --dry-run                    # print the sbatch line only
submit.sh projection --partition debug --time 00:30:00 --tasks 8
```

Any `KEY=VALUE` argument is exported into the job, which is how a run is
configured (knobs in `antarctica/README.md`). One argument holding several
pairs is refused, since that is how zsh passes an unquoted `$VAR`: spell the
pairs out, or write `${=VAR}` in zsh. A value that really contains a space
followed by `NAME=` can be exported in the calling shell instead. `--tasks`,
`--mem`, `--time`, `--partition`, `--constraint`, `--account` and `--name`
override site defaults for one submission, in either `--opt value` or
`--opt=value` form.

`submit.sh` sets every `KEY=VALUE`, with `ISMIP7_SITE` and `ISMIP7_REPO`, in
sbatch's own environment and passes a bare `--export=ALL`, which carries that
environment whole into the job and its chain links. The printed line is the
command that runs:
`env ISMIP7_SITE=... ISMIP7_REPO=... KEY=VALUE ... sbatch ... --export=ALL <script>`.
An `--export` list has two faults: sbatch splits it on commas
(`ISMIP7_SUBCYCLES=1,4,16,64` would arrive as `1`), and any list sets
`SLURM_GET_USER_ENV=1`, under which slurmd rebuilds the login environment when
the job starts and requeues and holds the job when that fails ("user env
retrieval failed requeued held" on IU Quartz). `ismip7_activate` clears that
variable and `SLURM_EXPORT_ENV`, the list itself, which srun reads as its own
`--export`, so a job submitted with a list passes neither on to its chain
links.

A job therefore starts with the submitting shell's environment alone. A shell
that never read the login scripts, such as the one `ssh host command` starts,
holds no module system, so a job it submits finds no `module`.
`ismip7_activate` then reads `ISMIP7_MODULE_INIT` (default `/etc/profile`, which
a site file can change) and stops with exit status 2 if that defines none.

`submit.sh script PATH` submits any job script that sources `site_core.sh`. It
is how `antarctica/Makefile` and `manage_timing_campaign.py` submit, and how to
run `budd_map_census.script` by hand:

```bash
submit.sh script scripts/batch_runners/budd_map_census.script --cd antarctica \
    --queue debug --tasks 16 --mem 64G --time 00:45:00 ISMIP7_MAP=$PWD/mesh/<map>.h5
```

`check_melt_bound.script` is the serial OCX tripwire core 11 waits on
(issue #11). Each site runs it for itself, and again whenever the OCX ocean on
its mirror changes. The job ends with the check's status, 1 when a basin or a
256 km block is flagged, and its log names the OCX version the reader opened:

```bash
submit.sh script scripts/batch_runners/check_melt_bound.script --cd antarctica \
    --queue debug --tasks 1 --mem 16G --time 00:30:00 \
    ISMIP7_LC=1000 ISMIP7_INV_H5=$PWD/mesh/<map or forward state>.h5
```

`scalar_processing.script` runs the organisers' `ismip7-scalars` over one
experiment of a tree written by `write_ismip7_output.py`, then
`compare_scalars.py` against the model's own scalars (issue #13). The tool
lives in a Python 3.11+ venv of its own, named by `ISMIP7_TOOLS_VENV` (a
per-user path, so `sites/local.env`; `antarctica/README.md` has the install).
The tree carries the model's densities in `AIS/<source_id>/<ism_id>/params.nc`,
written by `ismip7-scalars-set-params`, and the job reads that file and stops
without it. The job ends with the comparison's status. A run whose native scalars
integrate over true area, as a current forward writes them, adds
`ISMIP7_SCALAR_NATIVE_AF2=1`. A 286-year experiment on the 8 km grid took
3 min and at most 3.7 GiB on one core (IU Quartz, September 2026):

```bash
submit.sh script scripts/batch_runners/scalar_processing.script --cd antarctica \
    --queue short --tasks 1 --mem 16G --time 01:00:00 \
    ISMIP7_TOOLS_VENV=<venv> ISMIP7_SCALAR_TREE=<the writer's --out-dir> \
    ISMIP7_SCALAR_EXPERIMENT=ssp585 ISMIP7_SCALAR_CONFIGID=C007 \
    ISMIP7_SCALAR_REFYEAR=2016 ISMIP7_SCALAR_OUT=<dir>
```

`--queue short|long|debug` names the site's partition by class, `--cd DIR`
submits from `ISMIP7_REPO/DIR` (the timing scripts run from `antarctica/`),
`--dependency` and `--wait` pass through. Standard output is sbatch's alone, so
a caller can read the job id; the composed command goes to standard error.

## Your first day on a cluster

The same six steps for everyone. Nothing tracked is edited unless the cluster
itself is new.

1. `antarctica/scripts/batch_runners/site_recon.sh` on a login node. It only
   reads, and prints accounts, partitions and their limits, node features and
   memory.
2. If `submit.sh smoke --dry-run` says no site matches this host, the cluster
   is new: see the next section. Otherwise the site file already exists.
3. `cp sites/local.env.example sites/local.env` and keep the lines that are
   yours: `ISMIP7_ACCOUNT`, `ISMIP7_WORK`, a private `ISMIP7_FIREDRAKE` or
   `ISMIP7_CONTAINER`, and `ISMIP7_SBATCH_EXTRA` for mail. Write each as
   `VAR="${VAR:-value}"` so a value given for one submission still wins.
4. `submit.sh verify`: four ranks, minutes. Firedrake imports, MPI works across
   ranks, the partitioners partition; on a container site, the image starts in
   the checkout and sees the job's environment.
5. `submit.sh smoke`: a short real inversion on the debug partition.
6. `make -C antarctica timing-dry-run`: every line the timing campaign would
   submit here, with `NOT RUNNABLE` for lanes this site's nodes cannot hold.
   Then `make -C antarctica timing`.

## Adding your cluster

```bash
cp sites/template.sh sites/mycluster.sh
$EDITOR sites/mycluster.sh          # required fields are marked
ISMIP7_SITE=mycluster submit.sh smoke --dry-run
```

`site_recon.sh` on a login node prints what the file needs: associations and
accounts, partitions and their limits, node features and memory, queue depth.
It only reads. Add your hostname pattern to `ISMIP7_SITE_MATCH` and
`ISMIP7_SITE` becomes unnecessary.

Required: `ISMIP7_FIREDRAKE`, `ISMIP7_PART_LONG`, `ISMIP7_PART_SHORT`,
`ISMIP7_PART_DEBUG`, `ISMIP7_REPO`, `ISMIP7_WORK`. Everything else defaults.
`ISMIP7_REPO` defaults to `ISMIP7_REPO_SELF`, the checkout the command was run
from, in every site file. A site whose gitignored artifacts (the forcing tree,
the meshes, the MAPs, the observational rasters) do not live beside the code --
because the cluster holds more than one checkout, as Rice does -- names
`ISMIP7_DATA_ROOT` and `ISMIP7_OBS_DATA_ROOT` in its site file; the mesh and
the MAP follow `ISMIP7_REPO`, so a run from a second checkout states
`ISMIP7_MESH` on the submit line, with `ISMIP7_MAP_OUT` for an inversion (where
the MAP is written) or `ISMIP7_INVERSION` for a forward (which MAP to read).
A site that sets `ISMIP7_CONTAINER` is not asked for `ISMIP7_FIREDRAKE`, and
neither is any `--dry-run`.
A missing value is reported at submission with the file and variable named.
`submit.sh build` asks for the other five only, since it creates the venv that
`ISMIP7_FIREDRAKE` names.

Job sizes come in pairs so inversions and forwards can differ:
`ISMIP7_TASKS_INV` and `ISMIP7_MEM_INV`, `ISMIP7_TASKS_FWD` and
`ISMIP7_MEM_FWD`, each falling back to `ISMIP7_TASKS` and `ISMIP7_MEM`. The
node constraint splits the same way (`ISMIP7_CONSTRAINT_INV`, `_FWD`, both
falling back to `ISMIP7_CONSTRAINT`).

Optional, for every submission: `ISMIP7_SBATCH_EXTRA`, extra sbatch flags split
on spaces (a QOS the partition insists on; your own `--mail-type`/`--mail-user`
in `sites/local.env`). Chain successors are given it again.

Optional, for the timing campaign: `ISMIP7_CONSTRAINT_TIMING` (the node feature
`submit.sh script` asks for; defaults to `ISMIP7_CONSTRAINT_FWD`),
`ISMIP7_CORES_PER_NODE` and `ISMIP7_MEM_PER_NODE` (one such node; physical
cores, and memory as sbatch spells it), and `ISMIP7_TIMING_JIT_CACHE`. The
matrix's lanes are single-node at every site so that the matrices compare. A
lane one node cannot hold is refused by `submit.sh` with exit status 3 and
recorded `not_runnable reason=exceeds_site_cores|exceeds_site_mem`; it shows as
such in `TIMING_MATRIX.md`. A `--constraint` given by hand (`SLURM_CONSTRAINT=`
in the Makefile) skips the check, since the limits describe
`ISMIP7_CONSTRAINT_TIMING`'s nodes and not the ones you named. A lane's
`seconds_per_step` has no warm-up excluded, so the timing scripts keep a kernel
cache that persists between jobs (the one the site's modules or venv already
name, as IU's firedrake modulefile does; else Firedrake's default location; or
`ISMIP7_TIMING_JIT_CACHE`) rather than the private per-job one
`ismip7_activate` gives the runners. Only the kernel cache goes back that way:
loopy's persistent dict is switched off in every job (`LOOPY_NO_CACHE=1`),
because two jobs compiling the same kernel seconds apart race on a shared one,
which is what `make timing-scout` launches, and because the ranks of one job
writing a per-job one on a networked filesystem stalled forever on its sqlite
lock (`site_core.sh` has both incidents). Each record's `host` block says which
site and node measured it and whether that cache started empty.

### A container site

Set `ISMIP7_CONTAINER` to an Apptainer or Singularity image and that is the
whole difference. `ismip7_activate` puts `container_bin/` first on `PATH`; its
`python` runs `apptainer exec <binds> <image> python3 "$@"`, so every `python`
line in the job scripts runs in the image unchanged. `ismip7_mpirun` starts the
ranks with the image's own `mpiexec` inside one `exec`: every job here is one
node, so the MPI in the image never has to agree with the host's Slurm about
PMI. The checkout, `ISMIP7_WORK`, the data roots and the kernel cache are bound
already. `ISMIP7_CONTAINER_ARGS` adds binds, `ISMIP7_CONTAINER_MPIEXEC` gives
that launcher flags (`mpiexec --mca plm isolated` if it objects to the Slurm
allocation around it), `ISMIP7_CONTAINER_RUNTIME` names `singularity` where
that is what is installed. `--cleanenv` must never be passed: a job's
configuration is its `ISMIP7_*` environment.

### The sites that ship

| file | state |
|---|---|
| `sites/rice_nots.sh` | complete, in production. Details below. |
| `sites/iu_quartz.sh` | complete, from the IU Quartz runners on the upstream `timing_matrix` branch: partition `general` (`debug` for tests), account `r00905`, the IU module stack (`module use /N/u/dlilien/Quartz/modulefiles`, then gnu, openmpi, python, zlib, hdf5, openblas, patchelf, petsc, firedrake), 16 ranks for an inversion and 64 for a forward (the fastest production lane of `antarctica/TIMING_MATRIX_QUARTZ_SCPC_GAMG.md`), 128 cores and 515700 MB per node. A second IU user points `ISMIP7_FIREDRAKE` at their own build and `ISMIP7_ACCOUNT` at their own allocation in `sites/local.env`. |
| `sites/uchicago_midway.sh` | a container site, not yet run end to end. It describes the image `icepack2_midway3_source.def` builds and names a persistent kernel cache; the image path, partitions, account and per-node limits are empty, so the first submission refuses until they are filled (the file's header says where to read each) (issue #47). |
| `sites/local.sh` | no scheduler: every setting comes from the environment. For debugging a job script on a workstation, and what the chain tests use. |

## Rice NOTS, in detail

Measured on the login node in September 2026. `sites/rice_nots.sh` is the file
it produced.

The default account is `commons`, so `-A` is unnecessary. Four partitions are
granted:

| partition | wall limit | cap | nodes | cascadelake nodes |
|---|---|---|---|---|
| commons | 1 day | 1536 cpu / 2304 GB | 142 | 72 |
| long | 3 days | 384 cpu / 1 TB | 61 | 1 |
| debug | 30 min | 1 job | 2 | |
| scavenge | 1 hour | preemptible | 154 | many |

Cascade Lake is reachable today on `commons` and `scavenge`. Pinning it on
`long` would queue behind a single machine, so `sites/rice_nots.sh` sets
`ISMIP7_CONSTRAINT_INV=sapphirerapids` (`long` is where the 2 km inversions
run) and `ISMIP7_CONSTRAINT_FWD=cascadelake`, the generation the timings below
were measured on and the one the Firedrake build matches.

**deepsC is refused.** `sbatch --test-only -p deepsC` returns "Invalid account
or account/partition combination". It needs the separate `deepsc` account,
described in Slurm as Prof. Alan Levander's, with 76 users already in it and no
Slurm coordinators, so adding a user goes through the help desk. It is worth
having because it has no wall-time limit (a 2500 m inversion runs 20 to 26
hours against the one-day `commons` limit) and is less contended (182 pending
against 624). It does not solve memory: all 52 deepsC nodes are 187 GB, so the
2 km / 5 km-interior inversion at about 255 GB needs `long` and its 515 GB
Sapphire Rapids nodes.

Once granted, switch on the wrapper. Clear the constraint as well, since every
deepsC node is Cascade Lake while inversions default to Sapphire Rapids:

```bash
submit.sh inversion --partition deepsC --constraint '' --time 2-00:00:00
```

deepsC nodes are uniform: 40 physical cores (2 x 20), 80 logical
(`ThreadsPerCore=2`), 187135 MB, features `cascadelake,opath`, no wall limit.
The big-memory nodes live elsewhere (`commons` and `long` reach 1.5 TB, `as143`
748 GB). The solver is memory-bandwidth bound, so `submit.sh` passes
`--hint=nomultithread` and ranks land on physical cores. The CPU matches the
workstation the timings came from (Xeon Gold 5218R, 2.10 GHz), so they transfer
one to one and are conservative, the workstation having been oversubscribed.

## Building Firedrake on a cluster

The lessons generalise. The module names do not.
`build_firedrake_rice.sbatch` is the worked example and runs at Rice only.
Another site copies it to `build_firedrake_<site>.sbatch` and substitutes its
own stack; `submit.sh build` says so when the site is not `rice_nots`.

Five submissions were needed at Rice, each failure a gap between login and
compute nodes:

| what broke | why | fix |
|---|---|---|
| `firedrake-configure` refused | RHEL 9.4 is unsupported | `--os unknown`, PETSc builds deps from source |
| `git: command not found` | compute nodes have no git | fetch the PETSc tarball with `curl` |
| Bison would not build | no `flex`, `bison`, `m4` | load `M4/1.4.19 flex/2.6.4 Bison/3.8.2` |
| PnetCDF configure failed | OpenMPI adds `-levent_core`, libevent off `LIBRARY_PATH` | load `libevent/2.1.12` |
| NetCDF link failed | `ld: cannot find -lzstd` | load `zstd/1.5.5 Szip/2.1.1 libaec/1.0.6` |
| modules "loaded" but gcc was 11 | a bare name like `flex` resolves against the wrong GCCcore and drops the set back to `/usr/bin` | pin every version, abort if `gcc` is not under `/opt/apps` |

The last one is silent: `module load` returns 0 and you get gcc 11, Python 3.9
and no `mpicc`. The script fails loudly now. While debugging, never pipe
`module load` through `head` or `grep`; SIGPIPE breaks the load and you measure
an environment the job never sees.

On partitioners: upstream's own option set carries `--download-ptscotch`, so a
stock build works. The workstation build was configured `--with-ptscotch`
against a system scotch that supplied nothing, and PETSc fell back to `simple`,
which cuts by point index. That is harmless on a structured mesh and
catastrophic on a gmsh-ordered one (ghost/owned 16.0 at 4 ranks, 287 at 32).
The build adds `--download-parmetis` so two partitioners can be cross-checked.

Rice has no Firedrake module, so the build uses EasyBuild `foss/2023b` (GCC
13.2.0, OpenMPI 4.1.6, OpenBLAS, ScaLAPACK, FFTW) with Python 3.11.5 on the
same GCCcore, CMake 3.27.6, and the M4, flex, Bison, libevent, zstd, Szip and
libaec modules the gaps above require. `ISMIP7_MODULES` must name that build
set and `libGLU/9.0.3` on top of it, which the `gmsh` wheel dlopens at import.
A successful build prints the whole list at the end.

### After the build

`install_deps.sh` (once, on a login node) pip-installs the data stack and the
four editable packages: `icepack`, `icepack2`, `tlm_adjoint`, `icepack_tools`.
It takes the venv and work filesystem from the site file and expects the
sources under `$ISMIP7_WORK/sw/src` (`FD_PREFIX` moves that; at Rice it is
`/projects/ah301/sw/src`). `icepack`, `tlm_adjoint` and `icepack_tools` are
rsynced from a workstation, `icepack_tools` because it is a local project with
no remote at all.

`icepack2` is **cloned at a pin**, and is the one source a site must not rsync.
The inversion needs two lines that Firedrake 2026 forces on it
(`Mesh.geometric_dimension` became an attribute, and `viscous_power` and
`flow_law` call it). Those lines used to be uncommitted edits in one
workstation checkout, so two sites could differ with nothing to read (issue
#46); they are now a commit, on the branch of the open pull request
icepack/icepack2#3. `install_deps.sh` names that branch's head by SHA and
installs exactly it, which is icepack2 `main` (`40e848b`) plus the two files, so
every site runs one known tree and `--check-icepack2` prints which:

```
bash antarctica/scripts/batch_runners/install_deps.sh --check-icepack2
icepack2: e0a46c9ce3e95dc916660a200e695f0417aa3037 (fix/firedrake-2026-geometric-dimension from https://github.com/hoffmaao/icepack2.git)
```

A checkout that is not a clone, or that carries local changes, is reported and
left alone rather than reset: at a site still holding the rsynced copy, move it
aside (`mv icepack2 icepack2.rsynced`) and run the script again. When the pull
request merges, override `ICEPACK2_REMOTE`, `ICEPACK2_REF` and `ICEPACK2_SHA` to
follow `main` and change the defaults in `install_deps.sh` in the same pass;
`tests/test_icepack2_pin.py` pins them, so that edit is a deliberate one. As of
22 September 2026 the workstation and Rice carried these edits and Quartz ran a
clean `40e848b`, which is the disagreement the pin ends.

Two more gaps surfaced here: `/tmp` is not writable on the login
nodes (the script sets `TMPDIR`), and the `gmsh` wheel dlopens `libGLU.so.1`.

`verify.sbatch` proves the build works across ranks: four tasks under `srun`,
each partitioner on a unit square, then the real 2500 m mesh under `ptscotch`
and `simple`. The build job's own check ran `srun -n 4` inside a one-task
allocation, which fails for every type and proves nothing.

`smoke.sbatch` runs the 32 km inversion for two iterates on four ranks, a few
minutes on `scavenge`, so a multi-day job cannot die on a missing file hours
after queueing. Run it after any change to the stack. No mesh is tracked, so
at a site that has none the job first builds `antarctica_320000_32000_buffered0.msh`
and its boundary-id sidecar (a minute of serial gmsh, once per checkout); the
unsuffixed `antarctica_320000_32000.msh` an older generator left is used where
it exists. A mesh named through `ISMIP7_MESH` is never built.

## The job scripts

### `map_check_score.script` and `map_check_audit.script`

`make -C antarctica map-check` (`antarctica/MAP_CHECK.md`) takes one released
MAP through its checks with `timing_redistribute.script`,
`timing_prepare.script`, `timing_cache_audit.script`,
`timing_transient.script` and `projection.sbatch`, plus these two. The score
script runs `scripts/score_map.py --json` on the MAP's own mesh or, with
`ISMIP7_MAP_CHECK_RESTART`, on a prepared map-check cache, and under Budd the
`check_budd_map.py` shelf-gate census. The audit script runs
`check_ismip6_track.py`, `compare_runs.py` and `region_budget.py` over the two
ten-year controls and collects every exit code and output into one JSON. Both
source `site_core.sh` alone and walk their status file running to finished or
failed, as the timing scripts do.

### `partition_probe.sbatch`, run first after a build

Distributes the mesh at 1 to 32 ranks and reports the ghost-to-owned dof ratio.
Hundredths mean a working parallel partitioner. Tens or hundreds mean PETSc
fell back to `simple` and every scaling number would measure that. The
workstation build reports 2.97 at 2 ranks and 287 at 32, which is why no
scaling curve is quoted from it.

### `inversion.sbatch`, self-resuming

Defaults write `inversion_icepack2_rc_n3_dg0_logvelnet_<ISMIP7_LC>.h5` (1000
on the production mesh) under the settings the 2500 m result came from:
sigma-normalised velocity misfit with ISSM's logarithmic term, the pointwise
dH/dt term, and the integrated net mass-balance constraint that is off by
default in the repo.

`site_env.sh` defaults `ISMIP7_FRICTION` to `regularized_coulomb` everywhere.
Budd's shelf gate was a sign test on the roundoff residue of the effective
pressure, so Budd MAPs predating the fix need re-inverting; pass
`ISMIP7_FRICTION=budd` for those. The law picks the filename tag, so a Budd
re-inversion writes its own `_budd` MAP.

`submit_inversions.sh [B|C|BC]` submits the two 2 km strategies (B: cell-mean
BedMachine sampling on the 20 km-interior mesh; C: the 5 km-interior mesh with
vertex sampling), each named so the converged 2500 m MAP is never touched, with
partition, constraint and rank count from the site file.

**The chain.** A 2 km inversion can outlast a wall limit. The inversion
checkpoints `ISMIP7_MAP_OUT` every 20 iterates and warm-starts from it, so each
job queues its successor first with `--dependency=afterany`, copying partition,
constraint, memory, time and task layout from `scontrol`. Every link exits at
once if `<map>.done` exists, meaning the MAP reached disk. The driver writes
that marker as soon as the checkpoint write returns, so a kill in the tail
(final solve, summary figure) cannot lose it; the runner's post-`srun` grep for
the driver's `Saved MAP:` line is the fallback. Depth is capped by
`ISMIP7_CHAIN_MAX` (4). A warm start whose mesh dof ordering differs from the
run's own is refused, since a rank-count change mid-chain would scramble theta
and phi silently. Regression test: `tests/test_inversion_chain.py`.

### `projection.sbatch`, self-chaining

A 285-year projection is about two days on the production mesh, so this
resubmits itself with `--dependency=afterok` until the run reaches its end
year, resuming through `ISMIP7_AUTO_RESUME=1`. Each driver owns its end year
and the runner does not default `ISMIP7_T_END`; the chain reads the value the
run used from the driver's `Time-stepping: <start>-><end>` line. At 1000 m /
10 km on 64 ranks under `scpc_gamg`, 24 h buys at most 140 simulated years, so
a full projection is about three links. At Rice's default of 32 ranks on one
Cascade Lake node, 31 min a year buys about 46 years a link, so a full
projection is about seven.

```bash
submit.sh projection ISMIP7_EXPERIMENT=control
```

| `ISMIP7_EXPERIMENT` | core | period |
|-------|------|----------------|
| `control` | 9 or 10 (by `ISMIP7_ESM`) | 2015-2300 |
| `ssp126_cesm_waccm` / `ssp126_mri_esm2` | 5 / 6 | 2015-2300 |
| `ssp370_cesm_waccm` / `ssp370_mri_esm2` | 3 / 4 | 2015-2100 |
| `ssp585_cesm_waccm` / `ssp585_mri_esm2` | 7 / 8 | 2015-2300 |
| `ocx` | 11 | 1979-2025 |
| `hist_cesm_waccm` / `hist_mri_esm2` | 1 / 2 | 1850-2014 |

The runner writes the submission's yearly fields and scalars by default (`ISMIP7_OUTPUT=1`), because every experiment it offers is a core experiment and a projection that reaches 2300 without them has to be run again. `ISMIP7_OUTPUT=0` turns that off for a pipeline exercise.

**The chain stops on a non-zero exit and never retries.** The July
grounding-line blow-up looked like a run that needed more time, and chaining
through it would have burned days.

Each job asks its own final checkpoint what happened, and that verdict owns the
exit code. `stalled=1` (rescue ladder and subcycles both exhausted) exits 1 and
submits nothing, so a stall is never mailed out as a completed job. A clean
exit at the end year finishes. A missing `Time-stepping:` line stops the chain,
since the job cannot tell where it got to. A clean exit short of the end year
resubmits, unless `ISMIP7_CHAIN=0`, or auto-resume is off (a successor would
cold-start and repeat the years), or the final checkpoint is unreadable, or the
job advanced no years at all (setup alone spent the budget).

A relaunch from a stalled state has cleared the diagnostic-Newton wall in every
observed case (see Known issues in `antarctica/README.md`). The chain leaves
that judgement to you: resubmit, and auto-resume picks the run up.

**MPI-IO goes through romio.** Every rank is started with
`OMPI_MCA_io=romio321` (`ismip7_mpirun` in `site_core.sh`), because Open MPI's
default ompio component writes a parallel HDF5 file to an NFS file system at
about 1.5 MB/s: on NOTS, 32 ranks on the 1000 m mesh, one 1.48 GB yearly
output file took 1042 s to `/scratch` under ompio and 52 s under romio321
(24 September 2026), on a file system that takes a single stream at 600 MB/s.
Before the change the yearly write was a third of a 45-minute model year.
`ISMIP7_MPI_IO` names the component; set it empty to leave the MPI's own
default. An MPI other than Open MPI ignores the variable.

The wall budget is derived per job. Each link reads its own partition's
`TimeLimit`, holds back 25 minutes and passes the rest as
`ISMIP7_WALL_STOP_MIN`, so the model stops a step early and writes a complete
final checkpoint for the successor.

## Measured costs

| run | ranks | memory | wall | core-hours |
|---|---|---|---|---|
| 2500 m inversion, 60-80 iterates | 12 | 80 GB | 1 day | 300 |
| 2500 m forward, per simulated year | 12 | 80 GB | 26 min | 5 |
| 2500 m projection, 285 years | 12 | 80 GB | 5 days | 1,470 |
| 2 km / 20 km interior inversion | 12 | 120 GB | 1.5 days | 450 |
| 2 km / 5 km interior inversion | 12 | 255 GB | 2-3 days | 800 |
| full 11-experiment set at 2500 m | | | | 15,000 |
| 1000 m / 10 km forward, per simulated year (Quartz, `scpc_gamg`) | 64 | under 70 GB | 10 min | 11 |
| 1000 m / 10 km projection, 285 years (Quartz, `scpc_gamg`) | 64 | under 70 GB | 2 days | 3,040 |
| 1000 m / 10 km forward, per simulated year (Rice, `scpc_gamg`) | 32 | under 180 GB | 31 min | 16 |

The two Quartz rows are the production configuration, from
`antarctica/TIMING_MATRIX_QUARTZ_SCPC_GAMG.md`: the transient loop of a
ten-step lane at `dt = 0.05` under the matrix's strict contract, extrapolated.
Setup, forcing updates and output are not in them, and the memory is 64 times
the largest rank's peak. The Rice row is a 1 km control from a transferred
2 km MAP on one Cascade Lake node (job 1592597), 92 s per `dt = 0.05` step.
The 2500 m and 2 km rows are whole runs on Cascade Lake under `full_mumps`.

Per iterate at 2500 m on 12 ranks: forward median 1081 s (p10 932, p90 1365),
adjoint 92 s, iterate 1174 s. The adjoint is 8% of the iterate, so the cost
sits in the forward's Newton continuation. Per forward step: 155 s median, of
which the level set is 13 s.

The eleven-experiment set is about nineteen node-days, and the experiments are
independent, so a handful of nodes finishes it inside a week.

## Open items

- Rice forwards take 32 ranks and 180 GB, one Cascade Lake node. The
  partition probe on the 1 km / 10 km mesh (job 1592757, 22 September 2026)
  reports ghost/owned 0.010 at 32 ranks (max 0.017, halo 1.0 % of owned),
  so the build partitions by locality and rank counts up to a node are
  meaningful there. The forward cost is the Rice row of the table above.
- An inversion factors the complete mixed Jacobian with MUMPS, which sets its
  memory; `tlm_adjoint` differentiates through that solve, so no setting
  changes it. Cluster forwards took the field split this item asked for:
  `projection.sbatch` defaults to `scpc_gamg`, which eliminates the cell-wise
  stress and traction blocks exactly and puts multigrid on the condensed
  velocity operator (section 7 of `antarctica/README.md`).
