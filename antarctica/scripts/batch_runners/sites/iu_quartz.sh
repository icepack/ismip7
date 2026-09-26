# Indiana University Quartz. Taken from the IU runners on the upstream
# `timing_matrix` branch (antarctica/scripts/batch_runners/timing_*.script),
# which is where these module names, the account and the venv path come from.
# The defaults below point at an existing IU build. A second IU user keeps the
# module stack and, in sites/local.env, points ISMIP7_FIREDRAKE at their own
# build (or at that one, if it is readable to them) and ISMIP7_ACCOUNT at
# their own allocation.

ISMIP7_SITE_NAME="iu_quartz"
# The login nodes are h1 and h2 and the compute nodes c<N>, all under
# quartz.uits.iu.edu (`hostname -f`); the bare forms cover a short hostname.
ISMIP7_SITE_MATCH="*.quartz.uits.iu.edu quartz* h[0-9]*.quartz*"

ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-/N/u/dlilien/Quartz/sw/firedrake/2026.04/firedrake/bin/activate}"
# `module use` first: the firedrake and petsc modulefiles are in a personal
# tree, not the system one.
ISMIP7_MODULES="${ISMIP7_MODULES:-gnu/9.3.0 openmpi/4.0.5 python/3.14.5 zlib/1.2.13 hdf5/1.12 openblas/0.3.13 patchelf/0.18.0 petsc/3.25.5 firedrake/2026.04}"
ISMIP7_MODULE_USE="${ISMIP7_MODULE_USE:-/N/u/dlilien/Quartz/modulefiles}"

ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-general}"
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-general}"
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-debug}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-r00905}"
ISMIP7_CONSTRAINT="${ISMIP7_CONSTRAINT:-}"

# The checkout the submission came from (site_core.sh works it out), so a
# second worktree runs its own tree without anything being edited.
ISMIP7_REPO="${ISMIP7_REPO:-$ISMIP7_REPO_SELF}"
ISMIP7_WORK="${ISMIP7_WORK:-$HOME}"

# Quartz does not export SCRATCH, so site_core.sh's per-job kernel cache,
# ${SCRATCH:-$HOME}/.pyop2_cache/$SLURM_JOB_ID, lands in a home directory that
# nothing purges, while /N/scratch/$USER purges at 30 days
# (00_SCRATCH_FILES_DELETED_AFTER_30_DAYS.txt sits in its root). Naming it here
# puts that cache on the filesystem the firedrake modulefile already points
# PYOP2_CACHE_DIR and FIREDRAKE_TSFC_KERNEL_CACHE_DIR at
# (/N/scratch/$USER/firedrake.cache), which is the warm cache
# ismip7_persistent_jit_cache hands a timing lane back to. Left unexported:
# site_core.sh reads it in this same shell, and a job has no business
# inheriting a SCRATCH the scheduler never set. The id -un fallback is for
# site_core.sh's set -u.
SCRATCH="${SCRATCH:-/N/scratch/${USER:-$(id -un)}}"

# The IU timing runs use 12 to 16 ranks per node, 128 to 240 GB, up to 48 h.
ISMIP7_TASKS="${ISMIP7_TASKS:-16}"
ISMIP7_MEM="${ISMIP7_MEM:-240G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-48:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-24:00:00}"
# Forwards take 64 ranks, the fastest lane the matrix measured on the
# production mesh: 10.0 min per simulated year at 1000 m / 10 km and dt 0.05
# under scpc_gamg, against 15.7 on 32 ranks and 32.3 on 16
# (antarctica/TIMING_MATRIX_QUARTZ_SCPC_GAMG.md); the production step, 0.025,
# takes twice the steps. That lane peaked at 1.0 GiB
# a rank (sacct MaxRSS, job 10524648), so the 240G above is ample. The inversion keeps the 16: its full
# mixed-Jacobian MUMPS solve is sized by memory, not by this. A bare
# `submit.sh script` takes the forward's size too; give a small job --tasks.
ISMIP7_TASKS_FWD="${ISMIP7_TASKS_FWD:-64}"

# One node, from `sinfo -p general,debug -N -o "%c %m"` (September 2026): all 90
# are 128 cores and 515700 MB, with no feature to choose between. `submit.sh
# script` refuses a timing lane that asks for more than this.
ISMIP7_CORES_PER_NODE="${ISMIP7_CORES_PER_NODE:-128}"
ISMIP7_MEM_PER_NODE="${ISMIP7_MEM_PER_NODE:-515700M}"
