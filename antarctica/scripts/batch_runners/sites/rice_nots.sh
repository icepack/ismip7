# Rice NOTS. Every number here was measured on the login node in September
# 2026; antarctica/scripts/batch_runners/readme.md says how.

ISMIP7_SITE_NAME="rice_nots"
ISMIP7_SITE_MATCH="nots* login*.nots.rice.edu bb*"

# Built by build_firedrake_rice.sbatch (job 1288484). libGLU is in the list because
# the gmsh wheel dlopens libGLU.so.1 at import and the compute nodes lack it.
ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-/projects/ah301/sw/venv-firedrake/bin/activate}"
ISMIP7_MODULES="${ISMIP7_MODULES:-foss/2023b Python/3.11.5 CMake/3.27.6 M4/1.4.19 flex/2.6.4 Bison/3.8.2 libevent/2.1.12 zstd/1.5.5 Szip/2.1.1 libaec/1.0.6 libGLU/9.0.3}"
# The tree those modules live in, named because only the login nodes' login
# shell puts it on MODULEPATH: a login shell on a compute node gets the
# system's /etc/modulefiles and /usr/share/modulefiles, and a job submitted
# from there carries that path, so foss/2023b is "unknown", libbz2 is missing
# and the import of PETSc fails (the 16 October 2026 dH/dt continuation,
# submitted by a chain's follow-on command through `env -i bash -l`, failed
# eight links in a row that way).
ISMIP7_MODULE_USE="${ISMIP7_MODULE_USE:-/opt/apps/modules/Core}"

# The default account (commons) is granted and needs no -A. It gives commons
# (1 day, 142 nodes), long (3 days, 61 nodes, the only one that fits a 2 km
# inversion), debug (30 min) and scavenge (1 hour, preemptible). The deepsC
# condo (52 uniform Cascade Lake nodes, 187 GB each, no wall limit) needs the
# separate `deepsc` account from the EEPS condo owner, which we do not hold.
ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-long}"
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-commons}"
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-scavenge}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"
# Sapphire Rapids (192 threads, 257-515 GB) is what the 2 km inversions need,
# because `long` holds exactly one Cascade Lake node. Everything else stays on
# cascadelake (80 threads, 187 GB): the measured timings and the chain depth
# come from that generation, and build_firedrake_rice.sbatch configures with
# -march=native, so a build has to land where the runs do.
ISMIP7_CONSTRAINT_INV="${ISMIP7_CONSTRAINT_INV:-sapphirerapids}"
ISMIP7_CONSTRAINT_FWD="${ISMIP7_CONSTRAINT_FWD:-cascadelake}"

# /home is a 10 TB NFS export and /scratch is purged, so the 20 TB /projects
# share is what ISMIP7_WORK names. ISMIP7_REPO follows the checkout the command
# was run from rather than naming one, as every other site file here does:
# submit.sh cds to ISMIP7_REPO, so a pinned path makes the same command, run
# from a second checkout, submit the FIRST tree's code -- silently, once both
# trees carry the same script names.
ISMIP7_REPO="${ISMIP7_REPO:-$ISMIP7_REPO_SELF}"
ISMIP7_WORK="${ISMIP7_WORK:-/projects/ah301}"

# The code root moves with the invocation; the shared artifacts cannot. Here
# they are one tree a second checkout has none of: the AIS forcing tree is
# ~313 GB, the 2 km meshes, the budd/RC MAPs and the observational rasters are
# hundreds of MB, and all of them are gitignored rather than copied per
# checkout. Spelled out rather than taken from ISMIP7_WORK, which
# sites/local.env lets each user point at their own project space: that would
# move the data with it, into a tree holding none. The meshes and MAPs sit
# beside them; a run from a second checkout names ISMIP7_MESH on the submit
# line, as the 2 km inversions do, with ISMIP7_MAP_OUT for an inversion (where
# the MAP is written) or ISMIP7_INVERSION for a forward (which MAP to read).
ISMIP7_DATA_ROOT="${ISMIP7_DATA_ROOT:-/projects/ah301/ismip7/ISMIP7/AIS}"
ISMIP7_OBS_DATA_ROOT="${ISMIP7_OBS_DATA_ROOT:-/projects/ah301/ismip7/antarctica/data}"
# One PyOP2 kernel cache for every job (site_core.sh, ISMIP7_SHARED_JIT_CACHE):
# a chain link starts from the kernels its predecessors compiled.
ISMIP7_SHARED_JIT_CACHE="${ISMIP7_SHARED_JIT_CACHE:-/projects/ah301/sw/pyop2-cache}"

ISMIP7_TASKS="${ISMIP7_TASKS:-32}"
ISMIP7_MEM="${ISMIP7_MEM:-240G}"
# Forwards take a whole Cascade Lake node. The partition probe (job 1592757,
# 22 September 2026, antarctica_10000_1000_buffered20000, 3.7 M cells) shows
# this build partitions by locality: ghost/owned 0.001 at 2 ranks, 0.006 at
# 16 and 0.010 at 32 (max 0.017), against the simple partitioner's 2.97 at
# 2 and 287 at 32; halo dofs are 1.0 % of owned at 32 ranks. A 1 km / 10 km
# control from a transferred 2 km MAP then ran 92 s per 0.05 yr step on 32
# ranks under scpc_gamg (job 1592597), within 180 GB. The older 12 ranks /
# 96 GB figure was a 2500 m measurement (readme.md's timing table).
ISMIP7_TASKS_FWD="${ISMIP7_TASKS_FWD:-32}"
ISMIP7_MEM_FWD="${ISMIP7_MEM_FWD:-180G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-3-00:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-1-00:00:00}"

# The timing campaign (`submit.sh script`, antarctica/Makefile) runs where the
# forwards do unless told otherwise, and a lane that one such node cannot hold
# is recorded not runnable rather than resized. The figures are the ones above:
# Cascade Lake is 80 threads, so 40 physical cores under --hint=nomultithread,
# and 187 GB; Sapphire Rapids is 192 threads and 257 to 515 GB depending on the
# node, so it gets a core limit and no single memory limit. On cascadelake the
# 64-rank lanes and the 240G (500 m) and 192G (1 km inversion) requests do not
# fit; ISMIP7_CONSTRAINT_TIMING=sapphirerapids in sites/local.env runs them.
ISMIP7_CONSTRAINT_TIMING="${ISMIP7_CONSTRAINT_TIMING:-$ISMIP7_CONSTRAINT_FWD}"
case "$ISMIP7_CONSTRAINT_TIMING" in
    cascadelake)
        ISMIP7_CORES_PER_NODE="${ISMIP7_CORES_PER_NODE:-40}"
        ISMIP7_MEM_PER_NODE="${ISMIP7_MEM_PER_NODE:-187G}" ;;
    sapphirerapids)
        ISMIP7_CORES_PER_NODE="${ISMIP7_CORES_PER_NODE:-96}" ;;
esac
