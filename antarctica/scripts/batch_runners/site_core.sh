#!/bin/bash
# The site half of the ISMIP7 batch environment: which cluster this is, how to
# get Firedrake on it, and what the scheduler wants. Nothing here submits
# anything, and nothing here chooses a model configuration.
#
# site_env.sh sources this and adds the model defaults (ISMIP7_LC, ISMIP7_MESH,
# ISMIP7_FRICTION, ...) that the inversion and projection runners share. The
# timing campaign's job scripts source this file alone: the campaign manager
# exports every model setting a lane needs, run_timing.py refuses a matrix lane
# that arrives with ISMIP7_MESH set, and a defaulted ISMIP7_FRICTION would
# change the law under any stage that did not export one.
#
# One file per cluster lives in sites/ and answers three questions: where is
# Firedrake, what does the scheduler want, and where does the data live. This
# file picks the right one, checks it is filled in, and turns it into an
# environment. Rice NOTS, IU Quartz and a UChicago Midway stub ship with the
# repository; sites/template.sh is the blank to copy for anywhere else.
#
#   ISMIP7_SITE=iu_quartz sbatch ...     # name it
#   sbatch ...                           # or let the hostname choose
#
# Every value a site file sets can be overridden per submission, because they
# are all written as ${VAR:-default}: `ISMIP7_PART_LONG=debug submit.sh ...`
# works without editing anything.
#
# What differs between two people on the same cluster (the account to charge,
# a private Firedrake build, where the forcing tree lives, a mail address) goes
# in sites/local.env, which git ignores. It is read before the site file, so
# its values win by the same ${VAR:-default} rule, and it should be written in
# that form too so that a per-submission override still beats it.
# sites/local.env.example is the blank. ISMIP7_LOCAL_ENV names a different file
# (the tests point it at /dev/null to stay independent of the checkout's own).
set -u

_ISMIP7_BR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The checkout this file belongs to, for a site file to default ISMIP7_REPO
# to. `cd` and `pwd` stay logical on purpose: a sandbox that reaches this tree
# through a symlink (the chain tests do) has to resolve to the sandbox.
ISMIP7_REPO_SELF="$(cd "$_ISMIP7_BR_DIR/../../.." && pwd)"

_ISMIP7_LOCAL_ENV="${ISMIP7_LOCAL_ENV:-$_ISMIP7_BR_DIR/sites/local.env}"
if [ -f "$_ISMIP7_LOCAL_ENV" ] && [ -r "$_ISMIP7_LOCAL_ENV" ]; then
    # shellcheck disable=SC1090
    . "$_ISMIP7_LOCAL_ENV"
fi

ismip7_site_file() {
    local host want f name match pat
    host="$(hostname -f 2>/dev/null || hostname)"
    want="${ISMIP7_SITE:-}"
    if [ -n "$want" ]; then
        f="$_ISMIP7_BR_DIR/sites/$want.sh"
        [ -r "$f" ] && { echo "$f"; return 0; }
        echo "ERROR: ISMIP7_SITE=$want but $f is not readable." >&2
        echo "       Available: $(cd "$_ISMIP7_BR_DIR/sites" && ls *.sh | sed 's/\.sh$//' | tr '\n' ' ')" >&2
        return 2
    fi
    for f in "$_ISMIP7_BR_DIR"/sites/*.sh; do
        name="$(basename "$f" .sh)"
        [ "$name" = "template" ] && continue
        match="$(sed -nE 's/^ISMIP7_SITE_MATCH="([^"]*)".*/\1/p' "$f" | head -1)"
        # The patterns are globs for `case` to match the hostname against, so
        # keep the shell from expanding them against the invoking directory.
        set -f
        for pat in $match; do
            # shellcheck disable=SC2254
            case "$host" in $pat) set +f; echo "$f"; return 0 ;; esac
        done
        set +f
    done
    echo "ERROR: no site definition matches host '$host'." >&2
    echo "       Copy $_ISMIP7_BR_DIR/sites/template.sh to sites/<name>.sh, fill it in," >&2
    echo "       then submit with ISMIP7_SITE=<name> (or add the hostname to its" >&2
    echo "       ISMIP7_SITE_MATCH so it is found automatically)." >&2
    return 2
}

_ISMIP7_SITE_FILE="$(ismip7_site_file)" || exit 2
# shellcheck disable=SC1090
. "$_ISMIP7_SITE_FILE"

# The site file must supply these six. Everything else defaults below.
ISMIP7_REQUIRED="ISMIP7_FIREDRAKE ISMIP7_PART_LONG ISMIP7_PART_SHORT ISMIP7_PART_DEBUG ISMIP7_REPO ISMIP7_WORK"
# build_firedrake_rice.sbatch creates the venv, so it runs where ISMIP7_FIREDRAKE
# has nothing to name yet. It asks for the rest.
ISMIP7_REQUIRED_BUILD="ISMIP7_PART_LONG ISMIP7_PART_SHORT ISMIP7_PART_DEBUG ISMIP7_REPO ISMIP7_WORK"

# A site can run Firedrake from a container image instead of a venv: set
# ISMIP7_CONTAINER to the image (.sif) and ISMIP7_FIREDRAKE is not asked for.
#   ISMIP7_CONTAINER_RUNTIME   apptainer (default) or singularity
#   ISMIP7_CONTAINER_ARGS      extra `exec` flags, e.g. more --bind pairs. The
#                              checkout, ISMIP7_WORK, ISMIP7_DATA_ROOT and the
#                              kernel cache are bound already. Never --cleanenv:
#                              a job's configuration is its ISMIP7_* environment.
#   ISMIP7_CONTAINER_MPIEXEC   the launcher INSIDE the image (default mpiexec),
#                              with any flags it needs there
# ISMIP7_MODULES is still loaded first, for whatever provides the runtime.
ISMIP7_CONTAINER="${ISMIP7_CONTAINER:-}"
ISMIP7_CONTAINER_RUNTIME="${ISMIP7_CONTAINER_RUNTIME:-apptainer}"
ISMIP7_CONTAINER_ARGS="${ISMIP7_CONTAINER_ARGS:-}"
ISMIP7_CONTAINER_MPIEXEC="${ISMIP7_CONTAINER_MPIEXEC:-mpiexec}"
[ -n "$ISMIP7_CONTAINER" ] && ISMIP7_REQUIRED="$ISMIP7_REQUIRED_BUILD"

ismip7_site_require() {
    local missing=""
    local v
    for v in ${1:-$ISMIP7_REQUIRED}; do
        [ -z "${!v:-}" ] && missing="$missing $v"
    done
    if [ -n "$missing" ]; then
        echo "ERROR: site '$ISMIP7_SITE_NAME' ($_ISMIP7_SITE_FILE) is missing:$missing" >&2
        echo "       Fill them in there, in sites/local.env, or export them for this submission." >&2
        exit 2
    fi
}

# The resource request of the running job, as sbatch flags, for a chain link
# that resubmits its own script. The job scripts carry no resource directives
# (submit.sh composes them from the site file), so a successor submitted from
# inside a job would otherwise land on the cluster's defaults: wrong partition,
# wrong wall limit, and often one task, which the runners refuse outright.
# --hint=nomultithread and the job name are not readable back from scontrol.
# sbatch takes its name default from SBATCH_JOB_NAME, which --export=ALL does
# not carry, so a successor would show up in squeue as the script filename.
# Both are restated here; submit.sh passes them for the first link.
# Both chains call this, so the rule lives here once.
ismip7_chain_resources() {
    local info feat tlim memn
    info="$(scontrol show job "${SLURM_JOB_ID:-}" 2>/dev/null)"
    feat="$(echo "$info" | grep -oE 'Features=[^ ]+' | cut -d= -f2)"
    [ "$feat" = "(null)" ] && feat=""
    tlim="$(echo "$info" | grep -oE 'TimeLimit=[^ ]+' | cut -d= -f2)"
    memn="$(echo "$info" | grep -oE 'MinMemoryNode=[^ ]+' | cut -d= -f2)"
    ISMIP7_CHAIN_RES=(-N "${SLURM_JOB_NUM_NODES:-1}" -n "${SLURM_NTASKS:-1}"
                      --cpus-per-task="${SLURM_CPUS_PER_TASK:-1}"
                      --hint=nomultithread)
    # The link's per-node layout, when it asked for one. Slurm sets
    # SLURM_NTASKS_PER_NODE only for a job that requested it, --export=ALL
    # hands the link's value to the successor, and srun reads it there as its
    # own --ntasks-per-node, so a successor allocated without the same request
    # is asked for a layout it may not hold: the successor of a 2 x 32 link
    # (NOTS 1614035) died at srun with "More processors requested than
    # permitted".
    [ -n "${SLURM_NTASKS_PER_NODE:-}" ] &&
        ISMIP7_CHAIN_RES+=(--ntasks-per-node="$SLURM_NTASKS_PER_NODE")
    [ -n "${SLURM_JOB_NAME:-}" ] && ISMIP7_CHAIN_RES+=(-J "$SLURM_JOB_NAME")
    [ -n "${SLURM_JOB_PARTITION:-}" ] && ISMIP7_CHAIN_RES+=(-p "$SLURM_JOB_PARTITION")
    [ -n "$feat" ] && ISMIP7_CHAIN_RES+=(-C "$feat")
    [ -n "$tlim" ] && ISMIP7_CHAIN_RES+=(--time="$tlim")
    [ -n "$memn" ] && ISMIP7_CHAIN_RES+=(--mem="$memn")
    [ -n "${SLURM_JOB_ACCOUNT:-}" ] && ISMIP7_CHAIN_RES+=(-A "$SLURM_JOB_ACCOUNT")
    # A flag the site needs on every submission (a QOS, say) is not something
    # scontrol hands back either, so the successor is given it again.
    if [ -n "${ISMIP7_SBATCH_EXTRA:-}" ]; then
        # shellcheck disable=SC2206
        ISMIP7_CHAIN_RES+=(${ISMIP7_SBATCH_EXTRA})
    fi
    return 0
}

# The file that defines `module` for a job that arrives without it. A job's
# environment is the submitting shell's alone (submit.sh passes a bare
# --export=ALL, so slurmd adds no login environment), and a shell that never
# read the login scripts, such as the one `ssh host command` starts, holds no
# module system. The IU venv's python then cannot find libpython.
# /etc/profile is the first file a login shell reads; a site whose module
# system is set up elsewhere names that file.
ISMIP7_MODULE_INIT="${ISMIP7_MODULE_INIT:-/etc/profile}"

ismip7_load_modules() {
    if ! command -v module >/dev/null 2>&1; then
        # A site that loads no modules needs no module command either.
        [ -n "${ISMIP7_MODULES:-}${ISMIP7_MODULE_USE:-}" ] || return 0
        if [ -r "$ISMIP7_MODULE_INIT" ]; then
            # Login scripts expect a shell without -e or -u, and /etc/profile
            # sets a umask; the job's own options and umask come back after.
            local flags="$-" mask
            mask="$(umask)"
            set +eu
            # shellcheck disable=SC1090
            . "$ISMIP7_MODULE_INIT" >/dev/null 2>&1
            case "$flags" in *e*) set -e ;; *) set +e ;; esac
            case "$flags" in *u*) set -u ;; *) set +u ;; esac
            umask "$mask"
        fi
        if ! command -v module >/dev/null 2>&1; then
            echo "ERROR: this job has no 'module' command, even after reading" >&2
            echo "       ISMIP7_MODULE_INIT=$ISMIP7_MODULE_INIT, so the site's modules cannot load." >&2
            echo "       Submit from a login shell, or name the file that defines it." >&2
            exit 2
        fi
    fi
    module purge 2>/dev/null || true
    if [ -n "${ISMIP7_MODULE_USE:-}" ]; then
        # shellcheck disable=SC2086
        module use ${ISMIP7_MODULE_USE}
    fi
    if [ -n "${ISMIP7_MODULES:-}" ]; then
        # shellcheck disable=SC2086
        module load ${ISMIP7_MODULES}
    fi
}

# Container sites: check the image and the runtime, and put container_bin/
# (a `python` that runs in the image) first on PATH. Once only, because a chain
# successor inherits this PATH through --export=ALL.
ismip7_activate_container() {
    if [ ! -r "$ISMIP7_CONTAINER" ]; then
        echo "ERROR: ISMIP7_CONTAINER is unreadable: '$ISMIP7_CONTAINER'" >&2
        echo "       Set it in $_ISMIP7_SITE_FILE or sites/local.env." >&2
        exit 2
    fi
    if ! command -v "$ISMIP7_CONTAINER_RUNTIME" >/dev/null 2>&1; then
        echo "ERROR: '$ISMIP7_CONTAINER_RUNTIME' is not on PATH after loading" >&2
        echo "       ISMIP7_MODULES='${ISMIP7_MODULES:-}'. Add the module that provides it." >&2
        exit 2
    fi
    case ":$PATH:" in
        *":$_ISMIP7_BR_DIR/container_bin:"*) ;;
        *) PATH="$_ISMIP7_BR_DIR/container_bin:$PATH" ;;
    esac
    export PATH ISMIP7_CONTAINER ISMIP7_CONTAINER_RUNTIME ISMIP7_CONTAINER_ARGS
}

# What the image has to see besides $HOME, /tmp and the working directory,
# which the runtime binds itself: a checkout, forcing tree or kernel cache on
# a project or scratch filesystem is invisible inside it otherwise, and the
# runtime then quietly starts in $HOME instead.
ismip7_container_binds() {
    local d seen=" "
    ISMIP7_CONTAINER_BINDS=""
    for d in "$ISMIP7_REPO" "$ISMIP7_WORK" \
             "${ISMIP7_DATA_ROOT:-}" "${ISMIP7_OBS_DATA_ROOT:-}" \
             "${PYOP2_CACHE_DIR:-}" "${XDG_CACHE_HOME:-}" \
             "${ISMIP7_TIMING_JIT_CACHE:-}"; do
        [ -n "$d" ] && [ -d "$d" ] || continue
        case "$seen" in *" $d "*) continue ;; esac
        seen="$seen$d "
        ISMIP7_CONTAINER_BINDS="$ISMIP7_CONTAINER_BINDS --bind $d"
    done
    export ISMIP7_CONTAINER_BINDS
}

ismip7_activate() {
    ismip7_site_require
    # A job submitted with an --export list carries two variables that must go
    # no further. SLURM_GET_USER_ENV=1 has slurmd rebuild the login
    # environment when a job starts, and requeue and hold the job when that
    # fails; a chain resubmit's --export=ALL would hand it to every successor.
    # SLURM_EXPORT_ENV holds the list, which srun takes as its own --export,
    # so a variable the list names comes back in the job's steps after the
    # script unsets it (projection.sbatch unsets ISMIP7_RESTART for its
    # successor).
    unset SLURM_GET_USER_ENV SLURM_EXPORT_ENV
    ismip7_load_modules
    if [ -n "$ISMIP7_CONTAINER" ]; then
        ismip7_activate_container
    else
        if [ ! -r "$ISMIP7_FIREDRAKE" ]; then
            echo "ERROR: ISMIP7_FIREDRAKE is unreadable: '$ISMIP7_FIREDRAKE'" >&2
            echo "       Set it in $_ISMIP7_SITE_FILE, or export it before submitting." >&2
            exit 2
        fi
        # shellcheck disable=SC1090
        . "$ISMIP7_FIREDRAKE"
    fi
    # What the site's own modules or venv say about kernel caches, before the
    # per-job one below replaces it: ismip7_persistent_jit_cache puts the two
    # kernel caches back. IU's firedrake modulefile points both at a scratch
    # directory, and that is the warm cache every IU timing lane has been
    # measured with.
    _ISMIP7_SITE_PYOP2_CACHE_DIR="${PYOP2_CACHE_DIR:-}"
    _ISMIP7_SITE_TSFC_CACHE_DIR="${FIREDRAKE_TSFC_KERNEL_CACHE_DIR:-}"
    _ISMIP7_SITE_XDG_CACHE_HOME="${XDG_CACHE_HOME:-}"
    export OMP_NUM_THREADS=1          # one thread per rank; the solver is MPI-parallel
    export OPENBLAS_NUM_THREADS=1     # likewise for the BLAS under PETSc and numpy
    # Each rank compiles UFL kernels; a shared cache on a networked filesystem
    # corrupts under concurrent writes (seen locally: "undefined symbol:
    # wrap_form0_cell_integral"). Every job gets its own, keyed on the job id,
    # so a chain link killed mid compile cannot hand its successor a truncated
    # object through --export=ALL.
    export PYOP2_CACHE_DIR="${SCRATCH:-$HOME}/.pyop2_cache/${SLURM_JOB_ID:-manual}"
    # XDG_CACHE_HOME is loopy's knob and also matplotlib's, so pin the font
    # cache where it is already warm and move only loopy.
    export MPLCONFIGDIR="${MPLCONFIGDIR:-${_ISMIP7_SITE_XDG_CACHE_HOME:-$HOME/.cache}/matplotlib}"
    # loopy keeps its own persistent dict under XDG_CACHE_HOME (pytools/), not
    # under PYOP2_CACHE_DIR; two jobs compiling the same kernel seconds apart
    # raced on it (Rice 1559476, NoSuchEntryError in preprocess_program). It
    # sits beside the per-job kernel cache rather than inside it, so that
    # ismip7_persistent_jit_cache can retire an unused kernel cache while this
    # one stays in use for the whole job.
    export XDG_CACHE_HOME="${SCRATCH:-$HOME}/.pyop2_cache/xdg/${SLURM_JOB_ID:-manual}"
    # That dict is a sqlite file every rank of the job writes at once, and on
    # a networked filesystem the write can stall on a lock it never gets:
    # pytools retries SQLITE_BUSY without limit (persistent_dict._exec_sql_fn),
    # so the whole job then sits in the first kernel compile until its wall
    # time (IU 10569250 and 10569252, 16 and 64 ranks on one node, 260
    # retries at five seconds each and counting; three sibling jobs saw 20
    # to 60 retries and got through). The dict starts empty in every job and
    # only memoizes loopy's own preprocessing, which costs seconds per
    # kernel, so it buys nothing here: switch it off. The compiled kernels
    # still cache under PYOP2_CACHE_DIR, which is file based.
    export LOOPY_NO_CACHE=1
    mkdir -p "$PYOP2_CACHE_DIR" "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
    [ -n "$ISMIP7_CONTAINER" ] && ismip7_container_binds
    return 0
}

# For a job whose wall time is the measurement. ismip7_activate gives every job
# a private, empty kernel cache, which is right for a chain link and wrong for a
# timing lane: seconds_per_step is the whole transient loop over its steps with
# no warm-up excluded, so a cold cache times the compiler. Call this after
# ismip7_activate to go back to a cache that persists between jobs, in order:
# ISMIP7_TIMING_JIT_CACHE when the site or sites/local.env names one; else
# whatever the site's modules or venv had set before ismip7_activate replaced
# it (IU's modulefile names a scratch directory); else Firedrake's own default.
#
# loopy's persistent dict stays off, and XDG_CACHE_HOME stays where
# ismip7_activate put it, one per job. The race that killed Rice 1559476 is
# between lanes launched together, and `make timing-scout` submits one lane
# per mesh at once, so sharing that dict back is the exact condition that
# failed. Doing without it costs seconds of loopy preprocessing per lane; the
# kernel compile the shared cache protects is the expensive part, and it
# comes back below.
ismip7_persistent_jit_cache() {
    rmdir "$PYOP2_CACHE_DIR" 2>/dev/null || true
    if [ -n "${ISMIP7_TIMING_JIT_CACHE:-}" ]; then
        export PYOP2_CACHE_DIR="$ISMIP7_TIMING_JIT_CACHE/pyop2"
        export FIREDRAKE_TSFC_KERNEL_CACHE_DIR="$ISMIP7_TIMING_JIT_CACHE/tsfc"
        mkdir -p "$PYOP2_CACHE_DIR" "$FIREDRAKE_TSFC_KERNEL_CACHE_DIR"
    elif [ -n "${_ISMIP7_SITE_PYOP2_CACHE_DIR:-}" ]; then
        export PYOP2_CACHE_DIR="$_ISMIP7_SITE_PYOP2_CACHE_DIR"
        if [ -n "${_ISMIP7_SITE_TSFC_CACHE_DIR:-}" ]; then
            export FIREDRAKE_TSFC_KERNEL_CACHE_DIR="$_ISMIP7_SITE_TSFC_CACHE_DIR"
        fi
    else
        unset PYOP2_CACHE_DIR
    fi
    [ -n "$ISMIP7_CONTAINER" ] && ismip7_container_binds
    return 0
}

# Launch an MPI program on N ranks inside the running allocation:
#   ismip7_mpirun N python -u script.py [args...]
# Every job script starts its ranks through this, so how a site launches MPI
# is decided here once rather than in each script.
#
# A venv site starts the ranks with srun. A container site starts them with the
# image's own mpiexec, inside one `exec`: every job here is a single node, and
# that way the MPI in the image never has to agree with the host's Slurm about
# PMI. The image has python3 and no python, so that one word is translated.
#
# Open MPI's default MPI-IO component (ompio) writes a parallel HDF5 file to an
# NFS file system at a crawl: on NOTS, measured 24 September 2026 with 32 ranks
# on the 1000 m mesh, one 1.48 GB yearly output file took 1042 s to /scratch
# (VAST) under ompio and 52 s under romio321, on a file system that takes a
# single stream at 600 MB/s. At 45 minutes per model year that write was a
# third of the run. ISMIP7_MPI_IO names the component (default romio321; set
# it empty to leave the MPI's own default); an MPI other than Open MPI ignores
# the variable.
ismip7_mpirun() {
    local n="$1"; shift
    if [ -n "${ISMIP7_MPI_IO-romio321}" ]; then
        export OMPI_MCA_io="${ISMIP7_MPI_IO-romio321}"
    fi
    if [ -z "$ISMIP7_CONTAINER" ]; then
        srun -n "$n" "$@"
        return
    fi
    if [ "${1:-}" = python ]; then
        shift; set -- python3 "$@"
    fi
    # shellcheck disable=SC2086
    "$ISMIP7_CONTAINER_RUNTIME" exec ${ISMIP7_CONTAINER_BINDS:-} $ISMIP7_CONTAINER_ARGS \
        "$ISMIP7_CONTAINER" $ISMIP7_CONTAINER_MPIEXEC -n "$n" "$@"
}

# --- job size ------------------------------------------------------------
# Only the fields ismip7_site_require checks have to come from the site file.
# Everything below carries a working default, so a site file written from that
# list alone still composes a complete submission.
ISMIP7_TASKS="${ISMIP7_TASKS:-16}"
ISMIP7_MEM="${ISMIP7_MEM:-120G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-2-00:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-1-00:00:00}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"
# Extra sbatch flags for every submission from this site or this user, split
# on spaces: a QOS the partition insists on, or --mail-type/--mail-user.
ISMIP7_SBATCH_EXTRA="${ISMIP7_SBATCH_EXTRA:-}"

# A site that needs one number sets ISMIP7_TASKS/ISMIP7_MEM and both kinds take
# it. A site with measured per-kind values sets the pair, so a forward runs at
# a rank count measured for forwards and never inherits the inversion's size
# unvalidated (sites/rice_nots.sh carries its measurements).
ISMIP7_TASKS_INV="${ISMIP7_TASKS_INV:-$ISMIP7_TASKS}"
ISMIP7_MEM_INV="${ISMIP7_MEM_INV:-$ISMIP7_MEM}"
ISMIP7_TASKS_FWD="${ISMIP7_TASKS_FWD:-$ISMIP7_TASKS}"
ISMIP7_MEM_FWD="${ISMIP7_MEM_FWD:-$ISMIP7_MEM}"

# The node feature splits the same way. An inversion needs the partition with
# the memory, and the rest of the kinds go wherever the measured timings and
# the -march=native build came from.
ISMIP7_CONSTRAINT_INV="${ISMIP7_CONSTRAINT_INV:-${ISMIP7_CONSTRAINT:-}}"
ISMIP7_CONSTRAINT_FWD="${ISMIP7_CONSTRAINT_FWD:-${ISMIP7_CONSTRAINT:-}}"

# `submit.sh script`, which is how the timing campaign submits, takes its node
# feature from ISMIP7_CONSTRAINT_TIMING and refuses a job that cannot fit one
# node of that kind: ISMIP7_CORES_PER_NODE physical cores, ISMIP7_MEM_PER_NODE
# of requestable memory (240G, 187000M, ...). Empty means no limit is known,
# and nothing is refused. The campaign's lanes are single-node by design, so a
# lane this site's nodes cannot hold is recorded as not runnable here rather
# than sized differently and compared with the other sites' as if it were not.
ISMIP7_CONSTRAINT_TIMING="${ISMIP7_CONSTRAINT_TIMING:-$ISMIP7_CONSTRAINT_FWD}"
ISMIP7_CORES_PER_NODE="${ISMIP7_CORES_PER_NODE:-}"
ISMIP7_MEM_PER_NODE="${ISMIP7_MEM_PER_NODE:-}"

# The scheduler half of the job banner. site_env.sh defines
# ismip7_banner_model for the runners that carry a model configuration; a job
# that sourced this file alone prints the lines below and nothing else.
ismip7_banner() {
    echo "=== $(date -Is)  job ${SLURM_JOB_ID:-none} on $(hostname) ==="
    echo "    site    $ISMIP7_SITE_NAME ($(basename "$_ISMIP7_SITE_FILE"))"
    echo "    sched   partition ${SLURM_JOB_PARTITION:-?}  nodes ${SLURM_JOB_NUM_NODES:-?}" \
         " ntasks ${SLURM_NTASKS:-?}  mem ${SLURM_MEM_PER_NODE:-?}M" \
         "${ISMIP7_ACCOUNT:+ account $ISMIP7_ACCOUNT}"
    echo "    node    $(scontrol show node "$(hostname -s)" 2>/dev/null \
                      | tr ' ' '\n' | grep -E '^(CPUTot|RealMemory|ActiveFeatures)=' \
                      | paste -sd' ' || echo 'unknown')"
    echo "    repo    $ISMIP7_REPO"
    if declare -F ismip7_banner_model >/dev/null; then
        ismip7_banner_model
    fi
}
