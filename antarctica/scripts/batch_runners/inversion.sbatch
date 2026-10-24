#!/bin/bash
# ISMIP7 icepack2 inversion on any cluster, with a self-resuming chain.
#
# Submit through submit_inversions.sh, which sets the mesh, sampling and
# output name; or through submit.sh, which takes the resources from this
# cluster's site file:
#   antarctica/scripts/batch_runners/submit.sh inversion \
#     ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_BUFFER_M=0 \
#     ISMIP7_MESH=... ISMIP7_RASTER_SAMPLE=vertex_front ISMIP7_MAP_OUT=...
#
# CHAIN. A 2 km inversion can outlast a wall limit. The inversion checkpoints
# ISMIP7_MAP_OUT every ISMIP7_CHECKPOINT_EVERY_IT iterates (default 1) and can
# warm-start from it, so this script queues its own successor FIRST
# (--dependency=afterany, so it runs whether this job converged, timed out or
# lost its node), and each link exits at once
# if the "<map>.done" marker exists. Depth is capped by ISMIP7_CHAIN_MAX. The
# successor copies this job's partition, constraint, memory, time and task
# layout from scontrol, so a chain never changes rank count mid-run. The
# inversion reads a warm start by point location (interpolate_with_fill), so
# one written on another rank count or on another mesh loads correctly too.
#
# Defaults write inversion_icepack2_rc_n3_dg0_logvelnet_<ISMIP7_LC>.h5 (1000
# on the production mesh site_env.sh names) under the settings the 2500 m
# result came from: the sigma-normalised velocity chi^2
# plus ISSM's logarithmic misfit, the pointwise dH/dt term, and the integrated
# net mass-balance term (off by default in the repo) that takes grounding-line
# discharge from 1628 to 2211 Gt/yr. The friction law is the one deliberate
# difference: site_env.sh defaults ISMIP7_FRICTION to regularized_coulomb.

# Resources (partition, tasks, memory, wall time, account, constraint) come
# from submit.sh and this cluster's sites/<name>.sh, never from directives
# here: a directive cannot read a site file, and a header that disagrees
# with the command line is rejected outright by Slurm.
#SBATCH -o logs/ismip7_inv_%j.out
#SBATCH -e logs/ismip7_inv_%j.err
#SBATCH --mail-type=END,FAIL

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p logs

# Resources are the wrapper's job (submit.sh reads sites/<cluster>.sh). A job
# that arrives from the scheduler with one task was submitted by hand without
# them, and at 2 km that dies on memory hours into the queue rather than at
# submission. Off a scheduler there is no allocation to judge, so no check.
if [ -n "${SLURM_JOB_ID:-}" ] && [ "${SLURM_NTASKS:-1}" = 1 ]; then
    echo "ERROR: this job has SLURM_NTASKS=1, so no resources were requested." >&2
    echo "       Submit with: antarctica/scripts/batch_runners/submit.sh inversion [KEY=VALUE ...]" >&2
    exit 2
fi

. antarctica/scripts/batch_runners/site_env.sh
ismip7_activate
ismip7_banner

export ISMIP7_MISFIT_NORM="${ISMIP7_MISFIT_NORM:-sigma}"
export ISMIP7_SIGMA_U_FLOOR="${ISMIP7_SIGMA_U_FLOOR:-3}"
export ISMIP7_LOG_VEL_WEIGHT="${ISMIP7_LOG_VEL_WEIGHT:-auto}"
# The dH/dt term is OFF by default (Andrew, 26 Sep 2026): the inversion's
# prognostic step carries no forcing, so the term asks the velocity to absorb
# whatever SMB and melt are wrong by.
export ISMIP7_DHDT_WEIGHT="${ISMIP7_DHDT_WEIGHT:-0}"
export ISMIP7_DHDT_NET_SIGMA="${ISMIP7_DHDT_NET_SIGMA:-10}"
export ISMIP7_GAMMA_THETA="${ISMIP7_GAMMA_THETA:-1e5}"
export ISMIP7_GAMMA_PHI="${ISMIP7_GAMMA_PHI:-1e5}"
export ISMIP7_MAXITER="${ISMIP7_MAXITER:-200}"
# A chain link must continue the SAME objective from the accepted iterate
# (icepack2_tools.handoff); a link that would change it refuses.
export ISMIP7_WARM_START_STRICT="${ISMIP7_WARM_START_STRICT:-1}"
export ISMIP7_RASTER_SAMPLE="${ISMIP7_RASTER_SAMPLE:-vertex_front}"
# Name the output explicitly: without it every run, including a smoke test,
# writes the production filename and can silently replace a converged map.
# ISMIP7_MAP_DEFAULT is named from the law, so ISMIP7_FRICTION=budd writes a
# _budd file rather than resuming from and then overwriting the RC one. It is
# the same path projection.sbatch loads.
export ISMIP7_MAP_OUT="${ISMIP7_MAP_OUT:-$ISMIP7_MAP_DEFAULT}"
DONE="${ISMIP7_MAP_OUT}.done"
# Off a scheduler (sites/local.sh, debugging on a workstation) there is no
# job id, no allocation and nothing to chain to.
LOG="logs/ismip7_inv_${SLURM_JOB_ID:-manual}.out"

echo "    map out  $ISMIP7_MAP_OUT"
echo "    sampling $ISMIP7_RASTER_SAMPLE  net sigma $ISMIP7_DHDT_NET_SIGMA  log-vel $ISMIP7_LOG_VEL_WEIGHT  maxiter $ISMIP7_MAXITER"
# The driver prints the modes it resolves ("Inversion linear solver:"); an
# unset knob here means the driver's default, which this file does not repeat.
echo "    solver   ${ISMIP7_INVERSION_LINEAR_SOLVER:-default}  ramp ${ISMIP7_DIAGNOSTIC_LINEAR_SOLVER:-default}"

if [ -f "$DONE" ]; then
    echo "chain: $DONE exists, inversion already finished; nothing to do"
    exit 0
fi

# --- the iteration budget of the whole chain --------------------------------
# ISMIP7_MAXITER is one process's budget, and a link that resumes the
# checkpoint counts from zero again, so a chain of one-hour links spends
# ISMIP7_MAXITER per link until its depth runs out. ISMIP7_CHAIN_ITER_BUDGET
# caps the accepted iterations across the chain instead: each link adds the
# count its predecessor left in the checkpoint (objective_iteration, that
# process's count at the saved iterate) to the sidecar <map>.iters, keyed by
# the checkpoint's modification time and size so a predecessor that saved
# nothing is not counted twice, and runs the remainder as its ISMIP7_MAXITER.
# A chain whose budget is spent marks the MAP finished and runs the follow-on
# command. The staged initialization routine gives every stage its budget
# this way.
if [ -n "${ISMIP7_CHAIN_ITER_BUDGET:-}" ]; then
    remaining=$(python - "$ISMIP7_MAP_OUT" "$ISMIP7_CHAIN_ITER_BUDGET" <<'PY'
import os
import sys

map_fn, budget = sys.argv[1], int(sys.argv[2])
iters_fn = map_fn + ".iters"
done, stamp = 0, ""
if os.path.exists(iters_fn):
    parts = open(iters_fn).read().split()
    done = int(parts[0])
    stamp = parts[1] if len(parts) > 1 else ""
if os.path.exists(map_fn):
    st = os.stat(map_fn)
    this = f"{st.st_mtime_ns}:{st.st_size}"
    if this != stamp:
        import h5py
        with h5py.File(map_fn, "r") as h:
            done += int(round(float(h["/"].attrs.get("objective_iteration", 0))))
        with open(iters_fn, "w") as fh:
            fh.write(f"{done} {this}\n")
print(budget - done)
PY
)
    if [ -z "$remaining" ]; then
        echo "chain: WARNING could not account the iteration budget;" \
             "this link runs ISMIP7_MAXITER=$ISMIP7_MAXITER"
    elif [ "$remaining" -le 0 ]; then
        touch "$DONE"
        echo "chain: the iteration budget ISMIP7_CHAIN_ITER_BUDGET=$ISMIP7_CHAIN_ITER_BUDGET" \
             "is spent; marker $DONE"
        ismip7_chain_then_cmd
        exit 0
    else
        export ISMIP7_MAXITER="$remaining"
        echo "chain: iteration budget $ISMIP7_CHAIN_ITER_BUDGET," \
             "$remaining left for this link (ISMIP7_MAXITER)"
    fi
fi

# --- queue the successor before doing anything slow -----------------------
DEPTH="${ISMIP7_CHAIN_DEPTH:-0}"
MAXD="${ISMIP7_CHAIN_MAX:-4}"
if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo "chain: no scheduler, so no successor to queue"
elif [ "$DEPTH" -lt "$MAXD" ]; then
    ismip7_chain_resources
    # The depth goes in sbatch's own environment, which a bare --export=ALL
    # carries. The list form --export=ALL,ISMIP7_CHAIN_DEPTH=N has Slurm
    # rebuild the login environment when the successor starts, and IU Quartz
    # held successor 10555984 on 21 September 2026 when that failed.
    NEXT=$(ISMIP7_CHAIN_DEPTH=$((DEPTH + 1)) sbatch --parsable \
        --dependency=afterany:"$SLURM_JOB_ID" \
        "${ISMIP7_CHAIN_RES[@]}" \
        --export=ALL \
        antarctica/scripts/batch_runners/inversion.sbatch) \
      && echo "chain: successor $NEXT queued (depth $((DEPTH + 1)) of $MAXD)" \
      || echo "chain: WARNING could not queue a successor"
else
    echo "chain: depth $DEPTH reached ISMIP7_CHAIN_MAX=$MAXD, no successor"
fi

# --- resume if a checkpoint exists ------------------------------------------
# The driver takes theta and phi from the checkpoint and, under
# ISMIP7_LOG_VEL_WEIGHT=auto, the log-velocity weight it records, so every
# link minimises the objective the first one set (issue 68).
if [ -f "$ISMIP7_MAP_OUT" ]; then
    export ISMIP7_WARM_START="$ISMIP7_MAP_OUT"
    # Resuming this chain's own checkpoint is always strict. A first link may
    # set ISMIP7_WARM_START_STRICT=0 to change the objective of the MAP it
    # starts from on purpose (another prior), and the successor it queued
    # inherits that 0 through --export=ALL; here it is held to the one
    # objective the chain has been minimising.
    export ISMIP7_WARM_START_STRICT=1
    # A first link's fluidity and state from another MAP
    # (ISMIP7_WARM_START_FLUIDITY, ISMIP7_WARM_START_STATE) seed the chain
    # once; a resumed link takes the fluidity, geometry and state of its own
    # checkpoint.
    unset ISMIP7_WARM_START_FLUIDITY ISMIP7_WARM_START_STATE
    echo "chain: warm-starting theta/phi from the periodic checkpoint $ISMIP7_MAP_OUT (strict handoff)"
fi

# Read the warm start from node-local disk. Every rank reading the MAP over
# NFS through ROMIO (OMPI_MCA_io=romio321, chosen for the write speed) sits in
# NFS lock round trips: on NOTS /scratch, 30 Sep 2026, 32 ranks took 45+ min
# to load a 580 MB 2 km MAP, ~12 KB per rank per 20 s, so a one-hour link
# never reached its first iteration. A local copy takes a second. Same
# basename, so the provenance the driver records is unchanged; one node only,
# since other nodes' ranks cannot see this node's disk.
_WS_LOCAL=""
if [ -n "${ISMIP7_WARM_START:-}" ] && [ -f "$ISMIP7_WARM_START" ] \
        && [ "${SLURM_JOB_NUM_NODES:-1}" = 1 ]; then
    _WS_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ismip7_ws.XXXXXX")" \
      && cp "$ISMIP7_WARM_START" "$_WS_DIR/" \
      && _WS_LOCAL="$_WS_DIR/$(basename "$ISMIP7_WARM_START")" \
      && echo "warm start: local copy $_WS_LOCAL" \
      && export ISMIP7_WARM_START="$_WS_LOCAL"
fi

ismip7_mpirun "${SLURM_NTASKS:-1}" python -u antarctica/scripts/inversion_icepack2.py
rc=$?
echo "=== $(date -Is)  inversion exit ${rc} ==="
[ -n "$_WS_LOCAL" ] && rm -rf "$(dirname "$_WS_LOCAL")"

# Finished for good whenever L-BFGS-B returned AND the MAP reached disk:
# ISMIP7_MAXITER is the budget of a run, so ITERATIONS REACHED LIMIT is the
# normal way a production inversion ends. The successor exists for the other
# case, a job the wall clock killed before the MAP was written (no "Saved MAP:"
# in the log, which the driver prints only once the checkpoint write has
# returned); it warm-starts from the periodic checkpoint and finishes the
# budget. (The previous rule kept chaining after a completed budget: job
# 1339328 finished its 200 iterations and 1339342 started re-inverting the
# finished MAP, Sep 13 2026.) The driver writes this marker itself as soon as
# it has saved the MAP, so the rule below is the fallback for a driver that
# could not write it; touch is idempotent. A driver whose optimizer stopped
# right after failed trial evaluations prints "MAP not final:" before saving
# and writes no marker: scipy reports such a stop as convergence at a point
# whose gradient is not small (issue 153, job 11461566), so the successor
# resumes from the MAP instead.
if grep -q 'MAP not final:' "$LOG"; then
    echo "chain: the optimizer stopped on failed trials; no marker, the successor resumes from the MAP"
elif grep -q 'Saved MAP:' "$LOG"; then
    touch "$DONE"
    echo "chain: finished ($(grep -o 'Saved MAP:.*' "$LOG" | tail -1)); marker $DONE"
fi
# The link that finished the MAP hands the chain on (ISMIP7_CHAIN_THEN_CMD).
[ -f "$DONE" ] && ismip7_chain_then_cmd
exit $rc
