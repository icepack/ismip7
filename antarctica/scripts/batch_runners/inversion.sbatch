#!/bin/bash
# ISMIP7 icepack2 inversion on any cluster, with a self-resuming chain.
#
# Submit through submit_inversions.sh, which sets the mesh, sampling and
# output name; or through submit.sh, which takes the resources from this
# cluster's site file:
#   antarctica/scripts/batch_runners/submit.sh inversion \
#     ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_BUFFER_M=0 \
#     ISMIP7_MESH=... ISMIP7_RASTER_SAMPLE=vertex ISMIP7_MAP_OUT=...
#
# CHAIN. A 2 km inversion can outlast a wall limit. The inversion checkpoints
# ISMIP7_MAP_OUT every 20 iterates and can warm-start from it, so this script
# queues its own successor FIRST (--dependency=afterany, so it runs whether
# this job converged, timed out or lost its node), and each link exits at once
# if the "<map>.done" marker exists. Depth is capped by ISMIP7_CHAIN_MAX. The
# successor copies this job's partition, constraint, memory, time and task
# layout from scontrol, so a chain never changes rank count mid-run -- the
# inversion now refuses a warm start whose dof ordering differs anyway.
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
export ISMIP7_DHDT_WEIGHT="${ISMIP7_DHDT_WEIGHT:-1}"
export ISMIP7_DHDT_NET_SIGMA="${ISMIP7_DHDT_NET_SIGMA:-10}"
export ISMIP7_GAMMA_THETA="${ISMIP7_GAMMA_THETA:-1e5}"
export ISMIP7_GAMMA_PHI="${ISMIP7_GAMMA_PHI:-1e5}"
export ISMIP7_MAXITER="${ISMIP7_MAXITER:-200}"
export ISMIP7_RASTER_SAMPLE="${ISMIP7_RASTER_SAMPLE:-vertex}"
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

if [ -f "$DONE" ]; then
    echo "chain: $DONE exists, inversion already finished; nothing to do"
    exit 0
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
    echo "chain: warm-starting theta/phi from the periodic checkpoint $ISMIP7_MAP_OUT"
fi

ismip7_mpirun "${SLURM_NTASKS:-1}" python -u antarctica/scripts/inversion_icepack2.py
rc=$?
echo "=== $(date -Is)  inversion exit ${rc} ==="

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
# could not write it; touch is idempotent.
if grep -q 'Saved MAP:' "$LOG"; then
    touch "$DONE"
    echo "chain: finished ($(grep -o 'Saved MAP:.*' "$LOG" | tail -1)); marker $DONE"
fi
exit $rc
