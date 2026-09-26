#!/bin/bash
# ISMIP7 forward run (control, historical or projection) on any cluster.
# SELF-CHAINING: a 285-year projection is about four days on the production
# mesh and will not fit an ordinary queue, so this script resubmits itself
# until the run reaches its driver's own end year, resuming from its own
# newest checkpoint each time. The driver owns that year and this script does
# not default ISMIP7_T_END; the chain reads the end back out of the driver's
# "Time-stepping: <start>-><end>" line. Setting ISMIP7_T_END overrides every
# experiment, so leave it unset unless you mean to.
#
#   antarctica/scripts/batch_runners/submit.sh projection \
#          ISMIP7_EXPERIMENT=control
#
# ISMIP7_EXPERIMENT selects the driver:
#   control            antarctica/scripts/control/run.py
#   ssp126_cesm_waccm  antarctica/scripts/projections/ssp126_cesm_waccm.py
#   ssp126_mri_esm2, ssp370_cesm_waccm, ssp370_mri_esm2, ssp585_cesm_waccm,
#   ssp585_mri_esm2    the other core projections, same layout
#   ocx                antarctica/scripts/projections/ocx.py
#   hist_cesm_waccm    antarctica/scripts/historical/cesm_waccm.py
#   hist_mri_esm2      antarctica/scripts/historical/mri_esm2.py
#
# ISMIP7_CHAIN_THEN (space-separated experiments from the same list) queues
# follow-on runs once THIS run has reached its end year, and only then: a
# control or projection branches from the historical's endpoint, so it must
# not start before the historical has finished, and no earlier point in the
# chain knows that it has. Each follow-on starts a chain of its own with this
# run's environment and allocation, less ISMIP7_CHAIN_THEN itself and the
# period (ISMIP7_T_START, ISMIP7_T_END), which each follow-on's driver owns:
#
#   ISMIP7_EXPERIMENT=hist_mri_esm2 ISMIP7_ESM=MRI-ESM2-0 \
#       ISMIP7_CHAIN_THEN="control ssp585_mri_esm2" submit.sh projection
#
# Sizing at the production step, 0.025 yr (issue 20). Rice measured 22.5 min
# per simulated year at 1000 m / 10 km on 32 Cascade Lake ranks under
# scpc_gamg, 62 simulated years per 24 h job. The timing matrix on Quartz
# (TIMING_MATRIX_QUARTZ_SCPC_GAMG.md) measured 30 s per 0.05 yr step on 64
# ranks, 10.0 min per simulated year; at 0.025 that is about 20 min and 21
# core-hours per simulated year if a step costs the same, which is not yet
# measured. That is the transient loop alone, under the matrix's strict
# contract: setup, forcing updates and output come on top. 24 h of wall time
# then buys about 70 simulated years, so a 285-year projection is about five
# chained jobs and a historical from 2003 fits in one. At 0.05 the matrix's
# year cost 15.7 min on 32 ranks and 32.3 on 16; scpc_mumps does not get
# below 23 min on any rank count (TIMING_MATRIX_QUARTZ_SCPC_MUMPS.md).
#
# The older measurement, on Cascade Lake under full_mumps: 155 s per 0.1 yr
# step at 2500 m / 64 km on 12 ranks, ~26 min per simulated year, ~5
# core-hours per year, about six chained 24 h jobs for 285 years.
#
# NOTE (measured 2026-09-09): deepsC has NO wall-time limit, so on deepsC the
# 24 h below is a self-imposed courtesy on a shared condo and the chaining is
# optional. The chaining is what makes this runnable on `commons` (1-day limit)
# and `long` (3-day limit), and it is also the safety net that turns a node
# failure into a resume rather than a lost week. Keep it on.

# Resources (partition, tasks, memory, wall time, account, constraint) come
# from submit.sh and this cluster's sites/<name>.sh, never from directives
# here: a directive cannot read a site file, and a header that disagrees
# with the command line is rejected outright by Slurm.
#SBATCH -o logs/ismip7_fwd_%j.out
#SBATCH -e logs/ismip7_fwd_%j.err
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
    echo "       Submit with: antarctica/scripts/batch_runners/submit.sh projection [KEY=VALUE ...]" >&2
    exit 2
fi

. antarctica/scripts/batch_runners/site_env.sh
ismip7_activate
ismip7_banner

# Off a scheduler (sites/local.sh, debugging a job script on a workstation)
# there is no job id, no allocation and nothing to chain to.
JOB="${SLURM_JOB_ID:-manual}"
LOG="logs/ismip7_fwd_${JOB}.out"

# The driver of an experiment name, or nothing for a name it does not know.
driver_for() {
    case "$1" in
        control)           echo antarctica/scripts/control/run.py ;;
        ssp126_cesm_waccm) echo antarctica/scripts/projections/ssp126_cesm_waccm.py ;;
        ssp126_mri_esm2)   echo antarctica/scripts/projections/ssp126_mri_esm2.py ;;
        ssp370_cesm_waccm) echo antarctica/scripts/projections/ssp370_cesm_waccm.py ;;
        ssp370_mri_esm2)   echo antarctica/scripts/projections/ssp370_mri_esm2.py ;;
        ssp585_cesm_waccm) echo antarctica/scripts/projections/ssp585_cesm_waccm.py ;;
        ssp585_mri_esm2)   echo antarctica/scripts/projections/ssp585_mri_esm2.py ;;
        ocx)               echo antarctica/scripts/projections/ocx.py ;;
        hist_cesm_waccm)   echo antarctica/scripts/historical/cesm_waccm.py ;;
        hist_mri_esm2)     echo antarctica/scripts/historical/mri_esm2.py ;;
    esac
}
EXPERIMENT="${ISMIP7_EXPERIMENT:-control}"
DRIVER=$(driver_for "$EXPERIMENT")
if [ -z "$DRIVER" ]; then
    echo "ERROR: unknown ISMIP7_EXPERIMENT='$EXPERIMENT'" >&2; exit 2
fi
# A misspelt follow-on is found now, not days from now when this run ends.
for next in ${ISMIP7_CHAIN_THEN:-}; do
    if [ -z "$(driver_for "$next")" ]; then
        echo "ERROR: unknown experiment '$next' in ISMIP7_CHAIN_THEN" >&2; exit 2
    fi
done

# --- forward configuration ---------------------------------------------
# Point the forward at the MAP the inversion chain wrote. Without this it would
# derive the name through icepack2_tools.naming.map_basename, which has no
# `logvelnet` element, and abort on a missing file after the queue wait.
export ISMIP7_INVERSION="${ISMIP7_INVERSION:-$ISMIP7_MAP_DEFAULT}"
# The production solver: the stress fields eliminated exactly cell by cell and
# GAMG on the condensed velocity system. On the production mesh it takes the
# Newton path scpc_mumps takes (129 iterations over the matrix's ten steps on
# 64 ranks, outflux equal to 9e-13) in 43 % of the time. Named here and
# nowhere more general: site_env.sh and solverconfig's own default are read by
# the inversion too, whose linear solve is the full mixed-Jacobian MUMPS by
# construction and which stamps the solver it finds in the environment on its
# MAP. ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps restores the reference.
export ISMIP7_DIAGNOSTIC_LINEAR_SOLVER="${ISMIP7_DIAGNOSTIC_LINEAR_SOLVER:-scpc_gamg}"
# The production step (issue 20), equal to runconfig.DT_DEFAULT, which the
# drivers fall back to; a test holds the two equal. The matrix's rule
# (0.125 yr x lc / 2500) gives 0.05 at 1000 m, where a 1 km control from a
# transferred 2 km Budd MAP diverged at the Lambert confluence; at 0.025 the
# same control ran five years, and Rice's 1 km historicals ran 78 and 66 years
# with no rescue step. The stable step does not grow with the mesh (0.25 yr
# ran away at 2500 m and at 5000 m), so a coarser ISMIP7_LC wants at most
# 0.125. A chain inherits the step with the rest of the environment.
export ISMIP7_DT="${ISMIP7_DT:-0.025}"
# ISMIP7_T_END is deliberately NOT defaulted here. Every experiment this
# runner offers has its own end (historical 2015, the ssp370 pair 2101, the
# other projections and the control 2301, OCX 2026), and the driver owns it;
# one value exported over all of them ran ssp370 past its series and the
# historical a year past its handoff. An explicit ISMIP7_T_END still wins,
# and the chain reads the end the run actually used out of its own log.
export ISMIP7_CHECKPOINT_EVERY_YR="${ISMIP7_CHECKPOINT_EVERY_YR:-5}"
export ISMIP7_KEEP_CHECKPOINTS="${ISMIP7_KEEP_CHECKPOINTS:-3}"
# The protocol's control holds calving at end-of-2014 conditions, and
# ISMIP7_FIXED_FRONT holds the front here. The apparent mass-balance reference
# is this repository's choice: the protocol leaves initial conditions to each
# group, and whether the production runs keep the reference is a group
# decision (issue #104). Projections that should calve freely must override
# ISMIP7_CALVING and unset ISMIP7_FIXED_FRONT.
#
# These two default only when UNSET (${VAR-default}, no colon). An explicitly
# empty value is one of the documented off spellings, and `submit.sh
# projection VAR=` delivers one: VAR= goes into sbatch's environment, and ALL
# carries it into the job and its successors. :- would substitute the default
# and make that off switch unreachable.
#   ISMIP7_APPARENT_MB  off: 0, off, none, empty. On: 1 or balance (balanced
#                       control), div (divergence only).
#   ISMIP7_AUTO_RESUME  off: 0, empty. On: any other integer.
#
# Both values go to the drivers exactly as given. runconfig owns each value
# set and rejects anything else at startup, so normalising here would only
# hide a typo: an ISMIP7_AUTO_RESUME=off meant to disable resume would become
# a 1 and silently continue an old trajectory.
export ISMIP7_APPARENT_MB="${ISMIP7_APPARENT_MB-1}"
export ISMIP7_FIXED_FRONT="${ISMIP7_FIXED_FRONT:-1}"
# Resume from this experiment's own newest checkpoint. This is what makes the
# chain below safe to run unattended.
export ISMIP7_AUTO_RESUME="${ISMIP7_AUTO_RESUME-1}"
# The submission's own output: the yearly ISMIP7 fields and scalars
# (icepack2_tools/ismip7_output.py), which write_ismip7_output.py regrids and
# the compliance checker reads. Every experiment this runner offers is a core
# experiment, and off it produces a timeseries CSV and checkpoints alone, so
# a forgotten flag is found after the days of wall time it takes a projection
# to finish, and the whole run has to be repeated because the per-year fields
# are gone. Defaulted on for that reason; `ISMIP7_OUTPUT=0` turns it off for a
# pipeline exercise, and the empty spelling `ISMIP7_OUTPUT=` does too (the
# accessor's documented off values), which is why this is ${VAR-default} and
# not ${VAR:-default}. The timing campaign has its own runner and is unaffected.
export ISMIP7_OUTPUT="${ISMIP7_OUTPUT-1}"

echo "    experiment $EXPERIMENT -> $DRIVER"
echo "    map in  $ISMIP7_INVERSION"
[ -f "$ISMIP7_INVERSION" ] \
    || echo "    WARNING: that MAP does not exist; run inversion.sbatch first" \
            "or set ISMIP7_INVERSION"
echo "    solver=$ISMIP7_DIAGNOSTIC_LINEAR_SOLVER dt=$ISMIP7_DT" \
     "t_end=${ISMIP7_T_END:-driver default (2301 projections, 2101 ssp370, 2015 historical, 2026 OCX)}" \
     "auto_resume=$ISMIP7_AUTO_RESUME" \
     "apparent_mb=$ISMIP7_APPARENT_MB fixed_front=$ISMIP7_FIXED_FRONT"

# --- run, leaving time to resubmit cleanly ------------------------------
# Stop the model a little before the wall clock so the final checkpoint is
# written and the next job in the chain starts from a complete state. Slurm
# gives no warning at the limit, it just kills the step, so the budget has to
# be handed to the driver: ISMIP7_WALL_STOP_MIN is minutes of run time.
#
# Derived from THIS job's own TimeLimit, and deliberately not exported: the
# resubmit below uses --export=ALL, so an exported value would pin every later
# link in the chain to the first job's limit even when it lands on a partition
# with a shorter one.
wall_stop_min="${ISMIP7_WALL_STOP_MIN:-}"
tl=$(scontrol show job "$JOB" 2>/dev/null \
     | grep -oE 'TimeLimit=[^ ]+' | cut -d= -f2)
if [ -z "$wall_stop_min" ]; then
    wall_stop_min=$(python - "$tl" <<'PY'
import sys
tl = sys.argv[1]
# Minutes held back for the final checkpoint write and the resubmit below.
MARGIN_MIN = 25.0
# Slurm renders a limit as [days-]HH:MM:SS or, under one hour, MM:SS.
days, _, clock = tl.rpartition("-")
parts = clock.split(":")
try:
    h, m, s = ([0.0] * (3 - len(parts)) + [float(x) for x in parts])
    total = (float(days) if days else 0.0) * 1440.0 + h * 60.0 + m + s / 60.0
except ValueError:
    print(0)                      # UNLIMITED or unparseable: no budget
else:
    print(f"{max(total - MARGIN_MIN, 0.0):.1f}")
PY
)
fi
echo "    wall budget ${wall_stop_min} min (limit ${tl:-?})"
ISMIP7_WALL_STOP_MIN="$wall_stop_min" ismip7_mpirun "${SLURM_NTASKS:-1}" python -u "$DRIVER"
rc=$?
echo "=== $(date -Is)  forward exit ${rc} ==="

# --- chain --------------------------------------------------------------
# Two separate questions, in this order: what happened to this run, and should
# a successor be submitted. The verdict comes first and owns the exit code, so
# a stall is reported as a failure whether or not chaining was going to happen
# at all. A non-zero exit is a real failure and must not be retried blindly:
# the July grounding-line blow-up looked exactly like a run that "just needed
# more time", and chaining through it would have burned days.
CHAIN="${ISMIP7_CHAIN:-1}"
if [ "$rc" -ne 0 ]; then
    echo "chain: stopping, driver exited ${rc}"
    exit $rc
fi

# Ask the driver's own final checkpoint, which is both the state the next job
# would resume from and the record of why this one stopped: `t_yr` is the year
# actually reached and `stalled` is the driver's verdict that the solver gave
# up. It must be THIS job's checkpoint: the results directory is flat and
# shared, and the experiment name is built inside the driver, so neither the
# newest file nor a name rebuilt in the shell identifies this run. The driver
# prints "Saved: <final.h5>" on every exit, including an early stop, so take
# the path from our own Slurm log.
read -r LAST_YEAR STALLED START_YEAR <<EOF
$(NM_LOG="$LOG" python - <<'PY'
import os, re
log = os.environ["NM_LOG"]
chk_fn, start_year = None, ""
with open(log, errors="replace") as f:
    for line in f:
        if line.startswith("Saved: ") and line.rstrip().endswith("_final.h5"):
            chk_fn = line[len("Saved: "):].strip()
        else:
            m = re.search(r"Restart: evolved geometry .*t_yr=([-\d.eE+]+)", line)
            if m:
                start_year = m.group(1)
last, stalled = "", 0
if chk_fn and os.path.exists(chk_fn):
    import h5py
    with h5py.File(chk_fn, "r") as h:
        t = h["/"].attrs.get("t_yr")
        last = "" if t is None else repr(float(t))
        stalled = int(h["/"].attrs.get("stalled", 0))
print(last or "-", stalled, start_year or "-")
PY
)
EOF
if [ "$STALLED" = "1" ]; then
    echo "chain: the solver stalled (rescue ladder + subcycles exhausted) at ${LAST_YEAR}."
    echo "chain: NOT resubmitting. A stall is not a run that needs more time;"
    echo "chain: resuming would re-attempt the same years and stall again."
    exit 1
fi
if [ "$LAST_YEAR" = "-" ] || [ -z "$LAST_YEAR" ]; then
    echo "chain: could not read this job's final checkpoint from $LOG;" \
         "not resubmitting"
    exit 0
fi
# The end the run actually used, from the line run_simulation prints for
# every driver ("Time-stepping: <start>-><end>, dt=..."). Reading it back
# rather than re-declaring the per-experiment defaults here keeps the driver
# the single owner of its period.
T_END=$(sed -n 's/^Time-stepping: [^-]*->\([0-9.]*\),.*/\1/p' "$LOG" | tail -1)
if [ -z "$T_END" ]; then
    echo "chain: no 'Time-stepping:' line in $LOG, so the run's end year" \
         "is unknown; not resubmitting"
    exit 0
fi
# Finished before anything else: a completed run re-submitted under
# ISMIP7_AUTO_RESUME resumes at t_end, runs zero steps and ends at t_end, and
# that is success, not the no-progress case below.
done_yet=$(python -c "print(1 if float('$LAST_YEAR') >= float('$T_END') - 1e-6 else 0)")
if [ "$done_yet" = "1" ]; then
    echo "chain: reached t_end ${T_END} at ${LAST_YEAR}, finished"
    if [ -n "${ISMIP7_CHAIN_THEN:-}" ]; then
        if [ -z "${SLURM_JOB_ID:-}" ]; then
            echo "chain: no scheduler, so the follow-ons (${ISMIP7_CHAIN_THEN})" \
                 "are not queued; start them by hand"
            exit 0
        fi
        then_list="$ISMIP7_CHAIN_THEN"
        # Each follow-on is a fresh chain: it resumes its own checkpoints, not
        # this run's restart or wall budget, runs its driver's own period, not
        # this run's ISMIP7_T_START / ISMIP7_T_END, and queues nothing further.
        unset ISMIP7_CHAIN_THEN ISMIP7_WALL_STOP_MIN ISMIP7_RESTART \
              ISMIP7_T_START ISMIP7_T_END
        ismip7_chain_resources
        for next in $then_list; do
            echo "chain: queueing the follow-on ${next}"
            ISMIP7_EXPERIMENT="$next" \
                sbatch --dependency=afterok:"${SLURM_JOB_ID}" --export=ALL \
                       "${ISMIP7_CHAIN_RES[@]}" \
                       antarctica/scripts/batch_runners/projection.sbatch
        done
    fi
    exit 0
fi

# The verdict is settled and this run stopped short. Everything below only
# decides whether a successor is submitted.
if [ "$CHAIN" != "1" ]; then
    echo "chain: at ${LAST_YEAR} of ${T_END}, but chaining is disabled" \
         "(ISMIP7_CHAIN=$CHAIN); not resubmitting"
    exit 0
fi
# The chain is built on auto-resume: a successor inherits the environment and
# no ISMIP7_RESTART, so with resume off it would cold-start at t_start, end at
# the same year, and do it again forever. Ask runconfig, which owns the value
# set; the driver already ran, so the value is known good by now.
resume_on=$(python -c \
    "from icepack2_tools.runconfig import auto_resume; print(1 if auto_resume() else 0)")
if [ "$resume_on" != "1" ]; then
    echo "chain: at ${LAST_YEAR} of ${T_END}, but auto-resume is off" \
         "(ISMIP7_AUTO_RESUME=$ISMIP7_AUTO_RESUME), so a successor would"
    echo "chain: cold-start and repeat these years. Not resubmitting; set" \
         "ISMIP7_AUTO_RESUME=1 to chain, or continue by hand with ISMIP7_RESTART."
    exit 0
fi
# A clean exit short of t_end that also made no progress is not the wall-clock
# case the chain exists for: a solver give-up already exited above on the
# checkpoint's stalled attribute, and a job that died exited non-zero, so what
# is left is a run whose setup alone spent the whole budget. Another job under
# the same budget would do the same, so stop.
if [ "$START_YEAR" != "-" ]; then
    progressed=$(python -c "print(1 if float('$LAST_YEAR') > float('$START_YEAR') + 1e-6 else 0)")
    if [ "$progressed" = "0" ]; then
        echo "chain: no progress (resumed at ${START_YEAR}, ended at ${LAST_YEAR});" \
             "a successor under the same wall budget would not advance either. Stopping."
        exit 1
    fi
fi

echo "chain: at ${LAST_YEAR} of ${T_END}, resubmitting"
# Without the derived budget: the successor may land on a partition with a
# different limit and must derive its own. Without the restart either: an
# explicit one-off ISMIP7_RESTART wins over auto-resume in every driver, so
# inheriting it would send each link back to the same starting checkpoint and
# the chain would never advance past the first job's end year.
unset ISMIP7_WALL_STOP_MIN ISMIP7_RESTART
if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo "chain: no scheduler, so no successor to queue; resume by hand"
    exit 0
fi
# The successor also has to be told what to run on: this script carries no
# resource directives, so without the running job's own partition, constraint,
# wall limit, memory and task layout the link would land on the cluster's
# defaults, and a default of one task is refused by the guard above.
ismip7_chain_resources
sbatch --dependency=afterok:"${SLURM_JOB_ID}" --export=ALL \
       "${ISMIP7_CHAIN_RES[@]}" \
       antarctica/scripts/batch_runners/projection.sbatch
