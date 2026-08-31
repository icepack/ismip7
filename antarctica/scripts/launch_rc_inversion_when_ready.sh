#!/usr/bin/env bash
# Resource-gated launcher for the high-resolution regularized-Coulomb inversion.
#
# Polls free RAM (MemAvailable) and idle cores (nproc - load1); once BOTH clear
# their thresholds for STABLE consecutive checks, launches the RC inversion
# detached (nohup mpiexec) and exits. Single-instance via a lockfile; the lock is
# retained after launch as an "already fired" marker (remove it to re-arm).
#
# Run detached so it survives the shell/session:
#   setsid nohup antarctica/scripts/launch_rc_inversion_when_ready.sh >/dev/null 2>&1 &
# Tune via env (see defaults below). DRYRUN=1 prints one resource reading and exits.
#
# Status:  tail -f "$GATELOG"   (gate decisions) ;  tail -f the per-run inversion log
# Cancel:  kill the gate PID in the lockfile's first line; rm the lockfile to re-arm.
set -u

# ---- configuration (all env-overridable) ----
REPO="${REPO:-/home/andrew/projects/ismip7}"
PY="${PY:-$HOME/venv-firedrake/bin/python}"
# Default target: the 604 K-vertex / 1.2 M-cell 500 m mesh behind
# inversion_icepack2_500.h5 (the user's high-res precedent). Override ISMIP7_MESH
# for a different resolution. The inversion's auto path (antarctica_<lc*10>_<lc>.msh)
# may not match, so we pass ISMIP7_MESH explicitly.
MESH="${ISMIP7_MESH:-$REPO/antarctica/mesh/antarctica_64000_500_aniso.msh}"
LC="${ISMIP7_LC:-500}"
NRANKS="${NRANKS:-24}"                     # MPI ranks for the inversion
MIN_FREE_GB="${MIN_FREE_GB:-128}"          # require this much MemAvailable [GB] (~7M-DOF MUMPS, conservative)
MIN_FREE_CORES="${MIN_FREE_CORES:-28}"     # require this many idle cores (>= NRANKS)
STABLE="${STABLE:-3}"                      # consecutive passing checks before launch
INTERVAL="${INTERVAL:-120}"               # seconds between checks
MAXITER="${MAXITER:-500}"
# Set to 1 for BUFFERED meshes (h=0 at front -> no calving-terminus BC); leave
# empty for non-buffered meshes (keeps the calving back-pressure term).
NO_CALVING="${ISMIP7_NO_CALVING_TERMINUS:-}"
# Friction law for the inversion (budd -> inversion_icepack2_budd_<lc>.h5,
# regularized_coulomb -> _rc). Historically hardcoded to RC; the Budd
# production path (Jul 2026) needs per-resolution _budd MAPs.
FRICTION="${ISMIP7_FRICTION:-regularized_coulomb}"

# Logs under the repo so they survive reboots (the Jun 20 rc_500
# inversion log died with /tmp).
LOGDIR="$REPO/antarctica/results/logs"
mkdir -p "$LOGDIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
# n-tag so n=3 and n=4 inversion logs/locks don't collide on one machine
# (empty for n=4, matching the historical names; _n3 etc. otherwise).
NFLOW="${ISMIP7_N_FLOW:-3.0}"
NTAG="$(awk -v n="$NFLOW" 'BEGIN{ if ((n-4.0)^2 < 1e-12) print ""; else printf "_n%d", int(n+0.5) }')"
FTAG="$([ "$FRICTION" = "budd" ] && echo budd_inv || echo rc_inv)${NTAG}"
LOG="${LOG:-$LOGDIR/${FTAG}_${LC}_${STAMP}.log}"
GATELOG="${GATELOG:-$LOGDIR/${FTAG}_${LC}_gate.log}"
LOCK="${LOCK:-$LOGDIR/${FTAG}_${LC}.lock}"

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$GATELOG"; }

read_avail_gb(){ awk '/MemAvailable/{printf "%d", $2/1024/1024}' /proc/meminfo; }
read_idle_cores(){ awk -v n="$(nproc)" '{c=n-int($1+0.999); print (c<0)?0:c}' /proc/loadavg; }

if [ "${DRYRUN:-}" = "1" ]; then
  echo "mesh=$(basename "$MESH")  ranks=$NRANKS  need: RAM>=${MIN_FREE_GB}G cores>=${MIN_FREE_CORES}"
  echo "now: MemAvailable=$(read_avail_gb)G  idle_cores=$(read_idle_cores)  load1=$(awk '{print $1}' /proc/loadavg)"
  exit 0
fi

# single-instance guard (atomic: noclobber write fails if the lock exists)
if ! (set -C; echo "$$" > "$LOCK") 2>/dev/null; then
  log "lockfile $LOCK present (gate already running or inversion already launched). Exiting. rm it to re-arm."
  exit 0
fi
trap 'rm -f "$LOCK"' EXIT

[ -f "$MESH" ] || { log "ERROR: mesh not found: $MESH"; exit 1; }
[ -x "$PY" ]   || { log "ERROR: python not found: $PY"; exit 1; }
[ "$MIN_FREE_CORES" -ge "$NRANKS" ] || log "WARN: MIN_FREE_CORES ($MIN_FREE_CORES) < NRANKS ($NRANKS)"
log "ARMED: mesh=$(basename "$MESH") lc=$LC ranks=$NRANKS maxiter=$MAXITER no_calving='${NO_CALVING:-(off)}'"
log "  gate: RAM>=${MIN_FREE_GB}G AND idle_cores>=${MIN_FREE_CORES} for ${STABLE} checks @ ${INTERVAL}s; run log -> $LOG"

ok=0
while :; do
  gb="$(read_avail_gb)"; idle="$(read_idle_cores)"; load1="$(awk '{print $1}' /proc/loadavg)"
  if [ "$gb" -ge "$MIN_FREE_GB" ] && [ "$idle" -ge "$MIN_FREE_CORES" ]; then
    ok=$((ok+1))
    log "pass ${ok}/${STABLE}  (avail ${gb}G, idle ${idle} cores, load1 ${load1})"
    if [ "$ok" -ge "$STABLE" ]; then
      log "LAUNCHING RC inversion -> $LOG"
      cd "$REPO" || { log "ERROR: cd $REPO failed"; exit 1; }
      export OMP_NUM_THREADS=1 ISMIP7_LC="$LC" ISMIP7_MESH="$MESH" \
             ISMIP7_FRICTION="$FRICTION" ISMIP7_MAXITER="$MAXITER"
      [ -n "$NO_CALVING" ] && export ISMIP7_NO_CALVING_TERMINUS="$NO_CALVING"
      nohup mpiexec -n "$NRANKS" "$PY" antarctica/scripts/inversion_icepack2.py > "$LOG" 2>&1 &
      pid=$!
      echo "$pid" >> "$LOCK"
      log "launched: pid=$pid  (lock retained as run marker; rm $LOCK to re-arm)"
      trap - EXIT          # keep the lock so a re-run won't double-launch
      exit 0
    fi
  else
    [ "$ok" -gt 0 ] && log "reset (avail ${gb}G, idle ${idle} cores, load1 ${load1})"
    ok=0
  fi
  sleep "$INTERVAL"
done
