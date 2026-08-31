#!/usr/bin/env bash
# Resource-gated launcher for the Budd CTRL2015 validation run.
#
# The RC forward is prognostically unstable (velocity blows up ~2x/step -
# see the Jul 2026 diagnosis), so production forward runs use the proven-
# stable Budd residual friction (ISMIP7_FRICTION=budd: PISM-delta grounded
# floor via N_hat, exact-zero shelf drag, driven by the
# inversion_icepack2_budd MAP). This runs the 2500 m
# Budd CTRL with all the conservation work: RACMO SMB, OI-climatology
# ocean + per-basin K, exactly-conservative transport, fixed calving
# front, h_clamp_init=0 (true geometry), dt=0.1, 2015->T_END.
#
# Config note: h_clamp_init=0 + fixed front is NEW for Budd (May runs used
# h_clamp_init=10), so this is also the stability confirmation. A blown-up
# or non-converging run self-terminates within ~10 steps (continuation
# failure), so gating it is safe. Watch the budget columns in the
# timeseries CSV: outflux must stay ~1e3 Gt/yr (a runaway to 1e4-1e5 is
# the RC-style instability); resid must stay ~0.
#
#   setsid nohup antarctica/scripts/launch_ctrl_when_ready.sh >/dev/null 2>&1 &
# DRYRUN=1 prints one resource reading and exits.
set -u

REPO="${REPO:-/home/andrew/projects/ismip7}"
PY="${PY:-$HOME/venv-firedrake/bin/python}"
LC="${ISMIP7_LC:-2500}"
MESH="${ISMIP7_MESH:-$REPO/antarctica/mesh/antarctica_64000_2500.msh}"
T_END="${ISMIP7_T_END:-2025}"              # 10-yr validation by default
DT="${ISMIP7_DT:-0.1}"
K_NPZ="${ISMIP7_K_PER_BASIN_NPZ:-$REPO/antarctica/results/calibrated_K_per_basin_2500.npz}"
NRANKS="${NRANKS:-16}"
MIN_FREE_GB="${MIN_FREE_GB:-64}"           # 2500 m MUMPS is lighter than 500 m
MIN_FREE_CORES="${MIN_FREE_CORES:-20}"
STABLE="${STABLE:-5}"
INTERVAL="${INTERVAL:-180}"

LOGDIR="$REPO/antarctica/results/logs"
mkdir -p "$LOGDIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
# n-tag so n=3 and n=4 CTRL logs/locks don't collide and the gated MAP name
# matches setup_model (empty for n=4, matching legacy names; _n3 etc. otherwise).
NFLOW="${ISMIP7_N_FLOW:-3.0}"
NTAG="$(awk -v n="$NFLOW" 'BEGIN{ if ((n-4.0)^2 < 1e-12) print ""; else printf "_n%d", int(n+0.5) }')"
LOG="${LOG:-$LOGDIR/budd_ctrl${LC}${NTAG}_${STAMP}.log}"
GATELOG="${GATELOG:-$LOGDIR/budd_ctrl${LC}${NTAG}_gate.log}"
LOCK="${LOCK:-$LOGDIR/budd_ctrl${LC}${NTAG}.lock}"
BNDIDS="${ISMIP7_BNDIDS:-$REPO/antarctica/mesh/boundary_ids.json}"

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$GATELOG"; }
read_avail_gb(){ awk '/MemAvailable/{printf "%d", $2/1024/1024}' /proc/meminfo; }
read_idle_cores(){ awk -v n="$(nproc)" '{c=n-int($1+0.999); print (c<0)?0:c}' /proc/loadavg; }

if [ "${DRYRUN:-}" = "1" ]; then
  echo "friction=budd  mesh=$(basename "$MESH")  lc=$LC ranks=$NRANKS  need: RAM>=${MIN_FREE_GB}G cores>=${MIN_FREE_CORES}"
  echo "now: MemAvailable=$(read_avail_gb)G idle_cores=$(read_idle_cores) load1=$(awk '{print $1}' /proc/loadavg)"
  exit 0
fi

if ! (set -C; echo "$$" > "$LOCK") 2>/dev/null; then
  log "lockfile $LOCK present (gate running or run launched). Exiting. rm to re-arm."
  exit 0
fi
trap 'rm -f "$LOCK"' EXIT

[ -f "$MESH" ]   || { log "ERROR: mesh not found: $MESH"; exit 1; }
[ -x "$PY" ]     || { log "ERROR: python not found: $PY"; exit 1; }
[ -f "$K_NPZ" ]  || { log "ERROR: per-basin K npz not found: $K_NPZ"; exit 1; }
[ -f "$BNDIDS" ] || { log "ERROR: boundary ids not found: $BNDIDS"; exit 1; }
[ -f "$REPO/antarctica/mesh/inversion_icepack2_budd${NTAG}_${LC}.h5" ] \
  || { log "ERROR: Budd MAP inversion_icepack2_budd${NTAG}_${LC}.h5 not found"; exit 1; }
log "ARMED: Budd CTRL2015 lc=$LC dt=$DT t_end=$T_END ranks=$NRANKS mesh=$(basename "$MESH")"
log "  gate: RAM>=${MIN_FREE_GB}G AND idle>=${MIN_FREE_CORES} for ${STABLE}x${INTERVAL}s; run log -> $LOG"

ok=0
while :; do
  gb="$(read_avail_gb)"; idle="$(read_idle_cores)"; load1="$(awk '{print $1}' /proc/loadavg)"
  if [ "$gb" -ge "$MIN_FREE_GB" ] && [ "$idle" -ge "$MIN_FREE_CORES" ]; then
    ok=$((ok+1))
    log "pass ${ok}/${STABLE}  (avail ${gb}G, idle ${idle} cores, load1 ${load1})"
    if [ "$ok" -ge "$STABLE" ]; then
      log "LAUNCHING Budd CTRL -> $LOG"
      cd "$REPO" || { log "ERROR: cd $REPO failed"; exit 1; }
      export OMP_NUM_THREADS=1 ISMIP7_LC="$LC" ISMIP7_MESH="$MESH" \
             ISMIP7_FRICTION=budd ISMIP7_N_FLOW="$NFLOW" \
             ISMIP7_DT="$DT" ISMIP7_T_END="$T_END" \
             ISMIP7_K_PER_BASIN_NPZ="$K_NPZ" ISMIP7_BNDIDS="$BNDIDS" \
             ISMIP7_H_CLAMP_INIT=0 ISMIP7_FIXED_FRONT=1 \
             ISMIP7_CONTINUATION_STEPS="${ISMIP7_CONTINUATION_STEPS:-12}"
      nohup mpiexec -n "$NRANKS" "$PY" antarctica/scripts/control/run.py > "$LOG" 2>&1 &
      pid=$!
      echo "$pid" >> "$LOCK"
      log "launched: pid=$pid  (lock retained as run marker; rm $LOCK to re-arm)"
      trap - EXIT
      exit 0
    fi
  else
    [ "$ok" -gt 0 ] && log "reset (avail ${gb}G, idle ${idle} cores, load1 ${load1})"
    ok=0
  fi
  sleep "$INTERVAL"
done
