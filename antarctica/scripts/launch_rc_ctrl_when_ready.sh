#!/usr/bin/env bash
# Resource-gated launcher for the 500 m regularized-Coulomb CTRL2015 run.
#
# Runs antarctica/scripts/control/run.py with ISMIP7_FRICTION=regularized_coulomb
# against inversion_icepack2_rc_500.h5: constant 2015-2029 SMB climatology,
# OI-climatology ocean + per-basin calibrated K, dt=0.1, 2015->T_END.
# This is the clamp-free verification of the RC + h_clamp=0 pipeline: the
# initial state is the true BedMachine geometry and no thickness floor is
# ever applied, so total mass should track (SMB - melt - discharge) with no
# buffer-cell resurrection artifact.
#
# Polls free RAM and idle cores; once BOTH clear their thresholds for STABLE
# consecutive checks, launches detached (nohup mpiexec) and exits. Single-
# instance via a lockfile retained after launch ("already fired" marker —
# remove it to re-arm). Logs live under antarctica/results/logs/ so they
# survive reboots (the Jun 20 RC inversion log died with /tmp).
#
# Run detached so it survives the shell/session:
#   setsid nohup antarctica/scripts/launch_rc_ctrl_when_ready.sh >/dev/null 2>&1 &
# DRYRUN=1 prints one resource reading and exits.
set -u

# ---- configuration (all env-overridable) ----
REPO="${REPO:-/home/andrew/projects/ismip7}"
PY="${PY:-$HOME/venv-firedrake/bin/python}"
# The 500 m aniso mesh does not match the antarctica_<lc_coarse>_<lc>.msh
# default pattern, so pass it explicitly (must be the RC inversion's mesh;
# simulation.py verifies node ordering against the checkpoint).
MESH="${ISMIP7_MESH:-$REPO/antarctica/mesh/antarctica_64000_500_aniso.msh}"
LC="${ISMIP7_LC:-500}"
T_END="${ISMIP7_T_END:-2025}"              # 10-yr verification by default
DT="${ISMIP7_DT:-0.1}"
FIXED_FRONT="${ISMIP7_FIXED_FRONT:-1}"     # remove+tally ice past the initial extent
K_NPZ="${ISMIP7_K_PER_BASIN_NPZ:-$REPO/antarctica/results/calibrated_K_per_basin_2500.npz}"
NRANKS="${NRANKS:-24}"
MIN_FREE_GB="${MIN_FREE_GB:-128}"          # ~7M-DOF MUMPS, conservative
MIN_FREE_CORES="${MIN_FREE_CORES:-28}"     # >= NRANKS with headroom
STABLE="${STABLE:-5}"                      # consecutive passing checks
INTERVAL="${INTERVAL:-180}"                # seconds between checks

LOGDIR="$REPO/antarctica/results/logs"
mkdir -p "$LOGDIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
# n-tag so n=3 and n=4 CTRL logs/locks don't collide and the gated MAP name
# matches setup_model (empty for n=4, matching legacy names; _n3 etc. otherwise).
NFLOW="${ISMIP7_N_FLOW:-3.0}"
NTAG="$(awk -v n="$NFLOW" 'BEGIN{ if ((n-4.0)^2 < 1e-12) print ""; else printf "_n%d", int(n+0.5) }')"
# Geometry tag, matching icepack2_tools.naming.map_geom_tag: DG0 (the default)
# MAPs are _dg0, CG1 keeps the legacy untagged name. A MAP is only valid for
# the geometry space it was inverted under, so the gate must check the exact
# name the run will load.
GEOM="$(printf '%s' "${ISMIP7_GEOMETRY_SPACE:-dg0}" | tr 'A-Z' 'a-z')"
GTAG="$([ "$GEOM" = "cg1" ] && echo "" || echo "_dg0")"
LOG="${LOG:-$LOGDIR/rc_ctrl500${NTAG}_${STAMP}.log}"
GATELOG="${GATELOG:-$LOGDIR/rc_ctrl500${NTAG}_gate.log}"
LOCK="${LOCK:-$LOGDIR/rc_ctrl500${NTAG}.lock}"

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
  log "lockfile $LOCK present (gate already running or run already launched). Exiting. rm it to re-arm."
  exit 0
fi
trap 'rm -f "$LOCK"' EXIT

[ -f "$MESH" ]  || { log "ERROR: mesh not found: $MESH"; exit 1; }
[ -x "$PY" ]    || { log "ERROR: python not found: $PY"; exit 1; }
[ -f "$K_NPZ" ] || { log "ERROR: per-basin K npz not found: $K_NPZ"; exit 1; }
# Per-mesh sidecar preferred (icepack2_tools/boundary.py); the shared
# boundary_ids.json is the fallback and is overwritten by every mesh build.
MESH_STEM="$(basename "$MESH" .msh)"
BNDIDS="${ISMIP7_BNDIDS:-}"
if [ -z "$BNDIDS" ]; then
  if [ -f "$REPO/antarctica/mesh/boundary_ids_${MESH_STEM}.json" ]; then
    BNDIDS="$REPO/antarctica/mesh/boundary_ids_${MESH_STEM}.json"
  else
    BNDIDS="$REPO/antarctica/mesh/boundary_ids.json"
  fi
fi
[ -f "$BNDIDS" ] || { log "ERROR: boundary ids not found: $BNDIDS (untracked on this branch — restore the local copy or run make_boundary_ids.py)"; exit 1; }
[ -f "$REPO/antarctica/mesh/inversion_icepack2_rc${NTAG}${GTAG}_${LC}.h5" ] \
  || { log "ERROR: RC MAP not found: inversion_icepack2_rc${NTAG}${GTAG}_${LC}.h5 (geometry=$GEOM)"; exit 1; }
log "ARMED: RC CTRL2015 lc=$LC dt=$DT t_end=$T_END ranks=$NRANKS mesh=$(basename "$MESH")"
log "  gate: RAM>=${MIN_FREE_GB}G AND idle_cores>=${MIN_FREE_CORES} for ${STABLE} checks @ ${INTERVAL}s; run log -> $LOG"

ok=0
while :; do
  gb="$(read_avail_gb)"; idle="$(read_idle_cores)"; load1="$(awk '{print $1}' /proc/loadavg)"
  if [ "$gb" -ge "$MIN_FREE_GB" ] && [ "$idle" -ge "$MIN_FREE_CORES" ]; then
    ok=$((ok+1))
    log "pass ${ok}/${STABLE}  (avail ${gb}G, idle ${idle} cores, load1 ${load1})"
    if [ "$ok" -ge "$STABLE" ]; then
      log "LAUNCHING RC CTRL -> $LOG"
      cd "$REPO" || { log "ERROR: cd $REPO failed"; exit 1; }
      # Start the cold-start continuation fine (16 steps): on the 500 m
      # mesh each diagnostic solve is expensive, so starting finer is
      # cheaper in expectation than failing at 8 then re-ramping.
      export OMP_NUM_THREADS=1 ISMIP7_LC="$LC" ISMIP7_MESH="$MESH" \
             ISMIP7_FRICTION=regularized_coulomb ISMIP7_N_FLOW="$NFLOW" \
             ISMIP7_DT="$DT" ISMIP7_T_END="$T_END" \
             ISMIP7_K_PER_BASIN_NPZ="$K_NPZ" ISMIP7_BNDIDS="$BNDIDS" \
             ISMIP7_CONTINUATION_STEPS="${ISMIP7_CONTINUATION_STEPS:-16}"
      [ "$FIXED_FRONT" = "1" ] && export ISMIP7_FIXED_FRONT=1
      nohup mpiexec -n "$NRANKS" "$PY" antarctica/scripts/control/run.py > "$LOG" 2>&1 &
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
