#!/bin/bash
# run_core_matrix.sh - run the ISMIP7 AIS core-experiment matrix (cores 1-11)
# end to end, in protocol dependency order, with load gating and wall retry.
#
# Protocol order (cores 3-10 all branch from the historical endpoint, so the
# historicals must finish first; see experiment.py and control/run.py):
#
#   1,2   historical   CESM2-WACCM / MRI-ESM2-0        1850-2014
#   9,10  CTRL2015     both ESMs, constant 2015 climate   ->2300
#   3,4   ssp370       both ESMs                          ->2100
#   5,6   ssp126       both ESMs                          ->2300
#   7,8   ssp585       both ESMs                          ->2300
#   11    OCX          obs-constrained, cold start     1990-2025
#
# A core that branches from the historical is only launched once that ESM's
# historical endpoint is on disk, complete and post-fix. A missing, short or
# pre-fix historical would otherwise cold-start the dependent from the
# inversion geometry (or resume from a truncated endpoint), which breaks the
# projection-minus-CTRL cancellation the whole matrix rests on - and the
# drivers only WARN about it. This dependency check is provenance-only: it
# asks whether the endpoint is valid, never whether the caller asked to
# re-run, so FRESH=1 CORES=3,4,... still branches from a good historical.
#
# PROVENANCE: a completed timeseries is only reused when it postdates the
# forcing implementation ($PROV_REF). Results produced before the annual-mean
# atmosphere fix are invalid (January was applied as the whole year), so they
# are ARCHIVED under $R/archive_stale_<TS>/ - the timeseries, the final.h5 and
# the periodic checkpoints together - and re-run rather than silently skipped
# and fed to the audit. Moving them aside also stops a crashed relaunch from
# reading the superseded file back as its own progress, and a resume source is
# additionally required to postdate $PROV_REF so no run can continue from
# pre-fix geometry.
# Reuse precedence: FRESH=1 (re-run every selected core) beats REUSE=1 (reuse
# any completed timeseries unchecked) beats the provenance comparison. Neither
# flag affects the dependency check above or the resume-source check.
#
# WALL RETRY: the split-step diagnostic solve can exhaust its in-run rescue
# ladder late in a long projection and save-and-stop short of the target year.
# Relaunching from the saved state clears it - a fresh process re-runs the
# n=1->n continuation at the loaded geometry, which the in-run ladder cannot do
# (observed 3/3: ssp585-CESM 2096.7, CTRL-CESM 2268, CTRL-MRI 2250 all resumed
# and two reached 2300). This script therefore relaunches a short run from its
# own saved state, and keeps doing so while each attempt makes progress, up to
# MAX_ATTEMPTS. A run that stops advancing is left alone and reported. The
# resume source is the newest checkpoint whose recorded year is at or before
# the timeseries' last row, so a run is never continued from a state ahead of
# its own record. OCX is the exception: its driver cold-starts unconditionally
# and ignores ISMIP7_RESTART, so it runs at most once.
#
# Usage:
#   antarctica/scripts/run_core_matrix.sh                # full matrix
#   CORES=1,2,9 antarctica/scripts/run_core_matrix.sh    # a subset
#   ISMIP7_RUN_TAG=n3 NRANKS=8 antarctica/scripts/run_core_matrix.sh
#   FRESH=1 antarctica/scripts/run_core_matrix.sh        # ignore all prior output
#
# Env: ISMIP7_LC (32000), ISMIP7_N_FLOW (3), ISMIP7_RUN_TAG, ISMIP7_DT (0.1),
#      ISMIP7_FRICTION (budd), ISMIP7_FIXED_FRONT (1), ISMIP7_SUBCYCLES
#      (1,4,16,64), NRANKS (8), MAX_LOAD (cores-8, floored at 1),
#      MAX_ATTEMPTS (6), CORES (comma list, default all), FRESH, REUSE,
#      PROV_REF (icepack2_tools/forcing.py), ENS_HORIZON (2101),
#      VENV (~/venv-firedrake-2026).
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO"
R=antarctica/results
mkdir -p "$R/logs"

VENV="${VENV:-$HOME/venv-firedrake-2026}"
NRANKS="${NRANKS:-8}"          # 8 is MUMPS-safe at 32 km; raise for finer meshes
MAX_ATTEMPTS="${MAX_ATTEMPTS:-6}"
NCPU=$(python3 -c 'import os;print(os.cpu_count())')
MAX_LOAD="${MAX_LOAD:-$((NCPU - 8))}"
[ "$MAX_LOAD" -lt 1 ] && MAX_LOAD=1
CORES="${CORES:-1,2,9,10,3,4,5,6,7,8,11}"
FRESH="${FRESH:-}"
REUSE="${REUSE:-}"
PROV_REF="${PROV_REF:-icepack2_tools/forcing.py}"
ENS_HORIZON="${ENS_HORIZON:-2101}"   # ISMIP6 ensemble time base ends here

export OMP_NUM_THREADS=1
export ISMIP7_LC="${ISMIP7_LC:-32000}"
# ISMIP7_MESH is deliberately NOT defaulted here. The forward takes its mesh
# from the MAP checkpoint, which records its own mesh basename and uses that
# to pick the per-mesh boundary-id sidecar; a guessed name exported here would
# only be able to point a run at the wrong sidecar. Export it explicitly for a
# legacy MAP that predates the recorded attribute.
export ISMIP7_N_FLOW="${ISMIP7_N_FLOW:-3}"
export ISMIP7_FRICTION="${ISMIP7_FRICTION:-budd}"
export ISMIP7_DT="${ISMIP7_DT:-0.1}"
export ISMIP7_OUTPUT_INTERVAL="${ISMIP7_OUTPUT_INTERVAL:-10}"
export ISMIP7_APPARENT_MB="${ISMIP7_APPARENT_MB:-1}"
export ISMIP7_FIXED_FRONT="${ISMIP7_FIXED_FRONT:-1}"
export ISMIP7_SUBCYCLES="${ISMIP7_SUBCYCLES:-1,4,16,64}"
TAG="${ISMIP7_RUN_TAG:-}"
[ -n "$TAG" ] && export ISMIP7_RUN_TAG="$TAG"
SFX=${TAG:+_$TAG}
LC=$ISMIP7_LC
TS=$(date +%Y%m%d_%H%M%S)
START_EPOCH=$(date +%s)
PY="$VENV/bin/python"

PROV_MTIME=""
[ -e "$PROV_REF" ] && PROV_MTIME=$(stat -c %Y "$PROV_REF")
ARCHIVE="$R/archive_stale_$TS"

# core | label | driver | ESM | experiment-name stem | target year | kind
# kind: hist = historical, branch = branches from the historical endpoint,
#       cold = cold start with no restart support.
CORE_SPEC=(
  "1|hist_cesm|historical/cesm_waccm.py|CESM2-WACCM|hist_cesm2_waccm|2014|hist"
  "2|hist_mri|historical/mri_esm2.py|MRI-ESM2-0|hist_mri_esm2_0|2014|hist"
  "9|ctrl_cesm|control/run.py|CESM2-WACCM|ctrl2015_cesm2_waccm|2300|branch"
  "10|ctrl_mri|control/run.py|MRI-ESM2-0|ctrl2015_mri_esm2_0|2300|branch"
  "3|ssp370_cesm|projections/ssp370_cesm_waccm.py|CESM2-WACCM|ssp370_cesm2_waccm|2100|branch"
  "4|ssp370_mri|projections/ssp370_mri_esm2.py|MRI-ESM2-0|ssp370_mri_esm2_0|2100|branch"
  "5|ssp126_cesm|projections/ssp126_cesm_waccm.py|CESM2-WACCM|ssp126_cesm2_waccm|2300|branch"
  "6|ssp126_mri|projections/ssp126_mri_esm2.py|MRI-ESM2-0|ssp126_mri_esm2_0|2300|branch"
  "7|ssp585_cesm|projections/ssp585_cesm_waccm.py|CESM2-WACCM|ssp585_cesm2_waccm|2300|branch"
  "8|ssp585_mri|projections/ssp585_mri_esm2.py|MRI-ESM2-0|ssp585_mri_esm2_0|2300|branch"
  "11|ocx|projections/ocx.py|CESM2-WACCM|ocx|2025|cold"
)

declare -A HIST_STEM=() HIST_TARGET=() HIST_STATUS=() HIST_WHY=() \
           STATUS=() STEM_OF=() TARGET_OF=()
for spec in "${CORE_SPEC[@]}"; do
  IFS='|' read -r core label driver esm stem target kind <<< "$spec"
  STEM_OF[$core]=$stem
  TARGET_OF[$core]=$target
  if [ "$kind" = hist ]; then
    HIST_STEM[$esm]=$stem
    HIST_TARGET[$esm]=$target
  fi
done

wanted () { [[ ",$CORES," == *",$1,"* ]]; }
csv_of () { echo "$R/${1}${SFX}_${LC}_timeseries.csv"; }
h5_of ()  { echo "$R/${1}${SFX}_${LC}_final.h5"; }

last_year () {  # last simulated year of a timeseries, or empty
  local c="$1"
  [ -f "$c" ] || return 1
  awk -F, 'NR>1 && $1+0>0 {y=$1} END{if(y!="") print y}' "$c"
}

current_file () {  # file - true when it postdates $PROV_REF. A fact about the
                   # file only: FRESH/REUSE are reuse policy and must not
                   # change whether a dependency is judged valid.
  [ -n "$PROV_MTIME" ] || return 1
  [ -e "$1" ] || return 1
  [ "$(stat -c %Y "$1")" -ge "$PROV_MTIME" ]
}

reusable () {  # core label csv year - true when the existing output is current
  local core="$1" label="$2" csv="$3" year="$4"
  if [ -n "$FRESH" ]; then
    echo "[core $core] $label: FRESH=1, discarding existing timeseries at $year - re-running"
    return 1
  fi
  if [ -n "$REUSE" ]; then
    echo "[core $core] $label: REUSE=1, accepting existing timeseries at $year unchecked"
    return 0
  fi
  if [ -z "$PROV_MTIME" ]; then
    echo "[core $core] $label: provenance reference '$PROV_REF' not found - treating timeseries at $year as stale, re-running"
    return 1
  fi
  if current_file "$csv"; then
    echo "[core $core] $label: timeseries at $year postdates $PROV_REF - current"
    return 0
  fi
  echo "[core $core] $label: timeseries at $year PREDATES $PROV_REF (produced before the forcing fix) - re-running from scratch"
  return 1
}

archive_stale () {  # core label stem - move superseded output out of the way
  local core="$1" label="$2" stem="$3" f moved=() failed=() ckpts=0
  # The periodic checkpoints go too: pick_restart can resume from them, so a
  # pre-fix one left behind would become the resume state of a clean re-run.
  for f in "$(csv_of "$stem")" "$(h5_of "$stem")" \
           "$R/${stem}${SFX}_${LC}"_t*.h5; do
    [ -e "$f" ] || continue
    if mkdir -p "$ARCHIVE" && mv "$f" "$ARCHIVE/"; then
      case "$f" in
        *_"${LC}"_t*.h5) ckpts=$((ckpts + 1)) ;;
        *) moved+=("$(basename "$f")") ;;
      esac
    else
      failed+=("$f")
    fi
  done
  # The timeseries is what a later read would mistake for this run's result,
  # so assert it is gone rather than trusting mv's exit status alone.
  [ "${#failed[@]}" -eq 0 ] && [ -e "$(csv_of "$stem")" ] &&
    failed+=("$(csv_of "$stem")")
  if [ "${#failed[@]}" -gt 0 ]; then
    echo "[core $core] $label: CANNOT ARCHIVE ${failed[*]} - refusing to run." \
         "Superseded output left in place would be read back as this run's" \
         "result. Check permissions and free space on $R, then re-run."
    return 1
  fi
  if [ "${#moved[@]}" -gt 0 ] || [ "$ckpts" -gt 0 ]; then
    local what="${moved[*]:-}"
    [ "$ckpts" -gt 0 ] && what="${what:+$what + }$ckpts checkpoint(s)"
    echo "[core $core] $label: archived $what under $ARCHIVE/"
  fi
  return 0
}

pick_restart () {  # stem csv_last_year - "<t_yr>|<path>" of the newest
                   # checkpoint at or before the timeseries' last year, "|"
                   # when none qualifies, "?|" when the years cannot be read.
                   #
                   # final.h5 is only written when the time loop exits, so
                   # after a crash it still holds an EARLIER run's state and
                   # can sit ahead of the timeseries. Resuming from it would
                   # make run_simulation truncate the CSV at that later year
                   # and append from there, splicing two timelines together
                   # with a silent gap that reads back as progress.
                   #
                   # Candidates must also postdate $PROV_REF (or this run's
                   # launch when no reference exists): a pre-fix checkpoint is
                   # January-forcing geometry, and setup_model's restart guards
                   # check friction and a_ref_mb, not provenance.
  local out
  out=$("$PY" - "$R" "${1}${SFX}_${LC}" "${2:-0}" "${PROV_MTIME:-$START_EPOCH}" 2>/dev/null <<'PY'
import glob, os, sys
import h5py
d, base, want, floor = (sys.argv[1], sys.argv[2],
                        float(sys.argv[3]), float(sys.argv[4]))
cands = [os.path.join(d, base + "_final.h5")]
cands += sorted(glob.glob(os.path.join(d, base + "_t*.h5")))
best_t, best_fn = None, None
for fn in cands:
    if not os.path.exists(fn):
        continue
    try:
        if os.path.getmtime(fn) < floor:
            continue
        with h5py.File(fn, "r") as f:
            t = float(f["/"].attrs["t_yr"])
    except Exception:
        continue
    if t <= want + 0.05 and (best_t is None or t > best_t):
        best_t, best_fn = t, fn
print(f"{best_t:.4f}|{best_fn}" if best_fn else "|")
PY
) || out="?|"
  echo "${out:-?|}"
}

hist_ready () {  # esm - true when that ESM's historical endpoint is usable
  local esm="$1" stem target y h5
  if [ -n "${HIST_STATUS[$esm]:-}" ]; then
    [ "${HIST_STATUS[$esm]}" = ok ]
    return
  fi
  stem="${HIST_STEM[$esm]:-}"
  [ -n "$stem" ] || { HIST_STATUS[$esm]=ok; return 0; }
  target="${HIST_TARGET[$esm]}"
  h5=$(h5_of "$stem")
  y=$(last_year "$(csv_of "$stem")" || echo "")
  HIST_STATUS[$esm]=bad
  if [ ! -f "$h5" ]; then
    HIST_WHY[$esm]="no endpoint $(basename "$h5")"
  elif [ -z "$y" ] || ! awk "BEGIN{exit !($y >= $target)}"; then
    HIST_WHY[$esm]="endpoint stops at ${y:-nothing}, target $target"
  elif ! current_file "$h5"; then
    HIST_WHY[$esm]="endpoint predates $PROV_REF (pre-fix)"
  else
    HIST_STATUS[$esm]=ok
    return 0
  fi
  return 1
}

wait_for_load () {
  for _ in $(seq 1 720); do
    local l; l=$(python3 -c 'import os;print(int(os.getloadavg()[0]))')
    [ "$l" -le "$MAX_LOAD" ] && return 0
    sleep 30
  done
  echo "    (load never dropped below $MAX_LOAD; proceeding anyway)"
}

run_core () {  # core label driver esm stem target kind
  local core="$1" label="$2" driver="$3" esm="$4" stem="$5" target="$6" kind="$7"
  local csv h5 prev now attempt log attempts
  csv=$(csv_of "$stem"); h5=$(h5_of "$stem")

  now=$(last_year "$csv" || echo "")
  if [ -n "$now" ] && ! reusable "$core" "$label" "$csv" "$now"; then
    archive_stale "$core" "$label" "$stem" || return 1
    now=$(last_year "$csv" || echo "")
  fi
  if [ -n "$now" ] && awk "BEGIN{exit !($now >= $target)}"; then
    echo "[core $core] $label already at $now (target $target) - skipping"
    return 0
  fi

  # OCX cold-starts unconditionally, so a relaunch cannot resume it and would
  # truncate the previous timeseries instead of extending it.
  attempts="$MAX_ATTEMPTS"
  [ "$kind" = cold ] && attempts=1

  for attempt in $(seq 1 "$attempts"); do
    prev="$now"
    wait_for_load
    log="$R/logs/matrix_${TS}_core${core}_${label}_a${attempt}.log"
    # Attempt 1 uses the driver's own restart logic (projections and the CTRL
    # branch from the historical endpoint). Later attempts are wall retries and
    # must resume from this run's own saved state.
    local restart_env=() pick pick_t pick_fn
    if [ "$attempt" -gt 1 ]; then
      pick=$(pick_restart "$stem" "$now")
      pick_t=${pick%%|*}; pick_fn=${pick#*|}
      if [ "$pick_t" = "?" ]; then
        if [ ! -f "$h5" ] || [ "$csv" -nt "$h5" ]; then
          echo "[core $core] $label: cannot confirm a saved state consistent" \
               "with the timeseries at ${now:-nothing} - stopping rather than" \
               "resuming from a checkpoint that may be ahead of the record"
          return 1
        fi
        pick_fn="$h5"; pick_t="$now"
      elif [ -z "$pick_fn" ]; then
        echo "[core $core] $label: no checkpoint at or before the timeseries'" \
             "last year ${now:-nothing} - stopping rather than splicing two" \
             "timelines into one record"
        return 1
      fi
      restart_env=(ISMIP7_RESTART="$REPO/$pick_fn")
      echo "[core $core] $label retry $((attempt-1)) from $now" \
           "(resuming $(basename "$pick_fn") at t=$pick_t)"
    else
      echo "[core $core] $label starting (target $target)"
    fi
    local rc=0
    env "${restart_env[@]}" ISMIP7_ESM="$esm" \
      mpiexec -n "$NRANKS" "$PY" "antarctica/scripts/$driver" > "$log" 2>&1 || rc=$?
    now=$(last_year "$csv" || echo "")
    # A graceful wall stop saves and returns 0. A non-zero exit is a crash
    # (this machine is shared, so an OOM kill or MPI abort is plausible): it
    # costs one attempt and is never counted as reaching the target, but the
    # no-progress check below still ends a core that keeps dying in place.
    if [ "$rc" -ne 0 ]; then
      echo "[core $core] $label CRASHED (exit $rc) on attempt $attempt at" \
           "${now:-nothing} - see ${log##*/}"
      if [ -n "$prev" ] && [ -n "$now" ] &&
         awk "BEGIN{exit !($now <= $prev + 0.001)}"; then
        echo "[core $core] $label crashed without advancing past $prev - giving up"
        return 1
      fi
      continue
    fi
    echo "    -> reached ${now:-nothing}  (log ${log##*/})"
    [ -z "$now" ] && { echo "[core $core] produced no output - see log"; return 1; }
    awk "BEGIN{exit !($now >= $target)}" && { echo "[core $core] $label COMPLETE at $now"; return 0; }
    # only retry while the wall is actually moving forward
    if [ -n "$prev" ] && awk "BEGIN{exit !($now <= $prev + 0.001)}"; then
      echo "[core $core] $label stalled at $now (no progress) - giving up"
      return 1
    fi
  done
  if [ "$kind" = cold ]; then
    echo "[core $core] $label short of target at $now - not retried (its driver ignores ISMIP7_RESTART)"
  else
    echo "[core $core] $label short of target after $attempts attempts (at $now)"
  fi
  return 1
}

echo "ISMIP7 core matrix | lc=$LC n=$ISMIP7_N_FLOW tag='${TAG:-none}' dt=$ISMIP7_DT ranks=$NRANKS"
echo "cores: $CORES | max load $MAX_LOAD/$NCPU | wall retries up to $MAX_ATTEMPTS"
echo "fixed front=$ISMIP7_FIXED_FRONT subcycles=$ISMIP7_SUBCYCLES friction=$ISMIP7_FRICTION"
if [ -n "$FRESH" ]; then
  echo "reuse: FRESH=1, every selected core re-runs from scratch"
elif [ -n "$REUSE" ]; then
  echo "reuse: REUSE=1, any completed timeseries is accepted (provenance NOT checked)"
elif [ -n "$PROV_MTIME" ]; then
  echo "reuse: only output newer than $PROV_REF ($(date -d "@$PROV_MTIME" '+%Y-%m-%d %H:%M'))"
else
  echo "reuse: provenance reference '$PROV_REF' missing, all existing output treated as stale"
fi
echo

for spec in "${CORE_SPEC[@]}"; do
  IFS='|' read -r core label driver esm stem target kind <<< "$spec"
  wanted "$core" || continue
  if [ "$kind" = branch ] && ! hist_ready "$esm"; then
    echo "[core $core] $label BLOCKED: unusable $esm historical endpoint" \
         "(${HIST_WHY[$esm]:-unknown}). Running it now would cold-start (or" \
         "resume a truncated endpoint) and break the projection-minus-CTRL" \
         "cancellation. Run the $esm historical first."
    STATUS[$core]=blocked
    continue
  fi
  if run_core "$core" "$label" "$driver" "$esm" "$stem" "$target" "$kind"; then
    STATUS[$core]=ok
    [ "$kind" = hist ] && HIST_STATUS[$esm]=ok
  else
    STATUS[$core]=failed
    if [ "$kind" = hist ]; then
      HIST_STATUS[$esm]=bad
      HIST_WHY[$esm]="core $core did not reach $target in this run"
    fi
  fi
done

echo
echo "=== matrix summary ==="
for spec in "${CORE_SPEC[@]}"; do
  IFS='|' read -r core label driver esm stem target kind <<< "$spec"
  wanted "$core" || continue
  c=$(csv_of "$stem"); y=$(last_year "$c" || echo "-")
  printf "  core %-2s %-14s %8s / %s   %s\n" \
    "$core" "$label" "$y" "$target" "${STATUS[$core]:-not-run}"
done

# The audit only reads cores this run brought to their target; a blocked or
# failed core may still have an older timeseries on disk, which must not be
# reported as a matrix result.
echo
echo "=== observational audit ==="
for spec in "${CORE_SPEC[@]}"; do
  IFS='|' read -r core label driver esm stem target kind <<< "$spec"
  wanted "$core" || continue
  [ "${STATUS[$core]:-}" = ok ] || { printf "  %-14s (skipped: %s)\n" "$label" "${STATUS[$core]:-not-run}"; continue; }
  c=$(csv_of "$stem"); [ -f "$c" ] || continue
  printf "  %-14s %s\n" "$label" "$("$PY" antarctica/scripts/check_ismip6_track.py "$c" 2>&1 | tail -1)"
done

echo
# The CTRL only has to be a valid partner on disk - present, post-fix and long
# enough to cover the projection - not something this invocation re-ran, so a
# subset run such as CORES=5,6,7,8 still gets its comparison.
echo "=== ISMIP6 ensemble comparison (projection - same-ESM CTRL) ==="
for spec in "5|9" "6|10" "3|9" "4|10" "7|9" "8|10"; do
  IFS='|' read -r core ctrl_core <<< "$spec"
  wanted "$core" || continue
  p="${STEM_OF[$core]}"; pc=$(csv_of "$p")
  cc=$(csv_of "${STEM_OF[$ctrl_core]}")
  if [ "${STATUS[$core]:-}" != ok ]; then
    echo "--- core $core: $p - skipped (core $core ${STATUS[$core]:-not-run})"
    continue
  fi
  if [ ! -f "$cc" ]; then
    echo "--- core $core: $p - skipped (no CTRL timeseries $(basename "$cc"); run core $ctrl_core)"
    continue
  fi
  if [ "${STATUS[$ctrl_core]:-}" != ok ] && ! current_file "$cc"; then
    echo "--- core $core: $p - skipped (CTRL core $ctrl_core timeseries predates $PROV_REF)"
    continue
  fi
  # compare_ismip6.py clips to the projection-CTRL overlap and then to the
  # ISMIP6 ensemble time base, so the CTRL only has to reach the end of the
  # comparison window - a CTRL that stalled at 2268 still compares fine.
  proj_y=$(last_year "$pc" || echo "")
  ctrl_y=$(last_year "$cc" || echo "")
  need=$(awk -v a="${proj_y:-0}" -v b="$ENS_HORIZON" 'BEGIN{print (a<b)?a:b}')
  if [ -z "$ctrl_y" ] || ! awk "BEGIN{exit !($ctrl_y >= $need)}"; then
    echo "--- core $core: $p - skipped (CTRL core $ctrl_core reaches ${ctrl_y:-nothing}, needs $need to cover the comparison window)"
    continue
  fi
  echo "--- core $core: $p ---"
  "$PY" antarctica/scripts/compare_ismip6.py "$pc" "$cc" 2>&1 | tail -3
done
