#!/bin/bash
# The staged initialization, relaxation, re-inversion and forward routine of
# antarctica/INITIALIZATION_ROUTINE.md, on this repository's drivers and
# batch runners.
#
#   antarctica/scripts/initialization_routine.sh OUTDIR init
#   antarctica/scripts/initialization_routine.sh OUTDIR next
#   antarctica/scripts/initialization_routine.sh OUTDIR status
#   antarctica/scripts/initialization_routine.sh OUTDIR resume STAGE
#   antarctica/scripts/initialization_routine.sh OUTDIR watch
#
# `init` records the objective (every exported ISMIP7_* knob, the Python path
# and the MPI and kernel-cache settings of the calling shell, less the knobs
# the stages own) in OUTDIR/routine.env, and starts. From then on the routine
# is a state machine over the stage products in OUTDIR and the forwards'
# results directory: `next` runs or submits every stage whose inputs exist and
# which has not been started, and nothing else, so it can be called any number
# of times. `status` prints one line per stage. `resume STAGE` submits a stage
# again (a chain whose links all died); its chain resumes the stage's own
# checkpoint. `watch` (slurm mode) is the watchdog loop for a login node or a
# small allocation: it resubmits a started stage whose chain has no job left
# and cancels a forward link whose log has gone silent (a hung solve), with
# the node excluded from the stage's next submission. It also cancels a
# RUNNING forward link whose log holds ISMIP7_ROUTINE_FROZEN_MAX (5)
# zero-iteration momentum solves (`snes_its=0`, a frozen velocity) and parks
# its stage, in every routine; 0 disables the check, and a negative or
# non-numeric count is refused. ISMIP7_ROUTINE_STRICT=1 resubmits a chain only
# when the scheduler ended its newest link (TIMEOUT, PREEMPTED, NODE_FAIL,
# CANCELLED, BOOT_FAIL) and parks the stage on any other end; a link the
# accounting still holds live is left for the next pass. `resume STAGE`
# clears a park.
#
# Two launch modes (ISMIP7_ROUTINE_LAUNCH, recorded at init):
#   mpiexec   every stage runs in the foreground, in order, on
#             ISMIP7_ROUTINE_NRANKS ranks; `init` and `next` return when no
#             stage is ready (a workstation or a single long allocation).
#   slurm     every stage is a chain of batch_runners links (inversion.sbatch
#             with an iteration budget across its links, projection.sbatch
#             to the driver's end year), submitted with the resources in
#             ISMIP7_ROUTINE_SBATCH_INV / ISMIP7_ROUTINE_SBATCH_FWD. The link
#             that finishes a stage calls `next` (ISMIP7_CHAIN_THEN_CMD), so
#             the routine runs itself to the end without a login-node process.
#
# Stages, in order (ISMIP7_ROUTINE_STAGES selects a subset; a stage not
# listed is never started and its product, if present, still serves the ones
# after it):
#   shelf1 friction1 shelf2 friction2   the staged inversion: the floating
#             fluidity inside its bounds, then the friction, twice, the second
#             cycle at half the iterations (ISMIP7_ROUTINE_ITERS, 40)
#   relax     ISMIP7_ROUTINE_RELAX_YR (20) years of the control driver from
#             the last MAP with the apparent mass balance off and the front
#             pinned; its end state is stamped as a relaxation end state
#             (relaxed_state.h5, with the MAP's objective settings)
#   reinvert  = r_shelf1 r_friction1 r_shelf2 r_friction2, the four stages
#             again from the relaxed state; the result is map_forward.h5
#   hist      one historical per ESM (ISMIP7_ROUTINE_ESMS) from
#             ISMIP7_ROUTINE_HIST_START (2002) to 2015 with the apparent mass
#             balance off: its first year is the forced relaxation year of the
#             document (the jump from climatology to the ESM's forcing), a
#             spin-up year, and the years from 2003 are the core historical.
#             A start before 2003 takes the 2003 geometry
#             (ISMIP7_GEOMETRY_BACKDATE=12, 2015 - max(start, 2003)) as the
#             geometry at the start of the forced year, so the historical
#             proper starts in 2003 from the model's own state
#   control   the control per ESM from the historical's endpoint
#   proj      the projections per ESM (ISMIP7_ROUTINE_SCENARIOS, ssp126 ssp370
#             ssp585) from the same endpoint
#             Neither carries an apparent mass-balance reference: the
#             historical runs without one, and a reference can only be frozen
#             at the state a run starts balanced from, never built from an
#             evolved one (simulation.setup_model refuses), so the control and
#             the projections drift with the model's own imbalance and their
#             difference removes it.
#   ocx       the observed-forcing experiment from map_forward.h5
# map_forward.h5 is the product of the last inversion stage listed, so a run
# without `relax`/`reinvert` starts its forwards from the staged MAP.
#
# The front of every forward after the relaxation is a level-set calving law
# from icepack_tools.calving: ISMIP7_ROUTINE_CALVING (hfb), its parameters
# ISMIP7_ROUTINE_CALVING_PARAMS (sigma_max=0.15) and the advance rule
# ISMIP7_ROUTINE_FRONT_ADVANCE (none: retreat-only). The defaults are the
# submission front of ismip7_front_routine.sh; a campaign member
# (submit_calving_campaign.sh) is run as a full suite by naming its law,
# parameters and advance here. The collapse mask applies to the projections
# only; the relax stage pins the front.
# ISMIP7_ROUTINE_CALVING_FLOATING_ONLY (none) confines the law to floating ice
# (ISMIP7_CALVING_FLOATING_ONLY=1) in the projections (ssp) or every forward
# (all).
#
# An inversion stage moves one control alone, which only the TAO path does:
# `init` refuses a recorded ISMIP7_GRAD_PRECOND other than mass_consistent or
# prior when an inversion stage is listed.
set -euo pipefail

OUT=${1:?usage: initialization_routine.sh OUTDIR init|next|status|resume STAGE}
CMD=${2:-next}
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"
REPO=$(cd "$HERE/../.." && pwd)
ENVFILE="$OUT/routine.env"
LOGFILE="$OUT/routine.log"

log() { echo "[routine $(date '+%F %T')] $*" | tee -a "$LOGFILE" >&2; }
die() { log "ERROR: $*"; exit 1; }

# The knobs a stage owns: unset before every stage so one stage's setting
# never reaches another through the environment the hook inherits.
STAGE_KNOBS=(ISMIP7_INVERT ISMIP7_PHI_BOUNDS ISMIP7_FLUIDITY_CONTROL
    ISMIP7_DHDT_REWEIGHT ISMIP7_DHDT_REWEIGHT_SCALE ISMIP7_MAXITER
    ISMIP7_CHAIN_MAX ISMIP7_CHAIN_DEPTH ISMIP7_CHAIN_ITER_BUDGET ISMIP7_MAP_OUT
    ISMIP7_WARM_START ISMIP7_WARM_START_STRICT ISMIP7_WARM_START_STATE
    ISMIP7_WARM_START_FLUIDITY ISMIP7_WARM_START_THETA ISMIP7_WARM_START_PHI
    ISMIP7_WARM_START_GEOMETRY ISMIP7_WARM_START_PRIOR ISMIP7_INVERSION
    ISMIP7_RESTART ISMIP7_T_START ISMIP7_T_END ISMIP7_RUN_TAG ISMIP7_EXPERIMENT
    ISMIP7_ESM ISMIP7_APPARENT_MB ISMIP7_OUTPUT ISMIP7_FIXED_FRONT
    ISMIP7_CALVING ISMIP7_CALVING_PARAMS ISMIP7_FRONT_ADVANCE ISMIP7_FRACTURE
    ISMIP7_CALVING_FLOATING_ONLY
    ISMIP7_CHAIN_THEN ISMIP7_CHAIN_THEN_CMD ISMIP7_WALL_STOP_MIN
    ISMIP7_AUTO_RESUME ISMIP7_CHAIN ISMIP7_SNES_LOG ISMIP7_GEOMETRY_BACKDATE)

# ── the recorded environment ─────────────────────────────────────────────
record_env() {
    # The objective and the machine settings of the calling shell, less the
    # stage knobs and the routine's own resolved settings (appended below).
    local names pattern
    pattern='^declare -x (ISMIP7_|PYTHONPATH=|OMPI_|PYOP2_|FIREDRAKE_|OMP_NUM_THREADS=|PETSC_|TMPDIR=)'
    {
        echo "# recorded by initialization_routine.sh init, $(date -Is), $(hostname)"
        echo "# repository $REPO at $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo '?')"
        declare -px | grep -E "$pattern" \
            | grep -vE "^declare -x ($(IFS='|'; echo "${STAGE_KNOBS[*]}"))=" \
            | grep -v '^declare -x ISMIP7_ROUTINE_' \
            | sed 's/^declare -x /export /'
        echo "export ISMIP7_ROUTINE_LAUNCH=${ISMIP7_ROUTINE_LAUNCH:-mpiexec}"
        echo "export ISMIP7_ROUTINE_TAG=${ISMIP7_ROUTINE_TAG:-$(basename "$OUT")}"
        echo "export ISMIP7_ROUTINE_ITERS=${ISMIP7_ROUTINE_ITERS:-40}"
        echo "export ISMIP7_ROUTINE_RELAX_YR=${ISMIP7_ROUTINE_RELAX_YR:-20}"
        echo "export ISMIP7_ROUTINE_RELAX_START=${ISMIP7_ROUTINE_RELAX_START:-2015}"
        echo "export ISMIP7_ROUTINE_HIST_START=${ISMIP7_ROUTINE_HIST_START:-2002}"
        echo "export ISMIP7_ROUTINE_CALVING=${ISMIP7_ROUTINE_CALVING:-hfb}"
        echo "export ISMIP7_ROUTINE_CALVING_PARAMS=${ISMIP7_ROUTINE_CALVING_PARAMS-sigma_max=0.15}"
        echo "export ISMIP7_ROUTINE_FRONT_ADVANCE=${ISMIP7_ROUTINE_FRONT_ADVANCE:-none}"
        echo "export ISMIP7_ROUTINE_CALVING_FLOATING_ONLY=${ISMIP7_ROUTINE_CALVING_FLOATING_ONLY:-none}"
        echo "export ISMIP7_ROUTINE_FLUIDITY_CONTROL=${ISMIP7_ROUTINE_FLUIDITY_CONTROL:-floating}"
        echo "export ISMIP7_ROUTINE_PHI_BOUNDS=${ISMIP7_ROUTINE_PHI_BOUNDS:--2.3,2.3}"
        echo "export ISMIP7_ROUTINE_DHDT_REWEIGHT=${ISMIP7_ROUTINE_DHDT_REWEIGHT:-24}"
        echo "export ISMIP7_ROUTINE_NRANKS=${ISMIP7_ROUTINE_NRANKS:-1}"
        echo "export ISMIP7_ROUTINE_ESMS=\"${ISMIP7_ROUTINE_ESMS:-CESM2-WACCM MRI-ESM2-0}\""
        echo "export ISMIP7_ROUTINE_SCENARIOS=\"${ISMIP7_ROUTINE_SCENARIOS:-ssp126 ssp370 ssp585}\""
        echo "export ISMIP7_ROUTINE_STAGES=\"${ISMIP7_ROUTINE_STAGES:-shelf1 friction1 shelf2 friction2 relax reinvert hist control proj ocx}\""
        echo "export ISMIP7_ROUTINE_WARM_START=\"${ISMIP7_ROUTINE_WARM_START:-}\""
        echo "export ISMIP7_ROUTINE_SBATCH_INV=\"${ISMIP7_ROUTINE_SBATCH_INV:-}\""
        echo "export ISMIP7_ROUTINE_SBATCH_FWD=\"${ISMIP7_ROUTINE_SBATCH_FWD:-}\""
        echo "export ISMIP7_ROUTINE_CHAIN_MAX=${ISMIP7_ROUTINE_CHAIN_MAX:-200}"
        echo "export ISMIP7_ROUTINE_STRICT=${ISMIP7_ROUTINE_STRICT:-0}"
        echo "export ISMIP7_ROUTINE_FROZEN_MAX=${ISMIP7_ROUTINE_FROZEN_MAX:-5}"
    } > "$ENVFILE.tmp"
    mv "$ENVFILE.tmp" "$ENVFILE"
}

load_env() {
    # A clean slate: the stage knobs and the scheduler's variables of the
    # calling job go, then the recorded environment comes in.
    local v
    for v in "${STAGE_KNOBS[@]}"; do unset "$v"; done
    for v in $(compgen -v SLURM_ 2>/dev/null || true); do unset "$v"; done
    unset SBATCH_JOB_NAME SBATCH_PARTITION SBATCH_ACCOUNT
    [ -f "$ENVFILE" ] || die "$ENVFILE not found; run '$SELF $OUT init' first"
    # shellcheck disable=SC1090
    . "$ENVFILE"
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
    [ "$ISMIP7_ROUTINE_LAUNCH" = mpiexec ] || [ "$ISMIP7_ROUTINE_LAUNCH" = slurm ] \
        || die "ISMIP7_ROUTINE_LAUNCH must be mpiexec or slurm, got $ISMIP7_ROUTINE_LAUNCH"
    if [ "$ISMIP7_ROUTINE_LAUNCH" = slurm ]; then
        [ -n "$ISMIP7_ROUTINE_SBATCH_INV" ] || die "ISMIP7_ROUTINE_SBATCH_INV (sbatch resources of an inversion link) is empty"
        [ -n "$ISMIP7_ROUTINE_SBATCH_FWD" ] || die "ISMIP7_ROUTINE_SBATCH_FWD (sbatch resources of a forward link) is empty"
    fi
    # A routine recorded before the calving knobs ran the submission front;
    # without these its forwards would be handed an empty law.
    : "${ISMIP7_ROUTINE_CALVING:=hfb}"
    : "${ISMIP7_ROUTINE_FRONT_ADVANCE:=none}"
    [ -n "${ISMIP7_ROUTINE_CALVING_PARAMS+set}" ] || ISMIP7_ROUTINE_CALVING_PARAMS=sigma_max=0.15
    case "$ISMIP7_ROUTINE_FRONT_ADVANCE" in
        free|none) ;;
        *) die "ISMIP7_ROUTINE_FRONT_ADVANCE must be free or none, got $ISMIP7_ROUTINE_FRONT_ADVANCE" ;;
    esac
    # Which forwards confine the law to floating ice
    # (ISMIP7_CALVING_FLOATING_ONLY=1): none (the law's own rule), ssp (the
    # projections: the hfb front calved grounded cliffs under ssp370, 22 Oct
    # 2026) or all. A routine's running stages keep the rule their
    # checkpoints record, so change this only with the affected stages
    # restarted from scratch.
    : "${ISMIP7_ROUTINE_CALVING_FLOATING_ONLY:=none}"
    case "$ISMIP7_ROUTINE_CALVING_FLOATING_ONLY" in
        none|ssp|all) ;;
        *) die "ISMIP7_ROUTINE_CALVING_FLOATING_ONLY must be none, ssp or all, got $ISMIP7_ROUTINE_CALVING_FLOATING_ONLY" ;;
    esac
    case "$ISMIP7_ROUTINE_CALVING ${ISMIP7_ROUTINE_CALVING_PARAMS:-}" in
        *" "*" "*) die "ISMIP7_ROUTINE_CALVING and ISMIP7_ROUTINE_CALVING_PARAMS take no spaces (parameters are comma separated)" ;;
    esac
    : "${ISMIP7_ROUTINE_STRICT:=0}"
    : "${ISMIP7_ROUTINE_FROZEN_MAX:=5}"
    case "$ISMIP7_ROUTINE_FROZEN_MAX" in
        *[!0-9]*) die "ISMIP7_ROUTINE_FROZEN_MAX must be a count of zero-iteration momentum solves (0 disables the check), got $ISMIP7_ROUTINE_FROZEN_MAX" ;;
    esac
    TAG=$ISMIP7_ROUTINE_TAG
    LC=${ISMIP7_LC:?ISMIP7_LC must be recorded}
    RESULTS="$REPO/antarctica/results"
    STAGE_LIST=$(expand_stages "$ISMIP7_ROUTINE_STAGES")
}

check_objective() {   # the recorded objective can run every listed stage
    local s precond=${ISMIP7_GRAD_PRECOND:-none}
    for s in $STAGE_LIST; do
        if is_inversion "$s" && [ "$precond" != mass_consistent ] && [ "$precond" != prior ]; then
            log "ERROR: stage $s inverts one control alone, which needs the TAO path; export ISMIP7_GRAD_PRECOND=mass_consistent or prior (got $precond)"
            return 1
        fi
    done
}

# ── stage names ──────────────────────────────────────────────────────────
INV_ORDER="shelf1 friction1 shelf2 friction2 r_shelf1 r_friction1 r_shelf2 r_friction2"

esm_key() {   # the projection.sbatch experiment suffix of an ESM
    case "$1" in
        CESM2-WACCM) echo cesm_waccm ;;
        MRI-ESM2-0)  echo mri_esm2 ;;
        *) die "no experiment key for ESM $1" ;;
    esac
}
esm_tag() {   # the experiment.py name element of an ESM
    echo "$1" | tr 'A-Z-' 'a-z_'
}
key_esm() {   # the ESM of a stage name's suffix
    case "$1" in
        *_cesm_waccm) echo CESM2-WACCM ;;
        *_mri_esm2)   echo MRI-ESM2-0 ;;
        *) die "no ESM in stage name $1" ;;
    esac
}

expand_stages() {
    local s esm scen out=""
    for s in $1; do
        case "$s" in
            reinvert) out="$out r_shelf1 r_friction1 r_shelf2 r_friction2" ;;
            hist|control)
                for esm in $ISMIP7_ROUTINE_ESMS; do out="$out ${s}_$(esm_key "$esm")"; done ;;
            proj)
                for esm in $ISMIP7_ROUTINE_ESMS; do
                    for scen in $ISMIP7_ROUTINE_SCENARIOS; do out="$out ${scen}_$(esm_key "$esm")"; done
                done ;;
            *) out="$out $s" ;;
        esac
    done
    echo "$out"
}

in_list() { [[ " $2 " == *" $1 "* ]]; }
is_inversion() { in_list "$1" "$INV_ORDER"; }

# The inversion stage before STAGE in the routine's order (its warm start).
inversion_before() {
    local s prev=""
    for s in $INV_ORDER; do
        [ "$s" = "$1" ] && { echo "$prev"; return; }
        in_list "$s" "$STAGE_LIST" && prev=$s
    done
}
# The last inversion stage of the first cycle: the MAP the relaxation starts
# from and whose name the re-inversion carries.
relax_source_stage() {
    local s last=""
    for s in shelf1 friction1 shelf2 friction2; do in_list "$s" "$STAGE_LIST" && last=$s; done
    echo "$last"
}
# The last inversion stage listed: the MAP every forward starts from.
last_inversion_stage() {
    local s last=""
    for s in $INV_ORDER; do in_list "$s" "$STAGE_LIST" && last=$s; done
    echo "$last"
}

map_path() {
    # A re-inversion from a relaxation end state has to write the name
    # icepack2_tools.relaxation.relaxed_map_name gives it (the source MAP's
    # stem with _relax<year>); the inversion refuses any other.
    if [ "$1" = r_shelf1 ] && in_list relax "$STAGE_LIST"; then
        echo "$OUT/map_$(relax_source_stage)_relax${ISMIP7_ROUTINE_RELAX_START}.h5"
    else
        echo "$OUT/map_$1.h5"
    fi
}

# The forward experiment key and the results basename stem of a stage.
fwd_experiment() {
    case "$1" in
        relax)     echo control ;;
        control_*) echo control ;;
        hist_*|ssp*_*|ocx) echo "$1" ;;
        *) die "not a forward stage: $1" ;;
    esac
}
fwd_esm() {
    case "$1" in
        relax) set -- $ISMIP7_ROUTINE_ESMS; echo "$1" ;;
        ocx) echo "" ;;
        *) key_esm "$1" ;;
    esac
}
fwd_final() {
    local esm tag
    esm=$(fwd_esm "$1"); [ -n "$esm" ] && tag=$(esm_tag "$esm")
    case "$1" in
        relax)     echo "$RESULTS/ctrl2015_${tag}_${TAG}_relax_${LC}_final.h5" ;;
        control_*) echo "$RESULTS/ctrl2015_${tag}_${TAG}_${LC}_final.h5" ;;
        hist_*)    echo "$RESULTS/hist_${tag}_${TAG}_${LC}_final.h5" ;;
        ssp*_*)    echo "$RESULTS/${1%%_*}_${tag}_${TAG}_${LC}_final.h5" ;;
        ocx)       echo "$RESULTS/ocx_${TAG}_${LC}_final.h5" ;;
    esac
}
fwd_t_end() {
    case "$1" in
        relax)   echo $((ISMIP7_ROUTINE_RELAX_START + ISMIP7_ROUTINE_RELAX_YR)) ;;
        hist_*)  echo 2015 ;;
        ssp370_*) echo "${ISMIP7_ROUTINE_SSP370_END:-2101}" ;;
        control_*|ssp*_*) echo "${ISMIP7_ROUTINE_PROJ_END:-2301}" ;;
        ocx)     echo "${ISMIP7_ROUTINE_OCX_END:-2026}" ;;
    esac
}

# ── products ─────────────────────────────────────────────────────────────
h5attr() {   # h5attr FILE ATTR -> the root attribute, or nothing
    python - "$1" "$2" <<'PY' 2>/dev/null || true
import sys
import h5py
with h5py.File(sys.argv[1], "r") as h:
    v = h["/"].attrs.get(sys.argv[2])
if v is not None:
    print(v.decode() if isinstance(v, bytes) else v)
PY
}
reached() {   # reached FILE YEAR: the checkpoint FILE has t_yr >= YEAR
    local t
    [ -f "$1" ] || return 1
    t=$(h5attr "$1" t_yr); [ -n "$t" ] || return 1
    python -c "import sys; sys.exit(0 if float('$t') >= float('$2') - 1e-6 else 1)"
}
is_done() {
    if is_inversion "$1"; then [ -f "$(map_path "$1").done" ]
    else reached "$(fwd_final "$1")" "$(fwd_t_end "$1")"; fi
}
is_started() { [ -f "$OUT/$1.jobs" ]; }

forward_map() { echo "$OUT/map_forward.h5"; }

# Side products that follow from finished stages: the relaxation end state
# stamped from the relax run's final checkpoint, and the forward MAP link.
derive_products() {
    local src last
    src=$(relax_source_stage)
    if in_list relax "$STAGE_LIST" && is_done relax && [ ! -f "$OUT/relaxed_state.h5" ]; then
        log "stamping $(basename "$(fwd_final relax)") as the relaxation end state of map_$src.h5"
        stamp_relaxed_state "$(map_path "$src")" "$(fwd_final relax)" "$OUT/relaxed_state.h5" \
            "$ISMIP7_ROUTINE_RELAX_START" "$(fwd_t_end relax)" 2>&1 | tee -a "$LOGFILE" >&2
        [ -f "$OUT/relaxed_state.h5" ] || die "the relaxation end state was not stamped"
    fi
    last=$(last_inversion_stage)
    if [ -e "$(forward_map)" ]; then
        :
    elif [ -n "$last" ]; then
        if is_done "$last"; then
            ln -s "$(basename "$(map_path "$last")")" "$(forward_map)"
            log "forward MAP: $(forward_map) -> $(basename "$(map_path "$last")")"
        fi
    elif [ -n "$ISMIP7_ROUTINE_WARM_START" ]; then
        # No inversion stage listed: the forwards run from the MAP given.
        ln -s "$ISMIP7_ROUTINE_WARM_START" "$(forward_map)"
        log "forward MAP: $(forward_map) -> $ISMIP7_ROUTINE_WARM_START (no inversion stage listed)"
    fi
}

stamp_relaxed_state() {
    # The control driver's final checkpoint carries the geometry, the
    # controls, the observations and the mixed state a re-inversion needs; the
    # relaxation driver's end state carries, in addition, the record
    # icepack2_tools.relaxation defines (the source MAP, its objective
    # settings, the window, the step, the forcing). Written here onto a copy,
    # so the re-inversion takes the designed path: the relaxed geometry is
    # recorded in the MAP and the forwards start from it.
    PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}" python - "$@" <<'PY'
import os
import shutil
import sys

import h5py

repo = os.environ["PYTHONPATH"].split(":")[0]
sys.path.insert(0, repo)
from icepack2_tools.relaxation import relaxation_attrs, relax_fssa_attrs
from icepack2_tools.runconfig import file_sha256

src_map, state_in, state_out, t_start, t_end = sys.argv[1:6]
with h5py.File(src_map, "r") as h:
    src_attrs = dict(h["/"].attrs)
tmp = state_out + ".tmp"
shutil.copyfile(state_in, tmp)
with h5py.File(tmp, "r+") as h:
    a = h["/"].attrs
    dt = float(a.get("dt_yr", 0.0))
    attrs = relaxation_attrs(
        source_map=src_map, source_sha256=file_sha256(src_map),
        source_attrs=src_attrs, t_start=float(t_start), t_end=float(t_end),
        dt=dt, backdate_years=0.0, forcing="control-climatology")
    theta = float(os.environ.get("ISMIP7_FSSA_THETA", "0") or 0)
    ref = a.get("fssa_reference", "none")
    ref = ref.decode() if isinstance(ref, bytes) else str(ref)
    attrs.update(relax_fssa_attrs(theta, ref))
    for key, value in attrs.items():
        a[key] = value
os.replace(tmp, state_out)
print(f"relaxed state {state_out}: {len(attrs)} attributes stamped, dt {dt:g} yr")
PY
}

# ── readiness ────────────────────────────────────────────────────────────
ready() {
    local s
    case "$1" in
        shelf1) return 0 ;;
        friction1|shelf2|friction2|r_friction1|r_shelf2|r_friction2)
            s=$(inversion_before "$1"); [ -n "$s" ] && is_done "$s" ;;
        r_shelf1)
            if in_list relax "$STAGE_LIST"; then [ -f "$OUT/relaxed_state.h5" ]
            else s=$(inversion_before "$1"); [ -n "$s" ] && is_done "$s"; fi ;;
        relax) s=$(relax_source_stage); [ -n "$s" ] && is_done "$s" ;;
        hist_*|ocx) [ -e "$(forward_map)" ] ;;
        control_*|ssp*_*) reached "$(fwd_final "hist_${1#*_}")" 2015 ;;
        *) die "unknown stage $1" ;;
    esac
}

# ── launching ────────────────────────────────────────────────────────────
# Every stage runs under the recorded environment plus its own knobs, in a
# subshell, so nothing leaks between stages.
stage_env_inversion() {   # prints `export` lines for an inversion stage
    local stage=$1 warm iters
    case "$stage" in
        shelf1|friction1|r_shelf1|r_friction1) iters=$ISMIP7_ROUTINE_ITERS ;;
        *) iters=$((ISMIP7_ROUTINE_ITERS / 2)) ;;
    esac
    case "$stage" in
        *shelf*)
            echo "export ISMIP7_INVERT=phi ISMIP7_PHI_BOUNDS=$ISMIP7_ROUTINE_PHI_BOUNDS" ;;
        *friction*)
            echo "export ISMIP7_INVERT=theta ISMIP7_DHDT_REWEIGHT=$ISMIP7_ROUTINE_DHDT_REWEIGHT" ;;
    esac
    echo "export ISMIP7_FLUIDITY_CONTROL=$ISMIP7_ROUTINE_FLUIDITY_CONTROL"
    echo "export ISMIP7_MAXITER=$iters ISMIP7_CHAIN_ITER_BUDGET=$iters ISMIP7_CHAIN_MAX=$ISMIP7_ROUTINE_CHAIN_MAX"
    echo "export ISMIP7_MAP_OUT=$(map_path "$stage")"
    echo "export ISMIP7_SNES_LOG=$OUT/snes_$stage.log"
    if [ "$stage" = shelf1 ]; then
        warm=$ISMIP7_ROUTINE_WARM_START
    elif [ "$stage" = r_shelf1 ] && in_list relax "$STAGE_LIST"; then
        warm=$OUT/relaxed_state.h5
    else
        warm=$(map_path "$(inversion_before "$stage")")
    fi
    if [ -n "$warm" ]; then
        # Another stage's objective (its controls, bounds, weights) or a
        # relaxation end state: the strict handoff is off for the first link
        # and the controls are kept as they are; the chain's own links are
        # strict again (inversion.sbatch).
        echo "export ISMIP7_WARM_START=$warm ISMIP7_WARM_START_STRICT=0 ISMIP7_WARM_START_THETA=1 ISMIP7_WARM_START_PHI=1"
    fi
}

stage_env_forward() {   # prints `export` lines for a forward stage
    local stage=$1 esm
    esm=$(fwd_esm "$stage")
    echo "export ISMIP7_EXPERIMENT=$(fwd_experiment "$stage") ISMIP7_AUTO_RESUME=1 ISMIP7_CHAIN=1"
    [ -n "$esm" ] && echo "export ISMIP7_ESM=$esm"
    case "$stage" in
        relax)
            echo "export ISMIP7_INVERSION=$(map_path "$(relax_source_stage)")"
            echo "export ISMIP7_RUN_TAG=${TAG}_relax ISMIP7_T_END=$(fwd_t_end relax)"
            echo "export ISMIP7_APPARENT_MB=0 ISMIP7_OUTPUT=0 ISMIP7_FIXED_FRONT=1 ISMIP7_CALVING=none"
            echo "unset ISMIP7_CALVING_PARAMS ISMIP7_FRONT_ADVANCE ISMIP7_FRACTURE ISMIP7_CALVING_FLOATING_ONLY" ;;
        *)
            echo "export ISMIP7_INVERSION=$(forward_map) ISMIP7_RUN_TAG=$TAG ISMIP7_OUTPUT=1"
            case "$stage" in
                hist_*)
                    echo "export ISMIP7_T_START=$ISMIP7_ROUTINE_HIST_START ISMIP7_APPARENT_MB=0"
                    [ "$ISMIP7_ROUTINE_HIST_START" -lt 2003 ] \
                        && echo "export ISMIP7_GEOMETRY_BACKDATE=$((2015 - 2003))" ;;
                ocx)    echo "export ISMIP7_APPARENT_MB=0" ;;
                control_*|ssp*_*) echo "export ISMIP7_APPARENT_MB=0" ;;
            esac
            # the routine's calving law (see the header); the legacy pinned
            # front stays off under a law
            echo "export ISMIP7_CALVING=$ISMIP7_ROUTINE_CALVING ISMIP7_FRONT_ADVANCE=$ISMIP7_ROUTINE_FRONT_ADVANCE"
            if [ -n "$ISMIP7_ROUTINE_CALVING_PARAMS" ]; then
                echo "export ISMIP7_CALVING_PARAMS=$ISMIP7_ROUTINE_CALVING_PARAMS"
            else
                echo "unset ISMIP7_CALVING_PARAMS"
            fi
            echo "unset ISMIP7_FIXED_FRONT"
            case "$stage" in
                ssp*_*) echo "export ISMIP7_FRACTURE=mask_front" ;;
                *)      echo "export ISMIP7_FRACTURE=none" ;;
            esac
            case "$ISMIP7_ROUTINE_CALVING_FLOATING_ONLY:$stage" in
                all:*|ssp:ssp*_*) echo "export ISMIP7_CALVING_FLOATING_ONLY=1" ;;
                *)                echo "unset ISMIP7_CALVING_FLOATING_ONLY" ;;
            esac ;;
    esac
}

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
        *) die "no driver for experiment $1" ;;
    esac
}

launch() {
    local stage=$1 kind script res jobname envlines
    if is_inversion "$stage"; then
        kind=inversion; envlines=$(stage_env_inversion "$stage")
        script=antarctica/scripts/batch_runners/inversion.sbatch; res=$ISMIP7_ROUTINE_SBATCH_INV
    else
        kind=forward; envlines=$(stage_env_forward "$stage")
        script=antarctica/scripts/batch_runners/projection.sbatch; res=$ISMIP7_ROUTINE_SBATCH_FWD
    fi
    jobname="rt_${TAG}_${stage}"
    log "stage $stage ($kind): $(echo "$envlines" | grep -oE 'ISMIP7_[A-Z_]+=[^ ]+' | tr '\n' ' ')"
    if [ "$ISMIP7_ROUTINE_LAUNCH" = slurm ]; then
        local jid
        # shellcheck disable=SC2086
        local excl=""
        # nodes the watchdog saw a link hang on (OUT/exclude_nodes, one per line)
        [ -s "$OUT/exclude_nodes" ] && excl="--exclude=$(sort -u "$OUT/exclude_nodes" | paste -sd,)"
        jid=$(cd "$REPO" && eval "$envlines" && export ISMIP7_CHAIN_THEN_CMD="bash $SELF $OUT next" \
              && sbatch --parsable $res $excl -J "$jobname" --export=ALL "$script") \
            || die "sbatch failed for stage $stage"
        echo "$jid $(date -Is)" >> "$OUT/$stage.jobs"
        log "stage $stage: job $jid ($jobname)"
    else
        local driver rc
        mkdir -p "$REPO/antarctica/results"
        if [ "$kind" = inversion ]; then driver=antarctica/scripts/inversion_icepack2.py
        else driver=$(driver_for "$(fwd_experiment "$stage")"); fi
        echo "mpiexec $(date -Is)" >> "$OUT/$stage.jobs"
        ( cd "$REPO" && eval "$envlines" && unset ISMIP7_CHAIN_ITER_BUDGET ISMIP7_CHAIN_MAX \
          && mpiexec -n "$ISMIP7_ROUTINE_NRANKS" python -u "$driver" ) > "$OUT/$stage.log" 2>&1 \
          && rc=0 || rc=$?
        if [ "$kind" = inversion ]; then
            grep -E "iter +[0-9]+:" "$OUT/$stage.log" | tail -1 | tee -a "$LOGFILE" >&2 || true
        else
            grep -E "^Saved: |Time-stepping:" "$OUT/$stage.log" | tail -2 | tee -a "$LOGFILE" >&2 || true
        fi
        [ "$rc" = 0 ] || die "stage $stage exited $rc; see $OUT/$stage.log"
        is_done "$stage" || die "stage $stage ran but left no finished product; see $OUT/$stage.log"
        log "stage $stage: done"
    fi
}

# ── the state machine ────────────────────────────────────────────────────
LOCKDIR=""
lock_is_stale() {
    # The holder's PID and host are in the lock. On this host the lock is
    # stale when its process is gone; another host's lock (or one whose
    # holder was not yet written) only when it is older than 20 minutes.
    local pid="" host=""
    read -r pid host < "$LOCKDIR/holder" 2>/dev/null || true
    if [ -n "$pid" ] && [ "$host" = "$(hostname)" ]; then
        ! kill -0 "$pid" 2>/dev/null
    else
        [ -n "$(find "$LOCKDIR" -maxdepth 0 -mmin +20 2>/dev/null)" ]
    fi
}
with_lock() {
    local waited=0
    LOCKDIR="$OUT/.routine.lock"
    until mkdir "$LOCKDIR" 2>/dev/null; do
        if lock_is_stale; then rm -rf "$LOCKDIR"; continue; fi
        sleep 5; waited=$((waited + 5))
        [ "$waited" -lt 1200 ] || die "could not take $LOCKDIR in 20 minutes"
    done
    echo "$$ $(hostname)" > "$LOCKDIR/holder"
    trap '[ -n "$LOCKDIR" ] && rm -rf "$LOCKDIR" || true' EXIT
}

next_pass() {   # start every ready, unstarted stage; prints how many
    local stage n=0
    derive_products
    for stage in $STAGE_LIST; do
        is_done "$stage" && continue
        is_started "$stage" && [ "$ISMIP7_ROUTINE_LAUNCH" = slurm ] && continue
        if ready "$stage"; then
            launch "$stage"; n=$((n + 1))
            # a stage that finished in the foreground may have unlocked others
            [ "$ISMIP7_ROUTINE_LAUNCH" = mpiexec ] && derive_products
        fi
    done
    echo "$n"
}

cmd_next() {
    with_lock
    local n
    if [ "$ISMIP7_ROUTINE_LAUNCH" = slurm ]; then
        n=$(next_pass)
        log "next: $n stage(s) submitted"
    else
        while :; do
            n=$(next_pass)
            [ "$n" = 0 ] && break
        done
        log "next: no stage ready"
    fi
}

cmd_status() {
    local stage state where
    echo "routine $TAG in $OUT ($ISMIP7_ROUTINE_LAUNCH); stages:$STAGE_LIST"
    for stage in $STAGE_LIST; do
        if is_inversion "$stage"; then
            where=$(map_path "$stage")
            if is_done "$stage"; then state="done ($(h5attr "$where" objective_iteration) its, total $(h5attr "$where" objective_total))"
            elif [ -f "$where" ]; then state="in progress ($(cat "$where.iters" 2>/dev/null | cut -d' ' -f1 || echo 0)+$(h5attr "$where" objective_iteration) its)"
            elif is_started "$stage"; then state="started, no checkpoint yet"
            elif ready "$stage"; then state="ready"
            else state="waiting"; fi
        else
            where=$(fwd_final "$stage")
            if is_done "$stage"; then state="done (t=$(h5attr "$where" t_yr))"
            elif [ -f "$where" ]; then state="in progress (t=$(h5attr "$where" t_yr) of $(fwd_t_end "$stage"))"
            elif is_started "$stage"; then state="started, no checkpoint yet"
            elif ready "$stage"; then state="ready"
            else state="waiting"; fi
        fi
        [ -f "$OUT/$stage.blocked" ] && state="BLOCKED ($(cat "$OUT/$stage.blocked")); fix, then resume $stage"
        if [ "$ISMIP7_ROUTINE_LAUNCH" = slurm ] && is_started "$stage" && ! is_done "$stage" && command -v squeue >/dev/null; then
            state="$state; queue: $(squeue -h -u "${USER:-$(id -un)}" -n "rt_${TAG}_${stage}" -o '%i:%T' 2>/dev/null | tr '\n' ' ')"
        fi
        printf '  %-20s %s\n' "$stage" "$state"
    done
    [ -f "$OUT/relaxed_state.h5" ] && echo "  relaxed_state.h5     present"
    [ -e "$(forward_map)" ] && echo "  map_forward.h5       -> $(readlink "$(forward_map)")"
    return 0
}

# The watchdog (slurm mode): a chain whose link was killed by the scheduler
# queues no successor, and a link that hangs in a solve sits silent until the
# wall limit. Every ISMIP7_ROUTINE_WATCH_SEC (120) it resubmits a started,
# unfinished stage with no job in the queue, and cancels a RUNNING forward
# link whose log has been silent for ISMIP7_ROUTINE_STALL_MIN (15) minutes,
# noting its node in OUT/exclude_nodes so the resubmission avoids it. The
# inversion links are left alone: an iteration can take minutes and a hung one
# is bounded by the wall limit, after which its queued successor resumes.
watch_pass() {
    local stage jobname rows jid state node elapsed logf age stall
    stall=${ISMIP7_ROUTINE_STALL_MIN:-15}
    # A scheduler that cannot be reached answers nothing, which must not read
    # as "no job": during the NOTS controller outage of 14 October 2026 every
    # stage was resubmitted on top of its live chain. One probe, then the pass.
    if ! squeue -h -u "${USER:-$(id -un)}" -o %i >/dev/null 2>&1; then
        log "watch: squeue failed; skipping this pass"
        return 0
    fi
    for stage in $STAGE_LIST; do
        is_started "$stage" || continue
        is_done "$stage" && continue
        jobname="rt_${TAG}_${stage}"
        rows=$(squeue -h -u "${USER:-$(id -un)}" -n "$jobname" -o "%i %T %N %M %r" 2>/dev/null) \
            || { log "watch: squeue failed; skipping the rest of this pass"; return 0; }
        # A successor that can never start (its predecessor ended TIMEOUT
        # after saving, under an afterok dependency) sits PENDING and made
        # the stage look alive for a day (seq1 control_mri_esm2, 19-20 Oct
        # 2026). Cancel it; the stage then reads as empty and is resubmitted.
        local dead
        dead=$(awk '$2 == "PENDING" && $5 == "DependencyNeverSatisfied" {print $1}' <<< "$rows")
        if [ -n "$dead" ]; then
            log "watch: cancelling $(echo $dead) ($stage): dependency never satisfied"
            # shellcheck disable=SC2086
            scancel $dead
            rows=$(awk '!($2 == "PENDING" && $5 == "DependencyNeverSatisfied")' <<< "$rows")
        fi
        if [ -z "$rows" ]; then
            # A chain ends without a successor in two ways: the scheduler took
            # the link (timeout, preemption, node failure, an outage), which a
            # resubmission resumes, or the driver exited with an error or the
            # solver stalled where the link resumed, which a resubmission
            # repeats. The second parks the stage until `resume`
            # (OUT/<stage>.blocked holds the reason).
            [ -f "$OUT/$stage.blocked" ] && continue
            local why
            why=$(driver_failure "$stage")
            if [ -n "$why" ]; then
                echo "$why" > "$OUT/$stage.blocked"
                log "watch: ERROR stage $stage stopped on a driver failure or a stall, not resubmitting: $why"
                continue
            fi
            # A strict routine (ISMIP7_ROUTINE_STRICT=1) resubmits a link only
            # when the scheduler took it; any other end parks the stage.
            if [ "$ISMIP7_ROUTINE_STRICT" = 1 ]; then
                if ! why=$(strict_block_reason "$stage"); then
                    log "watch: stage $stage has no job in the queue and its newest link has not ended in the accounting; waiting"
                    continue
                fi
                if [ -n "$why" ]; then
                    echo "$why" > "$OUT/$stage.blocked"
                    log "watch: ERROR stage $stage parked (strict): $why"
                    continue
                fi
            fi
            log "watch: stage $stage has no job in the queue; resubmitting"
            mv "$OUT/$stage.jobs" "$OUT/$stage.jobs.$(date +%s)"
            launch "$stage"
            continue
        fi
        # A chain's successor takes its node exclusions from its predecessor
        # (site_core.sh), so a chain submitted before the watchdog excluded a
        # node never learns of the exclusion: successors of such chains hung
        # on the excluded bb6u15g1 and bb6u13g1 on 22 October 2026. Stamp the
        # current list onto any pending link that carries none; the
        # scheduler lets a pending job be changed, a running one not.
        if [ -s "$OUT/exclude_nodes" ] && command -v scontrol >/dev/null; then
            local excl pend cur
            excl=$(sort -u "$OUT/exclude_nodes" | paste -sd,)
            while read -r pend _; do
                [ -n "$pend" ] || continue
                cur=$(scontrol show job "$pend" 2>/dev/null | grep -oE 'ExcNodeList=[^ ]+' | cut -d= -f2)
                [ "$cur" = "(null)" ] || continue
                if scontrol update JobId="$pend" ExcNodeList="$excl" 2>/dev/null; then
                    log "watch: job $pend ($stage) pending without node exclusions; now excluding $excl"
                fi
            done <<< "$(awk '$2 == "PENDING" {print $1}' <<< "$rows")"
        fi
        is_inversion "$stage" && continue
        while read -r jid state node elapsed _; do
            [ "$state" = RUNNING ] || continue
            logf="$REPO/logs/ismip7_fwd_$jid.out"
            [ -f "$logf" ] || continue
            # A momentum solve that "converges" in zero iterations leaves the
            # velocity where the restart put it (the frozen controls of
            # 19-21 October 2026); a healthy link logs none. Cancel the link
            # and park the stage before the frozen state is checkpointed on.
            local frozen
            frozen=$(grep -c "snes_its=0 " "$logf" 2>/dev/null || true)
            if [ "$ISMIP7_ROUTINE_FROZEN_MAX" -gt 0 ] && [ "${frozen:-0}" -ge "$ISMIP7_ROUTINE_FROZEN_MAX" ]; then
                echo "job $jid: $frozen zero-iteration momentum solves in $logf; the velocity is frozen" > "$OUT/$stage.blocked"
                log "watch: ERROR job $jid ($stage) on $node: $frozen zero-iteration momentum solves; cancelling, stage parked"
                scancel "$jid"
                continue
            fi
            age=$(( ($(date +%s) - $(stat -c %Y "$logf")) / 60 ))
            if [ "$age" -ge "$stall" ]; then
                log "watch: job $jid ($stage) on $node silent for $age min; cancelling, node excluded"
                echo "$node" >> "$OUT/exclude_nodes"
                scancel "$jid"
            fi
        done <<< "$rows"
    done
}

# The driver error or solver stall that ended a stage's newest link, or
# nothing when the link ended otherwise (the scheduler took it, it finished,
# or its log is not there). A forward or inversion link logs "chain: stopping,
# driver exited N" or "inversion exit N" with N != 0 when the driver itself
# failed. A forward whose solver stalled exits 0 with its state saved and
# logs "chain: the solver stalled", and its chain does not resubmit it.
# Whether a resubmission helps depends on where the link stalled. A link that
# advanced and then stalled is resubmitted: the next one starts the substep
# controller afresh at the stall, and the refinement it climbs through can
# take the step (on 16 October 2026 a restart of the seq0 ssp585 projection
# passed the step its link had stalled in at 2021.25). A link that stalled
# where it resumed made no progress, and another would stall there again (the
# seq0 ssp126 projection at 2019.95, relaunched link after link): parked.
# The newest link of a stage. OUT/<stage>.jobs records the links this script
# submitted; a chain's own successors (projection.sbatch resubmits itself under
# the same job name) are not in it, so the file's last line is the chain's
# first link, and on 21 October 2026 the watchdog parked seq0f's MRI historical
# naming a link that had ended cleanly two links before the one that aborted.
# The accounting knows every link by name; the file is the fallback where
# sacct is absent or answers nothing.
newest_link() {
    local first since jid
    first=$(head -n 1 "$OUT/$1.jobs")
    since=$(echo "$first" | cut -d' ' -f2 | cut -c1-19)
    jid=$(sacct -X -n -u "${USER:-$(id -un)}" ${since:+-S "$since"} --name="rt_${TAG}_$1" -o JobID 2>/dev/null \
          | awk 'NF {print $1}' | sort -n | tail -n 1)
    [ -n "$jid" ] || jid=$(tail -n 1 "$OUT/$1.jobs" | cut -d' ' -f1)
    echo "$jid"
}

driver_failure() {
    local jid kind logf err at from
    [ -f "$OUT/$1.jobs" ] || return 0
    jid=$(newest_link "$1")
    if is_inversion "$1"; then kind=inv; else kind=fwd; fi
    logf="$REPO/logs/ismip7_${kind}_$jid.out"
    [ -f "$logf" ] || return 0
    if grep -q "^chain: the solver stalled" "$logf"; then
        at=$(sed -n 's/^chain: the solver stalled.* at \([-0-9.eE+]*\)\.$/\1/p' "$logf" | tail -n 1)
        from=$(sed -n 's/.*Restart: evolved geometry.*(t_yr=\([-0-9.eE+]*\),.*/\1/p' "$logf" | tail -n 1)
        if [ -n "$at" ] && [ -n "$from" ] \
           && awk -v a="$at" -v b="$from" 'BEGIN { exit !(a - b < 1e-6 && b - a < 1e-6) }'; then
            echo "job $jid: the solver stalled at $at, where the link resumed"
        fi
        return 0
    fi
    if grep -qE "chain: stopping, driver exited [1-9]|inversion exit [1-9]" "$logf"; then
        err=$(grep -hE "^[A-Za-z]*Error: " "${logf%.out}.err" "$logf" 2>/dev/null | tail -n 1 | cut -c1-300)
        echo "job $jid: ${err:-driver exited non-zero; see $logf}"
    fi
    return 0
}

# The reason a strict routine parks a stage whose newest link ended without a
# successor: the solver stalled anywhere, or the link ended in a state other
# than one the scheduler imposes (TIMEOUT, PREEMPTED, NODE_FAIL, CANCELLED,
# BOOT_FAIL). Nothing when the scheduler took it, which is resubmitted. A
# failure status, and nothing, while the accounting still holds the link live
# (RUNNING, COMPLETING, PENDING, REQUEUED: squeue has dropped it before
# slurmdbd records its end, or a requeue is in flight), so the pass neither
# resubmits nor parks. The
# accounting prints a cancelled link as "CANCELLED+" (its "by <uid>" is cut
# off), which read as a driver end and parked dg1f's ssp126 on 22 October
# 2026 after the watchdog itself had cancelled the hung link.
strict_block_reason() {
    local jid kind logf st at
    [ -f "$OUT/$1.jobs" ] || return 0
    jid=$(newest_link "$1")
    if is_inversion "$1"; then kind=inv; else kind=fwd; fi
    logf="$REPO/logs/ismip7_${kind}_$jid.out"
    if [ -f "$logf" ] && grep -q "^chain: the solver stalled" "$logf"; then
        at=$(sed -n 's/^chain: the solver stalled.* at \([-0-9.eE+]*\)\.$/\1/p' "$logf" | tail -n 1)
        echo "job $jid: the solver stalled at ${at:-?}"
        return 0
    fi
    st=$(sacct -j "$jid" -X -n -o State 2>/dev/null | head -n 1 | awk '{print $1}')
    case "${st:-UNKNOWN}" in
        TIMEOUT|PREEMPTED|NODE_FAIL|CANCELLED*|BOOT_FAIL|UNKNOWN|"") return 0 ;;
        RUNNING|COMPLETING|PENDING|REQUEUED*) return 1 ;;
        *) echo "job $jid: ended $st without a successor; see $logf" ;;
    esac
}

cmd_watch() {
    [ "$ISMIP7_ROUTINE_LAUNCH" = slurm ] || die "watch is for slurm mode"
    local every=${ISMIP7_ROUTINE_WATCH_SEC:-120}
    log "watch: every $every s, stall after ${ISMIP7_ROUTINE_STALL_MIN:-15} min"
    while :; do
        ( with_lock; derive_products; watch_pass; n=$(next_pass); [ "$n" = 0 ] || log "watch: $n stage(s) submitted" ) || log "watch: pass failed ($?)"
        [ "${ISMIP7_ROUTINE_WATCH_ONCE:-0}" = 1 ] && break
        sleep "$every"
    done
}

cmd_resume() {
    local stage=${1:?resume STAGE}
    in_list "$stage" "$STAGE_LIST" || die "$stage is not a stage of this routine"
    is_done "$stage" && die "$stage is done"
    ready "$stage" || die "$stage is not ready"
    with_lock
    [ -f "$OUT/$stage.jobs" ] && mv "$OUT/$stage.jobs" "$OUT/$stage.jobs.$(date +%s)"
    rm -f "$OUT/$stage.blocked"
    launch "$stage"
}

case "$CMD" in
    init)
        [ -f "$ENVFILE" ] && die "$ENVFILE exists; use next, status or resume (remove it to record the environment again)"
        record_env
        load_env
        check_objective || { rm -f "$ENVFILE"; exit 1; }
        log "init: $ENVFILE recorded; stages:$STAGE_LIST"
        cmd_next ;;
    next)   load_env; cmd_next ;;
    status) load_env; cmd_status ;;
    resume) load_env; cmd_resume "${3:-}" ;;
    watch)  load_env; cmd_watch ;;
    *) die "unknown command $CMD (init|next|status|resume STAGE)" ;;
esac
