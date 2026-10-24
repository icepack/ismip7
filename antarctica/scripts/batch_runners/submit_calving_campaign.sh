#!/bin/bash
# A calving campaign for the forward suite: the same chain of experiments run
# once per calving-law member, each member its own run tag, as the Pine
# Island study ran its laws (pine-island-calving/scripts/campaign_projection.sh).
#
#   submit_calving_campaign.sh screen [MEMBER ...]   historical 2003-2015 per member,
#                                                     every yearly state kept, for
#                                                     score_front.py against the Greene record
#   submit_calving_campaign.sh suite  [MEMBER ...]   historical, then the control and
#                                                     the projections chained from it
#
# Run from the forward checkout with the production environment sourced
# (fwd_env_sep1.sh on NOTS); ISMIP7_INVERSION, the mesh and the solver knobs come
# from there. Members not named run all. The member table below is the one
# knob a campaign tunes: each line is name|law|params|front advance; the
# parameters are the Pine Island calibrations against the Greene record
# (coupled model, 2015-2021) where that study fitted them, the submission
# front's own value for hfb, and the published Antarctic values otherwise.
set -euo pipefail
STAGE=${1:?screen or suite}; shift || true
MEMBERS_WANTED="$*"
cd "$(dirname "$0")/../../.."
: "${ISMIP7_INVERSION:?source the production environment first}"
ESM=${ISMIP7_ESM:-CESM2-WACCM}
HIST=${ISMIP7_HIST_EXPERIMENT:-hist_cesm_waccm}
PROJ=${ISMIP7_PROJECTIONS:-"control ssp126_cesm_waccm ssp370_cesm_waccm ssp585_cesm_waccm"}
TAGBASE=${ISMIP7_CAMPAIGN_TAG:-calv}
PART=${ISMIP7_PARTITION:-scavenge}; ACCT=${ISMIP7_ACCOUNT:-commons}; CONSTR=${ISMIP7_CONSTRAINT:-sapphirerapids}
# A link's wall time; the chain's successors inherit it. Scavenge allows
# one hour at most. One step deep in 64x subcycling on 32 Icelake ranks
# runs past the default 8 min wall margin, and a link killed at the limit
# cannot resubmit itself (10 Oct: five of eleven screening members stopped
# that way), so set ISMIP7_WALL_MARGIN_MIN to 20 or so on such links.
WALL=${ISMIP7_WALL:-01:00:00}
MEMBERS="fixed|fixed||none
hfb015|hfb|sigma_max=0.15|none
hfb015free|hfb|sigma_max=0.15|free
hfb0|hfb|sigma_max=0.0|free
vm075|vonmises|sigma_max_fl=0.075,sigma_max_gr=1.0|free
vm150|vonmises|sigma_max_fl=0.15,sigma_max_gr=1.0|free
vm275|vonmises|sigma_max_fl=0.275,sigma_max_gr=1.0|free
vms099|vonmises_strain|sigma_max_fl=0.099,sigma_max_gr=1.0|free
thk450|thickness|hc=450|free
eig1e8|eigencalving|K=1.1e8|free
eig3e8|eigencalving|K=3.0e8|free
nocalv|nocalving||free
vm130|vonmises|sigma_max_fl=0.13,sigma_max_gr=1.0|free
vm130g|vonmises|sigma_max_fl=0.13,sigma_max_gr=0.13|free
hfb086|hfb|sigma_max=0.086|free
eig5e6|eigencalving|K=4.9e6|free
eig1e7|eigencalving|K=9.3e6|free
thk235|thickness|hc=235|free"
# vm130 .. thk235 are the 17 Oct fits of each law's own parameter to the
# flux-balance target (calvingMIP tune_greene.py --target balance on the
# calvdg1h fixed member's 2010 state: the model's stationary-front flux per
# Mouginot basin less the Greene area trend 2003-2015): von Mises 0.13 MPa
# (per-basin and continent-wide agree; vm130 keeps the 1 MPa grounded
# threshold, vm130g applies the fit to both), hfb 0.086 MPa per basin (0.062
# continent-wide), eigencalving K 4.9e6 m yr per basin (9.3e6 continent-wide),
# thickness hc 235 m. The 2004 state gives the same values within 10 %. The
# diagnostic flux at the fixed front ran ~1.7x the flux the calvdg1h forwards
# booked under the same law (an advancing front thins), so vm075 brackets
# the von Mises fit from below.
# hfb015free is the submission front's law let advance (24 Oct 2026 screen:
# does a free front at 0.15 MPa hold the observed fronts, on the extent
# anchor and on the advected front, with ISMIP7_CALVING_FLOATING_ONLY=1 set
# in the environment as the production projections run it).
# `nocalving` with a free front = the front carried by the ice and nothing
# removed: the no-calving baseline every law has to beat. `fixed` is the
# protocol's pinned front. The submission front is hfb015 with FRONT_ADVANCE
# none and the collapse mask (ismip7_front_routine.sh).
mkdir -p logs
echo "campaign $TAGBASE, stage $STAGE, ESM $ESM, MAP $(basename "$ISMIP7_INVERSION")"
git log --oneline -1
while IFS='|' read -r name law params advance; do
    [ -n "$name" ] || continue
    if [ -n "$MEMBERS_WANTED" ] && ! grep -qw "$name" <<< "$MEMBERS_WANTED"; then continue; fi
    tag=${TAGBASE}_${name}
    job=$(
        export ISMIP7_EXPERIMENT="$HIST" ISMIP7_ESM="$ESM" ISMIP7_RUN_TAG="$tag" \
               ISMIP7_CALVING="$law" ISMIP7_CALVING_PARAMS="$params" ISMIP7_FRONT_ADVANCE="$advance" \
               ISMIP7_FRACTURE=none ISMIP7_CHAIN=1 ISMIP7_AUTO_RESUME=1
        if [ "$STAGE" = screen ]; then
            export ISMIP7_KEEP_CHECKPOINTS=0 ISMIP7_CHECKPOINT_EVERY_YR=1 ISMIP7_OUTPUT=0
            unset ISMIP7_CHAIN_THEN
        else
            export ISMIP7_CHAIN_THEN="$PROJ" ISMIP7_OUTPUT=1
        fi
        sbatch --parsable -J "fwd_$tag" -p "$PART" -A "$ACCT" -C "$CONSTR" \
            --nodes=1 --ntasks-per-node=32 --cpus-per-task=1 --hint=nomultithread \
            --mem=240G --time="$WALL" --export=ALL \
            antarctica/scripts/batch_runners/projection.sbatch)
    printf '%-10s %-15s %-36s advance=%-4s -> %s\n' "$name" "$law" "${params:-(defaults)}" "$advance" "$job"
done <<EOF
$MEMBERS
EOF
