# Antarctica timing benchmark: resolution vs. core-count matrix.
#
# Run from antarctica/:
#   make timing
#
# Targets (each idempotent — skips already-built artifacts):
#   meshes        build 10 meshes (5 LC x {10,20}*LC, buffer=20000 m)
#   inversion     verify the imported improved DG0 logvelnet inversion
#   redistribute  rewrite the complete inversion checkpoint on 1 core
#   timing-prepare  prepare/repack one exact-mesh state for all 10 meshes
#   timing-inversion  OPTIONAL per-mesh logvelnet invert (TIMING_INVERT=1 to include it in
#                     `make timing`; TIMING_INVERSION_MAXITER, default 250); republish cache
#                     (FOLLOW_PREPARE=1 queues each invert behind its active prepare job;
#                     500 m meshes are never re-inverted)
#   timing-cache-audit  report the exact initial FV tendency of prepared caches
#   timing-probe    one-mesh contract probe (TIMING_CONTRACT=strict|divfront|front|div);
#                   timing-amb-probe is the legacy alias for TIMING_CONTRACT=div
#   timing-scout    submit the lowest retained core count for each mesh
#   timing-scale    submit only scaling lanes whose exact scout passed
#                   (TIMING_SOLVER=scpc_mumps|scpc_gamg names the solver the lanes time,
#                   MATRIX_STEPS / MATRIX_DT_2500 the interval; all are part of the tag;
#                   TIMING_INITIAL_STATE=invert runs lanes from invert-published caches
#                   under the _reinverted tag; FOLLOW_INVERT=1 queues them behind an
#                   active invert job)
#   transient       direct qualification/debug compatibility entry point
#   qualify       gate one solver through 2-step and 5-step 2.5 km probes
#   timestep-probe restart at 2.5 km for 5 steps at a selectable coarser dt
#   debug         one condensed-MUMPS job with full convergence monitors
#   reference     one full-Jacobian MUMPS reference job
#   solver-smoke  tiny two-rank setup test; no Antarctic data or time steps
#   matrix        write TIMING_MATRIX.md from results/timing/*.json
#                   (MATRIX_OUTPUT names another file, e.g. one per site and
#                   solver: TIMING_MATRIX_QUARTZ_SCPC_GAMG.md)
#   map-check     one released MAP through its checks (MAP_CHECK.md): fetch and
#                 verify, score at t=0 on its own mesh and transferred onto the
#                 production mesh, strict 10-step lanes on both, 10-year
#                 controls on both, audit. Prints the stage table and submits
#                 the first runnable stage (MAP_CHECK_STAGES=all|<stage>);
#                 MAP_CHECK_FRICTION=regularized_coulomb|budd picks the MAP.
#   map-check-dry-run  the stage table and every stage's composed sbatch line

LCS ?= 500 1000 2000 2500 5000
RATIOS ?= 20 10
CORES ?= 16 32 64
BUFFER ?= 20000
INV_LC ?= 2500
INV_LC_COARSE ?= 25000
FRICTION ?= budd
N_FLOW ?= 3.0
# The solver a campaign's lanes time: scpc_mumps (condensed MUMPS, the
# default) or scpc_gamg (the same condensation under GAMG). It leads the tag,
# so the two campaigns' records never mix, and it is the solver `make timing`
# wants qualified. The initial-state caches are scpc_mumps's under either
# (timing_campaign.CACHE_SOLVER_MODE): a GAMG campaign needs no prepare of its
# own and starts every lane where the MUMPS lane started. Must be one of
# timing_campaign.LANE_SOLVER_MODES.
TIMING_SOLVER ?= scpc_mumps
# The solver of the direct lanes (qualify, timestep-probe, transient-direct),
# which may be any mode solverconfig knows; it follows TIMING_SOLVER, so
# `make qualify TIMING_SOLVER=scpc_gamg` gates the solver the campaign times.
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER ?= $(TIMING_SOLVER)
ISMIP7_SOLVER_VIEW ?=
ISMIP7_RESCUE_ENABLED ?= 0
T_END ?= 2020
DT ?= 1.0
# 2026-09-16 dt ladder on 2500/25000 (strict contract, pristine prepare state):
# 5 x 0.25 ran away at step 2, 10 x 0.125 and 20 x 0.0625 completed 1.25 yr with
# flat speed and thickness. 10 x 0.125 is the cheapest rung that passes.
MATRIX_STEPS ?= 10
MATRIX_DT_2500 ?= 0.125
MATRIX_DT_2500_TAG := $(subst .,p,$(MATRIX_DT_2500))
# Absolute cap on the LC-scaled step (2026-09-17: the 5000 m scouts ran away at
# dt 0.25 exactly as 2500/25000 does; the stable dt does not grow with LC).
# Not a campaign parameter -- it is not in the tag -- and it never undercuts
# MATRIX_DT_2500 itself. Must match timing_campaign.MATRIX_DT_MAX.
MATRIX_DT_MAX := 0.125
TIMING_KIND ?= matrix
# Lane physics contract: strict (no apparent MB, no calving sink) or
# divfront / front / div (production closure and/or fixed 2015 calving front).
# The solver, steps, the 2.5 km dt and a non-strict contract are all encoded
# in the tag (scripts/timing_campaign.py campaign_tag builds the identical
# string).
# v4 (2026-09-16): the Budd shelf gate tests height above flotation, not the
# sign of roundoff N (ported from hoffmaao/antarctica e602705). v3 caches and
# records are the old-law archive; the manager never reuses them.
TIMING_CONTRACT ?= strict
TIMING_TAG ?= $(TIMING_SOLVER)_$(MATRIX_STEPS)step_dt$(MATRIX_DT_2500_TAG)at2500_dg0_logvelnet_cached_strict_v4$(if $(filter-out strict,$(TIMING_CONTRACT)),_$(TIMING_CONTRACT),)
TIMING_EXPERIMENT ?= timing
# Optional exact-input overrides for a single transient lane.  The resolution
# matrix leaves these empty and constructs each target mesh parametrically.
TIMING_LC_COARSE ?=
TIMING_MESH ?=
TIMING_BNDIDS ?=
# Campaign source MAP (promoted 2026-09-17): the 250-iteration re-inversion
# of the imported dg0_logvelnet MAP on the 2500/25000 campaign mesh under the
# HAF-gated Budd law (quartz job 10489636), i.e. the invert stage's own output
# for that mesh; the manager never re-inverts or re-prepares it. The
# .parallel.h5 is the 16-rank original, the .h5 its one-rank repack written
# by the same job, which then deletes the original: the repack on its own is
# a complete source (`inversion` passes, `redistribute` skips), and the
# original is asked for only where there is no repack yet. Must agree with
# timing_campaign.SOURCE_INVERSION_BASENAME. The imported old-gate MAP stays
# the solver-qualification input (SOLVER_TEST_INVERSION) and census control.
SOURCE_TIMING_INVERSION ?= $(CURDIR)/results/timing/inversion/inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.parallel.h5
TIMING_INVERSION ?= $(CURDIR)/results/timing/inversion/inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5
TIMING_CACHE_DIR ?= $(CURDIR)/results/timing/cache
TIMING_RESTART ?=
TIMING_ONLY_MESH ?=
TIMING_SCOUT_MONITOR ?=
TIMING_INVERSION_MAXITER ?= 250
TIMING_INVERSION_TIME ?= 24:00:00
FOLLOW_PREPARE ?=
FOLLOW_INVERT ?=
# prepare: campaign lanes from the transferred prepare cache (default).
# invert: lanes from each mesh's invert-published cache (_reinverted tag).
TIMING_INITIAL_STATE ?= prepare
# The per-mesh invert stage is optional; TIMING_INVERT=1 runs it in `make timing`.
TIMING_INVERT ?=
SBATCH_WAIT ?=
FORCE_TIMING ?=
# Every job goes through scripts/batch_runners/submit.sh, which reads this
# cluster's site file (scripts/batch_runners/sites/, chosen by hostname or
# ISMIP7_SITE, with per-user settings in sites/local.env), so these targets are
# the same command at every site. SLURM_QUEUE is the site's partition by class
# (short, long or debug); SLURM_PARTITION and SLURM_CONSTRAINT name one
# outright and are empty unless you set them.
SLURM_QUEUE ?= short
SLURM_PARTITION ?=
SLURM_CONSTRAINT ?=
SLURM_TIME ?= 12:00:00
SLURM_SELECT = --queue "$(SLURM_QUEUE)" \
  $(if $(SLURM_PARTITION),--partition "$(SLURM_PARTITION)",) \
  $(if $(SLURM_CONSTRAINT),--constraint "$(SLURM_CONSTRAINT)",)
# submit.sh submits the checkout this Makefile is in, from this directory.
SUBMIT = ISMIP7_REPO="$(abspath $(CURDIR)/..)" scripts/batch_runners/submit.sh script
SUBMIT_CD = --cd "$(notdir $(CURDIR))"
# `make sync-results` pulls a cluster's results/ into this checkout's:
# REMOTE_RESULTS is an rsync source, host:path/ with the trailing slash.
QUARTZ_RESULTS ?= quartz:/N/project/ice_rheology/ISMIP7/antarctica/results/
REMOTE_RESULTS ?= $(QUARTZ_RESULTS)
LEGACY_TIMING_TAG ?= scpc_mumps_5step_dt0p25at2500
MATRIX_OUTPUT ?= TIMING_MATRIX.md
# `make matrix` renders the campaign `make timing` last launched, unless the
# command line (or environment) names one: the tag itself or any parameter of
# it, e.g. `make matrix TIMING_SOLVER=scpc_gamg` after a scouts-only start.
MATRIX_TAG_NAMED := $(filter-out file,$(foreach v,TIMING_TAG TIMING_SOLVER MATRIX_STEPS MATRIX_DT_2500 TIMING_CONTRACT,$(origin $(v))))

# One released MAP through its checks (scripts/manage_map_check.py). The
# 22 September 2 km snapshots are the defaults; MAP_CHECK_FRICTION picks the
# law and with it the file. The target is the production mesh; on Quartz name
# MAP_CHECK_TARGET_MESH explicitly from a scratch clone. The controls melt with
# the tracked calibration; MAP_CHECK_K_NPZ names a legacy per-basin K instead.
MAP_CHECK_FRICTION ?= regularized_coulomb
MAP_CHECK_RELEASE ?= maps-2km-snap-2026-09-22
MAP_CHECK_SNAP ?= snap20260922_0948
MAP_CHECK_MAP_BASENAME ?= inversion_icepack2_$(if $(filter regularized_coulomb,$(MAP_CHECK_FRICTION)),rc,budd)_n3_dg0_logvelnet_2000_int5000_bilap_$(MAP_CHECK_SNAP).h5
MAP_CHECK_URL ?= https://github.com/icepack/ismip7/releases/download/$(MAP_CHECK_RELEASE)/$(MAP_CHECK_MAP_BASENAME)
MAP_CHECK_MD5 ?= auto
MAP_CHECK_DIR ?= $(CURDIR)/results/map_check
MAP_CHECK_MAP ?= $(MAP_CHECK_DIR)/maps/$(MAP_CHECK_MAP_BASENAME)
MAP_CHECK_LC ?= 2000
MAP_CHECK_LC_COARSE ?= 5000
MAP_CHECK_BUFFER ?= 0
MAP_CHECK_TARGET_LC ?= 1000
MAP_CHECK_TARGET_LC_COARSE ?= 10000
MAP_CHECK_TARGET_BUFFER ?= 20000
MAP_CHECK_TARGET_MESH ?= $(CURDIR)/mesh/antarctica_$(MAP_CHECK_TARGET_LC_COARSE)_$(MAP_CHECK_TARGET_LC)_buffered$(MAP_CHECK_TARGET_BUFFER).msh
MAP_CHECK_STAGES ?= next
MAP_CHECK_SOLVER ?= scpc_gamg
MAP_CHECK_CORES ?= 64
MAP_CHECK_DIAG_CORES ?= 16
MAP_CHECK_T_END ?= 2025
MAP_CHECK_K_NPZ ?=
# MAP_CHECK_ATTRS=0 skips the h5 attribute check of the fetched file (tests);
# MAP_CHECK_REPACK=0 reads the fetched file as is. MAP_CHECK_CONTROLS_AFTER_FAILED_LANE=1
# lets a 10-year control run although the strict lane on its mesh failed: the
# control's production configuration (apparent mass balance, fixed front) is
# what a submission run uses, and its cost and drift are wanted even then.
MAP_CHECK_ATTRS ?= 1
MAP_CHECK_REPACK ?= 1
MAP_CHECK_CONTROLS_AFTER_FAILED_LANE ?=
# The controls' own request. The default is the production forward's; a
# smaller one backfills sooner on a busy partition, and a control that runs
# out of wall time chains itself, so a short request costs a restart and no
# result.
MAP_CHECK_CONTROL_MEM ?= 240G
MAP_CHECK_CONTROL_TIME ?= 12:00:00
MAP_CHECK_ARGS = run --root "$(CURDIR)" --map "$(MAP_CHECK_MAP)" --url "$(MAP_CHECK_URL)" \
  --md5 "$(MAP_CHECK_MD5)" --friction "$(MAP_CHECK_FRICTION)" \
  --lc $(MAP_CHECK_LC) --lc-coarse $(MAP_CHECK_LC_COARSE) --buffer-m $(MAP_CHECK_BUFFER) \
  --target-mesh "$(MAP_CHECK_TARGET_MESH)" --target-lc $(MAP_CHECK_TARGET_LC) \
  --target-lc-coarse $(MAP_CHECK_TARGET_LC_COARSE) --target-buffer-m $(MAP_CHECK_TARGET_BUFFER) \
  --root-dir "$(MAP_CHECK_DIR)" --stages "$(MAP_CHECK_STAGES)" --solver "$(MAP_CHECK_SOLVER)" \
  --cores $(MAP_CHECK_CORES) --diag-cores $(MAP_CHECK_DIAG_CORES) --t-end $(MAP_CHECK_T_END) \
  --check-map-attrs $(MAP_CHECK_ATTRS) --repack $(MAP_CHECK_REPACK) \
  $(if $(MAP_CHECK_K_NPZ),--k-npz "$(MAP_CHECK_K_NPZ)",) \
  $(if $(MAP_CHECK_CONTROLS_AFTER_FAILED_LANE),--controls-after-failed-lane,) \
  --control-mem "$(MAP_CHECK_CONTROL_MEM)" --control-time "$(MAP_CHECK_CONTROL_TIME)" \
  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" $(if $(FORCE_TIMING),--force,)

# The active solver qualification/reference lane uses the improved DG0
# log-velocity/net-balance MAP imported from the other branch and its exact
# source mesh.  Keep these paths together: MAPs and meshes are not
# interchangeable, and the per-mesh boundary sidecar is part of the pair.
SOLVER_TEST_LC ?= 2500
SOLVER_TEST_LC_COARSE ?= 64000
SOLVER_TEST_BUFFER ?= 20000
SOLVER_TEST_TAG ?= dg0_logvelnet
SOLVER_TEST_MESH ?= $(CURDIR)/mesh/antarctica_64000_2500.msh
SOLVER_TEST_BNDIDS ?= $(CURDIR)/mesh/boundary_ids_antarctica_64000_2500.json
# Pinned to the imported MAP on its own 64000/2500 mesh: these lanes pair
# mesh and MAP exactly, and the campaign source now lives on 25000/2500.
SOLVER_TEST_INVERSION ?= $(CURDIR)/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5
SOLVER_TEST_RESTART ?= $(CURDIR)/results/qualification_2step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_final.h5
TIMESTEP_PROBE_SOLVER ?= scpc_mumps
TIMESTEP_PROBE_DT ?= 0.25
TIMESTEP_PROBE_STEPS ?= 5
TIMESTEP_PROBE_END = $(shell awk 'BEGIN { printf "%.12g", 2015.2 + $(TIMESTEP_PROBE_STEPS) * $(TIMESTEP_PROBE_DT) }')
TIMESTEP_PROBE_DT_TAG = $(subst .,p,$(TIMESTEP_PROBE_DT))
TIMESTEP_PROBE_TAG = timestep_probe_$(TIMESTEP_PROBE_STEPS)step_dt$(TIMESTEP_PROBE_DT_TAG)_$(TIMESTEP_PROBE_SOLVER)_$(SOLVER_TEST_TAG)
TIMESTEP_PROBE_RESTART ?= $(CURDIR)/results/qualification_2step_$(TIMESTEP_PROBE_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_final.h5
ifeq ($(N_FLOW),4)
A4_FACTOR ?= 10.0
else ifeq ($(N_FLOW),4.0)
A4_FACTOR ?= 10.0
else
A4_FACTOR ?= 1.0
endif

ifeq ($(FRICTION),regularized_coulomb)
MAP_FRICTION_TAG := _rc
else ifeq ($(FRICTION),budd)
MAP_FRICTION_TAG := _budd
else
MAP_FRICTION_TAG :=
endif

ifeq ($(N_FLOW),4)
MAP_N_TAG :=
else ifeq ($(N_FLOW),4.0)
MAP_N_TAG :=
else ifeq ($(N_FLOW),3)
MAP_N_TAG := _n3
else ifeq ($(N_FLOW),3.0)
MAP_N_TAG := _n3
else
MAP_N_TAG := _n$(N_FLOW)
endif

INV_CHK := mesh/inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5
INV_CHK_MC := mesh/inversion_icepack2_budd_n3_dg0_logvelnet_2500.h5

empty :=
space := $(empty) $(empty)
colon := :
TIMING_LCS := $(subst $(space),$(colon),$(LCS))
TIMING_RATIOS := $(subst $(space),$(colon),$(RATIOS))
TIMING_CORES := $(subst $(space),$(colon),$(CORES))
DEBUG_STAMP := $(shell date +%Y%m%d_%H%M%S)

# Total-memory requests for timing jobs. These are keyed to the fine
# resolution: the 500 m case is the memory outlier, while the 2500 m debug
# case uses MUMPS and has a much smaller measured footprint. The ratio-10
# meshes are slightly larger than their ratio-20 counterparts, so one
# conservative budget per fine resolution covers both ratios.

.PHONY: timing meshes inversion redistribute transient transient-direct \
	qualify qualify-if-needed qualify-2step qualify-5step timestep-probe debug \
	reference solver-smoke matrix matrix-legacy sync-results timing-prepare timing-scout \
	timing-scale timing-dry-run timing-cache-smoke timing-cache-audit \
	timing-probe timing-amb-probe timing-inversion now now-check \
	map-check map-check-dry-run runlog runlog-check

timing: qualify-if-needed
	@mkdir -p results/timing
	@printf '%s\n' "$(TIMING_TAG)" > results/timing/latest_matrix_campaign.txt
	$(MAKE) timing-prepare
	$(if $(or $(TIMING_INVERT),$(filter invert,$(TIMING_INITIAL_STATE))),$(MAKE) timing-inversion,@echo "per-mesh invert stage skipped (TIMING_INVERT=1 runs it)")
	$(MAKE) timing-scout
	$(MAKE) timing-scale
	$(MAKE) matrix
	@echo "Timing state advanced. Re-run 'make timing' after the active stage settles."

meshes:
	@if command -v sbatch >/dev/null 2>&1; then \
	  $(SUBMIT) scripts/batch_runners/timing_meshes.script $(SUBMIT_CD) --wait \
	    --name ant_timing_mesh --tasks 1 --mem 32G --time 08:00:00 $(SLURM_SELECT) \
	    TIMING_LCS=$(TIMING_LCS) TIMING_RATIOS=$(TIMING_RATIOS) TIMING_BUFFER=$(BUFFER); \
	else \
	  for lc in $(LCS); do for r in $(RATIOS); do \
	    lcc=$$((lc * r)); fn=mesh/antarctica_$${lcc}_$${lc}_buffered$(BUFFER).msh; bnd=mesh/boundary_ids_antarctica_$${lcc}_$${lc}_buffered$(BUFFER).json; \
	    [ -f $$fn ] && [ -f $$bnd ] && { echo "skip $$fn"; continue; }; \
	    python scripts/mesh_antarctica.py --lc $$lc --lc-coarse $$lcc --buffer-m $(BUFFER); \
	  done; done; \
	fi

inversion:
	@if [ -f "$(SOURCE_TIMING_INVERSION)" ]; then \
	  echo "Using imported DG0 logvelnet inversion: $(SOURCE_TIMING_INVERSION)"; \
	elif [ -f "$(TIMING_INVERSION)" ]; then \
	  echo "Using one-rank timing inversion: $(TIMING_INVERSION) (its parallel original is gone)"; \
	else \
	  echo "missing timing inversion: neither $(TIMING_INVERSION)" >&2; \
	  echo "  nor its parallel original $(SOURCE_TIMING_INVERSION)" >&2; \
	  exit 2; \
	fi

# A repack whose original was deleted is up to date. Said outright rather than
# left to `-nt`, which bash calls true and dash false when the right-hand file
# does not exist.
redistribute: inversion
	@if [ -f "$(TIMING_INVERSION)" ] && { [ ! -f "$(SOURCE_TIMING_INVERSION)" ] || \
	    [ "$(TIMING_INVERSION)" -nt "$(SOURCE_TIMING_INVERSION)" ]; }; then \
	  echo "skip $(TIMING_INVERSION) (up to date)"; \
	elif command -v sbatch >/dev/null 2>&1; then \
	  $(SUBMIT) scripts/batch_runners/timing_redistribute.script $(SUBMIT_CD) --wait \
	    --name ant_timing_redist --tasks 1 --mem 16G --time 01:00:00 $(SLURM_SELECT) \
	    TIMING_REDISTRIBUTE_INPUT="$(SOURCE_TIMING_INVERSION)" \
	    TIMING_REDISTRIBUTE_OUTPUT="$(TIMING_INVERSION)"; \
	else \
	  python scripts/redistribute_checkpoint.py \
	    --input "$(SOURCE_TIMING_INVERSION)" --output "$(TIMING_INVERSION)"; \
	fi

transient:
	@if [ "$(TIMING_KIND)" = matrix ]; then \
	  $(MAKE) timing-scout; \
	  $(MAKE) timing-scale; \
	else \
	  $(MAKE) transient-direct TIMING_KIND="$(TIMING_KIND)"; \
	fi

transient-direct:
	@mkdir -p results/timing results/logs
	@if [ "$(TIMING_KIND)" = matrix ]; then \
	  printf '%s\n' "$(TIMING_TAG)" > results/timing/latest_matrix_campaign.txt; \
	fi
	@submit_fail=0; \
	for lc in $(LCS); do for r in $(RATIOS); do lcc=$$((lc * r)); for nc in $(CORES); do \
	  if [ -n "$(TIMING_LC_COARSE)" ]; then lcc=$(TIMING_LC_COARSE); fi; \
	  run_dt="$(DT)"; run_t_end="$(T_END)"; \
	  if [ "$(TIMING_KIND)" = matrix ]; then \
	    run_dt=$$(awk -v lc="$$lc" -v base="$(MATRIX_DT_2500)" -v cap="$(MATRIX_DT_MAX)" 'BEGIN { dt = base * lc / 2500; if (base + 0 > cap + 0) cap = base; if (dt > cap + 0) dt = cap; printf "%.12g", dt }'); \
	    run_t_end=$$(awk -v dt="$$run_dt" -v steps="$(MATRIX_STEPS)" 'BEGIN { printf "%.12g", 2015 + steps * dt }'); \
	  fi; \
	  run_experiment="$(TIMING_EXPERIMENT)"; \
	  if [ "$(TIMING_KIND)" = matrix ]; then \
	    run_experiment="$${run_experiment}_$(TIMING_TAG)_lcc$${lcc}_n$${nc}"; \
	  fi; \
	  out=$(CURDIR)/results/timing/timing_$(TIMING_TAG)_$${lc}_$${lcc}_$${nc}.json; \
	  status=$(CURDIR)/results/timing/status_$(TIMING_TAG)_$${lc}_$${lcc}_$${nc}.txt; \
	  status_state=; \
	  if [ -f "$$status" ]; then read -r status_state _ < "$$status"; fi; \
	  [ -z "$(FORCE_TIMING)" ] && { [ "$$status_state" = submitted ] || [ "$$status_state" = running ]; } && { echo "skip $$out ($$status_state stamp)"; continue; }; \
	  [ -z "$(FORCE_TIMING)" ] && [ -f "$$out" ] && [ "$$status_state" = finished ] && { echo "skip $$out (finished stamp)"; continue; }; \
	  case "$$lc" in \
	    500)  mem=240G;; \
	    1000) mem=96G;; \
	    2000) mem=80G;; \
	    2500) mem=64G;; \
	    5000) mem=32G;; \
	    *) echo "no timing memory budget for LC=$$lc" >&2; exit 2;; \
	  esac; \
	  mesh="$(TIMING_MESH)"; \
	  [ -n "$$mesh" ] || mesh=$(CURDIR)/mesh/antarctica_$${lcc}_$${lc}_buffered$(BUFFER).msh; \
	  bndids="$(TIMING_BNDIDS)"; \
	  [ -n "$$bndids" ] || bndids=$(CURDIR)/mesh/boundary_ids_antarctica_$${lcc}_$${lc}_buffered$(BUFFER).json; \
	  inversion="$(TIMING_INVERSION)"; \
	  restart="$(TIMING_RESTART)"; \
	  if [ ! -f "$$mesh" ]; then \
	    echo "mesh not found: $$mesh" >&2; printf 'not_runnable reason=missing_mesh path=%s\n' "$$mesh" > "$$status"; submit_fail=1; continue; \
	  fi; \
	  if [ ! -f "$$bndids" ]; then \
	    echo "boundary ids not found: $$bndids" >&2; printf 'not_runnable reason=missing_boundary_ids path=%s\n' "$$bndids" > "$$status"; submit_fail=1; continue; \
	  fi; \
	  if [ ! -f "$$inversion" ]; then \
	    echo "inversion not found: $$inversion" >&2; printf 'not_runnable reason=missing_inversion path=%s\n' "$$inversion" > "$$status"; submit_fail=1; continue; \
	  fi; \
	  if [ -n "$$restart" ] && [ ! -f "$$restart" ]; then \
	    echo "restart not found: $$restart" >&2; printf 'not_runnable reason=missing_restart path=%s\n' "$$restart" > "$$status"; submit_fail=1; continue; \
	  fi; \
	  if command -v sbatch >/dev/null 2>&1; then \
	    printf 'submitted timestamp=%s\n' "$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$$status"; \
	    submit_err=$$(mktemp); \
	    job_id=$$($(SUBMIT) scripts/batch_runners/timing_transient.script $(SUBMIT_CD) $(SBATCH_WAIT) \
	      --name timing_$${lc}_$${lcc}_$${nc} --tasks $$nc --mem $$mem --time $(SLURM_TIME) $(SLURM_SELECT) \
	      ISMIP7_LC=$$lc ISMIP7_LC_COARSE=$$lcc ISMIP7_BUFFER_M=$(BUFFER) "ISMIP7_MESH=$$mesh" "ISMIP7_BNDIDS=$$bndids" \
	      "ISMIP7_INVERSION=$$inversion" "ISMIP7_RESTART=$$restart" ISMIP7_FRICTION=$(FRICTION) ISMIP7_N_FLOW=$(N_FLOW) \
	      ISMIP7_A4_FACTOR=$(A4_FACTOR) ISMIP7_T_END=$$run_t_end ISMIP7_DT=$$run_dt \
	      ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER) ISMIP7_RESCUE_ENABLED=$(ISMIP7_RESCUE_ENABLED) \
	      ISMIP7_TIMING_KIND=$(TIMING_KIND) ISMIP7_TIMING_TAG=$(TIMING_TAG) ISMIP7_TIMING_EXPERIMENT=$$run_experiment \
	      "ISMIP7_TIMING_STATUS=$$status" ISMIP7_SNES_MONITOR=$(ISMIP7_SNES_MONITOR) "ISMIP7_SNES_LOG=$(ISMIP7_SNES_LOG)" \
	      ISMIP7_SOLVER_VIEW=$(ISMIP7_SOLVER_VIEW) 2> "$$submit_err"); \
	    submit_rc=$$?; cat "$$submit_err" >&2; \
	    if [ $$submit_rc -eq 0 ]; then \
	      job_id=$${job_id%%;*}; \
	      read -r status_state _ < "$$status"; \
	      if [ "$$status_state" = submitted ]; then \
	        printf 'submitted job_id=%s timestamp=%s\n' "$$job_id" "$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$$status"; \
	      fi; \
	      echo "submitted timing_$${lc}_$${lcc}_$${nc} as job $$job_id"; \
	    elif [ $$submit_rc -eq 3 ]; then \
	      reason=$$(sed -n 's/.* reason=\([^ ]*\).*/\1/p' "$$submit_err" | head -1); \
	      printf 'not_runnable reason=%s timestamp=%s\n' "$${reason:-exceeds_site_node}" "$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$$status"; \
	    else \
	      read -r status_state _ < "$$status"; \
	      if [ "$$status_state" = submitted ]; then \
	        printf 'submission_failed timestamp=%s\n' "$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$$status"; \
	      fi; \
	      submit_fail=1; \
	    fi; \
	    rm -f "$$submit_err"; \
	  else \
	    printf 'running timestamp=%s\n' "$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$$status"; \
	    if OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ISMIP7_LC=$$lc ISMIP7_LC_COARSE=$$lcc ISMIP7_BUFFER_M=$(BUFFER) ISMIP7_MESH=$$mesh ISMIP7_BNDIDS=$$bndids ISMIP7_INVERSION=$$inversion ISMIP7_RESTART=$$restart ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER) ISMIP7_RESCUE_ENABLED=$(ISMIP7_RESCUE_ENABLED) ISMIP7_TIMING_KIND=$(TIMING_KIND) ISMIP7_TIMING_TAG=$(TIMING_TAG) ISMIP7_TIMING_EXPERIMENT=$$run_experiment ISMIP7_SNES_MONITOR=$(ISMIP7_SNES_MONITOR) ISMIP7_SNES_LOG=$(ISMIP7_SNES_LOG) ISMIP7_SOLVER_VIEW=$(ISMIP7_SOLVER_VIEW) \
	      ISMIP7_FRICTION=$(FRICTION) ISMIP7_N_FLOW=$(N_FLOW) ISMIP7_A4_FACTOR=$(A4_FACTOR) ISMIP7_T_END=$$run_t_end ISMIP7_DT=$$run_dt \
	      mpiexec -n $$nc python scripts/run_timing.py; then \
	      printf 'finished exit_code=0 timestamp=%s\n' "$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$$status"; \
	    else \
	      rc=$$?; \
	      read -r status_state _ < "$$status" 2>/dev/null || status_state=; \
	      if [ "$$status_state" != failed ]; then \
	        printf 'failed category=external_termination exit_code=%s timestamp=%s\n' "$$rc" "$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$$status"; \
	      fi; \
	      submit_fail=1; \
	    fi; \
	  fi; \
	done; done; done; \
	exit $$submit_fail

qualify-if-needed: redistribute
	@if python scripts/check_solver_qualification.py \
	    --record results/timing/timing_qualification_5step_$(TIMING_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_$(SOLVER_TEST_LC_COARSE)_16.json \
	    --csv results/qualification_5step_$(TIMING_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_timeseries.csv \
	    --min-steps 5 >/dev/null 2>&1; then \
	  echo "$(TIMING_SOLVER) qualification already passes; skip"; \
	else \
	  $(MAKE) qualify ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(TIMING_SOLVER); \
	fi

timing-prepare:
	python scripts/manage_timing_campaign.py prepare \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  $(if $(FORCE_TIMING),--force,)

timing-inversion:
	python scripts/manage_timing_campaign.py invert \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(TIMING_INVERSION_TIME)" \
	  --maxiter "$(TIMING_INVERSION_MAXITER)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  $(if $(FORCE_TIMING),--force,) \
	  $(if $(FOLLOW_PREPARE),--follow-prepare,)

timing-cache-audit:
	python scripts/manage_timing_campaign.py audit \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  $(if $(FORCE_TIMING),--force,)

# One-mesh contract probe (kind cache_probe, tag <campaign>_probe): the
# experiment ladder. TIMING_CONTRACT selects strict/divfront/front/div.
timing-probe:
	python scripts/manage_timing_campaign.py probe \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  --matrix-steps "$(MATRIX_STEPS)" --matrix-dt-2500 "$(MATRIX_DT_2500)" \
	  --solver "$(TIMING_SOLVER)" --contract "$(TIMING_CONTRACT)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  $(if $(TIMING_SCOUT_MONITOR),--monitor,) \
	  $(if $(FORCE_TIMING),--force,)

# Legacy alias: the apparent-MB probe is the `div` contract probe.
timing-amb-probe:
	$(MAKE) timing-probe TIMING_CONTRACT=div

timing-scout:
	python scripts/manage_timing_campaign.py scout \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  --maxiter "$(TIMING_INVERSION_MAXITER)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  $(if $(TIMING_SCOUT_MONITOR),--monitor,) \
	  $(if $(FORCE_TIMING),--force,) \
	  $(if $(FOLLOW_INVERT),--follow-invert,) \
	  --initial-state "$(TIMING_INITIAL_STATE)" \
	  --matrix-steps "$(MATRIX_STEPS)" --matrix-dt-2500 "$(MATRIX_DT_2500)" \
	  --solver "$(TIMING_SOLVER)" --contract "$(TIMING_CONTRACT)"

timing-scale:
	python scripts/manage_timing_campaign.py scale \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  --maxiter "$(TIMING_INVERSION_MAXITER)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  $(if $(FORCE_TIMING),--force,) \
	  $(if $(FOLLOW_INVERT),--follow-invert,) \
	  --initial-state "$(TIMING_INITIAL_STATE)" \
	  --matrix-steps "$(MATRIX_STEPS)" --matrix-dt-2500 "$(MATRIX_DT_2500)" \
	  --solver "$(TIMING_SOLVER)" --contract "$(TIMING_CONTRACT)"

timing-dry-run:
	python scripts/manage_timing_campaign.py prepare --dry-run \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",)
	python scripts/manage_timing_campaign.py invert --dry-run --assume-valid-caches \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(TIMING_INVERSION_TIME)" \
	  --maxiter "$(TIMING_INVERSION_MAXITER)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",)
	python scripts/manage_timing_campaign.py scout --dry-run --assume-valid-caches \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  --maxiter "$(TIMING_INVERSION_MAXITER)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  $(if $(TIMING_SCOUT_MONITOR),--monitor,) \
	  --initial-state "$(TIMING_INITIAL_STATE)" \
	  --matrix-steps "$(MATRIX_STEPS)" --matrix-dt-2500 "$(MATRIX_DT_2500)" \
	  --solver "$(TIMING_SOLVER)" --contract "$(TIMING_CONTRACT)"
	python scripts/manage_timing_campaign.py scale --dry-run --assume-valid-caches \
	  --inversion "$(TIMING_INVERSION)" --cache-dir "$(TIMING_CACHE_DIR)" \
	  $(SLURM_SELECT) --walltime "$(SLURM_TIME)" \
	  --maxiter "$(TIMING_INVERSION_MAXITER)" \
	  $(if $(TIMING_ONLY_MESH),--only-mesh "$(TIMING_ONLY_MESH)",) \
	  --initial-state "$(TIMING_INITIAL_STATE)" \
	  --matrix-steps "$(MATRIX_STEPS)" --matrix-dt-2500 "$(MATRIX_DT_2500)" \
	  --solver "$(TIMING_SOLVER)" --contract "$(TIMING_CONTRACT)"

# Qualification is intentionally a single 2.5 km / 16-rank lane on the exact
# 64 km-coarse mesh used by the imported DG0 logvelnet MAP.  The checker gates
# requested versus completed step count, final year, convergence reasons and
# the mass-budget residual before the five-step probe starts.

qualify: redistribute
	$(MAKE) qualify-2step
	$(MAKE) qualify-5step

qualify-2step:
	$(MAKE) transient LCS=$(SOLVER_TEST_LC) RATIOS=1 CORES=16 \
	  BUFFER=$(SOLVER_TEST_BUFFER) TIMING_LC_COARSE=$(SOLVER_TEST_LC_COARSE) \
	  TIMING_MESH=$(SOLVER_TEST_MESH) TIMING_BNDIDS=$(SOLVER_TEST_BNDIDS) \
	  TIMING_INVERSION=$(SOLVER_TEST_INVERSION) \
	  ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER) \
	  T_END=2015.2 DT=0.1 TIMING_KIND=qualification \
	  TIMING_TAG=qualification_2step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG) \
	  TIMING_EXPERIMENT=qualification_2step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG) \
	  SBATCH_WAIT=--wait FORCE_TIMING=1
	python scripts/check_solver_qualification.py \
	  --record results/timing/timing_qualification_2step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_$(SOLVER_TEST_LC_COARSE)_16.json \
	  --csv results/qualification_2step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_timeseries.csv --min-steps 2

qualify-5step:
	$(MAKE) transient LCS=$(SOLVER_TEST_LC) RATIOS=1 CORES=16 \
	  BUFFER=$(SOLVER_TEST_BUFFER) TIMING_LC_COARSE=$(SOLVER_TEST_LC_COARSE) \
	  TIMING_MESH=$(SOLVER_TEST_MESH) TIMING_BNDIDS=$(SOLVER_TEST_BNDIDS) \
	  TIMING_INVERSION=$(SOLVER_TEST_INVERSION) \
	  TIMING_RESTART=$(SOLVER_TEST_RESTART) \
	  ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER) \
	  T_END=2015.7 DT=0.1 TIMING_KIND=qualification \
	  TIMING_TAG=qualification_5step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG) \
	  TIMING_EXPERIMENT=qualification_5step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG) \
	  SBATCH_WAIT=--wait FORCE_TIMING=1
	python scripts/check_solver_qualification.py \
	  --record results/timing/timing_qualification_5step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_$(SOLVER_TEST_LC_COARSE)_16.json \
	  --csv results/qualification_5step_$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_$(SOLVER_TEST_TAG)_$(SOLVER_TEST_LC)_timeseries.csv --min-steps 5

# Hold the exact 2.5 km state, solver, and number of steps fixed while varying
# only dt.  This is deliberately separate from `timing`: the matrix timestep
# stays unchanged until the largest passing value has been measured.
timestep-probe:
	$(MAKE) transient LCS=$(SOLVER_TEST_LC) RATIOS=1 CORES=16 \
	  BUFFER=$(SOLVER_TEST_BUFFER) TIMING_LC_COARSE=$(SOLVER_TEST_LC_COARSE) \
	  TIMING_MESH=$(SOLVER_TEST_MESH) TIMING_BNDIDS=$(SOLVER_TEST_BNDIDS) \
	  TIMING_INVERSION=$(SOLVER_TEST_INVERSION) \
	  TIMING_RESTART=$(TIMESTEP_PROBE_RESTART) \
	  ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(TIMESTEP_PROBE_SOLVER) \
	  ISMIP7_RESCUE_ENABLED=0 \
	  T_END=$(TIMESTEP_PROBE_END) DT=$(TIMESTEP_PROBE_DT) \
	  TIMING_KIND=qualification TIMING_TAG=$(TIMESTEP_PROBE_TAG) \
	  TIMING_EXPERIMENT=$(TIMESTEP_PROBE_TAG) \
	  SLURM_QUEUE=debug SLURM_TIME=01:00:00 \
	  SBATCH_WAIT=--wait FORCE_TIMING=1
	python scripts/check_solver_qualification.py \
	  --record results/timing/timing_$(TIMESTEP_PROBE_TAG)_$(SOLVER_TEST_LC)_$(SOLVER_TEST_LC_COARSE)_16.json \
	  --csv results/$(TIMESTEP_PROBE_TAG)_$(SOLVER_TEST_LC)_timeseries.csv \
	  --min-steps $(TIMESTEP_PROBE_STEPS)

debug: override LCS := $(SOLVER_TEST_LC)
debug: override RATIOS := 1
debug: CORES := 16
debug: override BUFFER := $(SOLVER_TEST_BUFFER)
debug: override TIMING_LC_COARSE := $(SOLVER_TEST_LC_COARSE)
debug: override TIMING_MESH := $(SOLVER_TEST_MESH)
debug: override TIMING_BNDIDS := $(SOLVER_TEST_BNDIDS)
debug: override TIMING_INVERSION := $(SOLVER_TEST_INVERSION)
debug: override TIMING_KIND := debug
debug: override TIMING_TAG := debug_scpc_mumps_$(SOLVER_TEST_TAG)
debug: override TIMING_EXPERIMENT := debug_scpc_mumps_$(SOLVER_TEST_TAG)
debug: override ISMIP7_DIAGNOSTIC_LINEAR_SOLVER := scpc_mumps
debug: override ISMIP7_SNES_MONITOR := 1
debug: override ISMIP7_SOLVER_VIEW := 1
debug: override ISMIP7_SNES_LOG = $(CURDIR)/results/logs/timing_scpc_mumps_debug_lc$(SOLVER_TEST_LC)_lcc$(SOLVER_TEST_LC_COARSE)_n$(CORES)_$(DEBUG_STAMP).log
debug: override SLURM_QUEUE = debug
debug: override SLURM_TIME = 01:00:00
debug: transient-direct

reference: override LCS := $(SOLVER_TEST_LC)
reference: override RATIOS := 1
reference: CORES := 16
reference: override BUFFER := $(SOLVER_TEST_BUFFER)
reference: override TIMING_LC_COARSE := $(SOLVER_TEST_LC_COARSE)
reference: override TIMING_MESH := $(SOLVER_TEST_MESH)
reference: override TIMING_BNDIDS := $(SOLVER_TEST_BNDIDS)
reference: override TIMING_INVERSION := $(SOLVER_TEST_INVERSION)
reference: override TIMING_KIND := reference
reference: override TIMING_TAG := reference_full_mumps_$(SOLVER_TEST_TAG)
reference: override TIMING_EXPERIMENT := reference_full_mumps_$(SOLVER_TEST_TAG)
reference: override ISMIP7_DIAGNOSTIC_LINEAR_SOLVER := full_mumps
reference: override ISMIP7_SNES_MONITOR := 1
reference: override ISMIP7_SOLVER_VIEW := 1
reference: override ISMIP7_SNES_LOG = $(CURDIR)/results/logs/timing_full_mumps_reference_lc$(SOLVER_TEST_LC)_lcc$(SOLVER_TEST_LC_COARSE)_n$(CORES)_$(DEBUG_STAMP).log
reference: override SLURM_QUEUE = debug
reference: override SLURM_TIME = 01:00:00
reference: transient-direct

solver-smoke:
	OMP_NUM_THREADS=1 mpiexec -n 2 python scripts/smoke_solver_modes.py

timing-cache-smoke:
	@set -e; \
	export ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_mumps; \
	export ISMIP7_RESCUE_ENABLED=0; \
	smoke_dir=$$(mktemp -d); \
	trap 'rm -rf "$$smoke_dir"' EXIT; \
	OMP_NUM_THREADS=1 mpiexec -n 2 python scripts/smoke_timing_cache.py \
	  write --output "$$smoke_dir/cache_parallel.h5"; \
	OMP_NUM_THREADS=1 python scripts/redistribute_checkpoint.py \
	  --input "$$smoke_dir/cache_parallel.h5" \
	  --output "$$smoke_dir/cache_1core.h5" \
	  --manifest "$$smoke_dir/cache_1core.json"; \
	OMP_NUM_THREADS=1 python scripts/smoke_timing_cache.py check \
	  --input "$$smoke_dir/cache_1core.h5" --signature "$$smoke_dir/n1.json"; \
	OMP_NUM_THREADS=1 mpiexec -n 2 python scripts/smoke_timing_cache.py check \
	  --input "$$smoke_dir/cache_1core.h5" --signature "$$smoke_dir/n2.json"; \
	OMP_NUM_THREADS=1 mpiexec -n 2 python scripts/audit_timing_cache.py \
	  --cache "$$smoke_dir/cache_1core.h5" \
	  --output "$$smoke_dir/cache_audit.json"; \
	python scripts/smoke_timing_cache.py compare \
	  "$$smoke_dir/n1.json" "$$smoke_dir/n2.json"

sync-results:
	@mkdir -p "$(CURDIR)/results/"
	rsync -a "$(REMOTE_RESULTS)" "$(CURDIR)/results/"
	@echo "Synced $(REMOTE_RESULTS) into canonical $(CURDIR)/results/"

matrix:
	@test ! -d "$(CURDIR)/antarctica/results" || { \
	  echo "nested results tree detected: $(CURDIR)/antarctica/results" >&2; \
	  echo "merge/remove it before building the matrix (see make sync-results)" >&2; \
	  exit 2; \
	}
	@tag="$(TIMING_TAG)"; \
	if [ -z "$(MATRIX_TAG_NAMED)" ] && [ -f results/timing/latest_matrix_campaign.txt ]; then \
	  read -r tag < results/timing/latest_matrix_campaign.txt; \
	fi; \
	python scripts/build_timing_matrix.py --tag "$$tag" --output "$(MATRIX_OUTPUT)"

map-check:
	python scripts/manage_map_check.py $(MAP_CHECK_ARGS)

map-check-dry-run:
	python scripts/manage_map_check.py $(MAP_CHECK_ARGS) --dry-run --assume-passed

now:
	python scripts/build_now.py --write --output "$(CURDIR)/../NOW.md"

now-check:
	python scripts/build_now.py --check --output "$(CURDIR)/../NOW.md"

# The run log: one JSON record per simulation in runlog/, rendered into
# reports/SIMULATIONS.md.
runlog:
	python scripts/build_runlog.py --write $(if $(RUNLOG_CSV),--csv "$(RUNLOG_CSV)",)

runlog-check:
	python scripts/build_runlog.py --check

matrix-legacy:
	@test ! -d "$(CURDIR)/antarctica/results" || { \
	  echo "nested results tree detected: $(CURDIR)/antarctica/results" >&2; \
	  exit 2; \
	}
	python scripts/build_timing_matrix.py --tag "$(LEGACY_TIMING_TAG)" \
	  --legacy-full-matrix --output TIMING_MATRIX_$(LEGACY_TIMING_TAG).md
