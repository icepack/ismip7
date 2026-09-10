# Antarctica timing benchmark: resolution vs. core-count matrix.
#
# Run from antarctica/:
#   make timing
#
# Targets (each idempotent — skips already-built artifacts):
#   meshes        build 10 meshes (5 LC x {10,20}*LC, buffer=20000 m)
#   inversion     single LC=2500 / LC_coarse=25000 inversion (Slurm or mpiexec)
#   redistribute  rewrite the complete inversion checkpoint on 1 core
#   transient     30 short runs (10 mesh combos x 16/32/64 cores),
#                 with a resolution-specific Slurm memory request
#   qualify       gate one solver through 2-step and 5-step 2.5 km probes
#   timestep-probe restart at 2.5 km for 5 steps at a selectable coarser dt
#   debug         one condensed-MUMPS job with full convergence monitors
#   reference     one full-Jacobian MUMPS reference job
#   solver-smoke  tiny two-rank setup test; no Antarctic data or time steps
#   matrix        write TIMING_MATRIX.md from results/timing/*.json

LCS ?= 500 1000 2000 2500 5000
RATIOS ?= 20 10
CORES ?= 16 32 64
BUFFER ?= 20000
INV_LC ?= 2500
INV_LC_COARSE ?= 25000
FRICTION ?= budd
N_FLOW ?= 3.0
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER ?= scpc_gamg
ISMIP7_SOLVER_VIEW ?=
T_END ?= 2020
DT ?= 1.0
TIMING_KIND ?= matrix
TIMING_TAG ?= $(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER)_t$(T_END)_dt$(DT)
TIMING_EXPERIMENT ?= timing
# Optional exact-input overrides for a single transient lane.  The resolution
# matrix leaves these empty and constructs each target mesh parametrically.
TIMING_LC_COARSE ?=
TIMING_MESH ?=
TIMING_BNDIDS ?=
TIMING_INVERSION ?= $(CURDIR)/$(INV_CHK)
TIMING_RESTART ?=
SBATCH_WAIT ?=
FORCE_TIMING ?=
SLURM_PARTITION ?= general
SLURM_TIME ?= 12:00:00

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
SOLVER_TEST_INVERSION ?= $(CURDIR)/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_2500.h5
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

INV_CHK := mesh/inversion_icepack2$(MAP_FRICTION_TAG)$(MAP_N_TAG)_$(INV_LC)_1core.h5
INV_CHK_MC := mesh/inversion_icepack2$(MAP_FRICTION_TAG)$(MAP_N_TAG)_$(INV_LC).h5

empty :=
space := $(empty) $(empty)
colon := :
TIMING_LCS := $(subst $(space),$(colon),$(LCS))
TIMING_RATIOS := $(subst $(space),$(colon),$(RATIOS))
DEBUG_STAMP := $(shell date +%Y%m%d_%H%M%S)

# Total-memory requests for timing jobs. These are keyed to the fine
# resolution: the 500 m case is the memory outlier, while the 2500 m debug
# case uses MUMPS and has a much smaller measured footprint. The ratio-10
# meshes are slightly larger than their ratio-20 counterparts, so one
# conservative budget per fine resolution covers both ratios.

.PHONY: timing meshes inversion redistribute transient qualify qualify-2step qualify-5step timestep-probe debug reference solver-smoke matrix

timing: qualify
	$(MAKE) matrix SBATCH_WAIT=--wait

meshes:
	@if command -v sbatch >/dev/null 2>&1; then \
	  sbatch --wait --job-name=ant_timing_mesh \
	    --export=ALL,TIMING_LCS=$(TIMING_LCS),TIMING_RATIOS=$(TIMING_RATIOS),TIMING_BUFFER=$(BUFFER) \
	    scripts/batch_runners/timing_meshes.script; \
	else \
	  for lc in $(LCS); do for r in $(RATIOS); do \
	    lcc=$$((lc * r)); fn=mesh/antarctica_$${lcc}_$${lc}_buffered$(BUFFER).msh; bnd=mesh/boundary_ids_antarctica_$${lcc}_$${lc}_buffered$(BUFFER).json; \
	    [ -f $$fn ] && [ -f $$bnd ] && { echo "skip $$fn"; continue; }; \
	    python scripts/mesh_antarctica.py --lc $$lc --lc-coarse $$lcc --buffer-m $(BUFFER); \
	  done; done; \
	fi

inversion:
	@[ -f $(INV_CHK_MC) ] && echo "skip inversion" || \
	if command -v sbatch >/dev/null 2>&1; then \
	  sbatch --wait --export=ALL,ISMIP7_LC=$(INV_LC),ISMIP7_LC_COARSE=$(INV_LC_COARSE),ISMIP7_BUFFER_M=$(BUFFER),ISMIP7_FRICTION=$(FRICTION),ISMIP7_N_FLOW=$(N_FLOW),ISMIP7_A4_FACTOR=$(A4_FACTOR) \
	    scripts/batch_runners/timing_inversion.script; \
	else ISMIP7_LC=$(INV_LC) ISMIP7_LC_COARSE=$(INV_LC_COARSE) ISMIP7_BUFFER_M=$(BUFFER) \
	  ISMIP7_FRICTION=$(FRICTION) ISMIP7_N_FLOW=$(N_FLOW) ISMIP7_A4_FACTOR=$(A4_FACTOR) \
	  mpiexec -n 12 python scripts/inversion_icepack2.py; fi

redistribute:
	@if command -v sbatch >/dev/null 2>&1; then \
	  sbatch --wait --job-name=ant_timing_redist --export=ALL,TIMING_INV_CHK=$(INV_CHK_MC) \
	    scripts/batch_runners/timing_redistribute.script; \
	else \
	  python scripts/redistribute_checkpoint.py --input $(INV_CHK_MC) --output $(INV_CHK); \
	fi

transient:
	@mkdir -p results/timing results/logs
	@for lc in $(LCS); do for r in $(RATIOS); do lcc=$$((lc * r)); for nc in $(CORES); do \
	  if [ -n "$(TIMING_LC_COARSE)" ]; then lcc=$(TIMING_LC_COARSE); fi; \
	  out=results/timing/timing_$(TIMING_TAG)_$${lc}_$${lcc}_$${nc}.json; \
	  [ -z "$(FORCE_TIMING)" ] && [ -f $$out ] && { echo "skip $$out"; continue; }; \
	  case "$$lc" in \
	    500)  mem=128G;; \
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
	  [ -f "$$mesh" ] || { echo "mesh not found: $$mesh" >&2; exit 2; }; \
	  [ -f "$$bndids" ] || { echo "boundary ids not found: $$bndids" >&2; exit 2; }; \
	  [ -f "$$inversion" ] || { echo "inversion not found: $$inversion" >&2; exit 2; }; \
	  [ -z "$$restart" ] || [ -f "$$restart" ] || { echo "restart not found: $$restart" >&2; exit 2; }; \
	  if command -v sbatch >/dev/null 2>&1; then \
	    sbatch $(SBATCH_WAIT) --job-name=timing_$${lc}_$${lcc}_$${nc} --ntasks-per-node=$$nc \
	      --partition=$(SLURM_PARTITION) --time=$(SLURM_TIME) --mem=$$mem \
	      --export=ALL,ISMIP7_LC=$$lc,ISMIP7_LC_COARSE=$$lcc,ISMIP7_BUFFER_M=$(BUFFER),ISMIP7_MESH=$$mesh,ISMIP7_BNDIDS=$$bndids,ISMIP7_INVERSION=$$inversion,ISMIP7_RESTART=$$restart,ISMIP7_FRICTION=$(FRICTION),ISMIP7_N_FLOW=$(N_FLOW),ISMIP7_A4_FACTOR=$(A4_FACTOR),ISMIP7_T_END=$(T_END),ISMIP7_DT=$(DT),ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER),ISMIP7_TIMING_KIND=$(TIMING_KIND),ISMIP7_TIMING_TAG=$(TIMING_TAG),ISMIP7_TIMING_EXPERIMENT=$(TIMING_EXPERIMENT),ISMIP7_SNES_MONITOR=$(ISMIP7_SNES_MONITOR),ISMIP7_SNES_LOG=$(ISMIP7_SNES_LOG),ISMIP7_SOLVER_VIEW=$(ISMIP7_SOLVER_VIEW) \
	      scripts/batch_runners/timing_transient.script; \
	  else \
	    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ISMIP7_LC=$$lc ISMIP7_LC_COARSE=$$lcc ISMIP7_BUFFER_M=$(BUFFER) ISMIP7_MESH=$$mesh ISMIP7_BNDIDS=$$bndids ISMIP7_INVERSION=$$inversion ISMIP7_RESTART=$$restart ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=$(ISMIP7_DIAGNOSTIC_LINEAR_SOLVER) ISMIP7_TIMING_KIND=$(TIMING_KIND) ISMIP7_TIMING_TAG=$(TIMING_TAG) ISMIP7_TIMING_EXPERIMENT=$(TIMING_EXPERIMENT) ISMIP7_SNES_MONITOR=$(ISMIP7_SNES_MONITOR) ISMIP7_SNES_LOG=$(ISMIP7_SNES_LOG) ISMIP7_SOLVER_VIEW=$(ISMIP7_SOLVER_VIEW) \
	      ISMIP7_FRICTION=$(FRICTION) ISMIP7_N_FLOW=$(N_FLOW) ISMIP7_A4_FACTOR=$(A4_FACTOR) ISMIP7_T_END=$(T_END) ISMIP7_DT=$(DT) \
	      mpiexec -n $$nc python scripts/run_timing.py; \
	  fi; \
	done; done; done

# Qualification is intentionally a single 2.5 km / 16-rank lane on the exact
# 64 km-coarse mesh used by the imported DG0 logvelnet MAP.  The checker gates
# requested versus completed step count, final year, convergence reasons and
# the mass-budget residual before the five-step probe starts.
qualify: qualify-2step
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
	  T_END=$(TIMESTEP_PROBE_END) DT=$(TIMESTEP_PROBE_DT) \
	  TIMING_KIND=qualification TIMING_TAG=$(TIMESTEP_PROBE_TAG) \
	  TIMING_EXPERIMENT=$(TIMESTEP_PROBE_TAG) \
	  SLURM_PARTITION=debug SLURM_TIME=01:00:00 \
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
debug: override SLURM_PARTITION = debug
debug: override SLURM_TIME = 01:00:00
debug: transient

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
reference: override SLURM_PARTITION = debug
reference: override SLURM_TIME = 01:00:00
reference: transient

solver-smoke:
	OMP_NUM_THREADS=1 mpiexec -n 2 python scripts/smoke_solver_modes.py

matrix: transient
	python scripts/build_timing_matrix.py
