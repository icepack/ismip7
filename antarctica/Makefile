# Antarctica timing benchmark: resolution vs. core-count matrix.
#
# Run from antarctica/:
#   make timing
#
# Stages (each idempotent — skips already-built artifacts):
#   meshes        build 10 meshes (5 LC x {10,20}*LC, buffer=20000 m)
#   inversion     single LC=2500 / LC_coarse=25000 inversion (Slurm or mpiexec)
#   redistribute  rewrite the complete inversion checkpoint on 1 core
#   transient     30 short runs (10 mesh combos x 16/32/64 cores)
#   matrix        write TIMING_MATRIX.md from results/timing/*.json

LCS ?= 500 1000 2000 2500 5000
RATIOS ?= 20 10
CORES ?= 16 32 64
BUFFER ?= 20000
INV_LC ?= 2500
INV_LC_COARSE ?= 25000
FRICTION ?= budd
N_FLOW ?= 3.0
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

.PHONY: timing meshes inversion redistribute transient matrix

timing: matrix

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
	@mkdir -p results/timing
	@for lc in $(LCS); do for r in $(RATIOS); do lcc=$$((lc * r)); for nc in $(CORES); do \
	  out=results/timing/timing_$${lc}_$${lcc}_$${nc}.json; \
	  [ -f $$out ] && { echo "skip $$out"; continue; }; \
	  mesh=$(CURDIR)/mesh/antarctica_$${lcc}_$${lc}_buffered$(BUFFER).msh; \
	  bndids=$(CURDIR)/mesh/boundary_ids_antarctica_$${lcc}_$${lc}_buffered$(BUFFER).json; \
	  if command -v sbatch >/dev/null 2>&1; then \
	    sbatch --job-name=timing_$${lc}_$${lcc}_$${nc} --ntasks-per-node=$$nc \
	      --export=ALL,ISMIP7_LC=$$lc,ISMIP7_LC_COARSE=$$lcc,ISMIP7_BUFFER_M=$(BUFFER),ISMIP7_MESH=$$mesh,ISMIP7_BNDIDS=$$bndids,ISMIP7_INVERSION=$(CURDIR)/$(INV_CHK),ISMIP7_FRICTION=$(FRICTION),ISMIP7_N_FLOW=$(N_FLOW),ISMIP7_A4_FACTOR=$(A4_FACTOR),ISMIP7_T_END=2020,ISMIP7_DT=1.0 \
	      scripts/batch_runners/timing_transient.script; \
	  else \
	    ISMIP7_LC=$$lc ISMIP7_LC_COARSE=$$lcc ISMIP7_BUFFER_M=$(BUFFER) ISMIP7_MESH=$$mesh ISMIP7_BNDIDS=$$bndids ISMIP7_INVERSION=$(CURDIR)/$(INV_CHK) \
	      ISMIP7_FRICTION=$(FRICTION) ISMIP7_N_FLOW=$(N_FLOW) ISMIP7_A4_FACTOR=$(A4_FACTOR) ISMIP7_T_END=2020 ISMIP7_DT=1.0 \
	      mpiexec -n $$nc python scripts/run_timing.py; \
	  fi; \
	done; done; done

matrix: transient
	python scripts/build_timing_matrix.py
