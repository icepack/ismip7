#!/usr/bin/env python3
"""
Run a short transient simulation and record wall-clock timing.

Used by the `make timing` benchmark to measure resolution vs. core-count
performance. Loads the mesh/inversion for the current ISMIP7_LC /
ISMIP7_LC_COARSE / ISMIP7_BUFFER_M env vars (set by the Makefile or Slurm
job), runs 5 years of zero-forcing time stepping, and writes a JSON record
under results/timing/.

Usage:
    ISMIP7_LC=2500 ISMIP7_LC_COARSE=25000 ISMIP7_BUFFER_M=20000 \
    ISMIP7_INVERSION=mesh/inversion_icepack2_budd_n3_2500.h5 \
    mpiexec -n 16 python scripts/run_timing.py
"""

import json
import os
import sys
from time import perf_counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from firedrake import COMM_WORLD
from firedrake.petsc import PETSc

from simulation import setup_model, run_simulation, lc

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIMING_DIR = os.path.join(_ROOT, "results", "timing")

T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2020"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "5"))
LINEAR_SOLVER = os.environ.get(
    "ISMIP7_LINEAR_SOLVER", "iterative"
).strip().lower()


def main():
    os.makedirs(TIMING_DIR, exist_ok=True)

    ncores = COMM_WORLD.size
    PETSc.Sys.Print(
        f"Timing run: lc={lc} "
        f"lc_coarse={os.environ.get('ISMIP7_LC_COARSE', 'unknown')} "
        f"buffer={os.environ.get('ISMIP7_BUFFER_M', 'unknown')} "
        f"ncores={ncores} solver={LINEAR_SOLVER} "
        f"t={T_START}->{T_END} dt={DT}"
    )

    t0 = perf_counter()
    ctx = setup_model()
    mesh = ctx["mesh"]
    target_lc_coarse = ctx["lc_coarse"]
    target_buffer_m = ctx["buffer_m"]
    nsteps = int(round((T_END - T_START) / DT))

    PETSc.Sys.Print(
        f"Timing compute mesh: {mesh.num_vertices()} vertices, "
        f"{mesh.num_cells()} cells"
    )

    run_simulation(
        ctx,
        experiment_name="timing",
        t_start=T_START,
        t_end=T_END,
        dt=DT,
        output_interval=OUTPUT_INTERVAL,
        checkpoint_interval=nsteps + 1,  # skip intermediate checkpoints
        forcing_callback=None,
    )
    run_seconds = perf_counter() - t0

    record = {
        "lc": lc,
        "lc_coarse": target_lc_coarse,
        "buffer_m": target_buffer_m,
        "ncores": ncores,
        "vertices": mesh.num_vertices(),
        "cells": mesh.num_cells(),
        "t_start": T_START,
        "t_end": T_END,
        "dt": DT,
        "nsteps": nsteps,
        "linear_solver": LINEAR_SOLVER,
        "run_seconds": run_seconds,
        "seconds_per_step": run_seconds / max(nsteps, 1),
    }

    if COMM_WORLD.rank == 0:
        out_fn = os.path.join(
            TIMING_DIR, f"timing_{lc}_{target_lc_coarse}_{ncores}.json"
        )
        with open(out_fn, "w") as f:
            json.dump(record, f, indent=2)
        PETSc.Sys.Print(
            f"Timing record: {run_seconds:.1f}s total "
            f"({record['seconds_per_step']:.2f}s/step) -> {out_fn}"
        )


if __name__ == "__main__":
    main()
