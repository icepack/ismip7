#!/usr/bin/env python3
"""Prepare a solved, zero-step transient state on one exact target mesh.

The parallel preparation job builds cell-averaged BedMachine geometry on the
requested mesh, transfers the continuous fields from the improved 2.5 km
inversion, performs the normal adaptive cold continuation, and writes a
complete mixed state. A separate one-rank repack publishes the final cache and
manifest.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.fspath(Path(__file__).resolve().parents[2]))

from firedrake import COMM_WORLD
from firedrake.petsc import PETSc

from icepack2_tools.solverconfig import (
    diagnostic_solver_mode,
    solver_provenance,
)
from simulation import save_model_state, setup_model
from timing_campaign import (
    CACHE_ROLE,
    CACHE_SCHEMA_VERSION,
    MATRIX_T_START,
    sha256_file,
    solver_configuration_fingerprint,
)
from icepack2_tools.runconfig import TARGET_MESH_GEOMETRY_METHOD


def _required(name):
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"{name} is required")
    return value


def main():
    inversion = os.path.realpath(_required("ISMIP7_INVERSION"))
    mesh_input = os.path.realpath(_required("ISMIP7_MESH"))
    output = os.path.realpath(_required("ISMIP7_TIMING_CACHE_RAW"))
    if diagnostic_solver_mode() != "scpc_mumps":
        raise RuntimeError("timing caches must be prepared with scpc_mumps")
    if os.path.realpath(inversion) == os.path.realpath(output):
        raise RuntimeError("cache output must differ from its source inversion")
    if not os.path.isfile(inversion):
        raise FileNotFoundError(inversion)
    if not os.path.isfile(mesh_input):
        raise FileNotFoundError(mesh_input)

    source_sha256 = sha256_file(inversion) if COMM_WORLD.rank == 0 else None
    source_sha256 = COMM_WORLD.bcast(source_sha256, root=0)
    mesh_sha256 = sha256_file(mesh_input) if COMM_WORLD.rank == 0 else None
    mesh_sha256 = COMM_WORLD.bcast(mesh_sha256, root=0)
    configuration = solver_provenance()
    fingerprint = solver_configuration_fingerprint(configuration)
    PETSc.Sys.Print(
        f"Preparing timing cache from {inversion} on {mesh_input} "
        f"with {COMM_WORLD.size} ranks"
    )

    ctx = setup_model(restart_from=None)
    geometry_method = ctx.get("geometry_source_method")
    if geometry_method != TARGET_MESH_GEOMETRY_METHOD:
        raise RuntimeError(
            "timing cache geometry was initialized with "
            f"{geometry_method!r}; expected {TARGET_MESH_GEOMETRY_METHOD!r}"
        )
    geometry_source = ctx.get("geometry_source")
    if not geometry_source:
        raise RuntimeError("timing cache geometry source was not recorded")
    attrs = {
        "timing_cache_schema_version": CACHE_SCHEMA_VERSION,
        "timing_cache_role": CACHE_ROLE,
        "source_inversion": inversion,
        "source_inversion_sha256": source_sha256,
        "source_mesh_sha256": mesh_sha256,
        "diagnostic_solver_mode": diagnostic_solver_mode(),
        "solver_configuration": json.dumps(configuration, sort_keys=True),
        "solver_configuration_fingerprint": fingerprint,
        "n_flow": float(ctx["n_flow_val"]),
        "a4_factor": float(os.environ.get("ISMIP7_A4_FACTOR", "1.0")),
        "geometry_source": os.path.realpath(geometry_source),
        "geometry_source_method": geometry_method,
    }
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    save_model_state(ctx, output, MATRIX_T_START, extra_attrs=attrs)
    PETSc.Sys.Print(
        f"Prepared zero-step full-state cache: {output} "
        f"(source sha256={source_sha256}, geometry={geometry_method}, "
        f"geometry source={os.path.basename(geometry_source)})"
    )


if __name__ == "__main__":
    main()
