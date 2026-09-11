#!/usr/bin/env python3
"""Tiny MPI-size portability smoke test for complete timing caches."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)

import firedrake as fd  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402
from mpi4py import MPI  # noqa: E402

from icepack2_tools.solverconfig import (  # noqa: E402
    diagnostic_solver_mode,
    diagnostic_solver_parameters,
    solver_provenance,
)
from simulation import save_model_state  # noqa: E402
from timing_campaign import (  # noqa: E402
    CACHE_ROLE,
    CACHE_SCHEMA_VERSION,
    atomic_write_json,
    solver_configuration_fingerprint,
)
from icepack2_tools.runconfig import TARGET_MESH_GEOMETRY_METHOD  # noqa: E402


def spaces(mesh):
    q = fd.FunctionSpace(mesh, "CG", 1)
    qg = fd.FunctionSpace(mesh, "DG", 0)
    velocity = fd.VectorFunctionSpace(mesh, "CG", 1)
    dg0 = fd.FiniteElement("DG", "triangle", 0)
    stress = fd.TensorFunctionSpace(mesh, dg0, symmetry=True)
    traction = fd.VectorFunctionSpace(mesh, dg0)
    return q, qg, velocity * stress * traction


def residual_form(state):
    mixed = state.function_space()
    u, membrane, basal = fd.split(state)
    v, q, r = fd.TestFunctions(mixed)
    residual = (
        fd.inner(fd.grad(u), fd.grad(v))
        + fd.inner(u, v)
        + fd.inner(membrane - fd.sym(fd.grad(u)), q - fd.sym(fd.grad(v)))
        + fd.inner(basal - u, r - v)
        - fd.inner(fd.as_vector((1.0, 1.0)), v)
    ) * fd.dx
    zero = fd.Constant(0.0)
    residual += fd.derivative(
        zero * membrane[0, 0] * basal[0] * fd.dx, state
    )
    return residual


def write_cache(path):
    if diagnostic_solver_mode() != "scpc_mumps":
        raise RuntimeError("smoke cache requires scpc_mumps")
    mesh = fd.UnitSquareMesh(4, 4)
    q, qg, mixed = spaces(mesh)
    state = fd.Function(mixed)
    residual = residual_form(state)
    problem = fd.NonlinearVariationalProblem(residual, state)
    solver = fd.NonlinearVariationalSolver(
        problem,
        solver_parameters=diagnostic_solver_parameters(),
        options_prefix="timing_cache_smoke_",
    )
    solver.solve()
    if solver.snes.getConvergedReason() <= 0:
        raise RuntimeError("tiny mixed solve did not converge")

    x, y = fd.SpatialCoordinate(mesh)
    theta = fd.Function(q).interpolate(0.1 * x)
    phi = fd.Function(q).interpolate(-0.2 * y)
    velocity_obs = fd.Function(state.subfunctions[0].function_space()).interpolate(
        fd.as_vector((1.0 + x, 2.0 - y))
    )
    bed = fd.Function(qg).interpolate(-100.0 - x)
    thickness = fd.Function(qg).interpolate(1000.0 + x + y)
    surface = fd.Function(qg).interpolate(bed + thickness)
    ctx = {
        "mesh": mesh,
        "Q": q,
        "Q_g": qg,
        "geom_dg": True,
        "mesh_basename": "tiny_timing_cache.msh",
        "z": state,
        "theta": theta,
        "phi": phi,
        "u_obs": velocity_obs,
        "b": bed,
        "h": thickness,
        "s": surface,
        "H_init": thickness.copy(deepcopy=True),
        "phi_eff": fd.Function(qg).assign(0.5),
        "C_w0": fd.Function(qg).assign(2.0),
        "N_ref": fd.Function(qg).assign(3.0),
        "A_prior": fd.Function(q).assign(4.0),
        "a_ref_mb": fd.Function(qg).assign(0.0),
        "friction": "budd",
        "lc": 1,
        "lc_coarse": 10,
        "buffer_m": 20000,
        "geometry_source": "/synthetic/BedMachine.nc",
        "geometry_source_method": TARGET_MESH_GEOMETRY_METHOD,
    }
    configuration = solver_provenance()
    attrs = {
        "timing_cache_schema_version": CACHE_SCHEMA_VERSION,
        "timing_cache_role": CACHE_ROLE,
        "source_inversion": "/synthetic/source.h5",
        "source_inversion_sha256": "synthetic-source",
        "source_mesh_sha256": "synthetic-mesh",
        "diagnostic_solver_mode": "scpc_mumps",
        "solver_configuration": json.dumps(configuration, sort_keys=True),
        "solver_configuration_fingerprint": (
            solver_configuration_fingerprint(configuration)
        ),
        "n_flow": 3.0,
        "a4_factor": 1.0,
        "geometry_source": "/synthetic/BedMachine.nc",
        "geometry_source_method": TARGET_MESH_GEOMETRY_METHOD,
    }
    save_model_state(ctx, path, 2015.0, extra_attrs=attrs)
    PETSc.Sys.Print(f"Wrote tiny cache on {mesh.comm.size} ranks: {path}")


def check_cache(path, signature_path):
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        q, _, mixed = spaces(mesh)
        fields = {}
        for name in (
            "log_friction",
            "log_fluidity",
            "velocity_obs",
            "bed",
            "thickness",
            "surface",
            "H_init",
            "phi_eff",
            "C_w0",
            "N_ref",
            "fluidity_prior",
            "a_ref_mb",
        ):
            fields[name] = chk.load_function(mesh, name=name)
        state = fd.Function(mixed)
        for target, name in zip(
            state.subfunctions,
            ("velocity", "membrane_stress", "basal_stress"),
            strict=True,
        ):
            target.assign(chk.load_function(mesh, name=name))
        attrs = {
            name: chk.get_attr("/", name)
            for name in (
                "timing_cache_schema_version",
                "timing_cache_role",
                "source_inversion_sha256",
                "source_mesh_sha256",
                "solver_configuration_fingerprint",
                "t_yr",
                "friction",
                "geometry_space",
                "lc",
                "lc_coarse",
                "geometry_source",
                "geometry_source_method",
            )
        }
    assembled = fd.assemble(residual_form(state))
    with assembled.dat.vec_ro as vector:
        residual_norm = vector.norm()
    signature = {
        "vertices": int(q.dim()),
        "cells": int(mesh.comm.allreduce(mesh.cell_set.size, op=MPI.SUM)),
        "field_norms": {
            name: float(fd.norm(field)) for name, field in fields.items()
        },
        "state_norm": float(fd.norm(state)),
        "residual_norm": float(residual_norm),
        "metadata": {
            name: value.item() if hasattr(value, "item") else value
            for name, value in attrs.items()
        },
    }
    if not math.isfinite(signature["residual_norm"]):
        raise RuntimeError("loaded full-state residual is not finite")
    if mesh.comm.rank == 0:
        atomic_write_json(signature_path, signature)
    PETSc.Sys.Print(
        f"Checked tiny cache on {mesh.comm.size} ranks: "
        f"residual={residual_norm:.3e}"
    )


def compare_signatures(paths):
    signatures = []
    for path in paths:
        with open(path) as stream:
            signatures.append(json.load(stream))
    reference = signatures[0]
    for index, signature in enumerate(signatures[1:], 2):
        if signature["vertices"] != reference["vertices"]:
            raise RuntimeError(f"signature {index} has different vertex count")
        if signature["cells"] != reference["cells"]:
            raise RuntimeError(f"signature {index} has different cell count")
        if signature["metadata"] != reference["metadata"]:
            raise RuntimeError(f"signature {index} has different metadata")
        for name in reference["field_norms"]:
            if not math.isclose(
                signature["field_norms"][name],
                reference["field_norms"][name],
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            ):
                raise RuntimeError(
                    f"signature {index} has different norm for {name}"
                )
        for name in ("state_norm", "residual_norm"):
            if not math.isclose(
                signature[name], reference[name],
                rel_tol=1.0e-11, abs_tol=1.0e-12,
            ):
                raise RuntimeError(f"signature {index} differs in {name}")
    print(f"PASS: {len(signatures)} MPI-size cache signatures are identical")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    write = subparsers.add_parser("write")
    write.add_argument("--output", required=True)
    check = subparsers.add_parser("check")
    check.add_argument("--input", required=True)
    check.add_argument("--signature", required=True)
    compare = subparsers.add_parser("compare")
    compare.add_argument("signatures", nargs="+")
    args = parser.parse_args()
    if args.command == "write":
        write_cache(args.output)
    elif args.command == "check":
        check_cache(args.input, args.signature)
    else:
        compare_signatures(args.signatures)


if __name__ == "__main__":
    main()
