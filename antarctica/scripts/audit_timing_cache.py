#!/usr/bin/env python3
r"""Audit the initial finite-volume flux divergence in a timing cache.

This is a read-only diagnostic: it loads the exact mesh and fields embedded in
one prepared cache and assembles the same DG0 upwind ``div(h * u)`` used to
construct ``a_ref_mb`` in ``simulation.run_simulation``.  It does not run a
diagnostic or transport solve.

The reported ``a_ref_equivalent`` is the source that would cancel the cached
state's no-forcing thickness tendency.  Consequently
``no_forcing_dhdt = -a_ref_equivalent``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

_SCRIPTS = Path(__file__).resolve().parent
_PROJECT = _SCRIPTS.parents[1]
sys.path.insert(0, os.fspath(_SCRIPTS))
sys.path.insert(0, os.fspath(_PROJECT))

import firedrake as fd  # noqa: E402
from firedrake import COMM_WORLD  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402
from mpi4py import MPI  # noqa: E402

from icepack2_tools.mpi_stats import global_count, global_size  # noqa: E402
from timing_campaign import (  # noqa: E402
    CACHE_ROLE,
    CACHE_SCHEMA_VERSION,
    atomic_write_json,
    validate_cache_manifest,
)


RHO_I = 917.0
RHO_W = 1024.0
RHO_GT = RHO_I / 1.0e12
AUDIT_SCHEMA_VERSION = 1


def _normal(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "item"):
        return value.item()
    return value


def _matches(actual, expected):
    if isinstance(expected, float):
        try:
            return math.isclose(float(actual), expected, abs_tol=1.0e-12)
        except (TypeError, ValueError):
            return False
    return actual == expected


def _read_manifest(path, cache_path):
    if path is None:
        return None
    with open(path) as stream:
        manifest = json.load(stream)
    valid, detail = validate_cache_manifest(
        manifest,
        lc=int(manifest["lc"]),
        lc_coarse=int(manifest["lc_coarse"]),
        cache_path=cache_path,
    )
    if not valid:
        raise RuntimeError(f"invalid timing cache manifest: {detail}")
    return manifest


def _load_cache(path, manifest):
    with fd.CheckpointFile(path, "r") as checkpoint:
        mesh = checkpoint.load_mesh()
        fields = {
            name: checkpoint.load_function(mesh, name=name)
            for name in (
                "velocity",
                "thickness",
                "bed",
                "surface",
                "H_init",
            )
        }
        attr_names = (
            "timing_cache_schema_version",
            "timing_cache_role",
            "source_inversion",
            "source_inversion_sha256",
            "source_mesh_sha256",
            "solver_configuration_fingerprint",
            "geometry_source",
            "geometry_source_method",
            "mesh_basename",
            "lc",
            "lc_coarse",
            "buffer_m",
            "t_yr",
            "friction",
            "geometry_space",
            "n_flow",
        )
        attrs = {
            name: _normal(checkpoint.get_attr("/", name))
            for name in attr_names
            if checkpoint.has_attr("/", name)
        }

    required_attrs = {
        "timing_cache_schema_version": CACHE_SCHEMA_VERSION,
        "timing_cache_role": CACHE_ROLE,
        "geometry_space": "dg0",
    }
    for name, expected in required_attrs.items():
        if not _matches(attrs.get(name), expected):
            raise RuntimeError(
                f"cache {name}={attrs.get(name)!r}; expected {expected!r}"
            )

    if manifest is not None:
        manifest_to_attr = {
            "cache_schema_version": "timing_cache_schema_version",
            "cache_role": "timing_cache_role",
            "source_inversion": "source_inversion",
            "source_inversion_sha256": "source_inversion_sha256",
            "source_mesh_sha256": "source_mesh_sha256",
            "solver_configuration_fingerprint": (
                "solver_configuration_fingerprint"
            ),
            "geometry_source": "geometry_source",
            "geometry_source_method": "geometry_source_method",
            "mesh_basename": "mesh_basename",
            "lc": "lc",
            "lc_coarse": "lc_coarse",
            "buffer_m": "buffer_m",
            "t_yr": "t_yr",
            "friction": "friction",
            "geometry_space": "geometry_space",
            "n_flow": "n_flow",
        }
        for manifest_name, attr_name in manifest_to_attr.items():
            expected = manifest[manifest_name]
            actual = attrs.get(attr_name)
            if not _matches(actual, expected):
                raise RuntimeError(
                    f"loaded cache {attr_name}={actual!r}; "
                    f"manifest records {expected!r}"
                )
    return mesh, fields, attrs


def _extreme_payload(field, supporting, coordinates, mode):
    values = np.asarray(field.dat.data_ro)
    xy = np.asarray(coordinates.dat.data_ro)
    if values.ndim == 2 and values.shape[1] == 1:
        values = values[:, 0]
    if values.ndim != 1 or xy.ndim != 2 or xy.shape[0] != values.shape[0]:
        raise ValueError("audit fields and coordinates have incompatible layouts")
    if mode not in {"min", "max", "absmax"}:
        raise ValueError("mode must be min, max, or absmax")

    comm = field.function_space().mesh().comm
    if values.size:
        if mode == "min":
            index = int(np.argmin(values))
            score = float(values[index])
        elif mode == "max":
            index = int(np.argmax(values))
            score = float(values[index])
        else:
            index = int(np.argmax(np.abs(values)))
            score = float(abs(values[index]))
        payload = {
            "a_ref_equivalent_m_per_yr": float(values[index]),
            "no_forcing_dhdt_m_per_yr": -float(values[index]),
            "xy_m": [float(value) for value in xy[index]],
            "rank": int(comm.rank),
        }
        for name, values_support in supporting.items():
            payload[name] = float(values_support[index])
        payload["grounded"] = bool(payload["height_above_flotation_m"] > 0.0)
        candidate = (score, int(comm.rank), payload)
    else:
        sentinel = math.inf if mode == "min" else -math.inf
        candidate = (sentinel, int(comm.rank), None)

    candidates = comm.allgather(candidate)
    if mode == "min":
        _score, _rank, winner = min(
            candidates, key=lambda item: (item[0], item[1])
        )
    else:
        _score, _rank, winner = max(
            candidates, key=lambda item: (item[0], -item[1])
        )
    if winner is None:
        raise RuntimeError("timing cache contains no owned DG0 cells")
    return winner


def _field_audit(field, supporting, coordinates):
    return {
        mode: _extreme_payload(field, supporting, coordinates, mode)
        for mode in ("min", "max", "absmax")
    }


def audit_cache(cache_path, manifest_path=None, front_hmin=1.0):
    cache_path = os.path.realpath(cache_path)
    manifest_path = (
        os.path.realpath(manifest_path) if manifest_path is not None else None
    )
    manifest = _read_manifest(manifest_path, cache_path)
    mesh, fields, attrs = _load_cache(cache_path, manifest)

    thickness = fields["thickness"]
    space = thickness.function_space()
    element = space.ufl_element()
    if element.degree() != 0 or thickness.ufl_shape:
        raise RuntimeError(
            "timing cache tendency audit requires scalar DG0 thickness"
        )
    for name in ("bed", "surface", "H_init"):
        if (
            fields[name].function_space().ufl_element() != element
            or fields[name].ufl_shape
        ):
            raise RuntimeError(f"cache field {name} is not in the thickness space")

    test = fd.TestFunction(space)
    normal = fd.FacetNormal(mesh)
    velocity = fields["velocity"]
    normal_velocity = fd.dot(velocity, normal)
    outward_velocity = (normal_velocity + abs(normal_velocity)) / 2.0
    flux = fd.assemble(
        (
            outward_velocity("+") * thickness("+")
            - outward_velocity("-") * thickness("-")
        )
        * fd.jump(test)
        * fd.dS
        + outward_velocity * thickness * test * fd.ds
    )
    area = fd.assemble(test * fd.dx)
    bad_area_local = bool(np.any(np.asarray(area.dat.data_ro) <= 0.0))
    if mesh.comm.allreduce(bad_area_local, op=MPI.LOR):
        raise RuntimeError("timing cache mesh contains a non-positive cell area")

    raw = fd.Function(space, name="a_ref_equivalent_raw")
    raw.dat.data[:] = flux.dat.data_ro / area.dat.data_ro
    bad_tendency_local = not bool(
        np.all(np.isfinite(np.asarray(raw.dat.data_ro)))
    )
    if mesh.comm.allreduce(bad_tendency_local, op=MPI.LOR):
        raise RuntimeError("initial FV tendency contains a non-finite value")

    coordinates = fd.Function(
        fd.VectorFunctionSpace(mesh, element), name="cell_coordinates"
    ).interpolate(fd.SpatialCoordinate(mesh))
    speed = fd.Function(space, name="cell_center_speed").interpolate(
        fd.sqrt(fd.dot(velocity, velocity))
    )
    h_data = np.asarray(thickness.dat.data_ro)
    b_data = np.asarray(fields["bed"].dat.data_ro)
    supporting = {
        "thickness_m": h_data,
        "initial_extent_thickness_m": np.asarray(
            fields["H_init"].dat.data_ro
        ),
        "bed_m": b_data,
        "surface_m": np.asarray(fields["surface"].dat.data_ro),
        "height_above_flotation_m": (
            h_data + np.minimum(b_data, 0.0) * (RHO_W / RHO_I)
        ),
        "cell_center_speed_m_per_yr": np.asarray(speed.dat.data_ro),
        "cell_area_m2": np.asarray(area.dat.data_ro),
    }

    raw_audit = _field_audit(raw, supporting, coordinates)
    raw_audit["net_gt_per_yr"] = float(fd.assemble(raw * fd.dx)) * RHO_GT

    initial_extent = np.asarray(fields["H_init"].dat.data_ro)
    beyond_front = initial_extent < front_hmin
    effective = raw.copy(deepcopy=True)
    effective.rename("a_ref_equivalent_fixed_front")
    effective.dat.data[beyond_front] = 0.0
    effective_audit = _field_audit(effective, supporting, coordinates)
    effective_audit["net_gt_per_yr"] = (
        float(fd.assemble(effective * fd.dx)) * RHO_GT
    )
    effective_audit["front_hmin_m"] = float(front_hmin)
    effective_audit["masked_cell_count"] = global_count(
        beyond_front, mesh.comm
    )

    metadata = {
        name: _normal(value) for name, value in attrs.items()
    }
    return {
        "audit_schema_version": AUDIT_SCHEMA_VERSION,
        "cache_path": cache_path,
        "manifest_path": manifest_path,
        "cells": global_size(thickness),
        "metadata": metadata,
        "sign_convention": (
            "a_ref_equivalent=FV div(h*u); no_forcing_dhdt=-a_ref_equivalent"
        ),
        "raw": raw_audit,
        "fixed_front_masked": effective_audit,
    }


def _print_payload(label, payload):
    xy = payload["xy_m"]
    PETSc.Sys.Print(
        f"  {label}: a_ref={payload['a_ref_equivalent_m_per_yr']:+.6e} "
        f"m/yr, no-forcing dH/dt={payload['no_forcing_dhdt_m_per_yr']:+.6e} "
        f"m/yr at ({xy[0]:.0f}, {xy[1]:.0f})"
    )
    PETSc.Sys.Print(
        "    "
        f"h={payload['thickness_m']:.3f} m, "
        f"H_init={payload['initial_extent_thickness_m']:.3f} m, "
        f"bed={payload['bed_m']:.3f} m, "
        f"surface={payload['surface_m']:.3f} m, "
        f"HAF={payload['height_above_flotation_m']:+.3f} m, "
        f"speed={payload['cell_center_speed_m_per_yr']:.3f} m/yr, "
        f"grounded={payload['grounded']}"
    )


def print_audit(record):
    metadata = record["metadata"]
    PETSc.Sys.Print("Timing cache initial-tendency audit")
    PETSc.Sys.Print(f"  cache: {record['cache_path']}")
    PETSc.Sys.Print(
        f"  mesh: {metadata.get('mesh_basename', 'unknown')} "
        f"({record['cells']} global DG0 cells)"
    )
    PETSc.Sys.Print(
        "  sign: a_ref_equivalent=FV div(h*u); "
        "no-forcing dH/dt=-a_ref_equivalent"
    )
    PETSc.Sys.Print(
        f"  raw net a_ref: {record['raw']['net_gt_per_yr']:+.6f} Gt/yr"
    )
    _print_payload("raw minimum", record["raw"]["min"])
    _print_payload("raw maximum", record["raw"]["max"])
    _print_payload("raw largest magnitude", record["raw"]["absmax"])

    masked = record["fixed_front_masked"]
    PETSc.Sys.Print(
        f"  fixed-front equivalent (H_init < {masked['front_hmin_m']:.3f} m "
        f"zeroed): {masked['masked_cell_count']} cells masked, "
        f"net a_ref={masked['net_gt_per_yr']:+.6f} Gt/yr"
    )
    _print_payload(
        "fixed-front largest magnitude", masked["absmax"]
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True)
    parser.add_argument(
        "--manifest",
        help="optional v3 manifest to validate before opening the cache",
    )
    parser.add_argument(
        "--output",
        help="optional JSON result path; written atomically by rank 0",
    )
    parser.add_argument(
        "--front-hmin",
        type=float,
        default=1.0,
        help="initial thickness threshold for the fixed-front equivalent",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    record = audit_cache(
        args.cache,
        manifest_path=args.manifest,
        front_hmin=args.front_hmin,
    )
    print_audit(record)
    if args.output and COMM_WORLD.rank == 0:
        atomic_write_json(args.output, record)
        PETSc.Sys.Print(f"Timing cache audit record -> {args.output}")
    COMM_WORLD.barrier()


if __name__ == "__main__":
    main()
