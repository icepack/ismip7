#!/usr/bin/env python3
"""
Shared mesh / boundary-id filename convention.

Every Antarctica mesh is built with a specific outline buffer
(`ISMIP7_BUFFER_M`, meters of ocean the ice outline is pushed into before
meshing - see `icepack2_tools/mesh.py`). The buffer changes the boundary
topology (gmsh physical-line count), so the mesh file and its boundary_ids
sidecar must be tagged with the exact buffer size used: a sidecar built for
one buffer size is not valid for a mesh built with a different one.

This module is the single place that defines those filenames, so the mesh
pipeline (`mesh_antarctica.py`, `make_boundary_ids.py`) and every solver
script that loads a mesh agree on where to find/write them.
"""

import os
import re

from icepack2_tools.runconfig import BUFFER_M_DEFAULT, buffer_m as _buffer_m

MESH_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "mesh")

# Owned by icepack2_tools.runconfig, which the outline extraction reads too.
DEFAULT_BUFFER_M = float(BUFFER_M_DEFAULT)

_ADAPT_RE = re.compile(r"_adapt(\d+)$")
_BUFFER_TAG_RE = re.compile(r"_buffered(\d+)(?=_|$)")


def adapt_lineage(name):
    """Split a mesh name into its unadapted root and its adaptation count.

    `adapt_mesh.py` writes each adapted mesh as `<root>_adapt<N>.msh`, so the
    filename itself records how many adaptations produced it, and the mesh on
    disk is the single source of truth for the count. Every writer of a
    checkpoint records `mesh_basename`, so the count is always recoverable;
    a separate counter attribute would only add a second record to disagree
    with, and a disagreement there lets an adaptation overwrite the mesh it
    was built from.

    Returns `(root, count)`, with `count == 0` for an unadapted mesh.
    """
    stem = os.path.splitext(os.path.basename(name))[0]
    m = _ADAPT_RE.search(stem)
    return (stem[: m.start()], int(m.group(1))) if m else (stem, 0)


def next_adapted_mesh_name(reference, experiment=None):
    """Basename (no extension) of the mesh one adaptation past `reference`.

    The counter comes from `reference`'s own name, so the result is always a
    name that mesh has not used: this is what keeps `adapt_mesh.py` from
    writing over the mesh it is reading.

    `experiment` is the run's own identity, the same string the driver puts in
    front of its checkpoint basenames (the experiment name with the run tag
    already applied). Every adapted mesh and its boundary_ids sidecar land in
    the one shared mesh directory, so without it a control and a projection
    adapting the same starting mesh write the same filenames and the second
    silently replaces the first's triangulation. The run tag alone will not do:
    it is a method-line suffix, so parallel experiments in one line share it.
    The identity is folded into the root once, so the lineage keeps extending
    it rather than repeating it.
    """
    root, count = adapt_lineage(reference)
    experiment = (experiment or "").strip()
    if experiment and not root.endswith(f"_{experiment}"):
        root = f"{root}_{experiment}"
    return f"{root}_adapt{count + 1}"


def buffer_tag(buffer_m):
    """Return the `_buffered<N>` suffix for a given outline buffer (meters)."""
    return f"_buffered{int(float(buffer_m))}"


def mesh_basename(lc_coarse, lc, buffer_m):
    """Basename (no extension) of the mesh built with the given resolution/buffer."""
    return f"antarctica_{lc_coarse}_{lc}{buffer_tag(buffer_m)}"


def mesh_filename(lc_coarse, lc, buffer_m):
    """Full path to the .msh file for the given resolution/buffer."""
    return os.path.join(MESH_DIR, mesh_basename(lc_coarse, lc, buffer_m) + ".msh")


def bndids_filename(lc_coarse, lc, buffer_m):
    """Full path to the boundary_ids sidecar matching the given mesh/buffer."""
    return os.path.join(MESH_DIR, f"boundary_ids_{mesh_basename(lc_coarse, lc, buffer_m)}.json")


def get_buffer_m():
    """Read ISMIP7_BUFFER_M from the environment with the shared default."""
    return _buffer_m()


def mesh_stem(name):
    """A mesh's name without directory or `.msh` extension, the form the
    checkpoints' `mesh_basename` attribute may or may not carry."""
    stem = os.path.basename(str(name))
    return stem[:-len(".msh")] if stem.endswith(".msh") else stem


def buffer_from_name(name):
    """The outline buffer (meters) a mesh name's `_buffered<N>` tag records, or
    None for a name without one (the untagged legacy names, and the adaptive
    preset's `antarctica_ua_*`). Directory and `.msh` extension are ignored."""
    m = _BUFFER_TAG_RE.search(mesh_stem(name))
    return float(m.group(1)) if m else None


def resolve_outline_buffer(first, name, last=None):
    """The outline buffer (meters) a new mesh is built with: ``first`` when
    given, else the `_buffered<N>` tag of ``name``, else ``last``; None when
    none of them says. No default is right for every legacy mesh, so a caller
    with None refuses rather than guessing."""
    for value in (first, buffer_from_name(name), last):
        if value is not None:
            return float(value)
    return None
