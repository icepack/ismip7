r"""Cross-mesh transfer of checkpoint fields with an explicit fill.

A forward that names a MAP with ``ISMIP7_INVERSION`` and a compute mesh with
``ISMIP7_MESH`` interpolates the MAP's continuous fields onto the compute
mesh. Firedrake's cross-mesh ``interpolate`` can only evaluate the source
where the source mesh exists; a target dof outside it is a "missing dof".
The stock behaviour writes ``0.0`` there, which is harmless for the log
controls (theta = phi = 0 is the prior) and fatal for the fluidity prior: the
2 km MAPs were inverted on a buffer-0 mesh and the production mesh carries a
20 km ocean buffer, so every dof in that ring is missing, and
``A_eff = A_prior * exp(phi) = 0`` there zeroes the dislocation term, the
``lin_reg`` regularizer and the ``alpha_gl`` collar in
``dual_friction.build_rc_residual`` at once: a singular membrane block.

``interpolate_with_fill`` makes the fill a stated choice and counts it.

A second artefact comes with the first. Firedrake locates a target point in a
source cell up to ``mesh.tolerance`` (0.5 of a reference cell by default)
outside that cell, and then evaluates the cell's own linear basis there:
a point in the ring within about half a source cell of the ice front is
extrapolated. Where the prior falls steeply towards the front that gave a
transferred fluidity prior spanning [-218.69, 1028.05] from a source spanning
[1.00, 783.69] (the 22 September 2 km MAPs onto the 1 km mesh), and
[-164.7, 921.1] on the 5 km production mesh (issue 81). The transfer
therefore locates strictly: the source mesh's tolerance is set to
``STRICT_TOLERANCE`` for the interpolation and restored afterwards, so a dof
the source does not contain is a missing dof and takes the fill, whatever
its distance from the outline. Behind that, inside a source cell linear
interpolation is a convex combination of the cell's vertex values, so a
value outside the source field's range can only come from extrapolation;
located dofs are still clamped to the source's range, component by
component, and the clamped dofs counted, as a guard on the tolerance.
"""
import contextlib
import numbers

import numpy as np
from firedrake import Function, FunctionSpace, SpatialCoordinate, dot
from mpi4py import MPI

from .mpi_stats import global_count, global_size

def meshes_match(mesh_a, mesh_b, comm=None):
    r"""Whether two meshes are the same mesh: the same global cell and vertex
    counts, and the same cells, compared through two sums over cells of their
    centroids ``c``: of ``c_x c_y`` and of ``|c|^2``. Counts alone take two
    triangulations of one vertex set (a flipped diagonal) for one mesh, and
    so does ``|c|^2`` alone, which is equal for both splits of a square.

    A warm start from a checkpoint on another mesh must not carry that mesh's
    cell-wise geometry across, so the inversion asks this before it decides
    what a warm start supplies; the forward asks it before it transfers a
    checkpoint onto a mesh file of the same name (two sites' builds of the
    production mesh differ, issue 20). Collective."""
    comm = comm if comm is not None else mesh_a.comm
    def counts(m):
        return (comm.allreduce(int(m.cell_set.size), op=MPI.SUM),
                comm.allreduce(int(m.coordinates.dat.data_ro.shape[0]), op=MPI.SUM))
    if counts(mesh_a) != counts(mesh_b):
        return False
    moments_a = _centroid_moments(mesh_a, comm)
    moments_b = _centroid_moments(mesh_b, comm)
    # Both measured against the |c|^2 sum, the size of the terms: on the polar
    # stereographic plane c_x c_y changes sign across the pole and its sum can
    # cancel far below them.
    scale = max(moments_a[1], moments_b[1])
    return all(abs(a - b) <= CENTROID_MOMENT_RTOL * scale
               for a, b in zip(moments_a, moments_b))


#: Relative agreement two identical meshes reach in :func:`meshes_match`'s
#: centroid sums under any partition: the owned cells are the same set, so
#: only the summation order differs, about 1e-13 over millions of cells.
CENTROID_MOMENT_RTOL = 1e-11


def _centroid_moments(mesh, comm):
    r"""Sums over the owned cells of every rank of ``c_x c_y`` and ``|c|^2``
    at the cell centroids ``c``, where a DG0 dof sits."""
    x = SpatialCoordinate(mesh)
    Q0 = FunctionSpace(mesh, "DG", 0)
    return tuple(
        comm.allreduce(float(np.sum(Function(Q0).interpolate(e).dat.data_ro)),
                       op=MPI.SUM)
        for e in (x[0] * x[1], dot(x, x)))


#: Relative point-location tolerance (on the reference cell) while a field is
#: transferred. Zero would let floating-point error push a dof that sits on
#: the shared outline outside; this keeps the outline and rejects anything a
#: source cell does not contain.
STRICT_TOLERANCE = 1e-8


@contextlib.contextmanager
def strict_location(mesh, tolerance=STRICT_TOLERANCE):
    r"""Locate points in ``mesh`` with ``tolerance`` inside the block and
    restore the mesh's own tolerance on exit (setting it clears the spatial
    index, which Firedrake rebuilds on the next location)."""
    previous = mesh.tolerance
    mesh.tolerance = tolerance
    try:
        yield
    finally:
        if isinstance(previous, numbers.Number):
            mesh.tolerance = previous
        else:
            mesh.clear_spatial_index()
            mesh._tolerance = previous


def _source_range(source, comm):
    """Per-component ``(lo, hi)`` of ``source`` over all ranks."""
    data = source.dat.data_ro
    data = data.reshape(data.shape[0], -1)
    width = data.shape[1]
    if data.shape[0]:
        lo, hi = data.min(axis=0), data.max(axis=0)
    else:
        lo, hi = np.full(width, np.inf), np.full(width, -np.inf)
    lo = np.array([comm.allreduce(float(v), op=MPI.MIN) for v in lo])
    hi = np.array([comm.allreduce(float(v), op=MPI.MAX) for v in hi])
    return lo, hi


def interpolate_with_fill(target, source, fill, comm=None):
    r"""Interpolate ``source`` into ``target`` across meshes; dofs of ``target``
    outside the source mesh take ``fill``.

    ``fill`` is a float, or a Function on ``target``'s space whose values are
    taken where the source has none (the raster-sampled velocity_obs, say).
    The source mesh locates strictly for the call (``STRICT_TOLERANCE``,
    restored afterwards), so every dof outside the source outline is a
    missing dof and takes the fill. Located dofs are clamped to the source
    field's own range (per component) behind that: a value beyond it can only
    be an extrapolation from a boundary cell, since interpolation inside a
    cell never leaves the range of its vertex values.
    Returns ``(n_missing, n_total, n_clamped)``, all reduced over ranks and
    counting dofs once (owned dofs only, whatever the value shape). A
    same-mesh call is a plain interpolate and reports nothing missing or
    clamped.
    """
    comm = comm if comm is not None else target.comm
    total = global_size(target, comm)
    if source.function_space().mesh() is target.function_space().mesh():
        target.interpolate(source)
        return 0, total, 0
    with strict_location(source.function_space().mesh()):
        target.interpolate(
            source, allow_missing_dofs=True, default_missing_val=np.nan
        )
    data = target.dat.data
    flat = data.reshape(data.shape[0], -1)
    missing = np.isnan(flat).any(axis=1)
    located = ~missing
    lo, hi = _source_range(source, comm)
    clamped = np.zeros(flat.shape[0], dtype=bool)
    if flat.shape[0]:
        # Count an excursion only past a roundoff margin: a source vertex
        # value reproduced at 1 ulp below the minimum is not an extrapolation.
        margin = 1e-9 * np.maximum(hi - lo, np.maximum(np.abs(lo), np.abs(hi)))
        margin = np.maximum(margin, 1e-12)
        below = (flat[located] < lo - margin).any(axis=1)
        above = (flat[located] > hi + margin).any(axis=1)
        clamped[located] = below | above
        flat[located] = np.clip(flat[located], lo, hi)
    if isinstance(fill, Function):
        if fill.function_space() != target.function_space():
            raise ValueError(
                "the fill Function must live on the target's function space"
            )
        flat[missing] = fill.dat.data_ro.reshape(flat.shape[0], -1)[missing]
    else:
        flat[missing] = float(fill)
    data[...] = flat.reshape(data.shape)
    n_missing = global_count(missing, comm)
    n_clamped = global_count(clamped, comm)
    if global_count(np.isnan(flat).any(axis=1), comm):
        raise RuntimeError(
            f"{target.name()}: NaN left after filling missing dofs; the fill "
            "itself carries NaN"
        )
    return n_missing, total, n_clamped
