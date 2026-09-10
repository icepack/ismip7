r"""Global (cross-rank) statistics of distributed data.

Anything read off ``Function.dat.data_ro`` covers only the dofs the calling
rank OWNS. ``PETSc.Sys.Print`` emits from rank 0 alone, so printing such a
statistic reports rank 0's slice: the number reads narrower than the field
really is and changes with the rank count. Worse, using one as a scalar in a
form (a normalization constant, say) gives each rank a different coefficient,
so the assembled residual and Jacobian disagree across ranks and the answer
depends on the partition -- ``u_c = Constant(u_speed.dat.data_ro.mean())`` did
exactly that in three drivers, and the eikonal solver rescaled the MESH
COORDINATES by a per-rank factor.

Reduce first, through these helpers. Every one is collective: call it on all
ranks or not at all. And the second rule is not optional: a collective is only
safe where the surrounding control flow is itself rank-invariant. Wrapping a
rank-local statistic in an allreduce inside a loop whose length differs per
rank converts a wrong number into a DEADLOCK -- which happened here, in a
``sorted(..., key=lambda n: global_sum(...))`` over a dict whose keys came
from rank-local masks. Agree the iteration set across ranks first, then
reduce in a deterministic order.

See AGENTS.md section 3.
"""

import numpy as np
from mpi4py import MPI


def _comm_of(f, comm):
    return f.function_space().mesh().comm if comm is None else comm


def _data(f):
    return np.asarray(f.dat.data_ro) if hasattr(f, "dat") else np.asarray(f)


def global_range(f, comm=None):
    r"""``(min, max)`` of a Function over all ranks."""
    d = _data(f)
    comm = _comm_of(f, comm)
    return (comm.allreduce(float(d.min()) if d.size else np.inf, op=MPI.MIN),
            comm.allreduce(float(d.max()) if d.size else -np.inf, op=MPI.MAX))


def global_extreme_location(f, coordinates, mode="max", comm=None):
    r"""Return a global scalar extreme and the coordinates where it occurs.

    ``f`` may be a scalar Firedrake Function, or a bare array when ``comm`` is
    supplied. ``coordinates`` must have one row per owned scalar dof. Only one
    candidate per rank is communicated; ties are resolved by rank so every
    rank receives the same deterministic result.
    """
    d = _data(f)
    xy = _data(coordinates)
    comm = _comm_of(f, comm)
    if d.ndim == 2 and d.shape[1] == 1:
        d = d[:, 0]
    if d.ndim != 1:
        raise ValueError("global_extreme_location requires scalar data")
    if xy.ndim != 2 or xy.shape[0] != d.shape[0]:
        raise ValueError(
            "coordinates must have one row per owned scalar dof"
        )
    if mode not in {"min", "max"}:
        raise ValueError("mode must be 'min' or 'max'")

    if d.size:
        index = int(np.argmin(d) if mode == "min" else np.argmax(d))
        candidate = (
            float(d[index]), int(comm.rank), tuple(float(v) for v in xy[index])
        )
    else:
        sentinel = np.inf if mode == "min" else -np.inf
        candidate = (sentinel, int(comm.rank), ())

    candidates = comm.allgather(candidate)
    if mode == "min":
        value, _rank, location = min(candidates, key=lambda item: (item[0], item[1]))
    else:
        value, _rank, location = max(candidates, key=lambda item: (item[0], -item[1]))
    return value, location


def global_max(f, comm=None):
    r"""Maximum of a Function over all ranks."""
    d = _data(f)
    return _comm_of(f, comm).allreduce(
        float(d.max()) if d.size else -np.inf, op=MPI.MAX)


def global_absmax(f, comm=None):
    r"""Largest magnitude over all ranks.

    The extent of a coordinate field is the usual case: a rank-local extent
    would non-dimensionalize each partition by a different factor. Accepts a
    Function, or a bare array with an explicit ``comm``.
    """
    d = _data(f)
    if comm is None and not hasattr(f, "dat"):
        raise TypeError("global_absmax on a bare array needs an explicit comm=")
    return _comm_of(f, comm).allreduce(
        float(np.abs(d).max()) if d.size else 0.0, op=MPI.MAX)


# Array-first alias kept for callers holding mesh.coordinates.dat.data.
def global_max_abs(arr, comm):
    r"""``max(|arr|)`` across ranks for a raw array (explicit comm required)."""
    return global_absmax(arr, comm)


def global_sum(f, comm=None):
    r"""Sum of entries over all ranks. Accepts a Function, or a bare array
    with an explicit ``comm`` (an ndarray carries no mesh to derive one)."""
    d = _data(f)
    if comm is None and not hasattr(f, "dat"):
        raise TypeError("global_sum on a bare array needs an explicit comm=")
    return _comm_of(f, comm).allreduce(float(d.sum()), op=MPI.SUM)


def global_mean(f, comm=None):
    r"""Mean over all ranks' ARRAY ENTRIES, matching ``ndarray.mean()`` on one
    rank: for a vector-valued Function that is the mean over every component,
    not a mean of magnitudes."""
    d = _data(f)
    comm = _comm_of(f, comm)
    total = comm.allreduce(float(d.sum()), op=MPI.SUM)
    n = comm.allreduce(int(d.size), op=MPI.SUM)
    return total / max(n, 1)


def global_size(f, comm=None):
    r"""Total number of owned DOFS across all ranks.

    Counts dofs, not array entries: a vector-valued Function stores
    ``(ndofs, dim)``, so summing ``.size`` would report ``dim`` times too
    many.
    """
    d = np.asarray(f.dat.data_ro)
    n = d.shape[0] if d.ndim else int(d.size)
    return int(_comm_of(f, comm).allreduce(int(n), op=MPI.SUM))


# Older name for the same quantity.
global_dof_count = global_size


def global_count(mask, comm):
    r"""Number of True entries across ranks. Takes an explicit comm because a
    bare boolean ndarray has no mesh to derive one from."""
    return int(comm.allreduce(int(np.count_nonzero(np.asarray(mask)))))
