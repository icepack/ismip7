r"""Hat-weighted patch means: rasters onto CG1 as the lumped L2 projection.

A CG1 nodal value stands for the field over the node's patch, the support of
its hat function, and a DG0 cell value for the cell mean. The vertex sampling
(``geometry.sample_to_geometry`` with ``method="vertex"``) reads one pixel
per vertex and so rests a 2 km cell on three of its sixteen BedMachine
pixels; the velocity observation does the same with one MEaSUReS pixel per
node. Neither conserves the raster's integral and the velocity carries a
sub-cell representativeness error that its formal error does not know about.

The patch mean puts a raster ``f`` on CG1 as

    f_i = ( \int phi_i f dx ) / ( \int phi_i dx ),

the lumped L2 projection, with both integrals taken by a fixed quadrature rule
whose points are where the raster is looked up (a Firedrake ``Quadrature``
element, so the assembly is the ordinary parallel one). It is as smooth as
the vertex sampling (neighbouring cells share two of three patches), and its
DG0 projection conserves the raster's mass exactly: the lumped masses are the
hat integrals, so

    sum_cells |c| (f_i + f_j + f_k) / 3  =  sum_i f_i \int phi_i  =  \int_q f.

Pixels without data (nodata, or outside the raster) are left out of both
integrals, so a node's value is the mean over the covered part of its patch
and ``coverage`` says how much of the patch that is.

For an observation the same operator gives a nodal datum consistent with what
the CG1 model field can represent: the patch mean of the observation, with
the patch mean of the formal variance as its sigma (no 1/N reduction, since
the MEaSUReS errors are correlated over far more than a cell) and the
sub-patch variance of the observation reported as the representativeness
error a point sample incurs.
"""
import os

import numpy as np
from mpi4py import MPI

from firedrake import (
    FiniteElement, Function, FunctionSpace, SpatialCoordinate, TestFunction,
    assemble, dx, interpolate,
)

# Quadrature degree of the patch integrals: the default rule puts about
# (d + 2)^2 / 4 points in a cell (25 at degree 10), of the order of the
# 500 m BedMachine pixels in a 2 km cell and a quarter of those in a 5 km one.
QDEGREE_DEFAULT = 10


def quadrature_points(mesh, degree=QDEGREE_DEFAULT):
    r"""The ``Quadrature`` space of ``degree`` on ``mesh`` and the physical
    coordinates of its points, ``(Qs, x, y)``, one dof per point."""
    Qse = FiniteElement("Quadrature", mesh.ufl_cell(), degree=degree,
                        quad_scheme="default")
    Qs = FunctionSpace(mesh, Qse)
    X = SpatialCoordinate(mesh)
    x = assemble(interpolate(X[0], Qs)).dat.data_ro.copy()
    y = assemble(interpolate(X[1], Qs)).dat.data_ro.copy()
    return Qs, x, y


def raster_at(dataset, x, y):
    r"""``dataset`` looked up by pixel containment at the points ``(x, y)``,
    NaN where the pixel is nodata or outside the raster. One window covering
    the rank's points is read (the whole raster under a scattered
    partition, as ``geometry.raster_cell_mean`` tolerates)."""
    import time
    from rasterio.windows import from_bounds

    out = np.full(x.shape, np.nan)
    if x.size == 0:
        return out
    t0 = time.perf_counter()
    bnd = dataset.bounds
    inside = (x >= bnd.left) & (x <= bnd.right) & (y >= bnd.bottom) & (y <= bnd.top)
    if not inside.any():
        return out
    rx, ry = abs(dataset.res[0]), abs(dataset.res[1])
    win = from_bounds(max(x[inside].min() - 2 * rx, bnd.left),
                      max(y[inside].min() - 2 * ry, bnd.bottom),
                      min(x[inside].max() + 2 * rx, bnd.right),
                      min(y[inside].max() + 2 * ry, bnd.top),
                      transform=dataset.transform)
    win = win.round_lengths(op="ceil").round_offsets(op="floor")
    arr = dataset.read(1, window=win, masked=True)
    arr = np.ma.filled(arr.astype("f8"), np.nan)
    t1 = time.perf_counter()
    # pixel containment from the affine transform directly (rasterio's
    # rowcol loops in Python in some versions; this is 46 M points at 2 km)
    T = dataset.transform
    inv = ~T
    xi, yi = x[inside], y[inside]
    fc = inv.a * xi + inv.b * yi + inv.c
    fr = inv.d * xi + inv.e * yi + inv.f
    r = np.floor(fr).astype(np.int64) - int(win.row_off)
    c = np.floor(fc).astype(np.int64) - int(win.col_off)
    ok = (r >= 0) & (r < arr.shape[0]) & (c >= 0) & (c < arr.shape[1])
    vals = np.full(r.shape, np.nan)
    vals[ok] = arr[r[ok], c[ok]]
    out[inside] = vals
    t2 = time.perf_counter()
    if os.environ.get("ISMIP7_PATCH_TIMING", "0") == "1":
        print(f"  [patchmean] {dataset.name.split(':')[-1]}: window {arr.shape} read "
              f"{t1 - t0:.1f}s, {xi.size} lookups {t2 - t1:.1f}s", flush=True)
    return out


def _measure(scheme):
    return dx(degree=scheme[0], scheme=scheme[1])


def patch_integrals(values, Qs, Q_cg, degree):
    r"""``(num, den)`` CG1 Functions: the hat-weighted integrals of
    ``values`` (NaN = no datum) and of its coverage, by the rule of ``Qs``."""
    v = TestFunction(Q_cg)
    f = Function(Qs)
    w = Function(Qs)
    valid = np.isfinite(values)
    f.dat.data[:] = np.where(valid, values, 0.0)
    w.dat.data[:] = valid.astype(float)
    meas = _measure((degree, "default"))
    num = assemble(f * v * meas)
    den = assemble(w * v * meas)
    return num, den


def patch_mean(values, Qs, Q_cg, degree, min_coverage=0.0):
    r"""The patch mean of ``values`` on ``Q_cg`` as ``(mean, coverage)``.

    ``mean`` is NaN where the covered fraction of the patch is at or below
    ``min_coverage`` (a caller fills it); ``coverage`` is that fraction.
    Collective."""
    num, den = patch_integrals(values, Qs, Q_cg, degree)
    m = lumped_mass(Q_cg)
    out = Function(Q_cg)
    cov = Function(Q_cg, name="coverage")
    d = den.dat.data_ro
    with np.errstate(invalid="ignore", divide="ignore"):
        cov.dat.data[:] = np.where(m.dat.data_ro > 0, d / m.dat.data_ro, 0.0)
        val = np.where(d > 0, num.dat.data_ro / np.where(d > 0, d, 1.0), np.nan)
    val[cov.dat.data_ro <= min_coverage] = np.nan
    out.dat.data[:] = val
    return out, cov


def lumped_mass(Q_cg):
    r"""``\int phi_i dx`` on ``Q_cg``, exact (the rule is at least degree 1)."""
    return assemble(TestFunction(Q_cg) * dx)


def raster_patch_mean(dataset, Qs, x, y, Q_cg, degree, min_coverage=0.0):
    r"""``patch_mean`` of a rasterio ``dataset`` looked up at the points of
    ``Qs``; also returns the point values for a second use."""
    vals = raster_at(dataset, x, y)
    mean, cov = patch_mean(vals, Qs, Q_cg, degree, min_coverage)
    return mean, cov, vals


def fill_uncovered(mean, fallback):
    r"""Replace NaN nodes of ``mean`` with ``fallback``'s values; returns the
    global count filled."""
    bad = ~np.isfinite(mean.dat.data_ro)
    n = int(mean.function_space().mesh().comm.allreduce(int(bad.sum()), op=MPI.SUM))
    if n:
        mean.dat.data[bad] = fallback.dat.data_ro[bad]
    return n


def patch_mean_velocity(vx_ds, vy_ds, ex_ds, ey_ds, V, Q_cg, degree=QDEGREE_DEFAULT,
                        min_coverage=0.5):
    r"""MEaSUReS onto the model's spaces as patch means.

    Returns ``u_obs`` (on ``V``, zero where unobserved), the per-component
    sigma ``(sig_x, sig_y)`` on ``Q_cg`` as the patch root-mean-square of the
    formal error (NaN where unobserved), the ``observed`` mask (covered
    fraction above ``min_coverage``, as a CG1 Function of 0/1), the
    ``coverage``, and ``(repr_x, repr_y)``: the sub-patch standard deviation
    of each velocity component, the representativeness error a point sample
    carries. A pixel counts as observed where MEaSUReS reports an error, as
    in the point sampling (ERR > 0 marks a measured pixel)."""
    mesh = V.mesh()
    Qs, x, y = quadrature_points(mesh, degree)
    vx = raster_at(vx_ds, x, y)
    vy = raster_at(vy_ds, x, y)
    ex = raster_at(ex_ds, x, y)
    ey = raster_at(ey_ds, x, y)
    emag = np.hypot(np.nan_to_num(ex), np.nan_to_num(ey))
    has = np.isfinite(vx) & np.isfinite(vy) & (emag > 0.0)
    nan = np.nan
    vx = np.where(has, vx, nan); vy = np.where(has, vy, nan)
    ex = np.where(has, ex, nan); ey = np.where(has, ey, nan)

    mx, cov = patch_mean(vx, Qs, Q_cg, degree, min_coverage)
    my, _ = patch_mean(vy, Qs, Q_cg, degree, min_coverage)
    mx2, _ = patch_mean(vx ** 2, Qs, Q_cg, degree, min_coverage)
    my2, _ = patch_mean(vy ** 2, Qs, Q_cg, degree, min_coverage)
    sx2, _ = patch_mean(ex ** 2, Qs, Q_cg, degree, min_coverage)
    sy2, _ = patch_mean(ey ** 2, Qs, Q_cg, degree, min_coverage)

    observed = Function(Q_cg, name="observed")
    obs = np.isfinite(mx.dat.data_ro) & (cov.dat.data_ro > min_coverage)
    observed.dat.data[:] = obs.astype(float)

    u_obs = Function(V, name="velocity_obs")
    u_obs.dat.data[:, 0] = np.where(obs, mx.dat.data_ro, 0.0)
    u_obs.dat.data[:, 1] = np.where(obs, my.dat.data_ro, 0.0)

    def _sqrt(f):
        g = Function(Q_cg)
        g.dat.data[:] = np.sqrt(np.maximum(f.dat.data_ro, 0.0))
        g.dat.data[~obs] = np.nan
        return g

    sig_x, sig_y = _sqrt(sx2), _sqrt(sy2)
    sig_x.rename("sigma_ux"); sig_y.rename("sigma_uy")
    vrx = Function(Q_cg); vry = Function(Q_cg)
    vrx.dat.data[:] = mx2.dat.data_ro - mx.dat.data_ro ** 2
    vry.dat.data[:] = my2.dat.data_ro - my.dat.data_ro ** 2
    repr_x, repr_y = _sqrt(vrx), _sqrt(vry)
    repr_x.rename("repr_ux"); repr_y.rename("repr_uy")
    return u_obs, (sig_x, sig_y), observed, cov, (repr_x, repr_y)
