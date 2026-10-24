r"""The hat-weighted patch mean (icepack2_tools.patchmean): a raster onto CG1
as its lumped L2 projection, mass kept through the DG0 projection, nodata
left out, and the velocity datum with its sigma and representativeness."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import icepack2_tools.dual_friction  # noqa: F401,E402  (icepack2 -> irksome first)

import numpy as np                                                    # noqa: E402
import pytest                                                         # noqa: E402
import rasterio                                                       # noqa: E402
from rasterio.io import MemoryFile                                    # noqa: E402
from rasterio.transform import from_origin                            # noqa: E402
from firedrake import (Function, FunctionSpace, RectangleMesh,        # noqa: E402
                       VectorFunctionSpace, assemble, dx)

from icepack2_tools import patchmean as P                             # noqa: E402
from icepack2_tools import geometry as G                              # noqa: E402

L = 20_000.0      # m, the domain
PIX = 250.0       # m, the raster pixel


PAD = 2           # pixels of raster beyond the mesh on every side


def _raster(arr, nodata=np.nan):
    r"""A float32 GTiff in memory covering [-PAD PIX, L + PAD PIX]^2 at PIX,
    rows from the top (the bilinear fallback needs the mesh inside it)."""
    prof = dict(driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
                dtype="float32", crs="EPSG:3031", nodata=nodata,
                transform=from_origin(-PAD * PIX, L + PAD * PIX, PIX, PIX), tiled=False)
    mf = MemoryFile()
    with mf.open(**prof) as d:
        d.write(arr.astype("float32"), 1)
    return mf.open()


@pytest.fixture
def mesh():
    return RectangleMesh(7, 9, L, L)


@pytest.fixture
def rough():
    n = int(L / PIX) + 2 * PAD
    yy, xx = np.mgrid[0:n, 0:n]
    x = (xx + 0.5 - PAD) * PIX
    y = L + PAD * PIX - (yy + 0.5) * PIX
    return 500.0 + 300.0 * np.sin(x / 700.0) * np.cos(y / 1100.0) + 40.0 * ((xx + yy) % 2)


def _mass_q(ds, mesh, degree):
    Qs, x, y = P.quadrature_points(mesh, degree)
    vals = P.raster_at(ds, x, y)
    f = Function(Qs)
    f.dat.data[:] = np.nan_to_num(vals)
    return float(assemble(f * P._measure((degree, "default"))))


def test_the_dg0_projection_keeps_the_raster_mass(mesh, rough):
    Q0 = FunctionSpace(mesh, "DG", 0)
    Qc = FunctionSpace(mesh, "CG", 1)
    ds = _raster(rough)
    degree = 10
    os.environ["ISMIP7_PATCH_QDEGREE"] = str(degree)
    try:
        H = G.sample_to_geometry(ds, Q0, Qc, method="patch")
    finally:
        del os.environ["ISMIP7_PATCH_QDEGREE"]
    mass_q = _mass_q(ds, mesh, degree)
    assert abs(float(assemble(H * dx)) - mass_q) <= 1e-9 * mass_q
    # the vertex sampling of the same raster does not
    Hv = G.sample_to_geometry(ds, Q0, Qc, method="vertex")
    assert abs(float(assemble(Hv * dx)) - mass_q) > 1e-4 * mass_q


def test_a_cg1_patch_mean_is_the_lumped_l2_projection(mesh, rough):
    Qc = FunctionSpace(mesh, "CG", 1)
    ds = _raster(rough)
    degree = 8
    Qs, x, y = P.quadrature_points(mesh, degree)
    mean, cov, vals = P.raster_patch_mean(ds, Qs, x, y, Qc, degree)
    assert np.allclose(cov.dat.data_ro, 1.0)
    num, den = P.patch_integrals(vals, Qs, Qc, degree)
    m = P.lumped_mass(Qc)
    assert np.allclose(den.dat.data_ro, m.dat.data_ro)
    assert np.allclose(mean.dat.data_ro, num.dat.data_ro / m.dat.data_ro)
    # a constant raster comes back as that constant at every node
    c, _, _ = P.raster_patch_mean(_raster(np.full_like(rough, 123.0)), Qs, x, y, Qc, degree)
    assert np.allclose(c.dat.data_ro, 123.0)


def test_nodata_is_left_out_and_reported_as_coverage(mesh, rough):
    Qc = FunctionSpace(mesh, "CG", 1)
    arr = rough.copy()
    arr[:, : arr.shape[1] // 2] = np.nan            # the left half unobserved
    ds = _raster(arr)
    degree = 8
    Qs, x, y = P.quadrature_points(mesh, degree)
    mean, cov, _ = P.raster_patch_mean(ds, Qs, x, y, Qc, degree, min_coverage=0.5)
    X = mesh.coordinates.dat.data_ro
    left = X[:, 0] < L / 2 - L / 7
    right = X[:, 0] > L / 2 + L / 7
    assert np.all(np.isnan(mean.dat.data_ro[left]))
    assert np.allclose(cov.dat.data_ro[left], 0.0)
    assert np.all(np.isfinite(mean.dat.data_ro[right]))
    assert np.allclose(cov.dat.data_ro[right], 1.0)
    # the covered half's values do not see the nodata
    full, _, _ = P.raster_patch_mean(_raster(rough), Qs, x, y, Qc, degree)
    assert np.allclose(mean.dat.data_ro[right], full.dat.data_ro[right])


def test_the_velocity_datum_carries_patch_sigma_and_representativeness(mesh, rough):
    V = VectorFunctionSpace(mesh, "CG", 1)
    Qc = FunctionSpace(mesh, "CG", 1)
    vx = rough - 500.0
    vy = np.full_like(rough, 7.0)
    ex = np.full_like(rough, 3.0)
    ey = np.full_like(rough, 2.0)
    ex[:, : ex.shape[1] // 2] = 0.0                  # left half: no measurement
    ey[:, : ey.shape[1] // 2] = 0.0                  # (MEaSUREs zeroes both)
    u, (sx, sy), observed, cov, (rx, ry) = P.patch_mean_velocity(
        _raster(vx), _raster(vy), _raster(ex), _raster(ey), V, Qc, degree=8)
    X = mesh.coordinates.dat.data_ro
    right = X[:, 0] > L / 2 + L / 7
    left = X[:, 0] < L / 2 - L / 7
    assert np.all(observed.dat.data_ro[right] == 1.0)
    assert np.all(observed.dat.data_ro[left] == 0.0)
    assert np.allclose(u.dat.data_ro[left], 0.0)
    assert np.allclose(u.dat.data_ro[right, 1], 7.0)
    assert np.allclose(sx.dat.data_ro[right], 3.0)
    assert np.allclose(sy.dat.data_ro[right], 2.0)
    assert np.all(np.isnan(sx.dat.data_ro[left]))
    # a constant component has no sub-patch spread; the rough one does
    assert np.allclose(ry.dat.data_ro[right], 0.0, atol=1e-6)
    assert np.median(rx.dat.data_ro[right]) > 20.0
