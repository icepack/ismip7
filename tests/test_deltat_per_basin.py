r"""The protocol's per-basin adjustment is a thermal-forcing offset at one
toolbox K, fitted on the forward's melt path and stamped onto any mesh."""
import importlib.util
import os

import numpy as np
import pytest
from mpi4py import MPI

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _calibrate_deltaT():
    spec = importlib.util.spec_from_file_location(
        "calibrate_deltaT",
        os.path.join(REPO, "antarctica", "scripts", "calibrate_deltaT.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _synthetic(K=8.5e-5, true_dT=(0.3, -0.5)):
    from icepack2_tools.forcing import quadratic_mixed_slope, _RHO_I
    rng = np.random.default_rng(0)
    n = 400
    tf = rng.uniform(0.5, 2.5, n)
    sal = np.full(n, 34.5)
    sin_a = np.full(n, 5.115e-3)
    area = np.full(n, 4.0e6)
    basin = np.repeat([3, 9], n // 2)
    floating = np.ones(n, bool)
    bids = np.array([3, 9])
    M_obs = np.zeros(2)
    for i, (bid, d) in enumerate(zip(bids, true_dT)):
        sel = basin == bid
        m = quadratic_mixed_slope(tf[sel] + d, sal[sel], sin_a[sel], K=K)
        M_obs[i] = float((m * area[sel]).sum()) * float(_RHO_I) / 1e12
    return dict(tf=tf, sal=sal, sin_a=sin_a, K=K, floating=floating,
                area=area, basin=basin, bids=bids, M_obs=M_obs)


def test_the_fit_recovers_a_known_offset():
    cd = _calibrate_deltaT()
    d = _synthetic()
    dT, M0, resid, sens, flagged = cd.fit_deltaT(
        d["tf"], d["sal"], d["sin_a"], d["K"], d["floating"], d["area"],
        d["basin"], d["bids"], d["M_obs"], MPI.COMM_WORLD)
    assert np.allclose(dT, (0.3, -0.5), atol=2e-4)
    assert np.all(M0 > 0) and np.all(np.sign(M0 - d["M_obs"]) == (-1, 1))
    assert np.all(np.abs(resid) < 1e-3)  # Gt/yr, the root tolerance in dT
    assert np.all(sens > 0)
    assert flagged == []


def test_an_unreachable_basin_takes_the_window_end_and_is_flagged():
    cd = _calibrate_deltaT()
    d = _synthetic()
    d["M_obs"][0] *= 50.0
    dT, M0, resid, sens, flagged = cd.fit_deltaT(
        d["tf"], d["sal"], d["sin_a"], d["K"], d["floating"], d["area"],
        d["basin"], d["bids"], d["M_obs"], MPI.COMM_WORLD)
    assert dT[0] == cd.DT_WINDOW[1]
    assert [b for b, _ in flagged] == [3]


class _TwinComm:
    r"""Two ranks holding identical cells: every reduction doubles."""
    def allreduce(self, x, op=MPI.SUM):
        return 2 * x if op == MPI.SUM else x


def test_the_fit_reduces_basin_totals_across_ranks():
    cd = _calibrate_deltaT()
    d = _synthetic()
    # the observation covers both ranks' cells; each rank sees only its own
    dT, M0, *_ = cd.fit_deltaT(
        d["tf"], d["sal"], d["sin_a"], d["K"], d["floating"], d["area"],
        d["basin"], d["bids"], 2 * d["M_obs"], _TwinComm())
    assert np.allclose(dT, (0.3, -0.5), atol=2e-4)


@pytest.fixture
def clean(monkeypatch):
    import icepack2_tools.forcing as forcing
    for k in ("ISMIP7_MELT_SLOPE", "ISMIP7_SIN_ALPHA_ANT", "ISMIP7_GEOMETRY_SPACE",
              "ISMIP7_K_SCALE", "ISMIP7_DELTAT_PER_BASIN_NPZ", "ISMIP7_FRACTURE",
              "ISMIP7_K_PER_BASIN_NPZ", "ISMIP7_K_MELT"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setattr(forcing, "_MELT_SLOPE_WARNED", False)
    monkeypatch.setattr(forcing, "_SLOPE_CAP_WARNED", False)
    monkeypatch.setattr(forcing, "_GEOMETRY_SPACE_WARNED", False)


def _offsets(tmp_path, **extra):
    r"""A basin grid with basin 3 on the left half of [0, 3]^2 and basin 9 on
    the right, and an offsets file for basins 3 and 9 pointing at it."""
    xr = pytest.importorskip("xarray")
    from icepack2_tools.runconfig import geometry_space
    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 1.0, 2.0, 3.0])
    bn = np.where(np.arange(4)[None, :] < 2, 3, 9).repeat(4, axis=0)
    imbie = tmp_path / "basins.nc"
    xr.Dataset({"basinNumber": (("y", "x"), bn)},
               coords={"x": x, "y": y}).to_netcdf(imbie)
    npz = tmp_path / "deltaT.npz"
    entries = dict(basin_ids=np.array([3, 9]), deltaT_basin=np.array([0.3, -0.5]),
                   K=8.5e-5, melt_slope="ant", sin_alpha_ant=5.115e-3,
                   geometry_space=geometry_space(), imbie2_nc=str(imbie))
    entries.update(extra)
    np.savez(npz, **entries)
    return str(npz)


def test_the_offset_is_stamped_by_basin(tmp_path, clean, capsys):
    from icepack2_tools.forcing import load_deltaT_per_basin
    npz = _offsets(tmp_path)
    mx = np.array([0.4, 2.6, 0.2, 2.9, 40.0])
    my = np.array([0.1, 0.1, 2.9, 2.9, 40.0])
    field, K = load_deltaT_per_basin(npz, mx, my)
    assert K == 8.5e-5
    assert np.allclose(field, [0.3, -0.5, 0.3, -0.5, 0.0])
    assert "WARNING" not in capsys.readouterr().out


def test_another_slope_constant_is_refused_with_the_deltaT_remedy(tmp_path, clean):
    r"""The forward has to apply the melt its offsets were fitted to, so an
    offsets file fitted under another slope constant stops the run."""
    from icepack2_tools.forcing import load_deltaT_per_basin
    npz = _offsets(tmp_path, sin_alpha_ant=5.7e-3)
    with pytest.raises(ValueError, match=r"sin\(alpha\) 0.0057") as e:
        load_deltaT_per_basin(npz, np.array([0.4]), np.array([0.1]))
    said = str(e.value)
    assert "calibrate_deltaT.py" in said and "ISMIP7_DELTAT_PER_BASIN_NPZ" in said
    assert "calibrate_melt.py" not in said


def test_another_geometry_is_refused_with_the_deltaT_remedy(tmp_path, clean):
    from icepack2_tools.forcing import load_deltaT_per_basin
    from icepack2_tools.runconfig import geometry_space
    other = "cg1" if geometry_space() == "dg0" else "dg0"
    npz = _offsets(tmp_path, geometry_space=other)
    with pytest.raises(ValueError, match=f"{other} geometry") as e:
        load_deltaT_per_basin(npz, np.array([0.4]), np.array([0.1]))
    assert "calibrate_deltaT.py" in str(e.value)


def test_another_slope_law_is_refused(tmp_path, clean):
    from icepack2_tools.forcing import load_deltaT_per_basin
    npz = _offsets(tmp_path, melt_slope="local")
    with pytest.raises(ValueError, match="ISMIP7_MELT_SLOPE local"):
        load_deltaT_per_basin(npz, np.array([0.4]), np.array([0.1]))


def test_a_missing_offsets_file_is_refused(tmp_path, clean, monkeypatch):
    from icepack2_tools.runconfig import deltat_per_basin_npz
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", str(tmp_path / "nope.npz"))
    with pytest.raises(FileNotFoundError, match="nope.npz"):
        deltat_per_basin_npz()


def test_a_scaled_K_is_refused_with_offsets(tmp_path, clean, monkeypatch):
    from icepack2_tools.runconfig import deltat_per_basin_npz
    npz = _offsets(tmp_path)
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", npz)
    monkeypatch.setenv("ISMIP7_K_SCALE", "1.26")
    with pytest.raises(ValueError, match="ISMIP7_K_SCALE=1.26") as e:
        deltat_per_basin_npz()
    assert "ISMIP7_DELTAT_PER_BASIN_NPZ" in str(e.value)
    monkeypatch.setenv("ISMIP7_K_SCALE", "1.0")
    assert deltat_per_basin_npz() == npz


class _Ocean:
    def get_thermal_forcing(self, yr, x, y, draft=None):
        return np.full(len(x), 1.5)

    def get_salinity(self, yr, x, y, draft=None):
        return np.full(len(x), 34.5)


def test_a_projection_melts_with_the_offsets_and_never_reads_the_K_file(
        tmp_path, clean, monkeypatch, capsys):
    fd = pytest.importorskip("firedrake")
    import icepack2_tools.forcing as forcing
    npz = _offsets(tmp_path)
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", npz)
    monkeypatch.setattr(forcing, "compute_sin_alpha",
                        lambda ctx: np.full(len(ctx["h"].dat.data_ro), 5.115e-3))
    # nearest basin cell: basin 3 for x < 1.5, basin 9 up to x = 3, then off
    # the grid, where the offset is zero and the ice still melts at K
    mesh = fd.RectangleMesh(4, 1, 6.0, 2.0)
    Q_g = fd.FunctionSpace(mesh, "DG", 0)
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    h = fd.Function(Q_g).assign(100.0)
    # the last cell is open ocean: it floats by the flotation test and holds
    # no ice, so the forward books no melt there, as the calibration fits none
    h.dat.data[-1] = 0.0
    ctx = {"mesh": mesh, "geom_xy": (xy[:, 0].copy(), xy[:, 1].copy()),
           "h": h, "b": fd.Function(Q_g).assign(-1000.0),
           "s": fd.Function(Q_g).assign(5.0),
           "ocean_melt": fd.Function(Q_g)}
    callback = forcing.make_forcing_callback(
        ocean=_Ocean(), K_per_basin_npz=str(tmp_path / "absent_K.npz"))
    callback(ctx, 2016.0)

    x = ctx["geom_xy"][0]
    dT = np.where(x < 1.5, 0.3, np.where(x <= 3.0, -0.5, 0.0))
    assert (dT == 0.0).any()
    expect = forcing.quadratic_mixed_slope(1.5 + dT, 34.5, 5.115e-3, K=8.5e-5)
    expect[-1] = 0.0
    assert np.allclose(ctx["ocean_melt"].dat.data_ro, expect)
    out = capsys.readouterr().out
    assert "outside the fitted basins" in out
