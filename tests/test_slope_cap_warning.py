r"""A K fitted against a capped draft slope announces itself once per run.

calibrate_melt.py records the cap it applied and load_K_per_basin compares
it with the forward's own slope, which applies none. Measured on the 2500 m
mesh, the uncapped cell slope integrates 3.7 times the capped melt at K = 1,
so the warning names that and the two consistent choices. It is printed
once, from rank 0, and not at all for a file fitted without a cap. It applies
to a legacy per-basin K named with ISMIP7_K_PER_BASIN_NPZ; an offsets file
fitted under another convention is refused instead (test_deltat_per_basin).
"""

import os

import numpy as np
import pytest

pytest.importorskip("firedrake")
xr = pytest.importorskip("xarray")

from icepack2_tools import forcing                                # noqa: E402

CAP = 5e-3


@pytest.fixture(autouse=True)
def _rearm(monkeypatch):
    r"""The warning fires once per process, so re-arm it around each test.
    The cap belongs to the local slope convention, so these run under it."""
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", "local")
    # The files here carry no geometry tag either; that warning has its own tests.
    monkeypatch.setenv("ISMIP7_GEOMETRY_SPACE", "cg1")
    monkeypatch.setattr(forcing, "_SLOPE_CAP_WARNED", False)
    monkeypatch.setattr(forcing, "_MELT_SLOPE_WARNED", False)
    monkeypatch.setattr(forcing, "_GEOMETRY_SPACE_WARNED", False)
    yield


@pytest.fixture
def basins(tmp_path, monkeypatch):
    r"""A 2 by 2 IMBIE2 grid with basins 1 and 2 side by side."""
    root = tmp_path / "ISMIP7" / "AIS"
    d = root / "parameterisations" / "ocean" / "imbie2"
    os.makedirs(d)
    xr.Dataset(
        {"basinNumber": (("y", "x"), np.array([[1, 2], [1, 2]]))},
        coords={"x": [0.0, 8000.0], "y": [0.0, 8000.0]},
    ).to_netcdf(d / "basin_numbers_ismip8km_v2.nc")
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(root))


def _npz(path, **extra):
    np.savez(path, basin_ids=np.array([1, 2]),
             K_basin=np.array([1e-4, 2e-4]), **extra)
    return str(path)


def _load(path):
    return forcing.load_K_per_basin(path, [0.0, 8000.0], [0.0, 0.0])


def test_a_capped_calibration_is_announced(basins, tmp_path, capsys):
    path = _npz(tmp_path / "calibrated_K_per_basin_2000.npz",
                sin_alpha_cap=CAP)
    K = _load(path)
    said = capsys.readouterr().out

    np.testing.assert_allclose(K, [1e-4, 2e-4])
    assert "calibrated_K_per_basin_2000.npz" in said
    assert "0.005" in said
    # The reader has to be able to find out what to do about it.
    assert "GEOMETRY_DISCRETIZATION" in said


def test_it_says_so_once(basins, tmp_path, capsys):
    path = _npz(tmp_path / "a.npz", sin_alpha_cap=CAP)
    _load(path)
    assert capsys.readouterr().out.strip()

    # A chained resume or a driver that reloads must not turn this into a
    # per-step banner.
    _load(path)
    assert capsys.readouterr().out == ""


def test_a_file_without_the_stamp_is_not_flagged(basins, tmp_path, capsys):
    r"""An older npz carries no cap, so there is nothing to compare and the
    loader stays quiet and still returns the K field."""
    path = _npz(tmp_path / "old.npz")
    K = _load(path)

    np.testing.assert_allclose(K, [1e-4, 2e-4])
    assert capsys.readouterr().out == ""
