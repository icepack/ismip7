r"""The tracked melt calibration every run reads (issue 26).

K50 of IU's rule-based toolbox selection, one K with a thermal-forcing offset
per IMBIE basin, fitted through the forward's own DG0 melt path on the
1000 m / 10 km production mesh. It is tracked, so a fresh clone melts with
it, and the forward refuses to apply it under settings it was not fitted
under.
"""
import importlib
import json
import os
import shutil
import sys

import numpy as np
import pytest

try:
    # dual_friction first: it pulls icepack2 -> irksome, which must be
    # imported before any UFL form is assembled.
    import icepack2_tools.dual_friction  # noqa: F401
except ImportError:
    pass
from icepack2_tools.forcing import (                                  # noqa: E402
    MELT_SLOPE_DEFAULT, SIN_ALPHA_ANT_DEFAULT, check_melt_contract,
    describe_melt_calibration, is_floating, load_deltaT_per_basin,
    melt_receiving,
)
from icepack2_tools.runconfig import (                                # noqa: E402
    GEOMETRY_SPACE_DEFAULT, LC_COARSE_DEFAULT, LC_DEFAULT,
    MELT_CALIBRATION_DEFAULT, RASTER_SAMPLE_DEFAULT, deltat_per_basin_npz,
    file_sha256, melt_calibration_contract, melt_calibration_sidecar,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NPZ = MELT_CALIBRATION_DEFAULT
PRODUCTION_MESH = f"antarctica_{LC_COARSE_DEFAULT}_{LC_DEFAULT}_buffered20000"


@pytest.fixture
def clean(monkeypatch):
    for k in ("ISMIP7_DELTAT_PER_BASIN_NPZ", "ISMIP7_K_PER_BASIN_NPZ",
              "ISMIP7_K_SCALE", "ISMIP7_K_MELT", "ISMIP7_MELT_SLOPE",
              "ISMIP7_SIN_ALPHA_ANT", "ISMIP7_GEOMETRY_SPACE", "ISMIP7_FRACTURE"):
        monkeypatch.delenv(k, raising=False)


def _basins(tmp_path):
    r"""An IMBIE2-style grid: basin 9 on the left half of [0, 3]^2, basin 14
    on the right."""
    xr = pytest.importorskip("xarray")
    x = np.array([0.0, 1.0, 2.0, 3.0])
    bn = np.where(np.arange(4)[None, :] < 2, 9, 14).repeat(4, axis=0)
    path = tmp_path / "basins.nc"
    xr.Dataset({"basinNumber": (("y", "x"), bn)},
               coords={"x": x, "y": x.copy()}).to_netcdf(path)
    return str(path)


def test_the_tracked_file_matches_its_sidecar():
    contract = melt_calibration_contract(NPZ)
    assert os.path.exists(NPZ)
    assert melt_calibration_sidecar(NPZ).endswith(
        "deltaT_per_basin_1000_K6.500e-05.source.json")
    assert contract["file"] == os.path.basename(NPZ)
    assert contract["sha256"] == file_sha256(NPZ)


def test_the_calibration_was_fitted_under_the_run_defaults():
    r"""A run that sets no knob melts under the settings the file was fitted
    under, so the forward applies the melt its calibration was fitted to."""
    contract = melt_calibration_contract(NPZ)
    with np.load(NPZ) as d:
        assert str(d["melt_slope"]) == MELT_SLOPE_DEFAULT == contract["melt_slope"]
        assert float(d["sin_alpha_ant"]) == SIN_ALPHA_ANT_DEFAULT == contract["sin_alpha_ant"]
        assert str(d["geometry_space"]) == GEOMETRY_SPACE_DEFAULT == contract["geometry_space"]
        assert float(d["K"]) == pytest.approx(6.5e-5, rel=1e-12)
        assert float(d["K"]) == pytest.approx(contract["K"], rel=1e-12)
        assert str(d["selected_as"]) == "K50" == contract["selected_as"]
        assert list(d["basin_ids"]) == list(range(16))
        assert float(np.sum(d["M_obs"])) == pytest.approx(1067.4, abs=0.05)
        assert np.all(np.abs(d["deltaT_basin"]) < 3.0)
    assert contract["raster_sample"] == RASTER_SAMPLE_DEFAULT
    assert contract["mesh"] == PRODUCTION_MESH
    assert contract["obs_table"]["total_gtyr"] == 1067.4


def test_a_run_with_nothing_set_reads_the_tracked_file(clean):
    assert deltat_per_basin_npz() == NPZ


def test_two_calibrations_are_refused(clean, monkeypatch, tmp_path):
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", NPZ)
    monkeypatch.setenv("ISMIP7_K_PER_BASIN_NPZ", str(tmp_path / "K.npz"))
    with pytest.raises(ValueError, match="both set"):
        deltat_per_basin_npz()


def test_a_scaled_K_is_refused_with_the_tracked_file(clean, monkeypatch):
    monkeypatch.setenv("ISMIP7_K_SCALE", "1.26")
    with pytest.raises(ValueError, match="ISMIP7_K_SCALE=1.26"):
        deltat_per_basin_npz()


def test_the_removed_scalar_K_knob_is_refused(clean, monkeypatch):
    monkeypatch.setenv("ISMIP7_K_MELT", "8.5e-5")
    with pytest.raises(ValueError, match="ISMIP7_K_MELT is no longer read"):
        deltat_per_basin_npz()


def test_a_file_that_differs_from_its_sidecar_is_refused(clean, monkeypatch, tmp_path):
    npz = tmp_path / os.path.basename(NPZ)
    shutil.copy(NPZ, npz)
    contract = melt_calibration_contract(NPZ)
    contract["sha256"] = "0" * 64
    with open(melt_calibration_sidecar(str(npz)), "w") as f:
        json.dump(contract, f)
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", str(npz))
    with pytest.raises(ValueError, match="does not match the sha256"):
        deltat_per_basin_npz()
    # a named file with no sidecar is read as given
    os.remove(melt_calibration_sidecar(str(npz)))
    assert deltat_per_basin_npz() == str(npz)


def test_the_offsets_load_under_the_run_defaults(clean, tmp_path):
    with np.load(NPZ) as d:
        want = dict(zip(d["basin_ids"].tolist(), d["deltaT_basin"].tolist()))
    field, K = load_deltaT_per_basin(NPZ, np.array([0.4, 2.6, 40.0]),
                                     np.array([0.1, 0.1, 40.0]),
                                     imbie2=_basins(tmp_path))
    assert K == pytest.approx(6.5e-5, rel=1e-12)
    assert np.allclose(field, [want[9], want[14], 0.0])


def test_another_slope_constant_refuses_the_tracked_file(clean, monkeypatch, tmp_path):
    monkeypatch.setenv("ISMIP7_SIN_ALPHA_ANT", "5.7e-3")
    with pytest.raises(ValueError, match="fitted under settings this run does not use"):
        load_deltaT_per_basin(NPZ, np.array([0.4]), np.array([0.1]),
                              imbie2=_basins(tmp_path))


def test_a_run_sampled_otherwise_is_refused():
    r"""The raster sampling moves the draft and the floating set the offsets
    were fitted on; a context that does not record it is not checked."""
    assert check_melt_contract(NPZ, {"raster_sample": "vertex"}) == PRODUCTION_MESH
    assert check_melt_contract(NPZ, {}) == PRODUCTION_MESH
    with pytest.raises(ValueError, match="raster_sample=vertex"):
        check_melt_contract(NPZ, {"raster_sample": "cell_mean"})


def test_the_provenance_line_names_the_file_its_hash_and_both_meshes(clean):
    same, = describe_melt_calibration(PRODUCTION_MESH + ".msh")
    other, = describe_melt_calibration("antarctica_320000_32000.msh")
    assert same.startswith("Forcing provenance: ocean melt calibration ")
    assert file_sha256(NPZ) in same and "K 6.500e-05 (K50)" in same
    assert f"fitted on {PRODUCTION_MESH}; this run's mesh is {PRODUCTION_MESH}" in same
    assert "differs" not in same
    assert "this run's mesh is antarctica_320000_32000, so its integrated melt" in other


def test_melt_falls_on_floating_cells_that_hold_ice():
    r"""Open ocean floats by the flotation test at draft 0; the melt set is
    the calibration's, floating cells holding ice."""
    b = np.array([-500.0, -500.0, 100.0])
    h = np.array([0.0, 300.0, 50.0])
    s = np.maximum(b + h, (1.0 - 917.0 / 1024.0) * h)
    assert is_floating(s, b).tolist() == [True, True, False]
    assert melt_receiving(s, b, h).tolist() == [False, True, False]


def test_the_climatology_callback_melts_no_ice_free_cell(clean, monkeypatch, tmp_path):
    fd = pytest.importorskip("firedrake")
    import icepack2_tools.forcing as forcing
    za = np.array([-1000.0, 0.0])
    monkeypatch.setattr(forcing, "build_oi_climatology_interpolators",
                        lambda data_root=None: {
                            "tf": (lambda pts: np.full(len(pts), 1.5), za),
                            "so": (lambda pts: np.full(len(pts), 34.5), za)})
    monkeypatch.setattr(forcing, "imbie2_basin_path",
                        lambda recorded=None: _basins(tmp_path))
    mesh = fd.RectangleMesh(4, 1, 3.0, 2.0)
    Q_g = fd.FunctionSpace(mesh, "DG", 0)
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    h = fd.Function(Q_g).assign(100.0)
    h.dat.data[0] = 0.0
    b = fd.Function(Q_g).assign(-1000.0)
    s = fd.Function(Q_g).interpolate(fd.max_value(b + h, (1.0 - 917.0 / 1024.0) * h))
    ctx = {"mesh": mesh, "Q": fd.FunctionSpace(mesh, "CG", 1),
           "V": fd.VectorFunctionSpace(mesh, "CG", 1), "Q_g": Q_g,
           "geom_xy": (xy[:, 0].copy(), xy[:, 1].copy()), "h": h, "b": b, "s": s,
           "ocean_melt": fd.Function(Q_g), "raster_sample": "vertex"}
    forcing.make_climatology_ocean_callback()(ctx, 0.0)
    melt = ctx["ocean_melt"].dat.data_ro
    assert melt[0] == 0.0
    assert (melt[1:] != 0.0).all()


def test_preflight_keeps_a_production_core_on_the_calibration_s_mesh(
        clean, monkeypatch, tmp_path):
    pytest.importorskip("firedrake")
    root = tmp_path / "ISMIP7" / "AIS"
    root.mkdir(parents=True)
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(root))
    monkeypatch.syspath_prepend(os.path.join(REPO, "antarctica", "scripts"))
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    monkeypatch.setattr(preflight, "imbie2_basin_path",
                        lambda recorded=None: _basins(tmp_path))
    assert preflight.melt_calibration_missing(f"/m/{PRODUCTION_MESH}.msh") == []
    miss, = preflight.melt_calibration_missing("/m/antarctica_320000_32000.msh")
    assert f"fitted on {PRODUCTION_MESH}" in miss and "ISMIP7_DELTAT_PER_BASIN_NPZ" in miss
    # naming the file, a refit or the tracked one, clears the mesh check
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", NPZ)
    assert preflight.melt_calibration_missing("/m/antarctica_320000_32000.msh") == []
    monkeypatch.setattr(preflight, "imbie2_basin_path",
                        lambda recorded=None: str(tmp_path / "absent.nc"))
    miss, = preflight.melt_calibration_missing(f"/m/{PRODUCTION_MESH}.msh")
    assert "IMBIE2 basin grid" in miss
    sys.modules.pop("preflight", None)
