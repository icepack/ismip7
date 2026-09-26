r"""The roots a second checkout does not carry.

A cluster can hold more than one checkout of this repository. The code moves
with the invocation; the forcing tree, the meshes, the MAPs and the
observational rasters do not - they are gitignored, so a fresh clone has none
of them. The melt calibration is the exception: it is tracked
(antarctica/calibration), so every clone melts with it. These check the
Python half of that split (site_env.sh and tests/test_site_core.py cover the
shell half): each root follows its own variable, and otherwise resolves into
the invoking checkout, so a site holding a single checkout is provably
unaffected.
"""
import os

import pytest

from icepack2_tools.runconfig import (
    MELT_CALIBRATION_DEFAULT, deltat_per_basin_npz, k_per_basin_npz,
    obs_data_root,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANT = os.path.join(REPO, "antarctica")


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for key in ("ISMIP7_OBS_DATA_ROOT", "ISMIP7_K_PER_BASIN_NPZ",
                "ISMIP7_DELTAT_PER_BASIN_NPZ", "ISMIP7_K_SCALE",
                "ISMIP7_DATA_ROOT", "ISMIP7_OBS_KIT"):
        monkeypatch.delenv(key, raising=False)


def test_the_obs_root_is_this_checkout_when_nothing_says_otherwise():
    assert obs_data_root() == os.path.join(ANT, "data")


def test_the_obs_root_follows_its_variable(monkeypatch):
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", "/projects/shared/antarctica/data")
    assert obs_data_root() == "/projects/shared/antarctica/data"


def test_the_melt_calibration_is_in_this_checkout():
    """A fresh clone melts with the tracked calibration: no results/ directory
    is searched, so a machine holding an old per-basin K file melts the same."""
    assert MELT_CALIBRATION_DEFAULT == os.path.join(
        ANT, "calibration", "deltaT_per_basin_1000_K6.500e-05.npz")
    assert deltat_per_basin_npz() == MELT_CALIBRATION_DEFAULT
    assert k_per_basin_npz() is None


def test_a_per_basin_K_is_read_only_when_named(tmp_path, monkeypatch):
    """The override names one file and replaces the tracked calibration; a name
    that does not exist is refused rather than skipped."""
    k_file = tmp_path / "mine.npz"
    k_file.write_bytes(b"")
    monkeypatch.setenv("ISMIP7_K_PER_BASIN_NPZ", str(k_file))
    assert k_per_basin_npz() == str(k_file)
    assert deltat_per_basin_npz() is None
    monkeypatch.setenv("ISMIP7_K_PER_BASIN_NPZ", str(tmp_path / "absent.npz"))
    with pytest.raises(FileNotFoundError, match="absent.npz"):
        k_per_basin_npz()


def test_the_inversion_driver_reads_the_obs_root(monkeypatch):
    """It is the driver that did NOT read it: submissions from a second Rice
    checkout died at 'Loading data...' looking for BedMachine while the shell
    layer had named the right tree."""
    import importlib
    import icepack2_tools.runconfig as rcmod
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", "/projects/shared/antarctica/data")
    importlib.reload(rcmod)
    assert rcmod.obs_data_root() == "/projects/shared/antarctica/data"


# ── The library half of the same split ──────────────────────────────────
# A checkout that carries its own data (or symlinks to one that does) resolves
# these whether or not they consult ISMIP7_OBS_DATA_ROOT, so a default built
# from __file__ looks correct everywhere it is usually run. It is wrong exactly
# once: a second clone with no data beside it, which is how the Rice 2 km
# inversions were staged. Reaching for the variable is what makes them agree.

def test_racmo_reads_the_obs_root_not_the_checkout(monkeypatch, tmp_path):
    from icepack2_tools.forcing import load_racmo_smb_climatology
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", str(tmp_path))
    # The path is resolved and opened before the function space is touched,
    # so the error names the file it went looking for.
    with pytest.raises(FileNotFoundError) as excinfo:
        load_racmo_smb_climatology(None)
    assert str(tmp_path) in str(excinfo.value)


def test_the_dhdt_cache_lands_under_the_obs_root(monkeypatch, tmp_path):
    from icepack2_tools.obs_dhdt import _cache_rasters
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", str(tmp_path))
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(tmp_path/"no-such-forcing-tree"))
    with pytest.raises(FileNotFoundError):
        _cache_rasters("dhdt_smith")
    assert (tmp_path/"dhdt_cache").is_dir()


def test_bedmachine_for_a_remesh_is_read_from_the_obs_root(monkeypatch, tmp_path):
    """The third library default, the one the adaptive remesh builds its
    outline from. The builder resolves and opens BedMachine before it touches
    gmsh, so the error names the directory it went looking in."""
    from icepack2_tools.adapt_mesh import antarctica_geometry_builder
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", str(tmp_path))
    with pytest.raises(FileNotFoundError) as excinfo:
        antarctica_geometry_builder(None)
    assert str(tmp_path) in str(excinfo.value)


def test_the_adapt_driver_reads_the_obs_root(monkeypatch):
    """The only caller of the builder above. It globs BedMachine and MEaSUReS
    out of its own module-level DATA_DIR, so a data-less clone died on an
    empty glob before the library default could help."""
    import importlib
    monkeypatch.syspath_prepend(os.path.join(REPO, "antarctica", "scripts"))
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", "/projects/shared/antarctica/data")
    driver = importlib.reload(importlib.import_module("adapt_mesh"))
    assert driver.DATA_DIR == "/projects/shared/antarctica/data"


def test_downloads_land_under_the_obs_root(monkeypatch, tmp_path):
    r"""download_data.py writes where every run reads: the obs root, not the
    checkout, once a site names one."""
    import importlib
    import sys
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", str(tmp_path / "shared"))
    monkeypatch.syspath_prepend(os.path.join(ANT, "scripts"))
    sys.modules.pop("download_data", None)
    mod = importlib.import_module("download_data")
    assert str(mod.DATA_DIR) == str(tmp_path / "shared")
    monkeypatch.delenv("ISMIP7_OBS_DATA_ROOT")
    sys.modules.pop("download_data")
    assert str(importlib.import_module("download_data").DATA_DIR) == os.path.join(ANT, "data")
