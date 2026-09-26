r"""What the preflight gate calls a missing input.

The reader bridges exactly one year past the end of a series (a CESM2-WACCM
atmosphere fetched while the empty 2300 files were withdrawn stops at 2299,
discussion #8), so a tree whose LAST year is absent runs correctly. The gate
has to agree: calling that a missing input reports cores 5 and 7 as BLOCKED
for such a tree. Any other hole stays an error, because the reader raises on
it.

``atm_years`` reads only filenames, so the synthetic trees here are empty
files in the layout ``atmosphere_path`` resolves. The shared mesh/MAP checks
are stubbed out so the status reflects the forcing coverage alone; everything
else, including the printed status line, is the shipped code.
"""
import importlib
import os
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

pytest.importorskip("firedrake")

ESM, SCENARIO = "CESM2-WACCM", "ssp585"
CORE_7 = "core  7"


def _tree(root, years, esm=ESM, scenario=SCENARIO):
    r"""Empty acabf-anomaly files for `years`, in the layout the reader uses."""
    for var in ("acabf-anomaly", "acabf"):
        d = os.path.join(root, esm, scenario, "SDBN1-8000m", var, "v2")
        os.makedirs(d, exist_ok=True)
        for y in years:
            head = f"{var}_AIS_{esm}_{scenario}_SDBN1-8000m_v2_{y}.nc"
            open(os.path.join(d, head), "wb").close()


def _run(monkeypatch, tmp_path, years, ocean_cover=lambda e, s: (2015, 2300)):
    r"""preflight.main() over a tree covering `years`, returning its output."""
    root = str(tmp_path / "ISMIP7" / "AIS")
    os.makedirs(root, exist_ok=True)
    _tree(root, years)
    monkeypatch.setenv("ISMIP7_DATA_ROOT", root)
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    # isolate the forcing-coverage verdict from the mesh/MAP/ocean checks
    monkeypatch.setattr(preflight, "shared_missing", lambda warn: [])
    monkeypatch.setattr(preflight, "racmo_ok", lambda: True)
    monkeypatch.setattr(preflight, "oi_ok", lambda: True)
    monkeypatch.setattr(preflight, "ocean_cover", ocean_cover)
    monkeypatch.setattr(preflight, "pool_status",
                        lambda *a, **k: ("ok", ""))
    preflight.main()
    return None


def _core_line(capsys, core=CORE_7):
    line = [ln for ln in capsys.readouterr().out.splitlines() if core in ln]
    assert len(line) == 1, line
    return line[0]


def test_an_absent_final_year_is_a_note_not_a_blocker(monkeypatch, tmp_path, capsys):
    r"""CESM2-WACCM ssp585 covers 2015-2300 and the tree stops at 2299: the
    reader persists 2299 for 2300, so the core is READY."""
    _run(monkeypatch, tmp_path, range(2015, 2300))       # 2015..2299
    line = _core_line(capsys)
    assert "READY" in line, line
    assert "BLOCKED" not in line
    assert "2300 is absent" in line and "persists 2299" in line


def test_a_hole_inside_the_series_still_blocks(monkeypatch, tmp_path, capsys):
    r"""The reader raises on an interior gap, so the gate must not clear it."""
    years = [y for y in range(2015, 2301) if y != 2100]
    _run(monkeypatch, tmp_path, years)
    line = _core_line(capsys)
    assert "BLOCKED" in line, line
    assert "2100..2100" in line


def test_a_short_tree_still_blocks(monkeypatch, tmp_path, capsys):
    r"""Two years short is not the one-year bridge; it is a short download."""
    _run(monkeypatch, tmp_path, range(2015, 2299))       # 2015..2298
    line = _core_line(capsys)
    assert "BLOCKED" in line, line


def test_a_complete_tree_is_ready_without_a_note(monkeypatch, tmp_path, capsys):
    _run(monkeypatch, tmp_path, range(2015, 2301))       # 2015..2300
    line = _core_line(capsys)
    assert "READY" in line, line
    assert "absent" not in line


# --- cores 9 and 10 melt under the ESM's own ctrl ocean ---------------------

def test_the_control_is_blocked_without_its_esms_ctrl_ocean(monkeypatch, tmp_path, capsys):
    r"""The control reads ``<ESM>/ctrl/ocean`` (icepack/ismip7#107), so that
    tree is what gates cores 9 and 10."""
    asked = []

    def cover(esm, scenario):
        asked.append((esm, scenario))
        return None if scenario == "ctrl" else (2015, 2300)

    _run(monkeypatch, tmp_path, range(2015, 2301), ocean_cover=cover)
    line = _core_line(capsys, "core  9")
    assert "BLOCKED" in line and "CESM2-WACCM/ctrl ocean tf/so" in line, line
    assert ("MRI-ESM2-0", "ctrl") in asked


def test_the_control_is_ready_on_its_ctrl_ocean(monkeypatch, tmp_path, capsys):
    _run(monkeypatch, tmp_path, range(2015, 2301))
    line = _core_line(capsys, "core  9")
    assert "READY" in line, line


# --- core 11 asks for the OCX product the driver will insist on -------------

def _ocx_tree(root, years):
    from icepack2_tools.forcing import OCX_ATMOSPHERE_SOURCE as SRC
    d = os.path.join(root, "OCX", SRC, "SDBN1-8000m", "acabf", "v1")
    os.makedirs(d, exist_ok=True)
    for y in years:
        open(os.path.join(d, f"acabf_AIS_{SRC}_OCX_SDBN1-8000m_v1_{y}.nc"), "wb").close()
    d = os.path.join(root, "OCX", "ocean", "main", "v1")
    os.makedirs(d, exist_ok=True)
    for var in ("tf", "so"):
        open(os.path.join(d, f"{var}_AIS_OCX_ocean_main_v1_1950-2025.nc"), "wb").close()


def _run_core_11(monkeypatch, tmp_path, years, forcing="protocol"):
    root = str(tmp_path / "ISMIP7" / "AIS")
    os.makedirs(root, exist_ok=True)
    if years is not None:
        _ocx_tree(root, years)
    monkeypatch.setenv("ISMIP7_DATA_ROOT", root)
    monkeypatch.setenv("ISMIP7_OCX_FORCING", forcing)
    monkeypatch.delenv("ISMIP7_OCX_OCEAN", raising=False)
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    monkeypatch.setattr(preflight, "shared_missing", lambda warn: [])
    monkeypatch.setattr(preflight, "racmo_ok", lambda: True)
    monkeypatch.setattr(preflight, "oi_ok", lambda: True)
    monkeypatch.setattr(preflight, "pool_status", lambda *a, **k: ("ok", ""))
    preflight.main()


def test_core_11_is_blocked_until_the_ocx_product_is_on_disk(monkeypatch, tmp_path, capsys):
    r"""RACMO and the OI climatology being present used to be enough, which is
    how the core ran on the stopgap with the product sitting unread."""
    _run_core_11(monkeypatch, tmp_path, None)
    line = _core_line(capsys, "core 11")
    assert "BLOCKED" in line and "OCX atmosphere acabf" in line and "absent" in line


def test_core_11_is_ready_on_the_product_and_points_at_the_tripwire(monkeypatch, tmp_path, capsys):
    _run_core_11(monkeypatch, tmp_path, range(1979, 2026))
    line = _core_line(capsys, "core 11")
    assert "READY" in line and "check_melt_bound.py --ocx" in line and "#48" in line


def test_core_11_on_the_stopgap_says_that_is_what_it_is(monkeypatch, tmp_path, capsys):
    _run_core_11(monkeypatch, tmp_path, None, forcing="stopgap")
    line = _core_line(capsys, "core 11")
    assert "READY" in line and "ISMIP7_OCX_FORCING=stopgap" in line


# ── the melt calibration a run reads ────────────────────────────────────────
def _shared(monkeypatch, tmp_path):
    r"""preflight.shared_missing under the environment the test set."""
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    return preflight.shared_missing([])


def test_a_deltat_file_is_the_melt_calibration(monkeypatch, tmp_path):
    r"""With ISMIP7_DELTAT_PER_BASIN_NPZ naming a file the drivers read no
    per-basin K, so the gate must not ask for one."""
    for knob in ("ISMIP7_K_PER_BASIN_NPZ", "ISMIP7_K_SCALE"):
        monkeypatch.delenv(knob, raising=False)
    dT = tmp_path / "deltaT_per_basin_1000_K8.500e-05.npz"
    np.savez(dT, K=8.5e-5)
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", str(dT))
    assert not [m for m in _shared(monkeypatch, tmp_path) if "per-basin" in m]


def test_a_named_deltat_file_that_is_absent_is_reported(monkeypatch, tmp_path):
    monkeypatch.delenv("ISMIP7_K_SCALE", raising=False)
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", str(tmp_path / "nope.npz"))
    missing = _shared(monkeypatch, tmp_path)
    assert any("ISMIP7_DELTAT_PER_BASIN_NPZ" in m and "does not exist" in m
               for m in missing)
