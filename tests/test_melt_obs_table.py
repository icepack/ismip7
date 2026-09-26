r"""The per-basin melt observation table: two published widths, one reader.

The ISMIP7 melt-calibration product ships ``Melt_Paolo_Davison_Adusumilli
_imbie2.csv`` with three columns (basin, melt, uncertainty). The table it
replaces, ``Melt_Paolo_Err_Adusumilli_imbie2_v3.csv``, carries area and two
per-area columns between the melt and its uncertainty. The replaced reader took
``row[3]``, which is the uncertainty in the old six-column table and past the
end of the new three-column one, where it raises IndexError. Any index chosen
for one width reads the wrong quantity or nothing at the other: ``row[2]``,
right for the new table, reads AREA as the uncertainty of every old-table basin,
values near 1e5 km^2 where the uncertainty is tens of Gt/yr. Columns are located
by header name instead, and this pins that.

The integrated targets differ by 23% (865.0 against 1067.4 Gt/yr), so which
table a calibration used is part of its provenance and belongs in the saved
npz. The second half pins the default search order, the ``[!]`` line printed
when that search ends on the older table, and ``ISMIP7_K_OUT``, which writes a
calibration away from the path every run reads.

Serial and needs no data files. Importing ``calibrate_melt`` loads Firedrake,
rasterio and icepack, so the tests skip where that stack is absent.
"""

import csv
import os
import sys

import numpy as np
import pytest

pytest.importorskip("firedrake")
pytest.importorskip("rasterio")
pytest.importorskip("icepack")

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SCRIPTS = os.path.join(_ROOT, "antarctica", "scripts")

OLD_HEADER = ("", "BMR (Gt/yr)", "Area (km^2)", "BMR uncert (Gt/yr)",
              "Average BMR (kg/m2/a)", "Average BMR uncert (kg/m2/a)")
NEW_HEADER = ("", "BMR (Gt/yr)", "BMR uncert (Gt/yr)")

# basin, melt, uncertainty, area; the area column exists only in the old layout
# and is deliberately far from the uncertainty, so reading it as the
# uncertainty fails the bound below.
ROWS = ((0, 39.85, 51.75, 124527.1), (1, 3.30, 3.50, 5455.0))


def _write(path, header, old_layout):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for bid, melt, sigma, area in ROWS:
            if old_layout:
                w.writerow([bid, melt, area, sigma, melt / area * 1e9,
                            sigma / area * 1e9])
            else:
                w.writerow([bid, melt, sigma])


def _load_obs_with(monkeypatch, path):
    r"""Import ``calibrate_melt`` against ``path`` and return its parse.

    The module resolves OBS_CSV at import, so the environment has to carry the
    path before the import and the module has to be dropped afterwards.
    """
    monkeypatch.setenv("ISMIP7_MELT_OBS_CSV", str(path))
    monkeypatch.syspath_prepend(_SCRIPTS)
    sys.modules.pop("calibrate_melt", None)
    try:
        import calibrate_melt
        return calibrate_melt._load_obs()
    finally:
        sys.modules.pop("calibrate_melt", None)


@pytest.mark.parametrize("old_layout", (True, False))
def test_both_published_layouts_read_the_same(tmp_path, monkeypatch,
                                              old_layout):
    path = tmp_path / ("old.csv" if old_layout else "new.csv")
    _write(path, OLD_HEADER if old_layout else NEW_HEADER, old_layout)

    bids, melt, sigma = _load_obs_with(monkeypatch, path)

    assert list(bids) == [row[0] for row in ROWS]
    assert np.allclose(melt, [row[1] for row in ROWS])
    # The uncertainty, never the area column that sits where it used to.
    assert np.allclose(sigma, [row[2] for row in ROWS])
    assert sigma.max() < 1e3


def test_a_table_without_the_named_columns_is_refused(tmp_path, monkeypatch):
    path = tmp_path / "wrong.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "melt", "error"])
        w.writerow([0, 39.85, 51.75])

    with pytest.raises(ValueError, match="bmr"):
        _load_obs_with(monkeypatch, path)


# Which table the default search opens, and where the result is written.
#
# IU Quartz staged the new table beside the old one, under
# parameterisations/ocean/meltobs/. A search with no candidate there ends on
# the old table, and the two integrate to 865.0 and 1067.4 Gt/yr, so the run
# has to say which one it opened. The default output is the file every forward
# and inversion in the checkout reads, so a calibration made as a check needs
# somewhere else to go.

NEW_NAME = "Melt_Paolo_Davison_Adusumilli_imbie2.csv"
OLD_NAME = "Melt_Paolo_Err_Adusumilli_imbie2_v3.csv"
BESIDE_OLD = os.path.join("parameterisations", "ocean", "meltobs")


@pytest.fixture
def calibrate(monkeypatch):
    r"""Import ``calibrate_melt`` afresh under a given environment.

    The knobs resolve at import, so each test sets its own and the module is
    dropped afterwards. Knobs a developer's shell may carry are cleared first.
    """
    def _import(**env):
        for key in ("ISMIP7_MELT_OBS_CSV", "ISMIP7_DATA_ROOT", "ISMIP7_K_OUT",
                    "ISMIP7_K_PER_BASIN_NPZ"):
            monkeypatch.delenv(key, raising=False)
        for key, value in env.items():
            monkeypatch.setenv(key, str(value))
        monkeypatch.syspath_prepend(_SCRIPTS)
        sys.modules.pop("calibrate_melt", None)
        import calibrate_melt
        return calibrate_melt

    yield _import
    sys.modules.pop("calibrate_melt", None)


def _stage(root, *relative):
    r"""Touch a table at ``root/relative`` and return its path."""
    path = os.path.join(str(root), *relative)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    _write(path, NEW_HEADER, old_layout=False)
    return path


def test_the_mirror_location_wins(tmp_path, calibrate, capsys):
    want = _stage(tmp_path, "meltobs", NEW_NAME)
    _stage(tmp_path, BESIDE_OLD, NEW_NAME)
    _stage(tmp_path, BESIDE_OLD, OLD_NAME)

    cm = calibrate(ISMIP7_DATA_ROOT=tmp_path)
    cm._announce_obs_table()
    said = capsys.readouterr().out

    assert cm.OBS_CSV == want
    assert want in said
    assert "[!]" not in said


def test_the_new_table_staged_beside_the_old_one_is_found(tmp_path, calibrate,
                                                         capsys):
    r"""The Quartz layout: both tables in one directory, no meltobs/ at the root."""
    want = _stage(tmp_path, BESIDE_OLD, NEW_NAME)
    _stage(tmp_path, BESIDE_OLD, OLD_NAME)

    cm = calibrate(ISMIP7_DATA_ROOT=tmp_path)
    cm._announce_obs_table()
    said = capsys.readouterr().out

    assert cm.OBS_CSV == want
    assert "[!]" not in said


def test_falling_through_to_the_old_table_is_announced(tmp_path, calibrate,
                                                       capsys):
    old = _stage(tmp_path, BESIDE_OLD, OLD_NAME)

    cm = calibrate(ISMIP7_DATA_ROOT=tmp_path)
    cm._announce_obs_table()
    said = capsys.readouterr().out

    assert cm.OBS_CSV == old
    assert "[!]" in said
    # The reader has to learn what was looked for, where, and the way out.
    assert NEW_NAME in said
    for searched in cm._OBS_CSV_CANDIDATES[:-1]:
        assert searched in said
    assert "ISMIP7_MELT_OBS_CSV" in said


def test_a_tree_with_neither_table_is_told_where_it_looked(tmp_path, calibrate,
                                                           capsys):
    r"""``main`` raises on the old table's path alone; this line precedes it."""
    cm = calibrate(ISMIP7_DATA_ROOT=tmp_path)
    cm._announce_obs_table()
    said = capsys.readouterr().out

    assert cm.OBS_CSV == cm._OBS_CSV_CANDIDATES[-1]
    assert "[!]" in said
    assert os.path.join(str(tmp_path), "meltobs", NEW_NAME) in said


def test_naming_the_old_table_is_a_choice_and_is_not_flagged(tmp_path,
                                                             calibrate, capsys):
    old = _stage(tmp_path, BESIDE_OLD, OLD_NAME)

    cm = calibrate(ISMIP7_DATA_ROOT=tmp_path, ISMIP7_MELT_OBS_CSV=old)
    cm._announce_obs_table()
    said = capsys.readouterr().out

    assert cm.OBS_CSV == old
    assert old in said
    assert "[!]" not in said


def test_the_default_output_is_the_mesh_s_file_under_results(calibrate):
    r"""Unset, ISMIP7_K_OUT leaves the calibration at the legacy path. No run
    searches for it: a forward reads a per-basin K only when
    ISMIP7_K_PER_BASIN_NPZ names it, and melts with the tracked calibration
    otherwise."""
    from icepack2_tools.runconfig import k_per_basin_npz

    cm = calibrate()
    results = os.path.join(_ROOT, "antarctica", "results")

    assert cm._k_out() == os.path.join(results,
                                       f"calibrated_K_per_basin_{cm.LC}.npz")
    assert k_per_basin_npz() is None


def test_the_output_can_be_named_away_from_the_default_path(tmp_path,
                                                            calibrate):
    want = str(tmp_path / "check" / "K_2000.npz")
    cm = calibrate(ISMIP7_K_OUT=want)

    assert cm._k_out() == want


def test_a_bare_output_name_resolves_under_results(calibrate):
    cm = calibrate(ISMIP7_K_OUT="K_check.npz")

    assert cm._k_out() == os.path.join(_ROOT, "antarctica", "results",
                                       "K_check.npz")


def test_the_printed_output_path_carries_the_suffix_savez_adds(tmp_path,
                                                               calibrate):
    cm = calibrate(ISMIP7_K_OUT=tmp_path / "K_check")

    assert cm._k_out() == str(tmp_path / "K_check.npz")
    np.savez(str(tmp_path / "K_check"), K_basin=np.zeros(2))
    assert os.path.exists(cm._k_out())
