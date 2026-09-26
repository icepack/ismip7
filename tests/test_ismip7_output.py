r"""The ISMIP7 writer's regridding policies and time encoding, on tiny
synthetic operators (no mesh, no data).

The time stamps are the ones the discussion board settled on (#16, #20):
a state written for year 2015 is stamped 2016-01-01, day 60630 since
1850-01-01 on the standard calendar; a flux for 2015 is stamped 2015-07-01,
day 60446. The fill policies follow the request's csv: ``forbidden``
means the uncovered part of a pixel counts as zero (sums are conserved),
``outside_domain`` means the mean over the covered part, ``no_ice`` and
friends mean over the masked part. A flux is a whole-pixel mean under any
policy, which then decides only where it is fill.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "antarctica", "scripts"))
scipy_sparse = pytest.importorskip("scipy.sparse")
wio = pytest.importorskip("write_ismip7_output")


def _operator():
    # two pixels, three cells: pixel 0 fully covered by cells 0 and 1 (half
    # each), pixel 1 one quarter covered by cell 2
    A = wio.PIXEL_AREA
    return scipy_sparse.csr_matrix(np.array([[A / 2, A / 2, 0.0], [0.0, 0.0, A / 4]]))


def test_forbidden_policy_conserves_the_sum():
    W = _operator(); v = np.array([100.0, 300.0, 400.0])
    out = wio.regrid(W, v, "forbidden")
    assert np.allclose(out, [200.0, 100.0])                   # pixel 1: a quarter of 400
    assert np.isclose((out * wio.PIXEL_AREA).sum(), (W @ v).sum())


def test_outside_domain_policy_is_the_covered_mean():
    W = _operator(); v = np.array([100.0, 300.0, 400.0])
    out = wio.regrid(W, v, "outside_domain")
    assert np.allclose(out, [200.0, 400.0])                   # partial coverage does not dilute an elevation


def test_masked_policy_averages_over_the_mask():
    W = _operator(); v = np.array([10.0, 30.0, 50.0]); mask = np.array([True, False, True])
    out = wio.regrid(W, v, "no_ice", mask)
    assert np.allclose(out, [10.0, 50.0])
    assert np.isnan(wio.regrid(W, v, "no_ice", np.zeros(3, bool))).all()


def test_a_flux_is_a_whole_pixel_mean_and_the_policy_only_fills():
    r"""ismip7-scalars sums a flux times the pixel area (forum thread 50,
    issue #96), so a flux mean is taken over the whole pixel and sums to the
    model's integral. The policy still fills where its own mean would have no
    area: outside the domain for acabf, where no ice floats for libmassbffl.
    A cell that no longer floats keeps its melt in a pixel that still does."""
    A = wio.PIXEL_AREA
    # pixel 0 fully covered by cells 0 and 1, pixel 1 a quarter by cell 2,
    # pixel 2 not at all
    W = scipy_sparse.csr_matrix(np.array([[A / 2, A / 2, 0.0], [0.0, 0.0, A / 4],
                                          [0.0, 0.0, 0.0]]))
    v = np.array([100.0, 300.0, 400.0])
    out = wio.regrid(W, v, "outside_domain", whole_pixel=True)
    assert np.allclose(out[:2], [200.0, 100.0]) and np.isnan(out[2])
    assert np.isclose(np.nansum(out) * A, (W @ v).sum())
    melt = np.array([-10.0, -30.0, -50.0])
    afloat = np.array([True, False, False])               # cell 1 grounded during the year
    out = wio.regrid(W, melt, "no_floating_ice", afloat, whole_pixel=True)
    assert np.isclose(out[0], -20.0)                      # both cells' melt, over the pixel
    assert np.isnan(out[1]) and np.isnan(out[2])          # no floating ice left: fill
    assert np.isnan(wio.regrid(W, melt, "no_floating_ice", afloat)[1])


def test_the_request_types_choose_the_mean():
    r"""acabf and orog share ``outside_domain``. The flux is a whole-pixel
    mean; the elevation stays a covered-part mean, which the writer's
    ``base := orog - lithk`` rests on."""
    req = wio.request_table()
    W = _operator(); v = np.array([100.0, 300.0, 400.0])
    assert np.allclose(wio.pixel_values(W, v, req["orog"], {}), [200.0, 400.0])
    assert np.allclose(wio.pixel_values(W, v, req["acabf"], {}),
                       np.array([200.0, 100.0]) * wio.CONVERT["kg m-2 s-1"])


def test_the_flux_files_say_how_their_pixel_means_were_taken(tmp_path):
    import netCDF4
    req = wio.request_table()
    for var, flux in (("acabf", True), ("orog", False)):
        ds, _ = wio.create_2d(str(tmp_path / f"{var}.nc"), var, req[var], [2015], flux)
        ds.close()
    with netCDF4.Dataset(tmp_path / "acabf.nc") as ds:
        assert ds.getncattr(wio.FLUX_MEAN_ATTR) == wio.FLUX_MEAN
    with netCDF4.Dataset(tmp_path / "orog.nc") as ds:
        assert wio.FLUX_MEAN_ATTR not in ds.ncattrs()


def test_time_encoding_matches_the_board():
    assert wio.days_since_1850(2016, 1, 1) == 60630           # state written for year 2015
    assert wio.days_since_1850(2015, 7, 1) == 60446           # flux for year 2015


# --- the names a submission is filed under ----------------------------------

def test_the_core_counter_follows_from_the_forcing():
    r"""It names the experiment in every filename and in the directory, it
    was typed by hand, and C007 on a run that is C008 is a well-formed
    submission of the wrong experiment, which no checker can catch."""
    assert wio.set_counter("CESM2-WACCM", "ssp585", "CORE") == "C007"
    assert wio.set_counter("MRI-ESM2-0", "ssp585", "CORE", "C008") == "C008"
    assert wio.set_counter("MRI-ESM2-0", "ctrl", "CORE") == "C010"
    assert wio.set_counter("CESM2-WACCM", "ocx", "CORE") == "C011"
    with pytest.raises(ValueError, match="MRI-ESM2-0 ssp585 is C008"):
        wio.set_counter("MRI-ESM2-0", "ssp585", "CORE", "C007")
    with pytest.raises(ValueError, match="not a core experiment"):
        wio.set_counter("CESM2-WACCM", "ssp534-over", "CORE")
    # the ten ESM-forced cores, numbered as the conventions document has them
    assert sorted(wio.CORE_COUNTER.values()) == [f"C{n:03d}" for n in range(1, 11)]


def test_outside_the_core_set_the_counter_has_to_be_given():
    assert wio.set_counter("IPSL-CM6A-LR", "ssp370", "ESM", "E003") == "E003"
    with pytest.raises(ValueError, match="required for the PPE set"):
        wio.set_counter("CESM2-WACCM", "ssp585", "PPE")
    with pytest.raises(ValueError, match="must look like"):
        wio.set_counter("CESM2-WACCM", "ssp585", "CORE", "c7")


@pytest.mark.parametrize("bad", ["icepack2_dg0", "icepack2.1", "ice pack", ""])
def test_an_id_that_would_shift_the_filename_fields_is_refused(bad):
    r"""Discussion #17: ``ISSMv2026p2``, not ``ISSM_2026.2``."""
    assert wio.submission_id("ism-id", "icepack2") == "icepack2"
    assert wio.submission_id("source-id", "RICE-IU") == "RICE-IU"
    with pytest.raises(ValueError, match="no underscore, dot or space"):
        wio.submission_id("ism-id", bad)


def test_the_grid_is_the_official_8_km_description():
    r"""``isschecker/data/gdfs/gdf_ISMIP7_AIS_08000m.txt`` at 0.5.0: xsize and
    ysize 761, xfirst and yfirst -3040000, xinc and yinc 8000. The constants
    were a hand transcription nothing held to it."""
    from icepack2_tools import regrid
    x, y = regrid.ismip7_grid_coords()
    assert (regrid.ISMIP7_NX, regrid.ISMIP7_NY) == (len(x), len(y)) == (761, 761)
    assert x[0] == y[0] == -3040000.0
    assert np.allclose(np.diff(x), 8000.0) and np.allclose(np.diff(y), 8000.0)
    assert regrid.ISMIP7_DX == 8000.0 and x[-1] == y[-1] == 3040000.0


# af2_AIS_08000m_v1.nc, the area factor ismip7-scalars weights every pixel by
# (/ISMIP7/Output-Processing/Data/AIS on Globus, float32, "after Snyder
# (1987)", sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f),
# at the centres of eleven pixels (row j, column i), read on IU Quartz on
# 24 September 2026: the pole, its neighbours, 71 S, the grid's corners.
ORGANISERS_AF2 = {
    (380, 380): 1.056770205, (380, 381): 1.056769252, (381, 381): 1.056768417,
    (380, 500): 1.044315219, (100, 380): 0.991517127, (0, 0): 0.843348742,
    (760, 760): 0.843348742, (200, 600): 0.989615023, (380, 0): 0.941123068,
    (600, 150): 0.973586619, (450, 300): 1.046977997,
}


def test_the_area_factor_is_the_organisers_grid():
    r"""Issue #97: the native scalars weight each cell by the factor the
    organisers' tool weights each pixel by, to the file's float32 rounding."""
    from icepack2_tools.regrid import area_factor, ismip7_grid_coords
    x, y = ismip7_grid_coords()
    j, i = (np.array(ax) for ax in zip(*ORGANISERS_AF2))
    want = np.array(list(ORGANISERS_AF2.values()))
    assert np.allclose(area_factor(x[i], y[j]), want, rtol=1.2e-7, atol=0.0)


def test_the_area_factor_is_snyder_s_scale_factor():
    r"""(1/k)^2 of the polar stereographic projection on WGS84, true at 71 S
    (Snyder 1987, eqs. 21-33 to 21-35), written out apart from pyproj at
    points placed by latitude, the pole included."""
    from icepack2_tools.regrid import area_factor
    a, f = 6378137.0, 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    e = np.sqrt(e2)

    def t(p):
        return np.tan(np.pi / 4 - p / 2) / ((1 - e * np.sin(p)) / (1 + e * np.sin(p))) ** (e / 2)

    def m(p):
        return np.cos(p) / np.sqrt(1 - e2 * np.sin(p) ** 2)

    pc = np.radians(71.0)
    phi = np.radians([85.0, 71.0, 65.0, 55.0])        # degrees south
    rho = a * m(pc) * t(phi) / t(pc)
    k = rho / (a * m(phi))
    assert np.isclose(k[1], 1.0, rtol=1e-12, atol=0.0)
    assert np.allclose(area_factor(rho / np.sqrt(2), -rho / np.sqrt(2)), k ** -2, rtol=1e-9, atol=0.0)
    k_pole = m(pc) / (2 * t(pc)) * np.sqrt((1 + e) ** (1 + e) * (1 - e) ** (1 - e))
    assert np.isclose(area_factor(np.zeros(1), np.zeros(1))[0], k_pole ** -2, rtol=1e-9, atol=0.0)


def test_the_area_factor_matches_the_organisers_whole_grid():
    r"""Every pixel centre, where the file is at hand: ``ISMIP7_AF2_GRID``,
    or where ``download_forcing.py --scalar-processing`` puts it."""
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.environ.get("ISMIP7_AF2_GRID") or os.path.join(
        repo, "ISMIP7", "Output-Processing", "Data", "AIS", "af2_AIS_08000m_v1.nc")
    if not os.path.exists(path):
        pytest.skip(f"no {path}")
    import netCDF4
    from icepack2_tools.regrid import area_factor
    with netCDF4.Dataset(path) as ds:
        grid = np.asarray(ds.variables["af2"][:, :], dtype=np.float64)
        x = np.asarray(ds.variables["x"][:], dtype=np.float64)
        y = np.asarray(ds.variables["y"][:], dtype=np.float64)
    X, Y = np.meshgrid(x, y)
    got = area_factor(X.ravel(), Y.ravel()).reshape(grid.shape)
    assert np.max(np.abs(got / grid - 1.0)) <= 1.2e-7


def test_the_writer_tells_which_area_the_scalars_csv_integrates_over():
    r"""A forward with the area factor writes true-area scalars and one from
    before map-plane ones; the annual file's own areas under each convention
    tell which a CSV row holds."""
    sums = {wio.TRUE_AREA: (1.2357e13, 1.4595e12), wio.MAP_PLANE: (1.2077e13, 1.4261e12)}
    row = {"year": "2015", "iareagr": "1.235700e+13", "iareafl": "1.459500e+12"}
    assert wio.scalar_area(row, sums) == wio.TRUE_AREA
    row = dict(row, iareagr="1.207700e+13", iareafl="1.426100e+12")
    assert wio.scalar_area(row, sums) == wio.MAP_PLANE
    with pytest.raises(ValueError, match="match none"):
        wio.scalar_area(dict(row, iareagr="1.300000e+13"), sums)


def test_a_series_that_mixes_the_two_areas_is_refused():
    r"""A chained run whose links straddled the change to true-area scalars
    would submit both in one series."""
    assert wio.series_area({2015: wio.TRUE_AREA, 2016: wio.TRUE_AREA}) == wio.TRUE_AREA
    with pytest.raises(ValueError, match="map_plane in 2 years, 2015 to 2016; true_area"):
        wio.series_area({2015: wio.MAP_PLANE, 2016: wio.MAP_PLANE, 2017: wio.TRUE_AREA})


def test_the_bundled_request_says_where_it_came_from():
    r"""Upstream edits the table as the forum finds problems (#16, #22, #23,
    #46), so the copy records its tag, and the writer's columns are there."""
    import csv
    import json
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "icepack2_tools")
    with open(os.path.join(here, "ismip7_variable_request.source.json")) as f:
        src = json.load(f)
    assert src["repository"] == "ismip/ISM_SimulationChecker" and src["tag"] and len(src["commit"]) == 40
    with open(os.path.join(here, "ismip7_variable_request.csv"), newline="") as f:
        rows = {r["Variable Name"]: r for r in csv.DictReader(f)}
    for column in ("Type", "units", "standard_name", "fill_policy", "range_severity"):
        assert column in rows["lithk"]
    assert rows["topg"]["max_value_ais"] == "5500"            # #23
    assert float(rows["licalvf"]["max_value_ais"]) == 0.0     # #16: calving is a loss, so negative
    assert float(rows["ligroundf"]["min_value_ais"]) < 0.0    # #22: negative at ice rumples is allowed


def test_the_request_audit_reports_a_moved_cell():
    import audit_variable_request as avr
    ours = avr.rows("Variable Name,units,range_severity\nhfgeoubed,W m-2,error\nlithk,m,error\n")
    theirs = avr.rows("Variable Name,units,range_severity\nhfgeoubed,W m-2,warning\nnewvar,1,error\n")
    assert avr.differences(ours, theirs) == [
        ("hfgeoubed", "range_severity", "error", "warning"),
        ("lithk", None, "present", "absent"),
        ("newvar", None, "absent", "present"),
    ]
    assert "units" in avr.WRITER_COLUMNS and "range_severity" not in avr.WRITER_COLUMNS


def test_a_cell_at_flotation_to_rounding_is_written_grounded():
    r"""isschecker reads a wholly floating pixel less than 1 cm above
    ``topg`` as ice resting on the bed. A floating cell 5 mm above its bed
    goes into the grounded mask, one 2 cm above stays floating, and a grounded
    cell is left alone; the geometry does not move."""
    cells = {
        "lithk": np.array([100.0, 100.0, 300.0]),
        "orog": np.array([10.45, 10.02, 100.0]),       # bases -89.55, -89.98, -200
        "topg": np.array([-89.555, -90.0, -200.0]),
        "sftflf": np.array([1.0, 1.0, 0.0]),
        "sftgrf": np.array([0.0, 0.0, 1.0]),
    }
    before = {k: v.copy() for k, v in cells.items()}
    assert wio.ground_near_flotation(cells) == 1
    assert np.array_equal(cells["sftflf"], [0.0, 1.0, 0.0])
    assert np.array_equal(cells["sftgrf"], [1.0, 0.0, 1.0])
    for k in ("lithk", "orog", "topg"):
        assert np.array_equal(cells[k], before[k])


def test_front_melt_is_never_filled_where_the_basal_melt_would_be():
    r"""lifmassbf carries the melt of ice that flowed into cells holding no
    ice at either end of the year (issue #109) under the request's
    ``forbidden`` policy: a whole-pixel mean in every pixel, summing to the
    model's integral, where the same melt booked as libmassbffl would be fill
    for want of floating ice."""
    req = wio.request_table()
    assert req["lifmassbf"]["fill_policy"] == "forbidden"
    W = _operator(); v = np.array([0.0, -30.0, -40.0])          # m/yr of ice
    masks = {"no_floating_ice": np.zeros(3, dtype=bool)}           # nothing floats
    k = wio.CONVERT["kg m-2 s-1"]
    out = wio.pixel_values(W, v, req["lifmassbf"], masks)
    assert np.isfinite(out).all()
    assert np.allclose(out, np.array([-15.0, -10.0]) * k)
    assert np.isclose((out * wio.PIXEL_AREA).sum(), (W @ v).sum() * k)
    assert np.isnan(wio.pixel_values(W, v, req["libmassbffl"], masks)).all()


def test_the_melt_summary_reports_what_leaves_and_what_front_melt_carries():
    r"""Per year: the melt the fill leaves out of libmassbffl, its part in
    pixels only the near-flotation rule emptied, and lifmassbf's grid sum;
    then each at its largest and at the marker years and the series' last."""
    W = _operator()
    cells = {"libmassbffl": np.array([-10.0, -20.0, -40.0]),
             "lifmassbf": np.array([0.0, 0.0, -8.0])}
    # pixel 1 floated before the near-flotation rule and not after it
    got = wio.melt_booking(W, cells, np.array([1.0, 1.0]), np.array([1.0, 0.0]))
    gt = wio.RHO_I / 1e12 * wio.PIXEL_AREA / 4                     # cell 2's quarter of pixel 1
    assert got["left out"] == pytest.approx(-40.0 * gt)
    assert got["near flotation"] == pytest.approx(-40.0 * gt)
    assert got["front melt"] == pytest.approx(-8.0 * gt)

    per_year = {y: {"left out": -1.0, "near flotation": 0.0, "front melt": -2.0}
                for y in (2099, 2100, 2101)}
    per_year[2101] = {"left out": -3.0, "near flotation": -0.5, "front melt": -1.0}
    lines = wio.melt_summary(per_year)
    assert len(lines) == 3
    assert "libmassbffl leaves out" in lines[0] and "largest -3.0 Gt/yr (2101)" in lines[0]
    assert "-1.0 in 2100, -3.0 in 2101" in lines[0]
    assert "lifmassbf carries the front melt" in lines[2]
    assert "largest -2.0 Gt/yr (2099)" in lines[2] and "-2.0 in 2100, -1.0 in 2101" in lines[2]


def test_a_series_that_mixes_two_melt_bookings_is_refused():
    r"""A chained run whose links straddled the front-melt booking would
    submit lifmassbf as zero in some years and as the melt of the same cells
    in others."""
    from icepack2_tools.ismip7_output import FRONT_MELT
    assert wio.series_front_melt({2015: FRONT_MELT, 2016: FRONT_MELT}) == FRONT_MELT
    old = {2015: wio.UNSTAMPED_MELT, 2016: wio.UNSTAMPED_MELT}
    assert wio.series_front_melt(old) == wio.UNSTAMPED_MELT
    with pytest.raises(ValueError, match="mix melt bookings.*2 years, 2015 to 2016"):
        wio.series_front_melt({**old, 2017: FRONT_MELT})
