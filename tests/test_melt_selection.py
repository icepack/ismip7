r"""The toolbox objective on this model's cells: our aggregation reproduces the
toolbox's own term functions on a regular grid, the objective runs unchanged
on it, and the files a selection writes load in a forward."""
import hashlib
import json
import os

import numpy as np
import pytest
from mpi4py import MPI

xr = pytest.importorskip("xarray")
pytest.importorskip("pandas")

from icepack2_tools import melt_selection as ms  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOLS = os.path.join(REPO, "icepack2_tools")
MODELS = [label for label, _ in ms.NOTEBOOK_MODELS]
N = 8                      # an 8 x 8 grid of 8 km points
RESO = 8000.0


def _grid():
    r"""Labels and fields on an 8 x 8 grid, laid out so every rule of the
    toolbox's term functions is exercised: 16 basins of 2 x 2 points, all 10
    buttressing bins with some points unbinned, a floating mask narrower
    than the melt mask, coverage gaps in the ocean states and zero melt in
    one region-year and one model-basin."""
    rng = np.random.default_rng(1)
    iy, ix = np.mgrid[0:N, 0:N]
    x = np.arange(N) * RESO
    y = np.arange(N) * RESO
    basins = ((iy // 2) * 4 + ix // 2).astype(float)
    bins = ((iy * N + ix) % 10).astype(float)
    bins[0, 3] = bins[5, 6] = np.nan
    melts = np.ones((N, N), bool)
    melts[7, :3] = False                      # outside floating_frac > 0.5
    fmask = melts.astype(float)
    fmask[6, 6] = fmask[1, 2] = 0.0           # melts, but outside the floating mask
    regions = np.zeros((N, N), int)
    regions[2:4, 4:6] = 1                     # basin 6: Pine Island
    regions[4:6, 6:8] = 2                     # basin 11: Dotson

    def field(lo, hi, gaps=0):
        f = rng.uniform(lo, hi, (N, N))
        if gaps:
            f[rng.integers(0, N, gaps), rng.integers(0, N, gaps)] = np.nan
        return np.where(melts, f, np.nan)

    present = (field(0.3, 3.0), field(34.0, 34.8))
    states = {}
    for kind in ("cold", "warm"):
        for label in MODELS:
            tf, so = field(0.1, 4.0, gaps=3), field(34.0, 34.8)
            if kind == "warm" and label == "mathiot":
                tf[basins == 3] = 0.0             # a zero basin mean
            states[(kind, label)] = (tf, so)
    for year in ms.NOTEBOOK_OBS_YEARS:
        tf, so = field(0.1, 4.0, gaps=2), field(34.0, 34.8)
        if year == 2007:
            tf[regions == 2] = 0.0                # a zero region total
        states[("obs", year)] = (tf, so)
    return dict(x=x, y=y, basins=basins, bins=bins, melts=melts, fmask=fmask,
                regions=regions, present=present, states=states)


def _targets(tmp_path):
    rng = np.random.default_rng(2)
    table = tmp_path / "melt.csv"
    M_obs = rng.uniform(5, 80, 16)
    sigma = rng.uniform(1, 10, 16)
    with open(table, "w") as f:
        f.write(",BMR (Gt/yr),BMR uncert (Gt/yr)\n")
        for b in range(16):
            f.write(f"{b},{float(M_obs[b])!r},{float(sigma[b])!r}\n")
    models7 = MODELS + ["haid"]
    basins = np.arange(16.0)

    def term3():
        m = rng.uniform(1, 30, (7, 16))
        m[4:6, :] = np.nan
        m[4:6, [9, 14]] = rng.uniform(1, 30, (2, 2))
        m[1, :] = m[6, :] = np.nan
        m[1, 14], m[6, 14] = 12.0, 14.0
        return xr.Dataset({"melt_rate": (("model", "basins"), m),
                           "melt_rate_uncert": (("model", "basins"), 0.2 * m)},
                          coords={"model": models7, "basins": basins})

    years = np.array(ms.NOTEBOOK_OBS_YEARS)
    t4 = rng.uniform(10, 90, (2, len(years)))
    t4[0, [0, 3]] = np.nan
    t4[1, [1, 5]] = np.nan
    return {
        "bids": np.arange(16), "M_obs": M_obs, "sigma": sigma, "table": str(table),
        "term2": xr.Dataset({"melt_mean": ("BFRN_bins", rng.uniform(10, 300, 10)),
                             "melt_mean_err": ("BFRN_bins", rng.uniform(5, 50, 10))},
                            coords={"BFRN_bins": np.arange(10)}),
        "cold": term3(), "warm": term3(),
        "term4": xr.Dataset({"melt_rate": (("region", "year"), t4),
                             "melt_rate_uncert": (("region", "year"), 0.2 * t4)},
                            coords={"region": ["pig", "dotson"], "year": years}),
        "t2_weight": np.concatenate([[0.0], rng.uniform(0.1, 25, 9)]),
    }


def _both_routes(tmp_path, K=(1e-5, 3e-5, 6e-5, 1e-4, 2e-4)):
    r"""The terms by the notebook's route on the grid, and by ours on the
    melting points treated as cells of 64 km2."""
    g, targets = _grid(), _targets(tmp_path)
    K = np.asarray(K)
    rho, sin_a = ms.NOTEBOOK_RHO_I, ms.NOTEBOOK_SIN_ALPHA
    coords = {"y": g["y"], "x": g["x"]}

    def da(a):
        return xr.DataArray(a, dims=("y", "x"), coords=coords)

    theirs = ms.notebook_terms(
        K, rho, sin_a, tuple(da(f) for f in g["present"]),
        {key: (da(t), da(s)) for key, (t, s) in g["states"].items()},
        da(g["fmask"]), da(g["basins"]).rename("basins"),
        da(g["bins"]).to_dataset(name="BFRN_bins"), g["regions"], targets,
        targets["table"])

    jy, jx = np.nonzero(g["melts"])
    cells = ms.Cells(np.full(len(jy), RESO ** 2), ms.as_labels(g["basins"][jy, jx]),
                     ms.as_labels(g["bins"][jy, jx]), g["regions"][jy, jx], sin_a,
                     agg=g["fmask"][jy, jx] > 0.5)
    present = tuple(f[jy, jx] for f in g["present"])
    states = [(kind, label, g["states"][(kind, label)][0][jy, jx],
               g["states"][(kind, label)][1][jy, jx])
              for kind, label, _, _ in ms.state_files("/nowhere")]
    rows = [ms.aggregate(k, np.zeros(len(jy)), present, states, cells, rho) for k in K]
    agg = {name: np.array([r[name] for r in rows]) for name in rows[0]}
    ours = ms.toolbox_terms(K, agg, targets)
    ours.update(ms.notebook_weights(ours, targets["t2_weight"], targets["term4"]))
    return ours, theirs


def test_our_cell_aggregation_reproduces_the_toolbox_term_functions(tmp_path):
    ours, theirs = _both_routes(tmp_path)
    for name, v in ms.compare_terms(ours, theirs).items():
        assert v["same_nan_pattern"], name
        assert v["max_rel_diff"] < 1e-12, (name, v)
    # the rules the comparison has to cover really fire here
    t3, t4 = ours["t3_model"], ours["t4_model"]
    assert np.isnan(t3.sel(model="mathiot", basins=3.0)).all()        # zero mean
    assert np.isnan(t3.sel(model="timmermann", basins=0.0)).all()     # restricted
    assert np.isfinite(t3.sel(model="timmermann", basins=14.0)).all()
    assert np.isnan(t4.sel(region="dotson", year=2007)).all()         # zero total
    assert np.isnan(t4.sel(region="pig", year=1994)).all()            # no target
    for name in ("t1_weights", "t2_weights", "t3_weights", "t4_weights"):
        assert ours[name].equals(theirs[name]), name


def test_the_objective_picks_the_same_K_from_either_route(tmp_path):
    ours, theirs = _both_routes(tmp_path)
    a = ms.select(ours, 400, 400, seed=3)
    b = ms.select(theirs, 400, 400, seed=3)
    assert np.array_equal(a, b)
    assert len(np.unique(a)) > 1


def _planted(K=np.arange(1, 6) * 1e-5, best=2):
    r"""Aggregates that match every target at K[best] and scale with K."""
    f = K / K[best]
    basins = np.arange(16.0)
    bins = np.arange(10)
    models = ["a", "b"]
    region, year = ["pig", "dotson"], [2009, 2012]
    p2 = [1.0]
    t = {
        "t1_model": xr.DataArray(10.0 * f[None, :, None] * np.ones((1, 1, 16)),
                                 dims=("p2", "p1", "basins"),
                                 coords={"p2": p2, "p1": K, "basins": basins}),
        "t1_obs_mean": xr.DataArray(np.full(16, 10.0), dims=["basin"],
                                    coords={"basin": np.arange(16)}),
        "t2_model": xr.DataArray(5.0 * f[None, :, None] * np.ones((1, 1, 10)),
                                 dims=("p2", "p1", "BFRN_bins"),
                                 coords={"p2": p2, "p1": K, "BFRN_bins": bins}),
        "t2_obs_mean": xr.DataArray(np.full(10, 5.0), dims=["BFRN_bins"],
                                    coords={"BFRN_bins": bins}),
        "t3_model": xr.DataArray(2.0 * f[None, None, :, None] * np.ones((2, 1, 1, 16)),
                                 dims=("model", "p2", "p1", "basins"),
                                 coords={"model": models, "p2": p2, "p1": K,
                                         "basins": basins}),
        "t3_obs_mean": xr.DataArray(np.full((2, 16), 2.0), dims=["model", "basins"],
                                    coords={"model": models, "basins": basins}),
        "t4_model": xr.DataArray(30.0 * f[None, None, :, None] * np.ones((2, 1, 1, 2)),
                                 dims=("year", "p2", "p1", "region"),
                                 coords={"year": year, "p2": p2, "p1": K,
                                         "region": region}),
        "t4_obs_mean": xr.DataArray(np.full((2, 2), 30.0), dims=["region", "year"],
                                    coords={"region": region, "year": year}),
    }
    for n in (1, 2, 3, 4):
        t[f"t{n}_obs_sigma"] = 0.01 * t[f"t{n}_obs_mean"]
    t["t1_weights"] = xr.DataArray(t["t1_model"].isel(p1=0, p2=0) * 0 + 1,
                                   dims=["basins"], coords={"basins": basins})
    t["t2_weights"] = xr.DataArray(np.ones(10), dims=["BFRN_bins"],
                                   coords={"BFRN_bins": bins})
    t["t3_weights"] = xr.DataArray(np.ones((2, 16)), dims=["model", "basins"],
                                   coords={"model": models, "basins": basins})
    t["t4_weights"] = xr.DataArray(np.ones((2, 2)), dims=["year", "region"],
                                   coords={"year": year, "region": region})
    return t, K[best]


def test_a_planted_optimum_is_found_in_every_sample_and_chunks_pool():
    from icepack2_tools import ismip7_parameter_selection_toolbox as pst
    terms, best = _planted()
    np.random.seed(7)
    direct, _ = pst.calculate_objective_function(
        600, None, *(terms[n] for n in ms.TERM_ARGS))
    # one chunk of the whole sample is the direct call, draw for draw
    assert np.array_equal(ms.select(terms, 600, 600, seed=7), direct)
    chunked = ms.select(terms, 600, 150, seed=7)
    assert len(chunked) == 600 and np.all(chunked == best)


def test_inputs_that_would_drop_K_silently_are_refused():
    terms, _ = _planted()
    shifted = dict(terms)
    p1 = terms["t2_model"].p1.values
    shifted["t2_model"] = terms["t2_model"].assign_coords(
        p1=np.nextafter(p1, np.inf))
    with pytest.raises(ValueError, match="other p1 labels"):
        ms.check_inputs(shifted)
    holed = dict(terms)
    t1 = terms["t1_model"].copy()
    t1[0, 1, 4] = np.nan
    holed["t1_model"] = t1
    with pytest.raises(ValueError, match="t1_model is not finite"):
        ms.check_inputs(holed)


def test_the_summary_reads_like_the_notebook():
    K = np.arange(1, 11) * 1e-5
    s = np.repeat(K, [0, 5, 10, 30, 30, 10, 5, 5, 3, 2])
    unrooted = np.zeros(10, bool)
    unrooted[-1] = True
    out = ms.summarise(s, K, unrooted)
    assert out["p50"] == np.percentile(s, 50) and out["p50_on_grid"]
    assert out["mode"] == 4e-5
    assert out["counts"][3] == 30
    assert out["share_on_grid_edges"] == pytest.approx(2 / 100)
    assert out["share_on_unrooted_K"] == pytest.approx(2 / 100)


def test_the_fit_recovers_known_offsets_without_firedrake():
    from icepack2_tools.forcing import quadratic_mixed_slope, _RHO_I
    rng = np.random.default_rng(0)
    n, K = 400, 8.5e-5
    tf = rng.uniform(0.5, 2.5, n)
    sal = np.full(n, 34.5)
    sin_a = np.full(n, 5.115e-3)
    area = np.full(n, 4.0e6)
    basin = np.repeat([3, 9], n // 2)
    bids, true = np.array([3, 9]), (0.3, -0.5)
    M_obs = np.array([float((quadratic_mixed_slope(tf[basin == b] + d, sal[:n // 2],
                                                   sin_a[:n // 2], K=K)
                             * area[:n // 2]).sum()) * _RHO_I / 1e12
                      for b, d in zip(bids, true)])
    dT, *_ = ms.fit_deltaT(tf, sal, sin_a, K, np.ones(n, bool), area, basin,
                           bids, M_obs, MPI.COMM_SELF)
    assert np.allclose(dT, true, atol=2e-4)


def test_the_plausibility_statistics_shift_with_the_offset():
    tf = np.array([0.5, 1.0, 4.0, 2.0, 6.0])
    area = np.array([1.0, 1.0, 2.0, 1.0, 1.0])
    basin = np.array([0, 0, 0, 1, -1])
    p = ms.Plausibility(tf, area, basin)
    d = np.zeros(16)
    d[0], d[1] = -0.75, 3.5
    out = p.at(d)
    assert out["tf_min"][0] == pytest.approx(-0.25)
    assert out["tf_max"][0] == pytest.approx(3.25)
    assert out["tf_mean"][0] == pytest.approx((0.5 + 1.0 + 8.0) / 4 - 0.75)
    assert out["frac_below_0"][0] == pytest.approx(0.25)
    assert out["frac_above_5"][1] == pytest.approx(1.0)
    assert np.isnan(out["tf_mean"][2])
    assert "frac_below_floor" not in out
    ruled = p.at(d, ms.TFRule(floor_area=(0.5, 0.25), cap_area=(5.25, 0.25)))
    assert ruled["frac_below_floor"][0] == pytest.approx(0.5)   # -0.25 and 0.25 of 4
    assert ruled["frac_above_cap"][1] == pytest.approx(1.0)     # 5.5
    assert ruled["frac_above_cap"][0] == 0.0
    assert np.isnan(ruled["frac_below_floor"][2])


def _per(n_K):
    r"""Per-K, per-basin statistics that pass the default rule everywhere;
    basin 15 has no floating cells."""
    per = {"unrooted": np.zeros((n_K, 16)), "tf_mean": np.full((n_K, 16), 1.0),
           "tf_min": np.full((n_K, 16), 0.2), "tf_max": np.full((n_K, 16), 4.0),
           "frac_below_floor": np.zeros((n_K, 16)),
           "frac_above_cap": np.zeros((n_K, 16))}
    for name in per:
        per[name][:, 15] = np.nan
    per["unrooted"][:, 15] = 1.0          # fit_deltaT flags a basin with no cells
    return per


def test_the_rule_admits_a_K_only_when_every_test_passes():
    per = _per(7)
    per["unrooted"][1, 9] = 1.0           # Amundsen outside the window
    per["tf_min"][2, 4] = -1.81           # one cell below the floor
    per["frac_below_floor"][3, 6] = 0.26  # over a quarter of basin 6 below -1 degC
    per["frac_below_floor"][4, 6] = 0.25  # a quarter exactly passes
    per["tf_max"][5, 9] = 6.81
    per["frac_above_cap"][6, 9] = 0.3
    ok, tests = ms.TFRule(cap=6.8, cap_area=(6.0, 0.25)).admits(per)
    assert ok.tolist() == [True, False, False, False, True, False, False]
    assert tests["rooted"].tolist() == [True, False] + [True] * 5
    assert tests["floor"].tolist() == [True, True, False] + [True] * 4
    assert tests["floor_area"].tolist() == [True] * 3 + [False] + [True] * 3
    assert tests["cap"].tolist() == [True] * 5 + [False, True]
    assert tests["cap_area"].tolist() == [True] * 6 + [False]
    # each test switches off on its own
    ok, tests = ms.TFRule(cap=None, cap_area=None, rooted=False).admits(per)
    assert set(tests) == {"floor", "floor_area"}
    assert ok.tolist() == [True, True, False, False, True, True, True]
    # the defaults bound no single cell on the warm side
    rule = ms.TFRule()
    assert (rule.floor, rule.floor_area, rule.cap, rule.cap_area, rule.rooted) == (
        -1.8, (-1.0, 0.25), None, (5.5, 0.25), True)
    ok, tests = rule.admits(per)
    assert "cap" not in tests
    assert ok.tolist() == [True, False, False, False, True, True, False]


def test_a_basin_beyond_the_toolbox_window_is_fitted_inside_three_K():
    from icepack2_tools.forcing import quadratic_mixed_slope, _RHO_I
    rng = np.random.default_rng(3)
    n, K = 200, 3.0e-5
    tf = rng.uniform(1.0, 3.5, n)
    sal, sin_a, area = np.full(n, 34.5), np.full(n, 5.115e-3), np.full(n, 4.0e6)
    basin, bids = np.full(n, 9), np.array([9])
    M_obs = np.array([float((quadratic_mixed_slope(tf + 2.6, sal, sin_a, K=K)
                             * area).sum()) * _RHO_I / 1e12])
    fit = lambda **kw: ms.fit_deltaT(tf, sal, sin_a, K, np.ones(n, bool), area,  # noqa: E731
                                     basin, bids, M_obs, MPI.COMM_SELF, **kw)
    dT, _, _, _, flagged = fit()
    assert ms.DT_WINDOW == (-3.0, 3.0) and flagged == []
    assert dT[0] == pytest.approx(2.6, abs=2e-4)
    dT, _, resid, _, flagged = fit(window=ms.TOOLBOX_DT_WINDOW)
    assert dT[0] == 2.0 and [b for b, _ in flagged] == [9] and resid[0] < 0


def test_the_objective_on_the_admitted_K_stays_on_them(tmp_path):
    g, targets = _grid(), _targets(tmp_path)
    K = np.array([1e-5, 3e-5, 6e-5, 1e-4, 2e-4])
    jy, jx = np.nonzero(g["melts"])
    cells = ms.Cells(np.full(len(jy), RESO ** 2), ms.as_labels(g["basins"][jy, jx]),
                     ms.as_labels(g["bins"][jy, jx]), g["regions"][jy, jx],
                     ms.NOTEBOOK_SIN_ALPHA, agg=g["fmask"][jy, jx] > 0.5)
    present = tuple(f[jy, jx] for f in g["present"])
    states = [(kind, label, g["states"][(kind, label)][0][jy, jx],
               g["states"][(kind, label)][1][jy, jx])
              for kind, label, _, _ in ms.state_files("/nowhere")]
    rows = [ms.aggregate(k, np.zeros(len(jy)), present, states, cells,
                         ms.NOTEBOOK_RHO_I) for k in K]
    agg = {name: np.array([r[name] for r in rows]) for name in rows[0]}
    ok = np.array([False, True, True, True, False])
    full = ms.toolbox_terms(K, agg, targets)
    sub = ms.toolbox_terms(K[ok], {n: v[ok] for n, v in agg.items()}, targets)
    for name in ("t1_model", "t2_model", "t3_model", "t4_model"):
        assert sub[name].equals(full[name].sel(p1=K[ok])), name
    sub.update(ms.notebook_weights(sub, targets["t2_weight"], targets["term4"]))
    picks = ms.select(sub, 300, 100, seed=5)
    assert np.isin(picks, K[ok]).all()


def test_the_notebook_slope_recipe_on_a_plane():
    a, b = 3e-3, -4e-3
    x = np.arange(12) * 8000.0
    y = np.arange(10) * 8000.0
    draft = -500.0 + a * x[None, :] + b * y[:, None]
    draft[4, 5] = np.nan                     # one-sided differences around it
    floating = np.ones(draft.shape, bool)
    floating[4, 5] = False
    got = ms.notebook_mean_slope(draft, floating, 8000.0, 8000.0)
    assert got == pytest.approx(np.sin(np.arctan(np.hypot(a, b))), rel=1e-12)


def _nc(path, data, x, y, z=None):
    dims = ("z", "y", "x") if z is not None else ("y", "x")
    coords = {"x": x, "y": y, **({"z": z} if z is not None else {})}
    xr.Dataset({"v": (dims, data)}, coords=coords).to_netcdf(path)
    return str(path)


def test_the_sampler_keeps_gaps_as_nan_or_fills_them(tmp_path):
    x = np.array([0.0, 10.0, 20.0])
    y = np.array([20.0, 10.0, 0.0])              # descending, as some files store it
    z = np.array([-30.0, -90.0])
    data = np.arange(18, dtype=float).reshape(2, 3, 3)
    data[1, 0, 2] = np.nan                       # z -90, y 20, x 20
    path = _nc(tmp_path / "f.nc", data, x, y, z)
    px, py, pd_ = np.array([20.0, 1.0, 99.0]), np.array([19.0, 1.0, 0.0]), np.array([-100.0, -40.0, -30.0])
    kept = ms.sample_nearest(path, "v", px, py, pd_)
    filled = ms.sample_nearest(path, "v", px, py, pd_, fill=0.0)
    assert np.isnan(kept[0]) and filled[0] == 0.0
    assert kept[1] == filled[1] == data[0, 2, 0]     # y 0 is the last stored row
    assert np.isnan(kept[2]) and filled[2] == 0.0    # off the grid
    labels = ms.as_labels(ms.sample_nearest(_nc(tmp_path / "l.nc", np.array(
        [[1.0, np.nan], [3.0, 4.0]]), x[:2], y[1:]), "v",
        np.array([0.0, 10.0]), np.array([10.0, 10.0])))
    assert labels.tolist() == [1, -1]


def test_the_grid_interp_the_calibrations_use_is_the_zero_fill_sampler(tmp_path, monkeypatch):
    pytest.importorskip("firedrake")
    pytest.importorskip("rasterio")
    pytest.importorskip("icepack")
    import sys
    monkeypatch.syspath_prepend(os.path.join(REPO, "antarctica", "scripts"))
    sys.modules.pop("calibrate_melt", None)
    try:
        import calibrate_melt
        x, y, z = np.array([0.0, 10.0]), np.array([0.0, 10.0]), np.array([-30.0, -90.0])
        data = np.array([[[1.0, np.nan], [3.0, 4.0]], [[5.0, 6.0], [np.nan, 8.0]]])
        path = _nc(tmp_path / "f.nc", data, x, y, z)
        px, py, pz = np.array([9.0, 1.0, 50.0]), np.array([1.0, 9.0, 0.0]), np.array([-20.0, -95.0, -30.0])
        assert np.array_equal(calibrate_melt._grid_interp(path, "v", px, py, draft=pz),
                              ms.sample_nearest(path, "v", px, py, pz, fill=0.0))
    finally:
        sys.modules.pop("calibrate_melt", None)


@pytest.fixture
def clean(monkeypatch):
    import icepack2_tools.forcing as forcing
    for k in ("ISMIP7_MELT_SLOPE", "ISMIP7_SIN_ALPHA_ANT", "ISMIP7_GEOMETRY_SPACE",
              "ISMIP7_K_SCALE", "ISMIP7_DELTAT_PER_BASIN_NPZ"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setattr(forcing, "_MELT_SLOPE_WARNED", False)
    monkeypatch.setattr(forcing, "_SLOPE_CAP_WARNED", False)
    monkeypatch.setattr(forcing, "_GEOMETRY_SPACE_WARNED", False)


def test_a_selected_offsets_file_loads_in_a_forward(tmp_path, clean, capsys):
    from icepack2_tools.forcing import load_deltaT_per_basin, sin_alpha_ant
    from icepack2_tools.runconfig import geometry_space
    x = np.array([0.0, 1.0, 2.0, 3.0])
    bn = np.where(np.arange(4)[None, :] < 2, 3, 9).repeat(4, axis=0).astype(float)
    imbie = _nc(tmp_path / "basins.nc", bn, x, x)
    ds = xr.open_dataset(imbie).rename({"v": "basinNumber"})
    ds.load().to_netcdf(tmp_path / "basinNumber.nc")
    imbie = str(tmp_path / "basinNumber.nc")
    bids, dT = np.array([3, 9]), np.array([0.25, -1.1])
    path = str(tmp_path / "deltaT_per_basin_2000_K4.750e-05.npz")
    ms.write_offsets(path, bids, dT, 4.75e-5, np.array([10.0, 20.0]),
                     np.array([8.0, 25.0]), np.zeros(2), np.array([3.0, 4.0]),
                     melt_slope="ant", sin_alpha_ant=sin_alpha_ant(),
                     sin_alpha_cap=float("inf"), geometry_space=geometry_space(),
                     imbie2_nc=imbie, selected_as="K05")
    mx, my = np.array([0.4, 2.6, 40.0]), np.array([0.1, 2.9, 40.0])
    field, K = load_deltaT_per_basin(path, mx, my)
    assert K == 4.75e-5
    basin = np.array([3, 9, -1])
    assert np.array_equal(field, ms.offsets_on_cells(dT, bids, basin))
    assert "WARNING" not in capsys.readouterr().out


def test_the_vendored_toolbox_is_the_recorded_upstream_file():
    with open(os.path.join(TOOLS, "ismip7_parameter_selection_toolbox.source.json")) as f:
        src = json.load(f)
    with open(os.path.join(TOOLS, "ismip7_parameter_selection_toolbox.py"), "rb") as f:
        body = f.read()
    assert hashlib.sha256(body).hexdigest() == src["sha256"]
    assert src["repository"] == "ismip/ismip7-antarctic-ocean-forcing"
    assert len(src["commit"]) == 40 and len(src["file_commit"]) == 40
    with open(os.path.join(TOOLS, "ismip7_parameter_selection_toolbox.LICENSE")) as f:
        assert f.readline().strip() == "MIT License"


def test_an_extended_K_grid_keeps_the_notebook_values():
    base = ms.notebook_K_grid()
    assert len(base) == 120 and base[-1] == pytest.approx(3.0e-4)
    assert np.array_equal(ms.notebook_K_grid(3.0e-4), base)
    wide = ms.notebook_K_grid(1.0e-3)
    assert len(wide) == 400 and wide[-1] == pytest.approx(1.0e-3)
    assert np.array_equal(wide[:120], base)
    assert np.allclose(np.diff(wide), 2.5e-6)
    with pytest.raises(ValueError, match="inside the notebook's grid"):
        ms.notebook_K_grid(2.0e-4)
