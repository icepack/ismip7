r"""Melt-parameter selection on this model's own cells (ISMIP7 ocean-forcing
protocol, section 2.3).

The ISMIP7 toolbox (``ismip7_parameter_selection_toolbox.py``, vendored
unchanged from ismip/ismip7-antarctic-ocean-forcing) draws K05, K50 and K95
from an objective over four aggregates of the parameterised melt at each K:

* term 1, the present-day total per IMBIE2 basin (Gt/yr);
* term 2, the present-day total per buttressing bin (Gt/yr);
* term 3, the warm minus cold basin mean under ocean-model states (kg/m2/yr);
* term 4, the Pine Island and Dotson totals under observed ocean states
  (Gt/yr).

Its ``calculate_term1..4`` build those on a regular grid alone. This module
builds them on any set of cells with areas, the forward's DG0 cells
included, in the layout ``calculate_objective_function`` reads, and that
function runs unchanged. ``fit_deltaT`` is the per-basin thermal-forcing
offset; fitted at every K before the objective runs, the order section 2.3
prefers, each K's offsets enter every term. ``TFRule`` says which K leave a
plausible present-day thermal forcing (section 2.1) and a fitted total in
every basin; the objective can then run on those K alone.

The defaults are the toolbox's quadratic example notebook
(``parameter_selection_quadratic_example.ipynb`` at 5f4d3df, cells 6 to 40):
the K grid, the six ocean-model states and their basin restrictions, the 13
observation years, the weights and 100000 samples.
"""
import csv
import os

import numpy as np

from .forcing import quadratic_mixed_slope, _RHO_I
from .mpi_stats import global_count, global_sum

N_BASINS = 16
N_BFRN_BINS = 10
NOTEBOOK_SAMPLE_SIZE = 100000
# The notebook converts m/yr of ice to kg/m2/yr with this density (cell 4).
NOTEBOOK_RHO_I = 918.0
# The constant slope behind the notebook's printed outputs: its gamma_T
# printout (cell 57) gives sin(alpha) = 5.1117e-3 within 7e-8 at each of the
# three percentiles, with the notebook's own constants.
NOTEBOOK_SIN_ALPHA = 5.1117e-3
# (label, file prefix) of the ocean-model states, in the notebook's order
# (cell 12); the prefix takes "cold" or "warm".
NOTEBOOK_MODELS = (
    ("mathiot", "Mathiot_NEMO_{state}_v3_"),
    ("timmermann", "Timmermann_FESOM_{state}_v3_"),
    ("naughten_ais_1", "Naughten_FESOM_ACCESS_{state}_v2_"),
    ("naughten_ais_2", "Naughten_FESOM_MMM_{state}_v2_"),
    ("jourdain_naughten", "Jourdain-Naughten_NEMO-MITgcm_{state}_"),
    ("naughten_naughten", "Naughten_MITamu-MITwed_{state}_"),
)
# Regional models count only in these basins (cell 13).
NOTEBOOK_MODEL_BASINS = {
    "timmermann": (14,),
    "jourdain_naughten": (9, 14),
    "naughten_naughten": (9, 14),
}
# Term 3 weighs these models one and the rest zero (cell 35).
NOTEBOOK_T3_MODELS = ("mathiot", "naughten_ais_1", "jourdain_naughten",
                      "naughten_naughten")
NOTEBOOK_OBS_YEARS = (1994, 2000, 2006, 2007, 2009, 2010, 2011, 2012, 2014,
                      2016, 2018, 2019, 2020)
# Term 4 weighs Pine Island in these years one and the rest zero (cell 36).
NOTEBOOK_T4_REGION = "pig"
NOTEBOOK_T4_YEARS = (2009, 2012)
# Shelf ids in shelf_mask_ismip8km.nc; Pine Island is its main trunk, east
# of x = -1.625e6 (cell 20).
PIG_ID, DOTSON_ID, PIG_X_MIN = 110, 97, -1.625e6
REGIONS = ("pig", "dotson")
# The deltaT search window, in K. The toolbox searches plus or minus 2 K
# (optimise_deltaT); the protocol sets no window and bounds the offsets
# through the thermal forcing they leave (TFRule), so this repository
# searches plus or minus 3 K.
TOOLBOX_DT_WINDOW = (-2.0, 2.0)
DT_WINDOW = (-3.0, 3.0)
# The order calculate_objective_function takes its arguments in.
TERM_ARGS = ("t1_model", "t1_obs_mean", "t1_obs_sigma", "t1_weights",
             "t2_model", "t2_obs_mean", "t2_obs_sigma", "t2_weights",
             "t3_model", "t3_obs_mean", "t3_obs_sigma", "t3_weights",
             "t4_model", "t4_obs_mean", "t4_obs_sigma", "t4_weights")


NOTEBOOK_K_MAX = 3.0e-4


def notebook_K_grid(k_max=None):
    r"""The notebook's K grid (cell 6): 2.5e-6 to 3.0e-4 in steps of 2.5e-6.

    ``k_max`` extends it upward on the same step, and its first 120 values
    stay the notebook's bit for bit. An extension changes more than the tail:
    the objective divides each term by its median over the whole grid, so
    the grid's extent sets the terms' relative weight.

    ``np.arange`` rounds differently across numpy builds, so a run makes the
    grid once and stores it, and every later step reads the stored values:
    the objective aligns its terms on these labels, and a K that differs in
    the last bit drops out of that inner join without a message."""
    grid = np.arange(0.25e-5, 3.025e-4, 0.25e-5)
    if k_max is None or np.isclose(k_max, NOTEBOOK_K_MAX, rtol=1e-9, atol=0):
        return grid
    if k_max < NOTEBOOK_K_MAX:
        raise ValueError(f"k_max {k_max:g} lies inside the notebook's grid, "
                         f"which ends at {NOTEBOOK_K_MAX:g}")
    wide = np.arange(0.25e-5, k_max + 0.25e-5, 0.25e-5)
    wide = wide[wide <= k_max * (1 + 1e-9)]
    if not np.array_equal(wide[:len(grid)], grid):
        raise ValueError("this numpy extends the grid with other values "
                         "than the notebook's")
    return wide


def toolbox_paths(root):
    r"""The toolbox's inputs under the ISMIP7 data root (``.../ISMIP7/AIS``)."""
    ocean = os.path.join(root, "parameterisations", "ocean")
    clim = os.path.join(root, "obs", "ocean", "climatology", "zhou_annual_06_nov")
    return {
        "term2": os.path.join(ocean, "meltobs", "melt_target_term2_v3.nc"),
        "term3_cold": os.path.join(ocean, "ocean_modelling_data",
                                   "melt_cold_target_term3_v2.nc"),
        "term3_warm": os.path.join(ocean, "ocean_modelling_data",
                                   "melt_warm_target_term3_v2.nc"),
        "term4": os.path.join(ocean, "ocean_observations_data",
                              "melt_observations_target_term4.nc"),
        "bfrn8": os.path.join(ocean, "bfrns", "BFRN_ismip8km_v2.nc"),
        "shelf": os.path.join(ocean, "shelfmask", "shelf_mask_ismip8km.nc"),
        "imbie2": os.path.join(ocean, "imbie2", "basin_numbers_ismip8km_v2.nc"),
        "floating8": os.path.join(ocean, "floatingmasks",
                                  "floatingmask_ismip8km.nc"),
        "bedmap3": os.path.join(root, "obs", "ocean", "topography", "bedmap3",
                                "v3", "bedmap3_AIS_obs_ocean_topography_v3.nc"),
        # The notebook's climatology, pinned by path (cell 7).
        "tf_06nov": os.path.join(
            clim, "tf", "v3",
            "tf_AIS_obs_ocean_climatology_zhou_annual_06_nov_v3_1972-2024.nc"),
        "so_06nov": os.path.join(
            clim, "so", "v4",
            "so_AIS_obs_ocean_climatology_zhou_annual_06_nov_v4_1972-2024.nc"),
    }


def bfrn_path(root, resolution_m):
    r"""The buttressing-bin file at a resolution in m (1000, 2000, ... 64000)."""
    km = int(round(float(resolution_m) / 1000.0))
    return os.path.join(root, "parameterisations", "ocean", "bfrns",
                        f"BFRN_ismip{km}km_v2.nc")


def state_files(root):
    r"""``(kind, label, tf path, so path)`` for every ocean state terms 3 and 4
    melt under: kind ``cold`` or ``warm`` with a model label, or ``obs`` with
    a year."""
    ocean = os.path.join(root, "parameterisations", "ocean")
    states = []
    for kind in ("cold", "warm"):
        for label, prefix in NOTEBOOK_MODELS:
            stem = os.path.join(ocean, "ocean_modelling_data",
                                prefix.format(state=kind))
            states.append((kind, label, stem + "TF.nc", stem + "S.nc"))
    for year in NOTEBOOK_OBS_YEARS:
        stem = os.path.join(ocean, "ocean_observations_data", f"Obs_{year}_")
        states.append(("obs", year, stem + "TF.nc", stem + "S.nc"))
    return states


def _ascending(axis, data, dim):
    r"""The axis ascending, with ``data`` flipped along ``dim`` to match."""
    axis = np.asarray(axis)
    if axis[0] > axis[-1]:
        return axis[::-1], np.flip(data, dim)
    return axis, data


def _nearest(axes, data, points, fill):
    from scipy.interpolate import RegularGridInterpolator
    if not np.isnan(fill):
        data = np.nan_to_num(data, nan=fill)
    interp = RegularGridInterpolator(axes, data, method="nearest",
                                     bounds_error=False, fill_value=fill)
    return interp(points)


def sample_nearest(nc_path, var, x, y, draft=None, fill=np.nan):
    r"""Nearest-neighbour sample of ``var`` in ``nc_path`` at points (x, y),
    and at ``draft`` (clipped into the depth range) when the field has a z
    axis. The field is loaded once, as float32, the forward's precision.

    ``fill`` stands in for NaN in the source and for points off the grid.
    The forward and the calibrations read the climatology with 0; NaN keeps a
    coverage gap visible, so a basin mean can skip it and a label field can
    say that a point has no label."""
    import xarray as xr

    ds = xr.open_dataset(nc_path)
    da = ds[var]
    zdim = next((d for d in da.dims if d.lower() in ("z", "depth", "lev")), None)
    if zdim is not None:
        data = da.transpose(zdim, "y", "x").values.astype(np.float32)
        z_arr, data = _ascending(ds[zdim].values, data, 0)
        y_arr, data = _ascending(ds["y"].values, data, 1)
        x_arr, data = _ascending(ds["x"].values, data, 2)
        z_clipped = np.clip(np.asarray(draft), z_arr[0], z_arr[-1])
        out = _nearest((z_arr, y_arr, x_arr), data,
                       np.column_stack([z_clipped, y, x]), fill)
    else:
        data = da.transpose("y", "x").values.astype(np.float32)
        y_arr, data = _ascending(ds["y"].values, data, 0)
        x_arr, data = _ascending(ds["x"].values, data, 1)
        out = _nearest((y_arr, x_arr), data, np.column_stack([y, x]), fill)
    ds.close()
    return out


def sample_grid(x_axis, y_axis, data, x, y, fill=np.nan):
    r"""Nearest-neighbour sample of an in-memory (y, x) field at (x, y)."""
    data = np.asarray(data, dtype=np.float64)
    y_arr, data = _ascending(y_axis, data, 0)
    x_arr, data = _ascending(x_axis, data, 1)
    return _nearest((y_arr, x_arr), data, np.column_stack([y, x]), fill)


def as_labels(values):
    r"""Integer labels from a sampled label field; -1 where it has none."""
    v = np.asarray(values, dtype=np.float64)
    out = np.full(v.shape, -1, dtype=np.int64)
    ok = np.isfinite(v)
    out[ok] = np.round(v[ok]).astype(np.int64)
    return out


def amundsen_regions(shelf_nc, x_axis=None, y_axis=None):
    r"""The term 4 regions on the shelf mask's 8 km grid: 1 on Pine Island's
    main trunk (shelf 110 east of x = -1.625e6), 2 on Dotson (shelf 97), 0
    elsewhere, built as notebook cell 20 builds them, x cut at the 8 km
    centres. Returns ``(x, y, labels)``; ``x_axis`` and ``y_axis`` stand in
    for a file that carries no coordinates."""
    import xarray as xr

    with xr.open_dataset(shelf_nc) as ds:
        shelves = ds["shelf_mask"]
        if "time" in shelves.dims:
            shelves = shelves.isel(time=0)
        s = shelves.transpose("y", "x").values
        xa = np.asarray(ds["x"].values) if "x" in ds.variables else np.asarray(x_axis)
        ya = np.asarray(ds["y"].values) if "y" in ds.variables else np.asarray(y_axis)
    pig = (s == PIG_ID) & (xa[None, :] > PIG_X_MIN)
    dotson = s == DOTSON_ID
    return xa, ya, pig.astype(np.int64) + 2 * dotson.astype(np.int64)


def read_melt_table(path):
    r"""Basin ids, observed melt and its uncertainty, both in Gt/yr.

    The two published tables differ in width: the Paolo+Adusumilli one carries
    area and per-area columns between melt and its uncertainty, the combined
    Paolo+Davison+Adusumilli one carries melt and uncertainty alone. Index 3 is
    the uncertainty in the first and past the end of the second, and any index
    chosen for one width reads the wrong quantity or nothing at the other.
    Columns are located by header name so the reader takes either.
    """
    bids, mobs, sobs = [], [], []
    with open(path) as f:
        r = csv.reader(f)
        header = next(r)

        def column(want):
            for i, name in enumerate(header):
                if name.strip().lower() == want:
                    return i
            raise ValueError(
                f"{path}: no {want!r} column in header {header}")

        i_m = column("bmr (gt/yr)")
        i_s = column("bmr uncert (gt/yr)")
        for row in r:
            if not row or not row[i_m]:
                continue
            bids.append(int(row[0]))
            mobs.append(float(row[i_m]))
            sobs.append(float(row[i_s]))
    return np.array(bids), np.array(mobs), np.array(sobs)


def fit_deltaT(tf, sal, sin_a, K, floating, area, basin, bids, M_obs, comm,
               window=None):
    r"""Per-basin offsets, with ``M(dT=0)``, the residual and ``dM/dT`` at
    the root, all in Gt/yr. The arrays are this rank's dofs; every basin
    total is reduced over ``comm``, so the root, and the file, are the same at
    any rank count. Collective: call on every rank.

    The root is searched in ``window`` (``DT_WINDOW`` when None); a basin
    with no root there takes the end point with the smaller residual and is
    flagged."""
    from scipy.optimize import brentq
    dT = np.full(len(bids), np.nan)
    M0 = np.zeros(len(bids))
    resid = np.full(len(bids), np.nan)
    sens = np.full(len(bids), np.nan)
    flagged = []
    for i, bid in enumerate(bids):
        sel = floating & (basin == bid)
        if global_count(sel, comm) == 0:
            flagged.append((bid, "no floating cells"))
            continue
        tf_b, s_b, a_b, A_b = tf[sel], sal[sel], sin_a[sel], area[sel]

        def M(d):
            m = quadratic_mixed_slope(tf_b + d, s_b, a_b, K=K)
            return global_sum(m * A_b, comm) * float(_RHO_I) / 1e12

        f = lambda d: M(d) - M_obs[i]  # noqa: E731
        M0[i] = M(0.0)
        lo, hi = DT_WINDOW if window is None else window
        if f(lo) * f(hi) < 0:
            dT[i] = brentq(f, lo, hi, xtol=1e-4)
        else:
            dT[i] = lo if abs(f(lo)) < abs(f(hi)) else hi
            flagged.append((bid, f"no root in [{lo:g}, {hi:g}] K, "
                                 f"took the end point"))
        resid[i] = f(dT[i])
        sens[i] = (M(dT[i] + 0.05) - M(dT[i] - 0.05)) / 0.1
    return dT, M0, resid, sens, flagged


def offsets_on_cells(dT, bids, basin):
    r"""Each cell's offset: its basin's, 0 for a cell in no fitted basin, as
    ``forcing.load_deltaT_per_basin`` stamps a file with ``fill`` 0."""
    field = np.zeros(len(basin))
    for bid, d in zip(bids, dT):
        if np.isfinite(d):
            field[basin == bid] = d
    return field


def notebook_mean_slope(draft, floating, dx, dy):
    r"""sin of the mean draft slope angle over ``floating``, by the notebook's
    recipe (cell 9): along each axis a centred difference, one-sided where a
    neighbour is missing and zero where both are; alpha = arctan |grad draft|,
    averaged over the floating points. ``draft`` is (y, x)."""
    d = np.asarray(draft, dtype=np.float64)

    def shift(a, n, axis):
        out = np.full_like(a, np.nan)
        src = [slice(None)] * a.ndim
        dst = [slice(None)] * a.ndim
        if n > 0:
            src[axis], dst[axis] = slice(None, -n), slice(n, None)
        else:
            src[axis], dst[axis] = slice(-n, None), slice(None, n)
        out[tuple(dst)] = a[tuple(src)]
        return out

    def along(axis, h):
        before, after = shift(d, 1, axis), shift(d, -1, axis)
        both = (after - before) / (2.0 * h)
        right = (d - before) / h
        left = (after - d) / h
        s = np.where(np.isfinite(both), both, np.where(np.isfinite(right), right, left))
        return np.where(np.isfinite(s), s, 0.0)

    alpha = np.arctan(np.sqrt(along(1, dx) ** 2 + along(0, dy) ** 2))
    sel = np.asarray(floating, bool) & np.isfinite(alpha)
    return float(np.sin(alpha[sel].mean()))


class Cells:
    r"""The cells a selection melts.

    ``basin``, ``bfrn`` and ``region`` are integer labels, -1 for none;
    region 1 is Pine Island and 2 Dotson. Terms 1 to 3 add up the cells in
    ``agg`` (the model's floating mask, all cells by default); term 4 adds up
    every cell with a region, as the toolbox's ``calculate_term4`` applies no
    floating mask."""

    def __init__(self, area, basin, bfrn, region, sin_alpha, agg=None):
        self.area = np.asarray(area, dtype=np.float64)
        n = len(self.area)
        self.basin = np.asarray(basin, dtype=np.int64)
        self.bfrn = np.asarray(bfrn, dtype=np.int64)
        self.region = np.asarray(region, dtype=np.int64)
        self.sin_alpha = np.broadcast_to(
            np.asarray(sin_alpha, dtype=np.float64), (n,))
        agg = np.ones(n, bool) if agg is None else np.asarray(agg, bool)
        self.basin_agg = np.where(agg, self.basin, -1)
        self.bfrn_agg = np.where(agg, self.bfrn, -1)


def label_totals(melt, area, labels, n):
    r"""Gt/yr per label 0..n-1 from melt in kg/m2/yr on cells of area m2,
    over the cells with finite melt; 0 where a label has none."""
    ok = np.isfinite(melt) & (labels >= 0) & (labels < n)
    w = np.asarray(melt, np.float64)[ok] * area[ok]
    return np.bincount(labels[ok], weights=w, minlength=n) / 1e12


def label_means(melt, area, labels, n):
    r"""Area-weighted mean per label 0..n-1 of melt over the cells with
    finite melt; NaN where a label has none."""
    ok = np.isfinite(melt) & (labels >= 0) & (labels < n)
    m = np.asarray(melt, np.float64)[ok]
    num = np.bincount(labels[ok], weights=m * area[ok], minlength=n)
    den = np.bincount(labels[ok], weights=area[ok], minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def melt_rate(tf, so, sin_alpha, K, rho):
    r"""kg/m2/yr by the forward's law; NaN wherever tf or so is."""
    return quadratic_mixed_slope(tf, so, sin_alpha, K=K) * rho


def aggregate(K, dT_cell, present, states, cells, rho):
    r"""The four aggregates at one K, before the toolbox's zero rules.

    ``present`` is ``(tf, so)`` of the present-day climatology on the cells,
    ``states`` is ``[(kind, label, tf, so), ...]`` as ``state_files`` lists
    them, sampled onto the cells, and ``dT_cell`` the offset every one of
    them melts with. Returns ``t1`` (16 basins), ``t2`` (10 bins), ``t3_cold``
    and ``t3_warm`` (6 models by 16 basins) and ``t4`` (2 regions by 13
    years)."""
    models = [label for label, _ in NOTEBOOK_MODELS]
    tf, so = present
    melt = melt_rate(tf + dT_cell, so, cells.sin_alpha, K, rho)
    out = {
        "t1": label_totals(melt, cells.area, cells.basin_agg, N_BASINS),
        "t2": label_totals(melt, cells.area, cells.bfrn_agg, N_BFRN_BINS),
        "t3_cold": np.full((len(models), N_BASINS), np.nan),
        "t3_warm": np.full((len(models), N_BASINS), np.nan),
        "t4": np.zeros((len(REGIONS), len(NOTEBOOK_OBS_YEARS))),
    }
    for kind, label, tf_s, so_s in states:
        m = melt_rate(tf_s + dT_cell, so_s, cells.sin_alpha, K, rho)
        if kind == "obs":
            j = NOTEBOOK_OBS_YEARS.index(label)
            out["t4"][:, j] = label_totals(m, cells.area, cells.region,
                                           len(REGIONS) + 1)[1:]
        else:
            out["t3_" + kind][models.index(label)] = label_means(
                m, cells.area, cells.basin_agg, N_BASINS)
    return out


class TFRule:
    r"""Which K the present-day thermal forcing admits once each basin's
    offset is applied.

    Section 2.1 asks the offsets to keep present-day thermal forcing "not
    significantly below 0 degC or above 5 degC". This repository reads
    significantly as a bound every floating cell must meet (``floor``,
    ``cap``) and a bound at most a share of each basin's area may cross
    (``floor_area``, ``cap_area``, each a (degC, share) pair). The defaults
    hold every cell at or above -1.8 degC, at most 25 % of any basin's area
    below -1.0 degC and at most 25 % above 5.5 degC, and bound no single cell
    on the warm side. ``rooted`` also refuses a K at which some basin cannot
    reach its observed total inside the offset window; the toolbox keeps such
    a K with a term 1 penalty. None switches a test off."""

    def __init__(self, floor=-1.8, floor_area=(-1.0, 0.25), cap=None,
                 cap_area=(5.5, 0.25), rooted=True):
        def pair(p):
            return None if p is None else (float(p[0]), float(p[1]))
        self.floor = None if floor is None else float(floor)
        self.cap = None if cap is None else float(cap)
        self.floor_area, self.cap_area = pair(floor_area), pair(cap_area)
        self.rooted = bool(rooted)

    def as_dict(self):
        return {"floor": self.floor, "floor_area": self.floor_area,
                "cap": self.cap, "cap_area": self.cap_area,
                "rooted": self.rooted}

    def admits(self, per):
        r"""A bool per K, true where every test passes, and each test's own
        bool per K. ``per`` holds the per-K, per-basin statistics
        (``unrooted``, ``tf_mean``, ``tf_min``, ``tf_max`` and, for the area
        tests, ``frac_below_floor`` and ``frac_above_cap`` from
        ``Plausibility.at`` with this rule). A basin with no floating cells
        has NaN statistics and takes part in no test."""
        has = np.isfinite(np.asarray(per["tf_mean"], np.float64))
        tests = {}
        with np.errstate(invalid="ignore"):
            if self.rooted:
                tests["rooted"] = ~((np.asarray(per["unrooted"]) > 0) & has).any(axis=1)
            if self.floor is not None:
                tests["floor"] = ~(np.asarray(per["tf_min"]) < self.floor).any(axis=1)
            if self.floor_area is not None:
                tests["floor_area"] = ~(np.asarray(per["frac_below_floor"])
                                        > self.floor_area[1]).any(axis=1)
            if self.cap is not None:
                tests["cap"] = ~(np.asarray(per["tf_max"]) > self.cap).any(axis=1)
            if self.cap_area is not None:
                tests["cap_area"] = ~(np.asarray(per["frac_above_cap"])
                                      > self.cap_area[1]).any(axis=1)
        ok = np.ones(has.shape[0], bool)
        for passed in tests.values():
            ok &= passed
        return ok, tests


class Plausibility:
    r"""Present-day thermal forcing with each basin's offset, against the
    protocol's plausible range of 0 to 5 degC (section 2.1): per basin the
    area-weighted mean, min and max, and the area fractions below 0 and above
    5, and beyond a ``TFRule``'s area thresholds when one is given. The
    offset is one number per basin, so min, max and mean shift with it and
    only the fractions are summed again at each K."""

    def __init__(self, tf, area, basin):
        tf = np.asarray(tf, np.float64)
        ok = np.isfinite(tf) & (basin >= 0) & (basin < N_BASINS)
        self.tf, self.area, self.basin = tf[ok], area[ok], basin[ok]
        self.total = np.bincount(self.basin, weights=self.area, minlength=N_BASINS)
        with np.errstate(invalid="ignore", divide="ignore"):
            self.mean = np.bincount(self.basin, weights=self.area * self.tf,
                                    minlength=N_BASINS) / self.total
        self.lo = np.full(N_BASINS, np.inf)
        self.hi = np.full(N_BASINS, -np.inf)
        np.minimum.at(self.lo, self.basin, self.tf)
        np.maximum.at(self.hi, self.basin, self.tf)
        empty = self.total <= 0
        self.mean[empty] = self.lo[empty] = self.hi[empty] = np.nan

    def at(self, dT_basin, rule=None):
        r"""The statistics per basin for offsets ``dT_basin`` (16): five, and
        with a ``rule`` the area fractions below its floor and above its cap
        area thresholds."""
        d = np.nan_to_num(np.asarray(dT_basin, np.float64), nan=0.0)
        shifted = self.tf + d[self.basin]

        def share(mask):
            with np.errstate(invalid="ignore", divide="ignore"):
                return np.bincount(self.basin, weights=self.area * mask,
                                   minlength=N_BASINS) / self.total
        out = {"tf_mean": self.mean + d, "tf_min": self.lo + d,
               "tf_max": self.hi + d, "frac_below_0": share(shifted < 0.0),
               "frac_above_5": share(shifted > 5.0)}
        if rule is not None and rule.floor_area is not None:
            out["frac_below_floor"] = share(shifted < rule.floor_area[0])
        if rule is not None and rule.cap_area is not None:
            out["frac_above_cap"] = share(shifted > rule.cap_area[0])
        return out


def load_targets(paths, melt_csv):
    r"""The observation targets and the term 2 weights, as the notebook reads
    them (cells 24 and 34)."""
    import xarray as xr

    bids, M_obs, sigma = read_melt_table(melt_csv)
    if not np.array_equal(bids, np.arange(N_BASINS)):
        raise ValueError(f"{melt_csv}: basins {bids.tolist()}, expected 0..15 "
                         f"in order (the toolbox indexes them by position)")
    bfrn8 = xr.load_dataset(paths["bfrn8"])
    return {
        "bids": bids, "M_obs": M_obs, "sigma": sigma,
        "term2": xr.load_dataset(paths["term2"]),
        "cold": xr.load_dataset(paths["term3_cold"]),
        "warm": xr.load_dataset(paths["term3_warm"]),
        "term4": xr.load_dataset(paths["term4"]),
        # the 8 km medians, whatever the model resolution (cell 34)
        "t2_weight": (bfrn8["BFRN_medians"] / bfrn8["BFRN_median"]).values,
    }


def toolbox_terms(K, agg, targets):
    r"""The model, observation mean and sigma DataArrays of the four terms, in
    the layout and with the rules of the toolbox's ``calculate_term1..4``.

    ``agg`` holds the aggregates stacked over K, first axis p1: ``t1``
    (p1, 16), ``t2`` (p1, 10 bins by label), ``t3_cold`` and ``t3_warm``
    (p1, 6 models, 16), ``t4`` (p1, 2 regions, 13 years). Labels come from
    the target files, so every inner join inside the objective keeps them:
    float basins 0..15, the term 2 target's bins, its region and year."""
    import xarray as xr

    K = np.asarray(K)
    p2 = np.ones(1)
    basins = np.asarray(targets["cold"]["basins"].values)
    if not np.array_equal(basins, np.arange(N_BASINS)):
        raise ValueError(f"term 3 target basins {basins}, expected 0..15")
    terms = {}

    t1 = xr.DataArray(np.asarray(agg["t1"])[None], dims=("p2", "p1", "basins"),
                      coords={"p2": p2, "p1": K, "basins": basins})
    terms["t1_model"] = t1.where(t1 != 0, np.nan)
    terms["t1_obs_mean"] = xr.DataArray(
        targets["M_obs"], name="melt_Gt_per_y", dims=["basin"],
        coords={"basin": range(N_BASINS)})
    terms["t1_obs_sigma"] = xr.DataArray(
        targets["sigma"], dims=["basin"], coords={"basin": range(N_BASINS)},
        name="melt_unc_Gt_per_y")

    term2 = targets["term2"]
    bins = term2["BFRN_bins"].values
    idx = np.round(np.asarray(bins, np.float64)).astype(int)
    if sorted(idx.tolist()) != list(range(N_BFRN_BINS)):
        raise ValueError(f"term 2 target bins {bins}, expected 0..9")
    t2 = xr.DataArray(np.asarray(agg["t2"])[:, idx][None],
                      dims=("p2", "p1", "BFRN_bins"),
                      coords={"p2": p2, "p1": K, "BFRN_bins": bins})
    terms["t2_model"] = t2.where(t2 != 0, np.nan)
    terms["t2_obs_mean"] = term2["melt_mean"]
    terms["t2_obs_sigma"] = term2["melt_mean_err"]

    models = [label for label, _ in NOTEBOOK_MODELS]
    missing = sorted(set(models) - set(targets["cold"]["model"].values.tolist()))
    if missing:
        raise ValueError(f"term 3 targets lack the models {missing}")
    means = {}
    for kind in ("cold", "warm"):
        m = np.array(agg["t3_" + kind], dtype=np.float64)
        for label, allowed in NOTEBOOK_MODEL_BASINS.items():
            keep = np.isin(np.arange(N_BASINS), allowed)
            m[:, models.index(label), ~keep] = np.nan
        means[kind] = np.where(m != 0, m, np.nan)
    diff = np.transpose(means["warm"] - means["cold"], (1, 0, 2))[:, None]
    terms["t3_model"] = xr.DataArray(
        diff, dims=("model", "p2", "p1", "basins"),
        coords={"model": models, "p2": p2, "p1": K, "basins": basins})
    warm, cold = targets["warm"], targets["cold"]
    terms["t3_obs_mean"] = warm.melt_rate - cold.melt_rate
    terms["t3_obs_sigma"] = np.sqrt(warm.melt_rate_uncert ** 2
                                    + cold.melt_rate_uncert ** 2)

    t4_obs = targets["term4"]
    regions = [str(r) for r in t4_obs["region"].values]
    years = [int(y) for y in t4_obs["year"].values]
    if sorted(regions) != sorted(REGIONS) or sorted(years) != sorted(NOTEBOOK_OBS_YEARS):
        raise ValueError(f"term 4 target covers {regions} x {years}")
    t4 = np.asarray(agg["t4"], np.float64)            # (p1, region, year)
    t4 = t4[:, [REGIONS.index(r) for r in regions]][:, :, [NOTEBOOK_OBS_YEARS.index(y) for y in years]]
    t4 = xr.DataArray(np.transpose(t4, (2, 0, 1))[:, None],
                      dims=("year", "p2", "p1", "region"),
                      coords={"year": t4_obs["year"].values, "p2": p2, "p1": K,
                              "region": t4_obs["region"].values})
    t4 = t4.where(t4 != 0, np.nan)
    t4 = t4.where(t4_obs.melt_rate.notnull())
    terms["t4_model"] = t4.reindex_like(t4_obs.melt_rate)
    terms["t4_obs_mean"] = t4_obs.melt_rate
    terms["t4_obs_sigma"] = t4_obs.melt_rate_uncert
    return terms


def notebook_weights(terms, t2_weight, t4_obs):
    r"""The weights of notebook cells 33 to 36, verbatim: one per basin in
    term 1, BFRN_medians / BFRN_median per bin in term 2, one for mathiot,
    naughten_ais_1, jourdain_naughten and naughten_naughten in term 3, and one
    for Pine Island in 2009 and 2012 in term 4, zero elsewhere."""
    import xarray as xr

    t1_model, t3_model = terms["t1_model"], terms["t3_model"]
    t1_weights = xr.DataArray(
        t1_model.isel(p1=0, p2=0) * 0 + 1,
        dims=["basins"],
        coords={"basins": t1_model.basins.values})
    t2_weights = xr.DataArray(
        t2_weight,
        dims=["BFRN_bins"],
        coords={"BFRN_bins": terms["t2_obs_mean"].BFRN_bins.values})
    t3_weights = xr.DataArray(
        np.ones(t3_model.isel(p1=0, p2=0).shape),
        dims=["model", "basins"],
        coords={"model": t3_model.model.values,
                "basins": t3_model.basins.values})
    keep = False
    for label in NOTEBOOK_T3_MODELS:
        keep = keep | (t3_weights.model == label)
    t3_weights = t3_weights.where(keep, other=0)
    t4_weights = xr.DataArray(
        np.ones(terms["t4_model"].isel(p1=0, p2=0).shape),
        dims=["year", "region"],
        coords={"year": t4_obs.year.values,
                "region": t4_obs.region.values})
    t4_weights = t4_weights.where(t4_weights.region == NOTEBOOK_T4_REGION, other=0)
    keep = False
    for year in NOTEBOOK_T4_YEARS:
        keep = keep | (t4_weights.year == year)
    t4_weights = t4_weights.where(keep, other=0)
    return {"t1_weights": t1_weights, "t2_weights": t2_weights,
            "t3_weights": t3_weights, "t4_weights": t4_weights}


def check_inputs(inputs):
    r"""Refuse inputs the objective would silently mishandle.

    Every term must carry the same p1 labels, since the objective's inner
    joins drop any K a term lacks. Terms 1 and 2 are averaged without
    skipping NaN, so a NaN basin or bin at one K removes that K, and a NaN
    term 1 weight removes every K and leaves idxmin nothing to return."""
    p1 = inputs["t1_model"].p1.values
    for name in ("t2_model", "t3_model", "t4_model"):
        if not np.array_equal(inputs[name].p1.values, p1):
            raise ValueError(f"{name} carries other p1 labels than t1_model")
    for name, label in (("t1_model", "basins"), ("t2_model", "BFRN_bins")):
        bad = ~np.isfinite(inputs[name].values)
        if bad.any():
            where = inputs[name].where(~np.isfinite(inputs[name]), drop=True)
            raise ValueError(
                f"{name} is not finite at {int(bad.sum())} entries "
                f"(p1 {where.p1.values[:4].tolist()}..., {label} "
                f"{where[label].values.tolist()}); the objective would drop "
                f"those K without a message")
    if not np.isfinite(inputs["t1_weights"].values).all():
        raise ValueError("t1_weights is not finite")
    if not np.array_equal(inputs["t2_model"].BFRN_bins.values,
                          inputs["t2_obs_mean"].BFRN_bins.values):
        raise ValueError("t2_model bins differ from the term 2 target's")


def select(inputs, sample_size=NOTEBOOK_SAMPLE_SIZE, chunk=10000, seed=0,
           reso=None):
    r"""``min_p1`` of ``sample_size`` draws of the toolbox's objective.

    The objective runs in chunks of ``chunk`` samples after one
    ``np.random.seed(seed)``. Every draw belongs to one sample and the median
    normalisation is taken over p1 and p2 within a sample, so chunks pool to
    the distribution one call draws from; a single chunk of ``sample_size``
    is that call, bit for bit. Chunks of 10000 hold term 3 near 0.1 GB, where
    one call of 100000 needs about 35 GB."""
    from . import ismip7_parameter_selection_toolbox as pst

    check_inputs(inputs)
    np.random.seed(seed)
    out, left = [], int(sample_size)
    while left > 0:
        n = min(int(chunk), left)
        p1, _ = pst.calculate_objective_function(
            n, reso, *(inputs[name] for name in TERM_ARGS))
        out.append(np.asarray(p1, dtype=np.float64))
        left -= n
    return np.concatenate(out)


def summarise(min_p1, K, unrooted=None):
    r"""The notebook's printout (cell 40) of a sample of ``min_p1`` on the K
    grid, and where the samples fell: the share on the grid's two edges and,
    given ``unrooted`` (a bool per K, true where some basin's offset hit the
    window), the share on such a K."""
    s = np.asarray(min_p1, dtype=np.float64)
    K = np.asarray(K, dtype=np.float64)
    edges = np.append(K[0] * 0.5, K + 1e-7)
    counts, _ = np.histogram(s, bins=edges)
    idx = np.abs(s[:, None] - K[None, :]).argmin(axis=1)
    out = {"n": int(len(s))}
    for q in (1, 5, 50, 95, 99):
        v = float(np.percentile(s, q))
        out[f"p{q}"] = v
        out[f"p{q}_on_grid"] = bool(np.isclose(v, K, rtol=1e-9, atol=0).any())
    out.update({
        "median": float(np.median(s)),
        "mode": float(K[np.argmax(counts)]),
        "min": float(s.min()), "max": float(s.max()),
        "counts": counts.tolist(),
        "share_on_grid_edges": float(np.mean((idx == 0) | (idx == len(K) - 1))),
    })
    if unrooted is not None:
        out["share_on_unrooted_K"] = float(np.mean(np.asarray(unrooted)[idx]))
    return out


def write_offsets(path, bids, dT, K, M_obs, M0, resid, sens, **provenance):
    r"""The offsets at one K in calibrate_deltaT.py's keys, so
    ``ISMIP7_DELTAT_PER_BASIN_NPZ`` can name the file; ``provenance`` adds
    the slope and geometry conventions and whatever else the caller
    records."""
    np.savez(path, basin_ids=bids, deltaT_basin=dT, K=K, M_obs=M_obs,
             M_dT0=M0, residual_gt=resid, sensitivity_gt_per_K=sens,
             **provenance)


def notebook_terms(K, rho, sin_alpha, present, states, mask_m, basins_m,
                   bfrn_m, region_map, targets, melt_table):
    r"""The four terms and weights by the notebook's own route (cells 10 to
    36): gridded melt ensembles through the vendored ``calculate_term1..4``.
    For checking ``toolbox_terms`` against on a regular grid.

    ``present`` is ``(tf, so)`` and ``states`` maps ``(kind, label)`` to
    ``(tf, so)``, all (y, x) DataArrays on the grid of ``mask_m`` with NaN
    where the notebook does not melt. ``basins_m`` is the basin field named
    ``basins``, ``bfrn_m`` a Dataset holding ``BFRN_bins``, ``region_map``
    the (y, x) labels of ``amundsen_regions`` and ``melt_table`` the per-basin
    observation CSV. Melt is computed as ``aggregate`` computes it, so both
    routes start from the same numbers."""
    import pandas as pd
    import xarray as xr
    from . import ismip7_parameter_selection_toolbox as pst

    K = np.asarray(K)
    coords = {"y": mask_m["y"].values, "x": mask_m["x"].values}
    reso = float(abs(coords["x"][1] - coords["x"][0]))
    cvt_m = reso ** 2 / 1e12

    def ensemble(tf, so):
        tf = np.asarray(tf.values, np.float64)
        so = np.asarray(so.values, np.float64)
        melt = np.stack([melt_rate(tf + 0.0, so, sin_alpha, k, rho) for k in K])
        return xr.DataArray(melt[None], dims=("p2", "p1", "y", "x"),
                            coords={"p2": np.ones(1), "p1": K, **coords})

    tb = {}
    pd_ensemble = ensemble(*present).to_dataset(name="melt_rate")
    melt_data = pd.read_csv(melt_table, index_col=0)
    tb["t1_model"], tb["t1_obs_mean"], tb["t1_obs_sigma"] = pst.calculate_term1(
        pd_ensemble, mask_m, basins_m, int(basins_m.max()), cvt_m, melt_data)
    tb["t2_model"], tb["t2_obs_mean"], tb["t2_obs_sigma"] = pst.calculate_term2(
        pd_ensemble, mask_m, bfrn_m, cvt_m, targets["term2"])
    del pd_ensemble

    models = [label for label, _ in NOTEBOOK_MODELS]
    by_kind = {}
    for kind in ("cold", "warm"):
        e = xr.concat([ensemble(*states[(kind, label)]).expand_dims({"model": [j]})
                       for j, label in enumerate(models)], dim="model")
        e = e.assign_coords(model=models).to_dataset(name="melt_rate")
        for label, allowed in NOTEBOOK_MODEL_BASINS.items():
            keep = e.model != label
            for b in allowed:
                keep = keep | (basins_m == b)
            e = e.where(keep, np.nan)
        by_kind[kind] = e
    tb["t3_model"], tb["t3_obs_mean"], tb["t3_obs_sigma"] = pst.calculate_term3(
        by_kind["cold"], by_kind["warm"], targets["cold"], targets["warm"],
        mask_m, basins_m)
    del by_kind

    obs = xr.concat([ensemble(*states[("obs", year)]).expand_dims({"year": [year]})
                     for year in NOTEBOOK_OBS_YEARS], dim="year")
    obs = obs.to_dataset(name="melt_rate")
    rlab = xr.DataArray(np.asarray(region_map), dims=("y", "x"), coords=coords)
    labels = np.vectorize({1: "pig", 2: "dotson"}.get, otypes=[object])(rlab.values)
    region_label_m = xr.DataArray(labels, dims=rlab.dims, coords=rlab.coords)
    tb["t4_model"], tb["t4_obs_mean"], tb["t4_obs_sigma"] = pst.calculate_term4(
        obs, region_label_m, targets["term4"], mask_m, cvt_m)
    del obs
    tb.update(notebook_weights(tb, targets["t2_weight"], targets["term4"]))
    return tb


def compare_terms(ours, theirs):
    r"""Per model term: whether the NaN patterns agree and the largest
    relative difference where both are finite, once both are in one layout
    (bins and basins sorted, regions and models in our order)."""
    layouts = {"t1_model": ("p2", "p1", "basins"),
               "t2_model": ("p2", "p1", "BFRN_bins"),
               "t3_model": ("model", "p2", "p1", "basins"),
               "t4_model": ("year", "p2", "p1", "region")}
    out = {}
    for name, dims in layouts.items():
        a = ours[name].transpose(*dims)
        b = theirs[name].transpose(*dims)
        if dims[-1] == "region":
            b = b.reindex(region=a.region.values)
        else:
            a, b = a.sortby(dims[-1]), b.sortby(dims[-1])
        if name == "t3_model":
            b = b.reindex(model=a.model.values)
        if a.shape != b.shape:
            out[name] = {"same_nan_pattern": False, "max_rel_diff": float("inf"),
                         "shapes": [list(a.shape), list(b.shape)]}
            continue
        av, bv = a.values, b.values
        ok = np.isfinite(av) & np.isfinite(bv)
        rel = (np.abs(av[ok] - bv[ok]) / np.maximum(np.abs(bv[ok]), 1e-300)
               if ok.any() else np.zeros(1))
        out[name] = {"same_nan_pattern": bool(np.array_equal(np.isnan(av),
                                                             np.isnan(bv))),
                     "max_rel_diff": float(rel.max()), "shape": list(av.shape)}
    return out
