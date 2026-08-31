r"""ISMIP7 forcing data reader for Antarctic simulations."""

import os
import numpy as np

_SEC_PER_YEAR = 31556926.0
_RHO_ICE = 917.0
_RHO_WATER = 1000.0

# Constants for the Burgard et al. 2022 quadratic-mixed-slope melt
# parameterization, taken verbatim from multimelt.constants
# (https://github.com/ClimateClara/multimelt). The ISMIP7 ocean-forcing
# pipeline calibrates K under exactly this decomposition.
_RHO_SW = 1028.0       # seawater, kg/m^3
_RHO_I = 917.0         # ice,      kg/m^3
_C_PO = 3974.0         # seawater specific heat, J/(kg K)
_L_I = 3.34e5          # latent heat of fusion of ice, J/kg
_BETA_S = 7.86e-4      # haline contraction coefficient (Lazeroms), 1/PSU
_G = 9.81              # gravity, m/s^2
_F_CORIOLIS = 1.4e-4   # representative Antarctic Coriolis parameter, 1/s

# melt_factor = (rho_sw * c_po) / (rho_i * L_i)   [1/K]
_MELT_FACTOR = (_RHO_SW * _C_PO) / (_RHO_I * _L_I)

# K50 median from Burgard 2022 calibration
# (parameter_selection_quadratic_example.ipynb).
_K_DEFAULT = 11.5e-5

_DEFAULT_DATA_ROOT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "ISMIP7", "AIS",
)


def smb_kgm2s_to_myr(smb_kgm2s):
    r"""Convert SMB from kg/m^2/s to m/yr ice equivalent."""
    return smb_kgm2s * _SEC_PER_YEAR / _RHO_ICE * (_RHO_WATER / _RHO_ICE)


def load_racmo_smb_climatology(Q, clim_start=2000, clim_end=2029, data_dir=None,
                               target_res=8000.0, rho_ice=_RHO_ICE):
    r"""RACMO2.4p1 mean-annual SMB (m/yr ice equiv) as a Function on Q's mesh.

    The RACMO ANT11 grid is rotated-pole, so this reprojects the climatology to
    an intermediate EPSG:3031 raster with rasterio, then samples it onto the mesh
    with icepack.interpolate -- the same path used for BedMachine -- avoiding any
    scattered-point interpolation. ``smbgl`` is a monthly mass sum (kg/m^2), so the
    annual SMB is the sum of the 12 months, averaged over the climatology window.
    """
    import xarray as xr
    import pyproj
    import icepack
    from affine import Affine
    from rasterio.io import MemoryFile
    from rasterio.warp import reproject, Resampling

    if data_dir is None:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data_dir = os.path.join(base, "antarctica", "data", "racmo")
    fn = os.path.join(
        data_dir, "smbgl_monthlyS_ANT11_RACMO2.4p1_ERA5_197901_202312.nc"
    )
    ds = xr.open_dataset(fn)

    smb = ds["smbgl"].squeeze("height")  # (time, rlat, rlon), kg/m^2 per month
    yrs = smb["time"].dt.year.values
    sel = (yrs >= clim_start) & (yrs <= clim_end)
    n_years = len(np.unique(yrs[sel]))
    src = (smb.isel(time=sel).sum("time").values / n_years / rho_ice).astype("float64")

    rlon, rlat = ds["rlon"].values, ds["rlat"].values
    d = float(rlon[1] - rlon[0])
    src_crs = pyproj.CRS.from_cf(ds["rotated_pole"].attrs)
    src_transform = Affine.translation(rlon[0] - d / 2, rlat[0] - d / 2) * Affine.scale(d, d)
    ds.close()

    # Reproject onto the standard ISMIP AIS grid (EPSG:3031, centers +/- 3040 km).
    half = 3040000.0
    N = int(round(2 * half / target_res)) + 1
    dst_transform = (
        Affine.translation(-half - target_res / 2, half + target_res / 2)
        * Affine.scale(target_res, -target_res)
    )
    dst = np.full((N, N), np.nan, dtype="float64")
    reproject(
        src, dst,
        src_transform=src_transform, src_crs=src_crs,
        dst_transform=dst_transform, dst_crs="EPSG:3031",
        src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.bilinear,
    )
    dst[~np.isfinite(dst)] = 0.0  # zero SMB over ocean / outside RACMO coverage

    with MemoryFile() as mf:
        with mf.open(
            driver="GTiff", height=N, width=N, count=1, dtype="float64",
            crs="EPSG:3031", transform=dst_transform,
        ) as out:
            out.write(dst, 1)
        with mf.open() as raster:
            return icepack.interpolate(raster, Q)


def _find_ismip7_data(data_root=None):
    if data_root is not None:
        return data_root
    env = os.environ.get("ISMIP7_DATA_ROOT")
    if env and os.path.isdir(env):
        return env
    if os.path.isdir(_DEFAULT_DATA_ROOT):
        return _DEFAULT_DATA_ROOT
    return None


def load_mean_annual_surface_temperature(Q, var="tas", data_root=None,
                                         fill_K=260.0):
    r"""Mean-annual surface temperature [K] on ``Q``, from the ISMIP7
    CESM2-WACCM SDBN1 climatology. ``var='tas'`` (near-surface air temp) is the
    standard englacial upper boundary condition for a thermal model; ``'ts'``
    (skin temp) is melt-capped at 273 K. Averages the 12 monthly slices and
    samples with a nearest-neighbour ``RegularGridInterpolator`` (xarray.interp
    blows up memory at mesh-node counts; same pattern as the OI ocean reader).
    NODATA/ocean (~0 K) is masked and filled with ``fill_K``."""
    import glob
    import xarray as xr
    import firedrake as fd
    from scipy.interpolate import RegularGridInterpolator
    root = _find_ismip7_data(data_root)
    if root is None:
        raise FileNotFoundError("ISMIP7 data root not found (set ISMIP7_DATA_ROOT)")
    hits = glob.glob(os.path.join(root, "CESM2-WACCM", "climatology", "SDBN1",
                                  var, "v1", f"{var}_*.nc"))
    if not hits:
        raise FileNotFoundError(f"no {var} climatology under {root}")
    ds = xr.open_dataset(hits[0])
    da = ds[var] if var in ds else ds[list(ds.data_vars)[0]]
    T = np.asarray(da.mean(dim="month").values, dtype=float)     # (y, x)
    x = np.asarray(ds["x"].values, dtype=float)
    y = np.asarray(ds["y"].values, dtype=float)
    T = np.where(T > 100.0, T, np.nan)                           # mask fill (~0 K)
    T_filled = np.where(np.isfinite(T), T, fill_K)
    interp = RegularGridInterpolator((y, x), T_filled, method="nearest",
                                     bounds_error=False, fill_value=fill_K)
    Tf = fd.Function(Q, name="T_srf")
    xy = Q.mesh().coordinates.dat.data_ro
    Tf.dat.data[:] = interp(np.column_stack([xy[:, 1], xy[:, 0]]))
    return Tf


def _version_subdirs(parent_dir):
    r"""(N, name) pairs for the v<N> subdirs of a directory, ascending."""
    versions = []
    if parent_dir is not None and os.path.isdir(parent_dir):
        for name in os.listdir(parent_dir):
            if name.startswith("v") and os.path.isdir(
                    os.path.join(parent_dir, name)):
                try:
                    versions.append((int(name[1:]), name))
                except ValueError:
                    pass
    return sorted(versions)


def _resolve_version(parent_dir, pinned):
    r"""The pinned version if that subdir exists, else the highest v<N>
    subdir present, else the pinned name unchanged (so missing trees keep
    producing the same non-existent path the callers already handle)."""
    if parent_dir is None or os.path.isdir(os.path.join(parent_dir, pinned)):
        return pinned
    versions = _version_subdirs(parent_dir)
    return versions[-1][1] if versions else pinned


def atmosphere_path(scenario, esm="CESM2-WACCM", variable="acabf-anomaly",
                    resolution="8000m", version="v2", data_root=None):
    root = _find_ismip7_data(data_root)
    if root is None:
        return None
    parent = os.path.join(root, esm, scenario, f"SDBN1-{resolution}", variable)
    return os.path.join(parent, _resolve_version(parent, version))


def ocean_path(scenario, esm="CESM2-WACCM", variable="tf",
               version="v3", data_root=None):
    root = _find_ismip7_data(data_root)
    if root is None:
        return None
    parent = os.path.join(root, esm, scenario, "ocean", variable)
    return os.path.join(parent, _resolve_version(parent, version))


class ISMIP7Atmosphere:
    r"""Read ISMIP7 downscaled atmosphere forcing for Antarctica."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585",
                 resolution="8000m", version="v2"):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self.resolution = resolution
        self.version = version
        self._cache = {}
        self._grid_x = None
        self._grid_y = None

    def _var_dir(self, variable):
        return atmosphere_path(
            self.scenario, self.esm, variable,
            self.resolution, self.version, self.data_root,
        )

    def _load_year(self, variable, year):
        import xarray as xr

        key = (variable, int(year))
        if key in self._cache:
            return self._cache[key]

        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return None

        version = os.path.basename(vdir)
        pattern = f"{variable}_AIS_{self.esm}_{self.scenario}_SDBN1-{self.resolution}_{version}_{int(year)}.nc"
        path = os.path.join(vdir, pattern)

        if not os.path.exists(path):
            return None

        ds = xr.open_dataset(path)

        if self._grid_x is None:
            for xname in ["x", "X", "lon"]:
                if xname in ds.coords or xname in ds.dims:
                    self._grid_x = ds[xname].values
                    break
            for yname in ["y", "Y", "lat"]:
                if yname in ds.coords or yname in ds.dims:
                    self._grid_y = ds[yname].values
                    break

        if variable in ds.data_vars:
            da = ds[variable]
        else:
            # Skip grid-mapping/bounds variables (crs, *_bnds, mapping):
            # this dataset family carries them and any of them ordered
            # first would silently become the forcing field.
            data_vars = [
                v for v in ds.data_vars
                if not v.endswith("_bnds")
                and v.lower() not in ("x", "y", "time", "crs", "mapping",
                                      "spatial_ref", "lat", "lon")
            ]
            da = ds[data_vars[0]] if data_vars else None
        if da is not None:
            if "time" in da.dims:
                da = da.isel(time=0)
            result = da.load()
            ds.close()
            self._cache[key] = result
            while len(self._cache) > 4:
                self._cache.pop(next(iter(self._cache)))
            return result
        ds.close()
        return None

    def available_years(self, variable="acabf-anomaly"):
        r"""List available years for a variable."""
        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return []
        years = []
        for fn in sorted(os.listdir(vdir)):
            if fn.endswith(".nc"):
                try:
                    yr = int(fn.rstrip(".nc").split("_")[-1])
                    years.append(yr)
                except ValueError:
                    pass
        return years

    def get_field(self, variable, year, mesh_x, mesh_y):
        r"""Get a forcing field interpolated to mesh coordinates."""
        import xarray as xr

        yr = int(round(year))
        da = self._load_year(variable, yr)
        if da is None:
            return np.zeros(len(mesh_x))

        mx = xr.DataArray(np.asarray(mesh_x), dims="node")
        my = xr.DataArray(np.asarray(mesh_y), dims="node")

        xdim = [d for d in da.dims if d.lower() == "x"]
        ydim = [d for d in da.dims if d.lower() == "y"]
        if xdim and ydim:
            vals = da.interp({xdim[0]: mx, ydim[0]: my}, method="nearest")
        else:
            dims = [d for d in da.dims if d not in ("time",)]
            if len(dims) >= 2:
                vals = da.interp({dims[-1]: mx, dims[-2]: my}, method="nearest")
            else:
                return np.zeros(len(mesh_x))

        return np.nan_to_num(vals.values.flatten(), nan=0.0)

    def get_smb(self, year, mesh_x, mesh_y, anomaly=True):
        r"""Get SMB field in m/yr ice equivalent."""
        var = "acabf-anomaly" if anomaly else "acabf"
        raw = self.get_field(var, year, mesh_x, mesh_y)
        return smb_kgm2s_to_myr(raw)

    def get_smb_gradient(self, year, mesh_x, mesh_y):
        r"""Get SMB elevation gradient (dacabfdz) for ice-elevation feedback."""
        return self.get_field("dacabfdz", year, mesh_x, mesh_y)

    def get_temperature(self, year, mesh_x, mesh_y, anomaly=True):
        r"""Get surface temperature (K or K anomaly)."""
        var = "ts-anomaly" if anomaly else "ts"
        return self.get_field(var, year, mesh_x, mesh_y)


class ISMIP7Ocean:
    r"""Read ISMIP7 ocean forcing for Antarctica."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585",
                 version="v3"):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self.version = version
        self._ds_cache = {}
        self._interp_cache = {}

    def _var_dir(self, variable):
        return ocean_path(
            self.scenario, self.esm, variable,
            self.version, self.data_root,
        )

    def _load_variable(self, variable):
        import xarray as xr

        if variable in self._ds_cache:
            return self._ds_cache[variable]

        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return None

        nc_files = sorted(
            os.path.join(vdir, f) for f in os.listdir(vdir)
            if f.endswith(".nc")
        )
        if not nc_files:
            return None

        ds = xr.open_mfdataset(nc_files, combine="by_coords")
        self._ds_cache[variable] = ds
        return ds

    def _year_field(self, variable, year, nan_fill):
        r"""(RegularGridInterpolator, ascending z array) for one variable-year.

        Loads the single year slice out of its decadal chunk file into an
        in-memory nearest-neighbour interpolator over (z, y, x) — one year
        of tf is ~66 MB and answers in milliseconds, where the previous
        pointwise xarray .interp over the open_mfdataset dask graph cost
        ~10 minutes per step even at 5k mesh nodes. Nearest-neighbour also
        matches the OI-climatology CTRL path and the per-basin K
        calibration. The last few (variable, year) fields stay cached, so
        sub-yearly time steps re-read nothing.
        """
        import re
        import xarray as xr
        from scipy.interpolate import RegularGridInterpolator

        yr = int(round(year))
        key = (variable, yr)
        if key in self._interp_cache:
            return self._interp_cache[key]

        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return None

        # The chunk file whose YYYY-YYYY range contains the year (nearest
        # range for years outside coverage).
        best, best_d = None, None
        for f in sorted(os.listdir(vdir)):
            m = re.search(r"_(\d{4})-(\d{4})\.nc$", f)
            if not m:
                continue
            y0, y1 = int(m.group(1)), int(m.group(2))
            d = 0 if y0 <= yr <= y1 else min(abs(yr - y0), abs(yr - y1))
            if best_d is None or d < best_d:
                best, best_d = os.path.join(vdir, f), d
        if best is None:
            return None

        ds = xr.open_dataset(best)
        da = None
        for name in ds.data_vars:
            if name.lower() in (variable.lower(), "thermal_forcing",
                                "thermalforcing", "salinity"):
                da = ds[name]
                break
        if da is None:
            cands = [v for v in ds.data_vars
                     if not v.endswith("_bnds") and v.lower() != "crs"]
            da = ds[cands[0]]

        if "time" in da.dims:
            times = ds["time"].values
            if hasattr(times[0], "year"):  # cftime calendars
                idx = min(range(len(times)),
                          key=lambda i: abs(times[i].year - yr))
            else:
                years = ds["time"].dt.year.values
                idx = int(np.argmin(np.abs(years - yr)))
            da = da.isel(time=idx)

        zdim = [d for d in da.dims if d.lower() in ("z", "depth", "lev")][0]
        za = ds[zdim].values.astype(float)
        ya = ds["y"].values.astype(float)
        xa = ds["x"].values.astype(float)
        data = da.transpose(zdim, "y", "x").values.astype(np.float32)
        ds.close()
        if za[0] > za[-1]:
            za = za[::-1]; data = data[::-1, :, :]
        if ya[0] > ya[-1]:
            ya = ya[::-1]; data = data[:, ::-1, :]
        if xa[0] > xa[-1]:
            xa = xa[::-1]; data = data[:, :, ::-1]
        data = np.nan_to_num(data, nan=nan_fill)
        interp = RegularGridInterpolator(
            (za, ya, xa), data,
            method="nearest", bounds_error=False, fill_value=float(nan_fill),
        )
        entry = (interp, za)
        self._interp_cache[key] = entry
        while len(self._interp_cache) > 4:
            self._interp_cache.pop(next(iter(self._interp_cache)))
        return entry

    def _lookup(self, variable, year, mesh_x, mesh_y, draft, nan_fill):
        entry = self._year_field(variable, year, nan_fill)
        if entry is None:
            return None
        interp, za = entry
        mx = np.asarray(mesh_x)
        my = np.asarray(mesh_y)
        if draft is None:
            d = np.full(len(mx), za[-1])   # shallowest level
        else:
            # z and draft are both negative depths below sea level
            d = np.clip(np.asarray(draft), za[0], za[-1])
        return interp(np.column_stack([d, my, mx])).astype(float)

    def get_thermal_forcing(self, year, mesh_x, mesh_y, draft=None):
        r"""Thermal forcing at the shelf base on mesh points [K]."""
        vals = self._lookup("tf", year, mesh_x, mesh_y, draft, nan_fill=0.0)
        if vals is None:
            return np.zeros(len(mesh_x))
        return vals

    def get_salinity(self, year, mesh_x, mesh_y, draft=None, fill=34.5):
        r"""Ambient salinity at ice draft depth on mesh points [PSU]."""
        vals = self._lookup("so", year, mesh_x, mesh_y, draft, nan_fill=fill)
        if vals is None:
            return np.full(len(mesh_x), fill)
        return vals

    def close(self):
        for ds in self._ds_cache.values():
            ds.close()
        self._ds_cache.clear()
        self._interp_cache.clear()


class ISMIP7Fracture:
    r"""Read ISMIP7 fracture / ice shelf collapse forcing."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585"):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self._collapse_mask = None
        self._excess_melt = None

    def _fracture_dir(self):
        root = _find_ismip7_data(self.data_root)
        if root is None:
            return None
        return os.path.join(root, self.esm, self.scenario, "fracture")

    def load(self):
        import xarray as xr

        fdir = self._fracture_dir()
        if fdir is None or not os.path.isdir(fdir):
            return self

        # Masks may sit flat in fracture/ (legacy mirror) or inside a
        # versioned fracture/v<N>/ subdir (the share's layout); the highest
        # version wins when both are present.
        scan_dirs = [fdir] + [
            os.path.join(fdir, name) for _, name in _version_subdirs(fdir)
        ]
        found = {}
        for d in scan_dirs:
            for fn in sorted(os.listdir(d)):
                path = os.path.join(d, fn)
                if not os.path.isfile(path):
                    continue
                if "collapse_mask" in fn:
                    found["collapse_mask"] = path
                elif "excess_melt" in fn:
                    found["excess_melt"] = path
        if "collapse_mask" in found:
            self._collapse_mask = xr.open_dataset(found["collapse_mask"])
        if "excess_melt" in found:
            self._excess_melt = xr.open_dataset(found["excess_melt"])

        return self

    def get_collapse_mask(self, year, mesh_x, mesh_y):
        r"""Get ice shelf collapse mask (0/1) at given year."""
        import xarray as xr

        if self._collapse_mask is None:
            return np.zeros(len(mesh_x))

        ds = self._collapse_mask
        var = list(ds.data_vars)[0]
        da = ds[var]

        if "time" in da.dims:
            da = da.sel(time=year, method="nearest")

        mx = xr.DataArray(np.asarray(mesh_x), dims="node")
        my = xr.DataArray(np.asarray(mesh_y), dims="node")

        xdim = [d for d in da.dims if d.lower() == "x"]
        ydim = [d for d in da.dims if d.lower() == "y"]
        if xdim and ydim:
            vals = da.interp({xdim[0]: mx, ydim[0]: my}, method="nearest")
        else:
            return np.zeros(len(mesh_x))

        return np.nan_to_num(vals.values.flatten(), nan=0.0)

    def close(self):
        if self._collapse_mask is not None:
            self._collapse_mask.close()
        if self._excess_melt is not None:
            self._excess_melt.close()


def quadratic_mixed_slope(tf, salinity, sin_alpha, K=_K_DEFAULT):
    r"""ISMIP7 / Burgard et al. 2022 quadratic-mixed-slope local melt.

    Matches multimelt.melt_functions.quadratic_mixed_slope with TF_avg = TF
    (local-quadratic variant):

        m = K * melt_factor * U_factor * TF * |TF| * sin(alpha)

    where

        melt_factor = (rho_sw * c_po) / (rho_i * L_i)             [1/K]
        U_factor    = (c_po / L_i) * beta_S * g/(2|f|) * S0       [m/s/K]

    Inputs:
        tf        : thermal forcing T - T_f at ice base, K (numpy array)
        salinity  : ambient salinity at ice draft, PSU (numpy array)
        sin_alpha : sin of local ice-draft slope, dimensionless (numpy array)
        K         : dimensionless tuning factor (scalar or per-node array).
                    Burgard K50 = 1.15e-4.

    Returns melt rate in m/yr ice equivalent (positive = melting).
    """
    U_factor = (_C_PO / _L_I) * _BETA_S * (_G / (2.0 * abs(_F_CORIOLIS))) \
        * salinity
    melt = K * _MELT_FACTOR * U_factor * tf * np.abs(tf) * sin_alpha
    return melt * _SEC_PER_YEAR


def load_K_per_basin(npz_path, mesh_x, mesh_y, fill=0.0):
    r"""Load a per-basin calibrated K and return a per-node array.

    `npz_path` should be the output of `antarctica/scripts/calibrate_melt.py`
    and is expected to contain `basin_ids` (int) and `K_basin` (float) plus
    the IMBIE2 basin file path (the IMBIE2 8 km grid is re-read here so
    that the K-field can be remapped to *any* mesh, not just the one used
    during calibration).

    Returns an array of shape (len(mesh_x),) of per-node K values, with
    `fill` outside the calibrated basin set or where K_basin is NaN.
    """
    import xarray as xr
    from scipy.interpolate import RegularGridInterpolator

    data = np.load(npz_path)
    bids = np.asarray(data["basin_ids"]).astype(int)
    Kbas = np.asarray(data["K_basin"]).astype(float)

    root = os.environ.get(
        "ISMIP7_DATA_ROOT",
        os.path.join(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))), "ISMIP7", "AIS"),
    )
    # v2 (calibration-era) first; the v3 release from the reorganized share
    # carries an IDENTICAL basinNumber field (verified Jul 2026), so it is a
    # safe fallback for fresh clones that only ran the new downloader.
    candidates = [
        os.path.join(root, "parameterisations", "ocean", "imbie2",
                     "basin_numbers_ismip8km_v2.nc"),
        os.path.join(root, "obs", "ocean", "IMBIE-basins", "v3",
                     "IMBIE-basins_AIS_obs_ocean_v3.nc"),
    ]
    imbie2 = next((p for p in candidates if os.path.exists(p)), candidates[0])
    ds = xr.open_dataset(imbie2)
    xa = ds["x"].values; ya = ds["y"].values
    bn = ds["basinNumber"].values
    if ya[0] > ya[-1]:
        ya = ya[::-1]; bn = bn[::-1, :]
    if xa[0] > xa[-1]:
        xa = xa[::-1]; bn = bn[:, ::-1]
    interp = RegularGridInterpolator(
        (ya, xa), bn.astype(np.float32),
        method="nearest", bounds_error=False, fill_value=-1.0,
    )
    pts = np.column_stack([np.asarray(mesh_y), np.asarray(mesh_x)])
    basin_node = np.round(interp(pts)).astype(int)
    ds.close()

    K_field = np.full(len(mesh_x), fill, dtype=float)
    for bid, kb in zip(bids, Kbas):
        if np.isfinite(kb):
            K_field[basin_node == bid] = kb
    return K_field


def compute_sin_alpha(ctx):
    r"""Return sin(alpha) of local ice-draft slope at CG1 nodes.

    Computes draft = s - h, projects grad(draft) into the vector CG1
    space, and returns sin(arctan(|grad|)) = |grad|/sqrt(1 + |grad|^2).
    """
    import firedrake as fd
    Q = ctx["Q"]
    V = ctx["V"]
    h = ctx["h"]
    s = ctx["s"]
    draft = fd.Function(Q).interpolate(s - h)
    grad_draft = fd.project(fd.grad(draft), V)
    g = grad_draft.dat.data_ro
    gmag = np.sqrt(g[:, 0] ** 2 + g[:, 1] ** 2)
    return gmag / np.sqrt(1.0 + gmag * gmag)


def _oi_climatology_path(root, var, version):
    r"""Path of one OI-climatology variable for a given release.

    `30_sep` is the 2025-09-30 release under meltMIP/ (the one the
    per-basin K was calibrated against — the default for that reason);
    `06_nov` is the 2026 re-release (1972-2024) mirrored from the
    reorganized GHub share under obs/ocean/climatology/.
    """
    if version == "30_sep":
        return os.path.join(
            root, "meltMIP", f"OI_Climatology_ismip8km_60m_{var}_extrap.nc"
        )
    return os.path.join(
        root, "obs", "ocean", "climatology", f"zhou_annual_{version}",
        var, "v3",
        f"{var}_AIS_obs_ocean_climatology_zhou_annual_{version}_v3_1972-2024.nc",
    )


def build_oi_climatology_interpolators(data_root=None, version=None):
    r"""Load the OI climatology TF and so into (z, y, x)
    RegularGridInterpolators (nearest, fill 0). Shared by the CTRL and any
    observationally-forced run (OCX). ISMIP7_OI_VERSION selects the
    release (default 30_sep, matching the per-basin K calibration;
    switching to 06_nov without recalibrating K shifts the melt)."""
    import xarray as xr
    from scipy.interpolate import RegularGridInterpolator

    root = _find_ismip7_data(data_root)
    version = version or os.environ.get("ISMIP7_OI_VERSION", "30_sep")
    interps = {}
    for var in ("tf", "so"):
        path = _oi_climatology_path(root, var, version)
        ds = xr.open_dataset(path)
        if var not in ds.data_vars:
            found = [v for v in ds.data_vars
                     if not v.endswith("_bnds") and v.lower() != "crs"]
            ds.close()
            raise KeyError(
                f"{path} has no '{var}' variable (found {found}). Known "
                f"upstream packaging bug (Jul 2026): the 06_nov release "
                f"ships the tf field inside the so/thetao files. Use "
                f"ISMIP7_OI_VERSION=30_sep until it is fixed."
            )
        da = ds[var]
        zdim = [d for d in da.dims if d.lower() in ("z", "depth", "lev")][0]
        za = ds[zdim].values
        ya = ds["y"].values
        xa = ds["x"].values
        data = da.transpose(zdim, "y", "x").values.astype(np.float32)
        if za[0] > za[-1]:
            za = za[::-1]; data = data[::-1, :, :]
        if ya[0] > ya[-1]:
            ya = ya[::-1]; data = data[:, ::-1, :]
        if xa[0] > xa[-1]:
            xa = xa[::-1]; data = data[:, :, ::-1]
        data = np.nan_to_num(data, nan=0.0)
        interps[var] = (
            RegularGridInterpolator(
                (za, ya, xa), data,
                method="nearest", bounds_error=False, fill_value=0.0,
            ),
            za,
        )
        ds.close()
    return interps


def make_climatology_ocean_callback(K_field, data_root=None):
    r"""Ocean-melt callback with CONSTANT OI-climatology TF/so and evolving
    geometry: the CTRL2015 / observationally-constrained ocean forcing.
    K_field is a scalar or per-node array (calibrated per-basin K)."""
    interps = build_oi_climatology_interpolators(data_root)

    def callback(ctx, t_yr):
        mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
        mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]
        h = ctx["h"].dat.data_ro
        b = ctx["b"].dat.data_ro
        s = ctx["s"].dat.data_ro
        draft = np.minimum(s - h, 0.0)

        tf_interp, za_tf = interps["tf"]
        so_interp, za_so = interps["so"]
        d_tf = np.clip(draft, za_tf[0], za_tf[-1])
        d_so = np.clip(draft, za_so[0], za_so[-1])
        tf = tf_interp(np.column_stack([d_tf, mesh_y, mesh_x]))
        sal = so_interp(np.column_stack([d_so, mesh_y, mesh_x]))
        sin_a = compute_sin_alpha(ctx)

        melt = quadratic_mixed_slope(tf, sal, sin_a, K=K_field)

        haf = s - (b + (_RHO_WATER / _RHO_ICE) * np.maximum(-b, 0.0))
        floating = haf <= 0
        ctx["ocean_melt"].dat.data[:] = np.where(floating, melt, 0.0)

    return callback


def make_forcing_callback(atm=None, ocean=None, fracture=None,
                          K=_K_DEFAULT, K_per_basin_npz=None,
                          smb_anomaly=True, smb_baseline=None):
    r"""Build a forcing callback for use with simulation.run_simulation().

    Ocean melt uses the ISMIP7 Burgard quadratic_mixed_slope formula
    (local-quadratic variant, TF_avg = TF).

    K can be:
      - a scalar (single dimensionless K applied everywhere), or
      - a per-node numpy array (same length as mesh CG1 dofs), or
      - left at default while `K_per_basin_npz` points to the output of
        `calibrate_melt.py`; the per-basin K is then looked up on the mesh
        on first call and reused for subsequent steps.

    smb_anomaly selects between `acabf-anomaly` (True) and the full
    `acabf` field (False). The anomaly files are referenced to the ESM's
    1960-1989 climatology, so anomaly-only SMB is NOT a usable total:
    either pass `smb_baseline` (a per-node array added to the anomaly
    every step — e.g. RACMO climatology minus the anomaly's mean over
    the control reference window) or set smb_anomaly=False to force
    with the full field.

    ISMIP7_K_SCALE multiplies whichever K is in effect (the calibrated
    per-basin K integrates 689 vs 865 Gt/yr observed on the 2500 m mesh,
    so 1.26 matches the observed total).
    """
    K_field_cache = {"arr": None}
    K_scale = float(os.environ.get("ISMIP7_K_SCALE", "1.0"))

    def callback(ctx, t_yr):
        mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
        mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]

        if atm is not None:
            smb = atm.get_smb(t_yr, mesh_x, mesh_y, anomaly=smb_anomaly)
            if smb_baseline is not None:
                smb = smb + smb_baseline
            ctx["accum"].dat.data[:] = smb

        if ocean is not None and "ocean_melt" in ctx:
            h = ctx["h"].dat.data_ro
            b = ctx["b"].dat.data_ro
            s = ctx["s"].dat.data_ro
            # Ice shelf draft (negative depth below sea level)
            draft = np.minimum(s - h, 0.0)

            tf = ocean.get_thermal_forcing(t_yr, mesh_x, mesh_y, draft=draft)
            sal = ocean.get_salinity(t_yr, mesh_x, mesh_y, draft=draft)
            sin_alpha = compute_sin_alpha(ctx)

            # Resolve K: per-basin npz takes precedence if supplied.
            if K_per_basin_npz is not None:
                if K_field_cache["arr"] is None:
                    K_field_cache["arr"] = load_K_per_basin(
                        K_per_basin_npz, mesh_x, mesh_y, fill=0.0
                    )
                K_use = K_field_cache["arr"]
            else:
                K_use = K

            melt = quadratic_mixed_slope(tf, sal, sin_alpha, K=K_use * K_scale)

            # Only apply melt where ice is floating (haf <= 0)
            haf = s - (b + (_RHO_WATER / _RHO_ICE) * np.maximum(-b, 0.0))
            floating = haf <= 0
            ctx["ocean_melt"].dat.data[:] = np.where(floating, melt, 0.0)

    return callback
