r"""ISMIP7 forcing data reader for Antarctic simulations."""

import os
import re
import warnings
import numpy as np

_SEC_PER_YEAR = 31556926.0
_RHO_ICE = 917.0
# Seawater for flotation, the value the forward builds its surface and its
# flotation surface with (simulation.py, rho_ratio = 917 / 1024). The melt
# callbacks used fresh water, 1000 kg/m^3, here until September 2026, which
# put the flotation surface too low and read every floating cell thicker than
# 78 percent of its flotation thickness as grounded, withholding its melt:
# 364 000 km2, 24 percent of BedMachine's shelf area, on the 2500 m mesh (the
# seawater test grounds 0.3 percent of it).
_RHO_SW_FLOTATION = 1024.0

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

# K50 of the ISMIP7 toolbox's standard sampling (parameter_selection_quadratic
# _example.ipynb, July 2026 update: K05 4.75e-5, K50 8.5e-5, K95 1.375e-4),
# sampled with the constant slope SIN_ALPHA_ANT_DEFAULT below. The earlier
# 1.15e-4 was the pre-update value. It is the notebook's reference value and
# the default argument of the melt law; a run takes its K from the melt
# calibration file (runconfig.deltat_per_basin_npz), 6.5e-5 in the tracked one.
_K_DEFAULT = 8.5e-5
_K_PERCENTILES = (4.75e-5, 8.5e-5, 1.375e-4)

_DEFAULT_DATA_ROOT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "ISMIP7", "AIS",
)


def _comm_rank():
    r"""World rank, or 0 without mpi4py: notices meant to be logged once must
    not repeat per rank."""
    try:
        from mpi4py import MPI
    except ImportError:
        return 0
    return MPI.COMM_WORLD.rank


def forcing_year(t_yr):
    r"""The calendar year a MODEL TIME lies in, from the time a step ENDS at.

    ``run_simulation`` hands the forcing callback the END of the step, and
    ``t = Y.0`` is 1 January of year Y, so the step from 2300.9 to 2301.0 lies
    in 2300 and the step from 2015.0 to 2015.1 lies in 2015. Rounding instead
    asks for the wrong year on every step past the half-year mark: a run that
    covers 2015-2300 would end by requesting 2301, two past CESM2-WACCM's
    2299, which the one-year end-of-series bridge cannot cover, so the run
    dies on its last step instead of banking 2300.

    This is also the year ``AnnualOutput`` is accumulating, so the forcing a
    step receives and the year its fluxes are booked into are the same.

    It belongs at the boundary where model time is known, which is the forcing
    callback. The readers below take a plain CALENDAR year: they are also
    called with years straight out of ``available_years`` (the climatology
    pools in ``smb_scheme`` and ``compute_climatology``), and shifting those
    would re-reference every pooled year by one and make the first year of a
    scenario read as preceding its own series.
    """
    import math
    return int(math.ceil(float(t_yr) - 1e-9)) - 1


def smb_kgm2s_to_myr(smb_kgm2s):
    r"""Convert SMB from kg/m^2/s to m/yr of ice.

    ``acabf`` and ``acabf-anomaly`` are mass fluxes, so the thickness rate the
    transport needs is the flux over the ice density; no water density enters.
    ``write_ismip7_output`` converts the model's ``acabf`` back with the same
    ice density. Until September 2026 this carried a further
    ``rho_water / rho_ice`` factor, which read every ESM SMB anomaly 9 percent
    too large (issue #29).
    """
    return smb_kgm2s * _SEC_PER_YEAR / _RHO_ICE


def _sample_raster(raster, Q):
    r"""Put an open raster onto ``Q``, as a CELL AVERAGE when ``Q`` is DG0.

    A DG0 dof sits at the cell centroid, so ``icepack.interpolate`` would take
    a one-point sample of the raster per cell. SMB sets the mass budget and the
    a_ref balance, so it goes through the same cell-averaging rule the geometry
    uses (see geometry.sample_to_geometry). CG1 is the nodal interpolant, as
    before.
    """
    import firedrake as fd
    import icepack

    if Q.ufl_element().degree() > 0:
        return icepack.interpolate(raster, Q)
    from .geometry import sample_to_geometry
    Q_cg = fd.FunctionSpace(Q.mesh(), "CG", 1)
    return sample_to_geometry(raster, Q, Q_cg)


def load_racmo_smb_climatology(Q, clim_start=2000, clim_end=2029, data_dir=None,
                               target_res=8000.0, rho_ice=_RHO_ICE):
    r"""RACMO2.4p1 mean-annual SMB (m/yr ice equiv) as a Function on Q's mesh.

    The RACMO ANT11 grid is rotated-pole, so this reprojects the climatology to
    an intermediate EPSG:3031 raster with rasterio, then samples it onto the mesh
    -- the same path used for BedMachine -- avoiding any scattered-point
    interpolation. On a DG0 ``Q`` the sample is a cell average rather than a
    centroid point sample (see :func:`_sample_raster`). ``smbgl`` is a monthly
    mass sum (kg/m^2), so the annual SMB is the sum of the 12 months, averaged
    over the climatology window.
    """
    import xarray as xr
    import pyproj
    from affine import Affine
    from rasterio.io import MemoryFile
    from rasterio.warp import reproject, Resampling

    if data_dir is None:
        from .runconfig import obs_data_root
        data_dir = os.path.join(obs_data_root(), "racmo")
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
            return _sample_raster(raster, Q)


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
    if var in ds:
        da = ds[var]
    else:
        # Skip grid-mapping/bounds variables, as _load_year does: any of them
        # ordered first would silently become the temperature field.
        cands = [
            v for v in ds.data_vars
            if not v.endswith("_bnds")
            and v.lower() not in ("x", "y", "time", "month", "crs", "mapping",
                                  "spatial_ref", "lat", "lon")
        ]
        if not cands:
            raise ValueError(f"no usable data variable in {hits[0]}")
        da = ds[cands[0]]
    da = da.mean(dim="month")
    if "y" in da.dims and "x" in da.dims:
        da = da.transpose("y", "x")
    T = np.asarray(da.values, dtype=float)                       # (y, x)
    x = np.asarray(ds["x"].values, dtype=float)
    y = np.asarray(ds["y"].values, dtype=float)
    # RegularGridInterpolator requires ascending axes (same flip as
    # _year_field / build_oi_climatology_interpolators / load_K_per_basin).
    if y[0] > y[-1]:
        y = y[::-1]
        T = T[::-1, :]
    if x[0] > x[-1]:
        x = x[::-1]
        T = T[:, ::-1]
    T = np.where(T > 100.0, T, np.nan)                           # mask fill (~0 K)
    T_filled = np.where(np.isfinite(T), T, fill_K)
    interp = RegularGridInterpolator((y, x), T_filled, method="nearest",
                                     bounds_error=False, fill_value=fill_K)
    Tf = fd.Function(Q, name="T_srf")
    xy = Q.mesh().coordinates.dat.data_ro
    Tf.dat.data[:] = interp(np.column_stack([xy[:, 1], xy[:, 0]]))
    return Tf


def version_key(name):
    r"""Sort key of a ``v<N>`` or ``v<N>.<M>`` directory name, or None if
    ``name`` is not one. Dotted versions are real (the fracture forcing went
    v2 -> v2.1 in September 2026), so the key is a tuple of integers, not
    ``int(name[1:])``, which raised on them and hid the directory."""
    if not re.fullmatch(r"v\d+(\.\d+)*", name):
        return None
    return tuple(int(x) for x in name[1:].split("."))


def _version_subdirs(parent_dir):
    r"""``[(key, name)]`` of the ``v<N>`` and ``v<N>.<M>`` subdirectories of
    ``parent_dir``, ascending."""
    versions = []
    if not os.path.isdir(parent_dir):
        return versions
    for name in os.listdir(parent_dir):
        if version_key(name) is not None and os.path.isdir(os.path.join(parent_dir, name)):
            versions.append((version_key(name), name))
    versions.sort()
    return versions


def _resolve_version(parent_dir, pinned):
    r"""The pinned version if that subdir exists, else the highest v<N>
    subdir present, else the pinned name unchanged (so missing trees keep
    producing the same non-existent path the callers already handle)."""
    if parent_dir is None or os.path.isdir(os.path.join(parent_dir, pinned)):
        return pinned
    versions = _version_subdirs(parent_dir)
    return versions[-1][1] if versions else pinned


def _newer_versions(parent_dir, version):
    r"""Version subdirs of ``parent_dir`` above ``version``: the reader stays
    on its pin so that a campaign does not change forcing under itself, and
    the run has to say so when something newer has landed beside it."""
    key = version_key(version)
    if key is None:
        return []
    return [name for k, name in _version_subdirs(parent_dir) if k > key]


# core_report.py lifts every line carrying this marker out of the run log, the
# way it lifts the climatology pool, so the committed report states which
# forcing a run opened. The focus groups ask for exactly that in the
# submission README (discussion #37), and an audit date is not it.
FORCING_PROVENANCE_MARKER = "Forcing provenance:"


def describe_forcing_provenance(*readers, variables=None):
    r"""One marker line per forcing variable the given readers resolve, plus
    one for every reader that resolves nothing, so an absent tree is a
    statement in the log and not a silence. ``variables`` narrows a reader
    kind to what the caller reads, ``{"atmosphere": ("acabf",)}``."""
    lines = []
    for reader in readers:
        if reader is None:
            continue
        wanted = (variables or {}).get(reader.kind)
        rows = reader.provenance(wanted) if wanted else reader.provenance()
        if not rows:
            lines.append(f"{FORCING_PROVENANCE_MARKER} {reader.kind} "
                         f"{reader.esm} {reader.scenario}: nothing on disk")
        for row in rows:
            newer = (f"  NEWER ON DISK, NOT READ: {' '.join(row['newer'])}"
                     if row["newer"] else "")
            lines.append(f"{FORCING_PROVENANCE_MARKER} {reader.kind} {row['variable']} "
                         f"{reader.esm} {reader.scenario} {row['product']} "
                         f"{row['version']}{newer}")
    return lines


def describe_observational_forcing(smb=None, ocean=False):
    r"""Marker lines for forcing that is not an ISMIP7 scenario tree: the
    RACMO SMB a control or the OCX stopgap runs on (``smb`` says which, in
    words) and the OI ocean climatology, whose release is
    ``ISMIP7_OI_VERSION``."""
    lines = []
    if smb:
        lines.append(f"{FORCING_PROVENANCE_MARKER} atmosphere {smb}")
    if ocean:
        release = os.environ.get("ISMIP7_OI_VERSION", "30_sep")
        lines.append(f"{FORCING_PROVENANCE_MARKER} ocean OI climatology tf+so, "
                     f"release {release}, constant in time")
    return lines


# The versions the readers ask for first. audit_forcing_versions.py reads
# these, so that "current" means what a run would open and not merely what is
# somewhere on disk.
ATMOSPHERE_VERSION = "v2"
OCEAN_VERSION = "v3"

# The oldest collapse-mask version the submission accepts, per ESM and
# scenario. The fracture reader takes the highest version on disk, so a tree
# synced from the mirror alone opens a mask the focus group has replaced when
# the fix is on Globus and not yet on the mirror: the MRI-ESM2-0 ssp585 v1
# flagged 593 cells of 8 km by 2250, the v2 of 22 September 2026 flags 25,149
# (discussion #30). A run under a mask mode refuses anything older
# (ISMIP7Fracture.check_min_version), and audit_forcing_versions.py fails on
# it, whatever the mirror publishes.
FRACTURE_MIN_VERSION = {
    ("CESM2-WACCM", "ssp126"): "v2.1",
    ("CESM2-WACCM", "ssp370"): "v2.1",
    ("CESM2-WACCM", "ssp585"): "v2.1",
    ("MRI-ESM2-0", "ssp585"): "v2",
}


# The observation-constrained experiment. It has no ESM and no scenario, and
# the focus groups decided it need not follow the ESM ordering, only be
# consistent with itself (discussion #41, item 6). So its tree is scenario
# first, ``OCX/<source>/<product>/<variable>/<version>/``, where every other
# tree is ``<ESM>/<scenario>/...``, while its FILENAMES keep source before
# OCX, which is the ``<esm>_<scenario>`` order the readers already build.
# Reading it as ``esm=OCX_ATMOSPHERE_SOURCE, scenario=OCX`` therefore needs
# the directory swapped and nothing else.
OCX = "OCX"
OCX_ATMOSPHERE_SOURCE = "RACMO2.3p2-ERA"
# The Antarctic OCX ocean cites no source: it is four expert-judgment
# scenarios (discussion #41), of which ``main`` is the core one.
OCX_OCEAN_VARIANTS = ("main", "cold", "warm", "vary")
OCX_OCEAN_SOURCE = "expert-judgment"


def _scenario_dir(root, esm, scenario):
    r"""``<root>/<esm>/<scenario>``, or ``<root>/OCX/<source>`` for OCX."""
    if scenario == OCX:
        return os.path.join(root, OCX, esm)
    return os.path.join(root, esm, scenario)


def atmosphere_path(scenario, esm="CESM2-WACCM", variable="acabf-anomaly",
                    resolution="8000m", version=ATMOSPHERE_VERSION, data_root=None):
    root = _find_ismip7_data(data_root)
    if root is None:
        return None
    parent = os.path.join(_scenario_dir(root, esm, scenario),
                          atmosphere_product(root, esm, scenario, resolution), variable)
    return os.path.join(parent, _resolve_version(parent, version))


ATMOSPHERE_PRODUCTS = ("SDBN1", "GEMB-SDBN1")


def atmosphere_product(root, esm, scenario, resolution="8000m"):
    r"""The downscaled-atmosphere directory name for this ESM and scenario.

    The core experiment uses ``SDBN1`` for CESM2-WACCM and ``GEMB-SDBN1`` for
    MRI-ESM2-0 (MRI's runoff needed an energy-balance step before the
    statistical downscaling; the directories were renamed in August 2026,
    discussion #37, data unchanged). Whichever exists on disk wins, ``SDBN1``
    first, so a tree fetched before the rename keeps working.
    """
    for product in ATMOSPHERE_PRODUCTS:
        if os.path.isdir(os.path.join(_scenario_dir(root, esm, scenario), f"{product}-{resolution}")):
            return f"{product}-{resolution}"
    return f"SDBN1-{resolution}"


def ocean_path(scenario, esm="CESM2-WACCM", variable="tf",
               version=OCEAN_VERSION, data_root=None, variant="main"):
    r"""``<root>/<esm>/<scenario>/ocean/<variable>/<version>``. The OCX ocean
    is ``<root>/OCX/ocean/<variant>/<version>`` instead: no source, one
    directory per expert-judgment scenario, and ``so``, ``tf`` and ``thetao``
    side by side in it with no directory of their own."""
    root = _find_ismip7_data(data_root)
    if root is None:
        return None
    if scenario == OCX:
        parent = os.path.join(root, OCX, "ocean", variant)
    else:
        parent = os.path.join(root, esm, scenario, "ocean", variable)
    return os.path.join(parent, _resolve_version(parent, version))


def _time_axis_days(values):
    r"""A time (or time-bounds) array as floating-point days.

    Only differences matter to the callers, so the origin is arbitrary and no
    calendar arithmetic is needed. Handles the three things xarray can hand
    back for a decoded CF time axis: plain numerics (undecoded), ``datetime64``
    (standard/proleptic_gregorian calendars) and object arrays of ``cftime``
    datetimes (noleap, 360_day, all_leap, ...), which have no ``__float__``
    and so cannot go through ``asarray(..., dtype=float)``.
    """
    arr = np.asarray(values)
    if arr.dtype.kind in "fiub":
        return arr.astype(float)
    if arr.dtype.kind in "Mm":
        unit = "datetime64[s]" if arr.dtype.kind == "M" else "timedelta64[s]"
        # NaT casts to a huge finite negative, which would sail through the
        # caller's finite/positive checks as a plausible weight; make it NaN.
        return np.where(
            np.isnat(arr), np.nan, arr.astype(unit).astype("float64") / 86400.0
        )
    flat = arr.reshape(-1)
    origin = flat[0]
    days = np.array(
        [(t - origin).total_seconds() for t in flat], dtype=float
    ) / 86400.0
    return days.reshape(arr.shape)


_WEIGHT_FALLBACK_WARNED = set()

# Real month lengths span 28-31 days (ratio 1.11); 360_day is uniform. Anything
# spanning more than 3x, or handing one month over half the annual weight, is a
# corrupt weight vector, not a calendar - and would quietly reinstate the
# single-month forcing this whole helper exists to eliminate.
_MAX_WEIGHT_RATIO = 3.0
_MAX_WEIGHT_SHARE = 0.5


def _warn_unweighted(da, ds, reason):
    r"""Warn once per (file, variable, reason) that month-length weighting was
    unavailable, so the degradation is visible in run logs instead of silent."""
    source = ds.encoding.get("source") if ds is not None else None
    source = source or "<unknown file>"
    name = getattr(da, "name", None) or "<unnamed variable>"
    key = (source, name, reason)
    if key in _WEIGHT_FALLBACK_WARNED:
        return
    _WEIGHT_FALLBACK_WARNED.add(key)
    warnings.warn(
        f"ISMIP7 forcing: {reason}; using the unweighted 12-month mean for "
        f"'{name}' in {source} (<=1% off a day-weighted annual mean)",
        RuntimeWarning,
        stacklevel=3,
    )


def _refuse_empty_time_axis(n, ds, what):
    r"""A forcing file with no time slices is a failed upload, not an empty
    year. The 2300 ``dacabfdz``/``dtsdz``/``dmrrodz`` files on the share were
    exactly that until 27 May 2026 (discussion #8), and a tree fetched before
    then still holds them. Left alone it surfaces as numpy's "zero-size array
    to reduction operation", which names neither the file nor the cause."""
    if n == 0:
        source = (ds.encoding.get("source") if ds is not None else None) or "<unknown file>"
        raise ValueError(
            f"ISMIP7 forcing: {what} in {source} has an empty time axis (0 "
            f"slices). The share held empty 2300 files until 2026-05-27 "
            f"(discussion #8); delete it and fetch it again."
        )


def _open_forcing(path):
    r"""``xarray.open_dataset`` for an ISMIP7 forcing file, its time axis
    decoded to ``cftime`` datetimes whatever the calendar or the year.

    The readers take only the year of a slice (:func:`_nearest_year_index`)
    and the day spacing of a monthly axis (:func:`_time_axis_days`), and both
    accept ``cftime`` objects. Left to its default, xarray decodes a
    standard-calendar axis to ``datetime64[ns]`` while the dates fit and falls
    back to ``cftime`` past 2262, with a ``SerializationWarning`` on every
    open: a 2015-2300 projection printed one per chunk it read past 2262, and
    the unit suite carried two. Asking for ``cftime`` up front makes the
    decode the same for every file, and silent. An axis with no CF units
    (plain years) is left numeric either way.
    """
    import xarray as xr
    try:
        coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    except AttributeError:          # xarray before 2025.01
        return xr.open_dataset(path, use_cftime=True)
    return xr.open_dataset(path, decode_times=coder)


def _nearest_year_index(time, year, ds=None, what="time"):
    r"""Index of the slice of ``time`` (a DataArray) nearest to ``year``.

    The axis may be decoded to ``datetime64`` (standard calendars), to
    ``cftime`` objects (noleap, or dates past 2262) or, where a file carries
    plain years and no CF units, left numeric. Only the year is ever read, so
    the calendars the ISMIP7 products disagree on (discussions #9 and #24)
    cannot shift a slice.
    """
    values = time.values
    _refuse_empty_time_axis(len(values), ds, what)
    if np.issubdtype(values.dtype, np.number):
        years = values
    elif hasattr(values[0], "year"):
        years = [t.year for t in values]
    else:
        years = time.dt.year.values
    return int(np.argmin(np.abs(np.asarray(years, dtype=float) - int(year))))


def _annual_mean_over_time(da, ds=None):
    r"""Collapse a per-year forcing file's time axis to the ANNUAL MEAN.

    The ISMIP7 SDBN1 atmosphere files carry 12 MONTHLY slices per year
    (``time`` = days since <year>-01-15, values 0, 31, 60, ...), so a single
    slice is one month, not the year. Taking ``isel(time=0)`` grabs JANUARY -
    peak austral summer, the maximum-ablation month - and applies it as the
    whole year's forcing. For MRI-ESM2-0 ssp585 2108 that is -16583 Gt/yr
    against an annual mean of -1276 Gt/yr: a 13x overestimate of ablation
    that grows with warming (the summer melt trend is far steeper than the
    annual one), which drove wildly negative post-2100 SMB, unphysical
    +/-6000 Gt/yr year-to-year swings, and a sea-level contribution above the
    ISMIP6 envelope.

    Months are weighted by their length (from ``time_bnds`` when present, else
    from the spacing of the time coordinate) so the result is a true annual
    mean rather than a 12-month unweighted average. A length-1 time axis
    collapses to that single value, so annual files are unaffected. Any
    failure to derive usable weights - an exotic calendar, a malformed bounds
    variable, an unmasked fill value - degrades to the unweighted mean (<=1%
    off) with a one-time warning rather than raising: a forcing read must not
    abort on a calendar variant, and must never let a degenerate weight vector
    concentrate the year onto one month.
    """
    import xarray as xr

    n = da.sizes["time"]
    _refuse_empty_time_axis(n, ds, f"'{da.name}'")
    if n == 1:
        return da.isel(time=0)

    w = None
    reason = None
    try:
        bnds_name = da.attrs.get("bounds") or (
            ds["time"].attrs.get("bounds")
            if ds is not None and "time" in ds else None
        )
        if ds is not None and bnds_name and bnds_name in ds:
            b = _time_axis_days(ds[bnds_name].values)
            if b.ndim == 2 and b.shape[0] == n:
                w = b[:, 1] - b[:, 0]
        if w is None and ds is not None and "time" in ds:
            t = _time_axis_days(ds["time"].values)
            if t.size == n:
                # month length = spacing to the next slice; the last month
                # reuses the previous spacing (the file ends at the year
                # boundary).
                d = np.diff(t)
                w = np.concatenate([d, d[-1:]])
    except Exception as exc:
        w = None
        reason = (f"month-length weights could not be read "
                  f"({type(exc).__name__}: {exc})")

    if w is None:
        reason = reason or "no time bounds and no usable time coordinate"
    else:
        w = np.asarray(w, dtype=float)
        if not np.all(np.isfinite(w)) or not np.all(w > 0):
            reason = "month-length weights are non-finite or non-positive"
        elif w.max() > _MAX_WEIGHT_RATIO * w.min():
            reason = (f"month lengths span {w.min():.4g}-{w.max():.4g} days, "
                      f"beyond any real calendar")
        elif w.max() / w.sum() > _MAX_WEIGHT_SHARE:
            reason = (f"one slice carries {100 * w.max() / w.sum():.1f}% of "
                      f"the annual weight")
        if reason is not None:
            w = None

    if w is None:
        _warn_unweighted(da, ds, reason)
        return da.mean("time")   # equal weights: <=1% off a day-weighted mean

    # weighted().mean() renormalizes by the weights of the non-NaN months, so
    # a partially masked node matches the unweighted mean instead of being
    # biased low, and an all-NaN node stays NaN rather than collapsing to 0.
    return da.weighted(xr.DataArray(w, dims="time")).mean("time")


class ISMIP7Atmosphere:
    r"""Read ISMIP7 downscaled atmosphere forcing for Antarctica."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585",
                 resolution="8000m", version=ATMOSPHERE_VERSION):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self.resolution = resolution
        self.version = version
        self._cache = {}
        self._persisted = set()          # variables already reported as persisted past the series end
        self._grid_x = None
        self._grid_y = None

    kind = "atmosphere"

    def _var_dir(self, variable):
        return atmosphere_path(
            self.scenario, self.esm, variable,
            self.resolution, self.version, self.data_root,
        )

    def provenance(self, variables=("acabf-anomaly", "acabf")):
        r"""``[{variable, product, version, dir, newer}]`` for the variables
        that are on disk: what ``_load_year`` would open, by its own rules."""
        rows = []
        for variable in variables:
            vdir = self._var_dir(variable)
            if vdir is None or not os.path.isdir(vdir):
                continue
            rows.append({
                "variable": variable,
                "product": os.path.basename(os.path.dirname(os.path.dirname(vdir))),
                "version": os.path.basename(vdir), "dir": vdir,
                "newer": _newer_versions(os.path.dirname(vdir), os.path.basename(vdir)),
            })
        return rows

    def _years_on_disk(self, vdir, variable, product, version):
        r"""Sorted years for which ``vdir`` holds a file ``_load_year`` would
        open: the full name is matched, so a file filed under the wrong
        scenario (the ssp585 anomalies named ``historical``, discussion #41)
        is not a year of this series."""
        head = f"{variable}_AIS_{self.esm}_{self.scenario}_{product}_{version}_"
        return sorted(int(m.group(1)) for f in os.listdir(vdir)
                      for m in [re.fullmatch(re.escape(head) + r"(\d{4})\.nc", f)] if m)

    def _year_span(self, vdir, variable, product, version):
        r"""``(first, last)`` year for which ``vdir`` holds a file, or None."""
        years = self._years_on_disk(vdir, variable, product, version)
        return (years[0], years[-1]) if years else None

    def _load_year(self, variable, year):
        key = (variable, int(year))
        if key in self._cache:
            return self._cache[key]

        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return None

        version = os.path.basename(vdir)
        product = os.path.basename(os.path.dirname(os.path.dirname(vdir)))   # SDBN1-8000m or GEMB-SDBN1-8000m
        pattern = f"{variable}_AIS_{self.esm}_{self.scenario}_{product}_{version}_{int(year)}.nc"
        path = os.path.join(vdir, pattern)

        if not os.path.exists(path):
            # End of the series: a tree whose last atmosphere year is 2299
            # while a 2015-2300 run needs 2300. Bridge exactly that one year,
            # once per variable in the log, rather than failing at the last
            # step. CESM2-WACCM was that case when the empty 2300 files were
            # removed (discussion #8); as of 22 September 2026 the mirror
            # carries 2300 again, padded by the atmosphere group with the
            # 2290-2299 mean, so those files are read as given and this branch
            # no longer fires for them. It still covers an older mirror copy,
            # and the ocean has no 2300 at all. On 23 September 2026 the
            # organisers answered that a single forcing year at the end makes
            # no significant difference and that a 2299 duplicate is
            # acceptable, so a run still ends in 2300 (discussion #49,
            # icepack/ismip7#78).
            # Anything further past the end is a short tree, not the end of
            # the series, and repeating one year of SMB for decades would be a
            # scientifically wrong run reported as a success, so it raises.
            span = self._year_span(vdir, variable, product, version)
            if span is None:
                # The variable has no files at all here: an optional product
                # (dacabfdz, ts-anomaly) the callers may legitimately run
                # without. Only a year missing from a series that exists is
                # an error.
                return None
            first, last = span
            if int(year) - last == 1:
                if variable not in self._persisted:
                    self._persisted.add(variable)
                    if _comm_rank() == 0:
                        print(f"  ISMIP7Atmosphere: {variable} has no year {int(year)}; "
                              f"persisting {last}, the last year on disk", flush=True)
                self._cache[key] = self._load_year(variable, last)
                return self._cache[key]
            # get_field would turn a None into a field of zeros and the run
            # would report a whole year of zero anomaly as a success, so the
            # reader refuses instead. A year BEFORE the series is its own
            # case: a projection asking for one means its timeline starts
            # earlier than the scenario does, not that the tree is short.
            where = ("precedes the series there" if int(year) < first
                     else "is missing from the series there")
            raise FileNotFoundError(
                f"ISMIP7Atmosphere: {variable} for {self.esm} {self.scenario} "
                f"has no year {int(year)} in {vdir}: it {where} "
                f"({first}-{last}; only the single year after the end is "
                f"bridged, 2300 after 2299)."
            )

        ds = _open_forcing(path)

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
                da = _annual_mean_over_time(da, ds)
            result = da.load()
            ds.close()
            self._cache[key] = result
            while len(self._cache) > 4:
                self._cache.pop(next(iter(self._cache)))
            return result
        ds.close()
        return None

    def available_years(self, variable="acabf-anomaly"):
        r"""Years ``_load_year`` can open for a variable. The same strict
        match as the loader: a looser one let a misnamed file pass the
        availability gate, after which ``get_field`` found nothing under the
        name it builds and returned a year of zeros."""
        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return []
        version = os.path.basename(vdir)
        product = os.path.basename(os.path.dirname(os.path.dirname(vdir)))
        return self._years_on_disk(vdir, variable, product, version)

    def get_field(self, variable, year, mesh_x, mesh_y):
        r"""Get a forcing field interpolated to mesh coordinates.

        The value is that year's ANNUAL MEAN, not a single slice: the SDBN1
        files are monthly, so ``_load_year`` collapses the time axis via
        ``_annual_mean_over_time`` (see that docstring for the weighting)."""
        import xarray as xr

        yr = int(year)
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
        r"""Get SMB elevation gradient (dacabfdz) for ice-elevation feedback.

        Nothing calls this: the model has no SMB-height feedback, and the
        submission README says so. Before anything does, three things from the
        forum. The protocol prefers the RUNOFF gradient ``dmrrodz``, since
        ``dacabfdz`` is dominated in places by precipitation patterns that
        have nothing to do with elevation; either is accepted if the README
        names it (discussion #36). Runoff is counted positive for mass LOSS,
        so the SMB correction is MINUS ``dmrrodz`` times the elevation change
        (#35), which a group found out from its results. And the AIS OCX
        ``dacabfdz`` was spatially shifted until it was replaced in place
        around 8 September 2026 (#45), so a copy fetched before then is wrong
        under the right name: ``download_mirror.py`` will say REPLACED.
        """
        return self.get_field("dacabfdz", year, mesh_x, mesh_y)

    def get_temperature(self, year, mesh_x, mesh_y, anomaly=True):
        r"""Get surface temperature (K or K anomaly)."""
        var = "ts-anomaly" if anomaly else "ts"
        return self.get_field(var, year, mesh_x, mesh_y)


class ISMIP7Ocean:
    r"""Read ISMIP7 ocean forcing for Antarctica."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585",
                 version=OCEAN_VERSION, variant="main"):
        self.data_root = _find_ismip7_data(data_root)
        if scenario == OCX:
            if variant not in OCX_OCEAN_VARIANTS:
                raise ValueError(f"the OCX ocean is one of {OCX_OCEAN_VARIANTS}, got {variant!r}")
            esm = OCX_OCEAN_SOURCE
        self.esm = esm
        self.scenario = scenario
        self.variant = variant
        self.version = version
        self._ds_cache = {}
        self._interp_cache = {}
        self._persisted = set()          # variables already reported as held past the series end

    kind = "ocean"

    def _var_dir(self, variable):
        return ocean_path(
            self.scenario, self.esm, variable,
            self.version, self.data_root, self.variant,
        )

    def provenance(self, variables=("tf", "so")):
        r"""``[{variable, product, version, dir, newer}]`` for the variables
        that are on disk: what ``_year_field`` would open."""
        rows = []
        for variable in variables:
            vdir = self._var_dir(variable)
            if vdir is None or not os.path.isdir(vdir) or not self.spans(variable):
                continue
            rows.append({
                "variable": variable,
                "product": f"ocean/{self.variant}" if self.scenario == OCX else "ocean",
                "version": os.path.basename(vdir), "dir": vdir,
                "newer": _newer_versions(os.path.dirname(vdir), os.path.basename(vdir)),
            })
        return rows

    def spans(self, variable="tf"):
        r"""Sorted ``[(first, last, path)]`` of the chunk files on disk."""
        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return []
        # by the variable's own name: the OCX ocean keeps so, tf and thetao
        # in one directory, and thetao's chunks are not tf's
        return sorted((int(m.group(1)), int(m.group(2)), os.path.join(vdir, f))
                      for f in os.listdir(vdir) if f.startswith(variable + "_")
                      for m in [re.search(r"_(\d{4})-(\d{4})\.nc$", f)] if m)

    def coverage(self, variable="tf"):
        r"""``(first, last)`` forcing year on disk, or None. The run gate and
        preflight both read this, so they agree with what ``_year_field``
        will serve."""
        spans = self.spans(variable)
        return (spans[0][0], max(sp[1] for sp in spans)) if spans else None

    def require_years(self, first, last, variables=("tf", "so")):
        r"""Refuse a run that needs years ``first`` to ``last`` from a tree
        that does not hold them, and return the ``(first, last)`` span every
        variable covers.

        A variable with no files at all reads as zero thermal forcing
        (``get_thermal_forcing``), so a run on an absent tree would melt
        nothing and report success. A driver calls this before its model
        setup. The single year after the series is held (``_chunk_for``), so
        ``last`` may lie one past the files.
        """
        covers = []
        for var in variables:
            cover = self.coverage(var)
            if cover is None:
                raise FileNotFoundError(
                    f"No ocean {var} data for {self.esm}/{self.scenario}. Download the "
                    f"ocean tree first."
                )
            if cover[0] > first or cover[1] + 1 < last:
                raise FileNotFoundError(
                    f"Ocean {var} for {self.esm}/{self.scenario} covers {cover[0]}-{cover[1]}, "
                    f"and this run needs {first}-{last} (one year past the end is "
                    f"held, no more). Download the rest of the ocean tree, or move "
                    f"ISMIP7_T_START / ISMIP7_T_END inside it."
                )
            covers.append(cover)
        return max(c[0] for c in covers), min(c[1] for c in covers)

    def _chunk_for(self, variable, yr):
        r"""The chunk file holding year ``yr``, or the last one for the single
        year after the series ends; None when the variable has no files.

        CESM2-WACCM stops at 2299 while a 2015-2300 run needs 2300 (discussion
        #8), so exactly that one year is held, and said once per variable in
        the log, the same rule as the atmosphere. Any other year outside the
        files used to be served from the nearest chunk without a word: a
        historical run starting before its ocean, or a tree with a chunk
        missing, ran on the wrong decade and reported success.
        """
        spans = self.spans(variable)
        if not spans:
            return None
        for first, last, path in spans:
            if first <= yr <= last:
                return path
        end = max(sp[1] for sp in spans)
        if yr == end + 1:
            if variable not in self._persisted:
                self._persisted.add(variable)
                if _comm_rank() == 0:
                    print(f"  ISMIP7Ocean: {variable} has no year {yr}; "
                          f"holding {end}, the last year on disk", flush=True)
            return max(spans, key=lambda sp: sp[1])[2]
        where = ("precedes the series there" if yr < spans[0][0]
                 else "is past the end of the series there" if yr > end
                 else "falls between the chunk files there")
        raise FileNotFoundError(
            f"ISMIP7Ocean: {variable} for {self.esm} {self.scenario} has no "
            f"year {yr} in {self._var_dir(variable)}: it {where} "
            f"({', '.join(f'{a}-{b}' for a, b, _ in spans)}; only the single "
            f"year after the end is held, 2300 after 2299)."
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
        from scipy.interpolate import RegularGridInterpolator

        yr = int(year)
        key = (variable, yr)
        if key in self._interp_cache:
            return self._interp_cache[key]

        best = self._chunk_for(variable, yr)
        if best is None:
            return None

        ds = _open_forcing(best)
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
            da = da.isel(time=_nearest_year_index(ds["time"], yr, ds, f"'{variable}'"))

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

    kind = "fracture"

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585"):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self._collapse_mask = None
        self._excess_melt = None
        self._collapse_mask_path = None

    def _fracture_dir(self):
        root = _find_ismip7_data(self.data_root)
        if root is None:
            return None
        return os.path.join(root, self.esm, self.scenario, "fracture")

    def fracture_dir(self):
        r"""Where the collapse mask is looked for, for error messages."""
        return self._fracture_dir() or f"<no ISMIP7 data root under {self.data_root}>"

    def has_collapse_mask(self):
        r"""True once ``load`` has found an ice-shelf collapse mask."""
        return self._collapse_mask is not None

    def load(self):
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
            self._collapse_mask = _open_forcing(found["collapse_mask"])
            self._collapse_mask_path = found["collapse_mask"]
        if "excess_melt" in found:
            self._excess_melt = _open_forcing(found["excess_melt"])

        return self

    def version(self):
        r"""The version of the collapse mask ``load`` opened, from its filename
        (``..._8km-v2.1.nc``), or None when none was opened or it names none."""
        if self._collapse_mask_path is None:
            return None
        m = re.search(r"[_-](v\d+(?:\.\d+)*)\.nc$", os.path.basename(self._collapse_mask_path))
        return m.group(1) if m else None

    def check_min_version(self):
        r"""Raise when the opened collapse mask is older than
        ``FRACTURE_MIN_VERSION`` allows for this ESM and scenario. A run under
        a mask mode calls this, so a tree still holding a replaced mask stops
        at startup even while the mirror serves nothing newer."""
        floor = FRACTURE_MIN_VERSION.get((self.esm, self.scenario))
        if floor is None or self._collapse_mask_path is None:
            return
        have = self.version()
        key = version_key(have) if have else None
        if key is not None and key >= version_key(floor):
            return
        raise RuntimeError(
            f"{self.esm} {self.scenario}: the collapse mask on disk is "
            f"{have or 'unversioned'} ({self._collapse_mask_path}), and the "
            f"submission needs {floor} or newer (FRACTURE_MIN_VERSION in "
            f"icepack2_tools/forcing.py). The mirror may still serve the older "
            f"one; {floor} is on Globus. Fetch it with `python "
            f"antarctica/scripts/download_forcing.py --scenarios --esm {self.esm} "
            f"--scenario {self.scenario}` (after `--login`), or copy it into "
            f"{self.fracture_dir()}/{floor}/, then rerun.")

    def provenance(self):
        r"""The collapse mask ``load`` opened, the only fracture product that
        is read: its version is in the filename (``..._8km-v2.1.nc``), and the
        highest version on disk is the one taken, so nothing newer is unread."""
        if self._collapse_mask_path is None:
            return []
        return [{"variable": "collapse_mask",
                 "product": os.path.basename(self._collapse_mask_path),
                 "version": self.version() or "unversioned",
                 "dir": os.path.dirname(self._collapse_mask_path), "newer": []}]

    def get_collapse_mask(self, year, mesh_x, mesh_y):
        r"""Get ice shelf collapse mask (0/1) at given year."""
        import xarray as xr

        if self._collapse_mask is None:
            return np.zeros(len(mesh_x))

        ds = self._collapse_mask
        # By name, then by shape, never by position: the files also carry the
        # scalar grid-mapping variable ``mapping`` as a data variable, and
        # 2-D ``lon``/``lat`` that are coordinates only because ``mask`` has
        # a ``coordinates`` attribute naming them.
        named = [v for v in ds.data_vars
                 if v == "mask" or ds[v].attrs.get("standard_name") == "ice_shelf_collapse_mask"]
        gridded = [v for v in ds.data_vars
                   if {d.lower() for d in ds[v].dims} >= {"x", "y"} and v.lower() not in ("lon", "lat")]
        if not (named or gridded):
            raise KeyError(f"{self._collapse_mask_path} has no gridded variable to read a "
                           f"collapse mask from (found {list(ds.data_vars)})")
        var = (named or gridded)[0]
        da = ds[var]

        if "time" in da.dims:
            # not sel(time=year): that only works on an axis of plain years,
            # and raises on one xarray has decoded to dates
            da = da.isel(time=_nearest_year_index(ds["time"], year, ds, f"'{var}'"))

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
                    ISMIP7 toolbox K50 = 8.5e-5 (K05 4.75e-5, K95 1.375e-4).

    Returns melt rate in m/yr ice equivalent (positive = melting).
    """
    U_factor = (_C_PO / _L_I) * _BETA_S * (_G / (2.0 * abs(_F_CORIOLIS))) \
        * salinity
    melt = K * _MELT_FACTOR * U_factor * tf * np.abs(tf) * sin_alpha
    return melt * _SEC_PER_YEAR


_SLOPE_CAP_WARNED = False


_GEOMETRY_SPACE_WARNED = False
_MELT_SLOPE_WARNED = False


def _warn_melt_slope(npz_path, fitted_under, running_under,
                     refit="calibrate_melt.py", knob="ISMIP7_K_PER_BASIN_NPZ"):
    r"""Say once that the file on disk was fitted under the other slope
    convention. Melt is linear in sin(alpha), so the fit does not transfer."""
    global _MELT_SLOPE_WARNED
    if _MELT_SLOPE_WARNED:
        return
    _MELT_SLOPE_WARNED = True
    if _comm_rank() == 0:
        print(
            f"  WARNING: {os.path.basename(npz_path)} was calibrated under "
            f"ISMIP7_MELT_SLOPE={fitted_under} and this run melts under "
            f"{running_under}; melt is linear in sin(alpha), so the fit does "
            f"not transfer. Refit with {refit} under this run's "
            f"ISMIP7_MELT_SLOPE and ISMIP7_SIN_ALPHA_ANT, or name a matching "
            f"file with {knob}.",
            flush=True,
        )


def _warn_geometry_space(npz_path, fitted_on, running_on,
                         refit="calibrate_melt.py", knob="ISMIP7_K_PER_BASIN_NPZ"):
    r"""Say once that the file on disk was fitted on a geometry other than
    the one this run melts with."""
    global _GEOMETRY_SPACE_WARNED
    if _GEOMETRY_SPACE_WARNED:
        return
    _GEOMETRY_SPACE_WARNED = True
    if _comm_rank() == 0:
        print(
            f"  WARNING: {os.path.basename(npz_path)} was calibrated on "
            f"{fitted_on} geometry and this run melts on {running_on}, so the "
            f"melt it applies may not be the melt it was fitted to. Refit "
            f"with ISMIP7_GEOMETRY_SPACE={running_on} "
            f"{refit}, or name a matching file with {knob}.",
            flush=True,
        )


def _warn_slope_cap(npz_path, cap):
    r"""Say once that the K on disk was fitted against a capped draft slope
    while the forward applies an uncapped one."""
    global _SLOPE_CAP_WARNED
    if _SLOPE_CAP_WARNED:
        return
    _SLOPE_CAP_WARNED = True
    if _comm_rank() == 0:
        print(
            f"  WARNING: {os.path.basename(npz_path)} was calibrated with the "
            f"draft slope capped at sin(alpha) = {cap:g}, and this forward "
            f"applies no cap, so it melts with a field the K was not fitted "
            f"against. Measured on the 2500 m mesh at the reference state "
            f"(GEOMETRY_DISCRETIZATION.md), the uncapped cell slope integrates "
            f"3.7 times the capped melt at K = 1, so this forward applies about "
            f"four times the total the K was fitted to. Capping the forward's "
            f"slope the same way or refitting with ISMIP7_SIN_ALPHA_CAP=inf "
            f"are the two consistent choices.",
            flush=True,
        )


def _check_melt_provenance(npz_path, data, refit, knob, strict=False):
    r"""Compare a melt calibration's recorded slope convention, slope cap and
    geometry with this run's, and say once per kind what differs.

    A fit is only valid for the draft slope it was fitted against, because
    melt is linear in sin(alpha). A file without `melt_slope` predates the
    knob and was fitted on the local slope. Under ant the constant is part of
    the convention, so a fit with another constant does not transfer either.
    Under local, the calibration records the cap it applied (none under dg0
    by default, 5e-3 under cg1) and compute_sin_alpha applies none
    (GEOMETRY_DISCRETIZATION.md). A fit is likewise only valid for the
    geometry it was fitted on: the calibration records the space it melted
    (cell by cell under dg0, on nodes under cg1), and a file without the entry
    predates the tag and was fitted on nodes. With the same slope cap the two
    fits agree within about 10 percent per basin, 22 percent in basin 7
    (GEOMETRY_DISCRETIZATION.md); the mismatch is still reported so a file's
    provenance is never silent.

    ``strict`` is for a thermal-forcing offsets file, the tracked calibration
    or one named with ISMIP7_DELTAT_PER_BASIN_NPZ: the forward has to apply
    the melt that file was fitted to, so any difference is refused. A legacy
    per-basin K file keeps the warnings."""
    from .runconfig import geometry_space
    differences = []

    def differs(sentence, warn, *args):
        differences.append(sentence)
        if not strict:
            warn(npz_path, *args)

    fitted_slope = str(data["melt_slope"]) if "melt_slope" in data else "local"
    if fitted_slope != melt_slope():
        differs(f"ISMIP7_MELT_SLOPE {fitted_slope} against this run's "
                f"{melt_slope()}",
                _warn_melt_slope, fitted_slope, melt_slope(), refit, knob)
    elif fitted_slope == "ant" and "sin_alpha_ant" in data:
        fitted_sin = float(data["sin_alpha_ant"])
        if np.isfinite(fitted_sin) and abs(fitted_sin / sin_alpha_ant() - 1.0) > 0.01:
            differs(f"sin(alpha) {fitted_sin:g} against this run's "
                    f"{sin_alpha_ant():g}",
                    _warn_melt_slope, f"ant with sin(alpha) = {fitted_sin:g}",
                    f"ant with sin(alpha) = {sin_alpha_ant():g}", refit, knob)
    if fitted_slope == "local" and melt_slope() == "local" and "sin_alpha_cap" in data:
        cap = float(data["sin_alpha_cap"])
        if np.isfinite(cap) and cap > 0.0:
            differs(f"the draft slope capped at {cap:g}, where this forward "
                    f"applies no cap", _warn_slope_cap, cap)
    fitted_on = str(data["geometry_space"]) if "geometry_space" in data else "cg1"
    if fitted_on != geometry_space():
        differs(f"{fitted_on} geometry against this run's {geometry_space()}",
                _warn_geometry_space, fitted_on, geometry_space(), refit, knob)
    if strict and differences:
        raise ValueError(
            f"{os.path.basename(npz_path)} was fitted under settings this run "
            f"does not use: {'; '.join(differences)}. The forward has to "
            f"apply the melt its calibration was fitted to. Run under the "
            f"file's settings, or refit with {refit} under this run's and "
            f"name the file with {knob}.")


def imbie2_basin_path(recorded=None):
    r"""The IMBIE2 8 km basin grid a per-basin calibration is stamped
    through: ``recorded`` (the path a calibration file names) where it exists
    on this machine, else the v2 calibration-era file under the data root,
    then the v3 release (an identical basinNumber field, verified July 2026).
    Returns the first candidate that exists, or the v2 path when none does,
    so a caller reports the file it looked for."""
    root = os.environ.get(
        "ISMIP7_DATA_ROOT",
        os.path.join(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))), "ISMIP7", "AIS"),
    )
    candidates = [
        os.path.join(root, "parameterisations", "ocean", "imbie2",
                     "basin_numbers_ismip8km_v2.nc"),
        os.path.join(root, "obs", "ocean", "IMBIE-basins", "v3",
                     "IMBIE-basins_AIS_obs_ocean_v3.nc"),
    ]
    if recorded:
        candidates.insert(0, recorded)
    return next((p for p in candidates if os.path.exists(p)),
                candidates[1 if recorded else 0])


def _basin_on_mesh(mesh_x, mesh_y, imbie2=None):
    r"""IMBIE2 basin number of each mesh dof (nearest 8 km cell, -1 outside).

    The 8 km basin grid is re-read here so a per-basin calibration can be
    stamped onto ANY mesh, not only the one it was fitted on. ``imbie2``
    names the file; otherwise :func:`imbie2_basin_path` finds it."""
    import xarray as xr
    from scipy.interpolate import RegularGridInterpolator

    if imbie2 is None:
        imbie2 = imbie2_basin_path()
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
    return basin_node


def load_deltaT_per_basin(npz_path, mesh_x, mesh_y, fill=0.0, imbie2=None):
    r"""Per-basin thermal-forcing offset deltaT_b [K] stamped onto the mesh.

    ``npz_path`` is the output of ``antarctica/scripts/calibrate_deltaT.py``
    (``basin_ids``, ``deltaT_basin``, ``K``): the ISMIP7 toolbox's
    ``optimise_deltaT`` run through this model's own melt path, i.e. the
    protocol's per-basin adjustment at one toolbox K. Returns
    ``(deltaT_field, K)``; a run that applies the offset melts with that K
    everywhere. ``fill`` is the offset outside the fitted basins: a floating
    dof in a basin the observation table does not carry, in a coverage gap
    (basin 0) or off the 8 km grid (basin -1) still melts at K with its
    unadjusted TF. This is intended, since the toolbox applies its one K
    everywhere, so a run's integrated melt exceeds the fitted basins' total
    by what those dofs carry; the ocean callbacks report that amount when the
    offsets load. The slope convention, slope cap and geometry are checked,
    and a difference from this run is refused, since the fit depends on all
    three and the forward has to apply the melt the file was fitted to."""
    data = np.load(npz_path)
    bids = np.asarray(data["basin_ids"]).astype(int)
    dT = np.asarray(data["deltaT_basin"]).astype(float)
    K = float(data["K"])
    _check_melt_provenance(npz_path, data, "calibrate_deltaT.py",
                           "ISMIP7_DELTAT_PER_BASIN_NPZ", strict=True)
    if imbie2 is None:
        imbie2 = imbie2_basin_path(
            str(data["imbie2_nc"]) if "imbie2_nc" in data else None)
    basin_node = _basin_on_mesh(mesh_x, mesh_y, imbie2)
    field = np.full(len(mesh_x), fill, dtype=float)
    for bid, d in zip(bids, dT):
        if np.isfinite(d):
            field[basin_node == bid] = d
    return field, K


def _deltaT_for_run(cache, npz, mesh_x, mesh_y, ctx=None):
    r"""``(deltaT_field, K)`` for the offsets file ``npz``, loaded once per
    callback. Dofs outside the fitted basins get a zero offset and are
    remembered for `_announce_deltaT`. With the run's ``ctx`` the file's
    contract is checked first (`check_melt_contract`)."""
    if "dT" not in cache:
        if ctx is not None:
            check_melt_contract(npz, ctx)
        field, K = load_deltaT_per_basin(npz, mesh_x, mesh_y, fill=np.nan)
        cache["fitted"] = np.isfinite(field)
        cache["dT"], cache["K"] = np.nan_to_num(field, nan=0.0), K
    return cache["dT"], cache["K"]


def check_melt_contract(npz_path, ctx):
    r"""Refuse a run whose initial geometry differs from the kind its melt
    calibration was fitted on, and return the mesh the calibration was fitted
    on (None for a file with no sidecar).

    The sidecar records the raster sampling that put BedMachine on the
    calibration's cells; ``setup_model`` records the one this run's geometry
    came from (``ctx["raster_sample"]``, read back from the MAP). Another
    sampling moves the draft and the floating set the offsets were fitted on.
    The calibration's cells also carry BedMachine's thickness unfloored, so
    an ice-free cell holds none and takes no melt (`melt_receiving`); a cold
    start that floors the initial thickness (``ctx["thickness_floor"]``,
    ISMIP7_H_CLAMP_INIT, 10 m under the legacy friction law) turns those
    cells into floating ice the offsets were never fitted over, and is
    refused. A context without these entries, such as the inversion's
    reference geometry, is not checked. The mesh is only reported: the
    offsets are per basin and stamp onto any mesh, so a coarse probe runs
    with them, and preflight.py keeps a production core on the calibration's
    mesh."""
    from .runconfig import melt_calibration_contract
    contract = melt_calibration_contract(npz_path)
    if contract is None:
        return None
    fitted = contract.get("raster_sample")
    running = ctx.get("raster_sample")
    if fitted and running and fitted != running:
        raise ValueError(
            f"{os.path.basename(npz_path)} was fitted on geometry sampled "
            f"from BedMachine with raster_sample={fitted}, and this run's "
            f"geometry was sampled with {running}. The forward has to apply "
            f"the melt its calibration was fitted to: refit the offsets on "
            f"this sampling (calibrate_deltaT.py) and name the file with "
            f"ISMIP7_DELTAT_PER_BASIN_NPZ.")
    floor = float(ctx.get("thickness_floor") or 0.0)
    if floor > 0.0:
        raise ValueError(
            f"This cold start floors the initial thickness at {floor:g} m "
            f"(ISMIP7_H_CLAMP_INIT), which makes every ice-free cell hold "
            f"floating ice, and {os.path.basename(npz_path)} was fitted over "
            f"cells holding BedMachine's own thickness. Run with "
            f"ISMIP7_H_CLAMP_INIT=0, the default of the budd and "
            f"regularized_coulomb laws.")
    return contract.get("mesh")


def describe_melt_calibration(dT_npz, K_npz=None, mesh_basename=None):
    r"""Marker lines naming the melt calibration a run melts with: the file,
    its sha256 and K, and the mesh it was fitted on beside this run's.

    ``dT_npz`` is the offsets file ``runconfig.deltat_per_basin_npz``
    resolved, or None on the legacy path, where ``K_npz`` is the per-basin K
    file. They are the paths the caller melts with, so the file the line
    names is the file the run melted with. ``core_report.py`` lifts the line
    into the run's record, which is how every submitted run can be shown to
    have read the same calibration."""
    from .runconfig import (MELT_CALIBRATION_DEFAULT, file_sha256,
                            melt_calibration_contract)
    if dT_npz is None:
        return [f"{FORCING_PROVENANCE_MARKER} ocean melt calibration "
                f"{os.path.basename(K_npz)} sha256 {file_sha256(K_npz)}: a "
                f"legacy per-basin K named with ISMIP7_K_PER_BASIN_NPZ, not "
                f"the tracked calibration"]
    contract = melt_calibration_contract(dT_npz) or {}
    with np.load(dT_npz) as data:
        K = float(data["K"])
        selected = (str(data["selected_as"]) if "selected_as" in data
                    else str(contract.get("selected_as", "")))
    fitted_on = contract.get("mesh", "a mesh its file does not record")
    build = (f" ({contract['mesh_build']})" if contract.get("mesh_build") else "")
    named = ("the tracked default"
             if os.path.abspath(dT_npz) == os.path.abspath(MELT_CALIBRATION_DEFAULT)
             else "named with ISMIP7_DELTAT_PER_BASIN_NPZ")
    line = (f"{FORCING_PROVENANCE_MARKER} ocean melt calibration "
            f"{os.path.basename(dT_npz)} sha256 {file_sha256(dT_npz)} ({named}): "
            f"K {K:.3e}{f' ({selected})' if selected else ''} with a "
            f"thermal-forcing offset per basin, fitted on {fitted_on}{build}")
    if mesh_basename:
        stem = os.path.splitext(os.path.basename(mesh_basename))[0]
        same = stem == fitted_on
        line += (f"; this run's mesh is {stem}"
                 + ("" if same else ", so its integrated melt differs from "
                                    "the fitted basin totals"))
    return [line]


def _announce_deltaT(cache, npz, ctx):
    r"""Say once, after the first melt is written, which offsets the run
    applies and how much of its melt falls outside the fitted basins.
    Collective: every rank enters the callback, so every rank reduces."""
    if cache.get("announced"):
        return
    cache["announced"] = True
    from firedrake import Function, PETSc, assemble, dx
    from .mpi_stats import global_count, global_range, global_size
    comm = ctx["mesh"].comm
    fitted = cache["fitted"]
    lo, hi = global_range(cache["dT"][fitted], comm)
    melt = ctx["ocean_melt"]
    unfitted = Function(melt.function_space())
    unfitted.dat.data[:] = np.where(fitted, 0.0, 1.0)
    gt = float(_RHO_I) / 1e12
    total = float(assemble(melt * dx)) * gt
    outside = float(assemble(melt * unfitted * dx)) * gt
    PETSc.Sys.Print(
        f"  Per-basin deltaT from {npz}: K={cache['K']:.3e} everywhere, "
        f"deltaT {lo:+.2f}..{hi:+.2f} K on "
        f"{global_count(fitted, comm)}/{global_size(fitted, comm)} dofs in "
        f"fitted basins. Melt {total:.1f} Gt/yr on floating cells holding "
        f"ice, of which {outside:.1f} outside the fitted basins at the "
        f"unadjusted TF.")


def load_K_per_basin(npz_path, mesh_x, mesh_y, fill=0.0):
    r"""Load a per-basin calibrated K and return a per-node array.

    `npz_path` should be the output of `antarctica/scripts/calibrate_melt.py`
    and is expected to contain `basin_ids` (int) and `K_basin` (float) plus
    the IMBIE2 basin file path (the IMBIE2 8 km grid is re-read here so
    that the K-field can be remapped to *any* mesh, not just the one used
    during calibration). Optional `melt_slope` and `sin_alpha_ant` record the
    slope convention the K was fitted under, `sin_alpha_cap` the cap on the
    local slope, and `geometry_space` the geometry; a mismatch with this run
    triggers a once-per-run warning.

    Returns an array of shape (len(mesh_x),) of per-node K values, with
    `fill` outside the calibrated basin set or where K_basin is NaN.
    """
    data = np.load(npz_path)
    bids = np.asarray(data["basin_ids"]).astype(int)
    Kbas = np.asarray(data["K_basin"]).astype(float)
    _check_melt_provenance(npz_path, data, "calibrate_melt.py",
                           "ISMIP7_K_PER_BASIN_NPZ")

    basin_node = _basin_on_mesh(mesh_x, mesh_y)

    K_field = np.full(len(mesh_x), fill, dtype=float)
    for bid, kb in zip(bids, Kbas):
        if np.isfinite(kb):
            K_field[basin_node == bid] = kb
    return K_field


def forcing_coords(ctx):
    r"""(x, y) of the dofs the forcing fields are written to.

    Forcing callbacks assign straight into `ctx["accum"].dat.data` etc., so
    they must sample the climatology at THOSE dofs. With CG1 geometry those
    coincide with the mesh vertices, which is what the callbacks used to
    assume; with DG0 geometry they are cell centroids and there are ~2x as
    many, so the vertex assumption would mismatch length outright or - worse
    on a mesh where the counts happened to be close - scramble the mapping.
    `ctx["geom_xy"]` is supplied by the driver; fall back to vertices for
    callers that predate it (all of which are CG1).
    """
    xy = ctx.get("geom_xy")
    if xy is not None:
        return xy
    coords = ctx["mesh"].coordinates.dat.data_ro
    return coords[:, 0], coords[:, 1]


def height_above_flotation(s, b, rho_water=_RHO_SW_FLOTATION, rho_ice=_RHO_ICE):
    r"""Height of the surface ``s`` above the flotation surface over bed ``b``,
    ``s - (b + (rho_water / rho_ice) max(-b, 0))``: zero or negative where
    the ice floats. The forward's flotation surface uses seawater, so the
    melt callbacks and the calibration test flotation with the same density
    the dynamics do; ``is_floating`` is the test itself."""
    b = np.asarray(b, dtype=float)
    return np.asarray(s, dtype=float) - (b + (rho_water / rho_ice) * np.maximum(-b, 0.0))


def is_floating(s, b):
    r"""The flotation test: ``height_above_flotation(s, b) <= 0``. An ice-free
    cell passes it too, open ocean at draft 0 and bare land at exactly 0; the
    melt law acts on :func:`melt_receiving`."""
    return height_above_flotation(s, b) <= 0.0


def melt_receiving(s, b, h):
    r"""The cells the melt law acts on: floating and holding ice (``h > 0``).

    This is the set the calibrations fit on (calibrate_melt.forward_geometry,
    select_melt_parameters.py), so the forward melts exactly the cells whose
    melt its calibration was fitted to. An ice-free cell also passes the
    flotation test, open ocean at draft 0 and bare land at a height above
    flotation of exactly 0. Melting it booked melt that the transport limiter
    then withheld, and a negative thermal forcing there made a source that
    grows ice wherever no front mask clears the cell."""
    return is_floating(s, b) & (np.asarray(h, dtype=float) > 0.0)


# The slope the quadratic law sees. The ISMIP7 reference example is "quadratic
# local with mean Antarctic slope (= no slope dependency)": one constant
# sin(alpha) for every shelf, the mean of the 8 km local draft slope over the
# shelves. The local slope is the notebook's other option, with its caveat
# that gridded slopes are bumpy. The toolbox's K05, K50 and K95 belong to the
# constant-slope law: its own gamma_T conversion gives the value they were
# sampled with, sin(alpha) = 5.1117e-3 for all three (K = gamma_T * 2|f| rho_sw
# / (rho_i g beta_S S0 sin(alpha) yr) with the notebook's constants), and the
# notebook's slope recipe on the ISMIP7 8 km BedMap3 v3 topography gives the
# same 5.1117e-3 (select_melt_parameters.py --geometry notebook8km, September
# 2026). Burgard et al. (2022) tuned against 2.9e-3. The default rounds the
# value the percentiles carry up by 0.065 percent, which moves melt by that
# fraction, so a K read against them means the same thing here.
MELT_SLOPES = ("ant", "local")
MELT_SLOPE_DEFAULT = "ant"
SIN_ALPHA_ANT_DEFAULT = 5.115e-3


def melt_slope():
    r"""``ISMIP7_MELT_SLOPE``: ``ant`` (one constant slope, the protocol's
    reference) or ``local`` (the slope of the draft on this mesh)."""
    value = os.environ.get("ISMIP7_MELT_SLOPE", MELT_SLOPE_DEFAULT).lower()
    if value not in MELT_SLOPES:
        raise ValueError(
            f"ISMIP7_MELT_SLOPE must be one of {MELT_SLOPES}, got {value!r}")
    return value


def sin_alpha_ant():
    r"""``ISMIP7_SIN_ALPHA_ANT``: the constant ``sin(alpha)`` under ``ant``."""
    return float(os.environ.get("ISMIP7_SIN_ALPHA_ANT", SIN_ALPHA_ANT_DEFAULT))


def compute_sin_alpha(ctx):
    r"""Return sin(alpha) of the ice-draft slope, on the geometry space.

    Under ``ISMIP7_MELT_SLOPE=ant`` (the default) this is one constant,
    :func:`sin_alpha_ant`, on every dof: the protocol's mean Antarctic slope.
    Under ``local`` it is the slope of this mesh's draft, computed as below.

    Computes draft = s - h and returns sin(arctan(|grad draft|)) =
    |grad|/sqrt(1 + |grad|^2), as a plain array aligned with the dofs of the
    geometry space (so it matches `ocean_melt` elementwise).

    A DG0 draft has an identically zero cell gradient - its slope lives in the
    inter-cell jumps - so reconstruct a CG1 draft first and differentiate that.
    Same device as the Weertman anchor uses via geometry.surface_slope, and
    legitimate for the same reason: this feeds a melt PARAMETERIZATION, not a
    force in the momentum residual.
    """
    import firedrake as fd
    from .geometry import cg1_lift
    Q = ctx["Q"]
    V = ctx["V"]
    h = ctx["h"]
    s = ctx["s"]
    Q_g = ctx.get("Q_g", Q)
    if melt_slope() == "ant":
        return np.full(h.dat.data_ro.shape[0], sin_alpha_ant())
    if Q_g.ufl_element().degree() == 0:
        draft = fd.Function(Q_g).interpolate(s - h)
        gd = fd.grad(cg1_lift(draft))
        gmag = fd.Function(Q_g).interpolate(
            fd.sqrt(fd.inner(gd, gd))
        ).dat.data_ro
    else:
        draft = fd.Function(Q).interpolate(s - h)
        grad_draft = fd.project(fd.grad(draft), V)
        g = grad_draft.dat.data_ro
        gmag = np.sqrt(g[:, 0] ** 2 + g[:, 1] ** 2)
    return gmag / np.sqrt(1.0 + gmag * gmag)


def _oi_climatology_path(root, var, version):
    r"""Path of one OI-climatology variable for a given release.

    `30_sep` is the 2025-09-30 release under meltMIP/ (the one the melt
    calibration was fitted against, and the default for that reason);
    `06_nov` is the 2026 re-release (1972-2024) mirrored from the
    reorganized GHub share under obs/ocean/climatology/.
    """
    if version == "30_sep":
        return os.path.join(
            root, "meltMIP", f"OI_Climatology_ismip8km_60m_{var}_extrap.nc"
        )
    # Versions are per variable here as everywhere (discussion #37): in
    # September 2026 the share holds tf at v3 and so and thetao at v4, the v4
    # being the fix for the July fault below. Take the highest on disk; v3 is
    # only the name of the path that is reported missing when there is none.
    parent = os.path.join(root, "obs", "ocean", "climatology",
                          f"zhou_annual_{version}", var)
    found = _version_subdirs(parent)
    v = found[-1][1] if found else "v3"
    return os.path.join(
        parent, v,
        f"{var}_AIS_obs_ocean_climatology_zhou_annual_{version}_{v}_1972-2024.nc",
    )


def build_oi_climatology_interpolators(data_root=None, version=None):
    r"""Load the OI climatology TF and so into (z, y, x)
    RegularGridInterpolators (nearest, fill 0). Shared by the observationally
    forced runs (the OCX stopgap) and the inversion's melt. The control reads
    its ESM's ``ctrl`` ocean through :class:`ISMIP7Ocean` instead
    (icepack/ismip7#107). ISMIP7_OI_VERSION selects the
    release (default 30_sep, the one the melt calibration was fitted
    against; on Quartz its files are byte-identical to 06_nov's tf v3 and
    so v4)."""
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
                f"upstream packaging fault (Jul 2026): the 06_nov release "
                f"shipped the tf field inside its v3 so/thetao files. The "
                f"share's v4 of those holds the right variable; fetch it, or "
                f"use ISMIP7_OI_VERSION=30_sep."
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


def make_climatology_ocean_callback(K_field=None, data_root=None):
    r"""Ocean-melt callback with CONSTANT OI-climatology TF/so and evolving
    geometry: the observationally constrained ocean forcing of the OCX
    stopgap, and the climatology every melt calibration is fitted against.

    The run melts with its melt calibration (runconfig.deltat_per_basin_npz):
    the file's one K and its per-basin TF offset (`load_deltaT_per_basin`),
    the tracked calibration unless another is named. ``K_field``, a scalar or
    per-dof array, is read only on the legacy per-basin K path
    (ISMIP7_K_PER_BASIN_NPZ, `load_K_per_basin`)."""
    from .runconfig import deltat_per_basin_npz
    interps = build_oi_climatology_interpolators(data_root)
    dT_npz = deltat_per_basin_npz()
    if dT_npz is None and K_field is None:
        raise ValueError(
            "ISMIP7_K_PER_BASIN_NPZ selects the legacy per-basin K path, and "
            "the caller passed no K for it (load_K_per_basin).")
    dT_cache = {}

    def callback(ctx, t_yr):
        mesh_x, mesh_y = forcing_coords(ctx)
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

        # The protocol's per-basin adjustment: a TF offset at one K.
        K_use = K_field
        if dT_npz is not None:
            dT, K_use = _deltaT_for_run(dT_cache, dT_npz, mesh_x, mesh_y, ctx)
            tf = tf + dT
        melt = quadratic_mixed_slope(tf, sal, sin_a, K=K_use)

        ctx["ocean_melt"].dat.data[:] = np.where(melt_receiving(s, b, h), melt, 0.0)
        if dT_npz is not None:
            _announce_deltaT(dT_cache, dT_npz, ctx)

    return callback


def reject_collapse_mask(what):
    r"""Refuse a mask mode of ``ISMIP7_FRACTURE`` where no collapse mask can
    be applied.

    ``run_simulation`` allocates ``ctx["collapse"]`` from the knob alone and
    announces the forcing, but only a forcing callback carrying an
    :class:`ISMIP7Fracture` ever fills it. A driver that has none would print
    the banner and apply nothing, so it says so at startup instead. The
    protocol defines no collapse mask for the control or the OCX experiment,
    which makes a mask mode a wrong request there rather than a no-op.
    """
    from .runconfig import FRACTURE_MASK_MODES, fracture as _fracture_mode
    if _fracture_mode() in FRACTURE_MASK_MODES:
        raise ValueError(
            f"ISMIP7_FRACTURE={_fracture_mode()} but {what} carries no ice-shelf collapse "
            f"forcing, so no mask can ever be applied. The protocol defines "
            f"no collapse mask for the control or the OCX experiment; run "
            f"them with ISMIP7_FRACTURE=none."
        )


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

    The run's melt calibration (runconfig.deltat_per_basin_npz), the
    tracked one unless another is named, supplies a TF offset per basin and
    its one K, which replace both. `K` and `K_per_basin_npz` are read only
    on the legacy per-basin K path (ISMIP7_K_PER_BASIN_NPZ).

    ISMIP7_K_SCALE multiplies the legacy per-basin K (calibrated against the
    older Paolo/Adusumilli table it integrates 689 vs 865 Gt/yr observed on
    the 2500 m mesh, so 1.26 matched that table's total; the factor belongs
    to that calibration). With an offsets file it is refused.
    """
    if fracture is None:
        reject_collapse_mask("this run's forcing callback")
    from .runconfig import deltat_per_basin_npz
    K_field_cache = {"arr": None}
    dT_npz = deltat_per_basin_npz()
    dT_cache = {}
    K_scale = float(os.environ.get("ISMIP7_K_SCALE", "1.0"))

    def callback(ctx, t_yr):
        mesh_x, mesh_y = forcing_coords(ctx)
        # t_yr is the END of the step; the readers want a calendar year
        yr = forcing_year(t_yr)

        if atm is not None:
            smb = atm.get_smb(yr, mesh_x, mesh_y, anomaly=smb_anomaly)
            if smb_baseline is not None:
                smb = smb + smb_baseline
            ctx["accum"].dat.data[:] = smb

        if ocean is not None and "ocean_melt" in ctx:
            h = ctx["h"].dat.data_ro
            b = ctx["b"].dat.data_ro
            s = ctx["s"].dat.data_ro
            # Ice shelf draft (negative depth below sea level)
            draft = np.minimum(s - h, 0.0)

            tf = ocean.get_thermal_forcing(yr, mesh_x, mesh_y, draft=draft)
            sal = ocean.get_salinity(yr, mesh_x, mesh_y, draft=draft)
            sin_alpha = compute_sin_alpha(ctx)

            # Resolve K: the protocol's per-basin adjustment (a TF offset at
            # one K) first, then a per-basin npz, then the scalar.
            if dT_npz is not None:
                dT, K_use = _deltaT_for_run(dT_cache, dT_npz, mesh_x, mesh_y, ctx)
                tf = tf + dT
            elif K_per_basin_npz is not None:
                if K_field_cache["arr"] is None:
                    K_field_cache["arr"] = load_K_per_basin(
                        K_per_basin_npz, mesh_x, mesh_y, fill=0.0
                    )
                K_use = K_field_cache["arr"]
            else:
                K_use = K

            melt = quadratic_mixed_slope(tf, sal, sin_alpha, K=K_use * K_scale)

            # Melt only floating cells holding ice, the set the calibration
            # was fitted on
            ctx["ocean_melt"].dat.data[:] = np.where(melt_receiving(s, b, h), melt, 0.0)
            if dT_npz is not None:
                _announce_deltaT(dT_cache, dT_npz, ctx)

        if fracture is not None and ctx.get("collapse") is not None:
            # The year's ice-shelf collapse mask on the geometry cells; the
            # transport removes floating cells it flags
            # (simulation.run_simulation, ISMIP7_FRACTURE=mask or mask_front).
            ctx["collapse"][:] = fracture.get_collapse_mask(yr, mesh_x, mesh_y) > 0.5
    return callback
