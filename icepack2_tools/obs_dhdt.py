r"""Observed mean thickness-change (dH/dt) for the transient inverse procedure.

The standard velocity-only inversion fits ``u`` to observations but never
constrains ``div(h u)``, so nothing stops the initial state from carrying a
flux divergence inconsistent with the observed geometry.  The forward then
either drifts or needs a large frozen apparent-MB correction to stand in for
it.  Constraining the modelled thickness tendency against an observed dH/dt
map closes that gap directly.

Source: ``AntarcticaObsISMIP7-v*.nc`` (the ISMIP7 observations MIPkit, Morlighem
on behalf of the observations focus group).  Two fields are available on the
1 km grid:

``dhdt_smith``
    Smith et al. 2020, ICESat + ICESat-2, **m/yr ice equivalent** (firn
    corrected).  A single mean map over the observation period -- the default.

    PROTOCOL NOTE (declared, Andrew Aug 2026): the Smith mean spans 2003-2019,
    so ~4 of its 16 years post-date the ISMIP7 2015 projection branch point,
    while the protocol dates assimilated initial states "prior to 2015"
    (cheat sheet Fig. 1; initial conditions "any time between 1850 to 2014").
    We keep Smith anyway -- it is the only firn-corrected product in the kit,
    and the alternative (a pre-2015 ``dhdt_cpom`` mean) trades a small
    protocol asterisk for an uncorrected firn signal that rivals the ice
    signal over much of East Antarctica.  The target is therefore interpreted
    as the best estimate of the mean thickness-change rate AROUND the
    initialization epoch, and any post-2015 model-vs-obs comparison (e.g.
    OCX) is mildly contaminated by construction.  Declare this in write-ups.
``dhdt_cpom``
    CPOM radar-altimetry time series (27 epochs), m/yr equivalent and
    explicitly **not** firn corrected, so it is NOT interchangeable with
    ``dhdt_smith``.  Averaged over epochs when selected.

Validation of ``dhdt_smith`` against the BedMachine mask on its own grid:
grounded ice integrates to **-85.6 Gt/yr** over 11 207 x10^3 km^2, which sits
in the IMBIE-3 range for the era, so the field is usable as-is.  Floating ice
integrates to +212 Gt/yr (mean +0.18 m/yr) -- altimetry over shelves is noisy
and confounded by firn/tide/ocean signals, and floating thickness change does
not affect VAF, so the constraint should be restricted to GROUNDED ice.

:func:`load_dhdt_obs` does NOT do that restriction itself: it returns a
**coverage** mask only (where the raster has observations), and the caller must
multiply that mask by its own grounded indicator.  A caller that skips this
step ships a constraint fitted to the +212 Gt/yr shelf signal above.

Sampling
--------
NODATA is *not* filled with zero: zero is a meaningful dH/dt, so a filled cell
would read as "observed to be in balance".  Only finite pixels contribute, and
cells whose observed coverage falls below ``min_coverage`` are masked out.

The pixels are **binned into mesh cells** rather than sampled through the
geometry space.  Routing a 1 km raster through ``sample_to_geometry`` evaluates
it at CG1 vertices, i.e. about three point samples standing in for the ~1000
pixels inside a 32 km cell, and that is badly biased for a field whose signal
is concentrated in narrow outlet regions: measured on the 32 km mesh it
recovered only +67.3 Gt/yr of the raster's true +127.0 Gt/yr integral, roughly
halving both the grounded and floating bands.  Binning every pixel and taking
the area-weighted mean is conservative by construction and reproduces the
integral.

Pixels are assigned to the nearest cell centroid rather than by exact
point-in-triangle location. On these meshes that is a close approximation and
it keeps the assignment a handful of vectorised KD-tree queries (walked in
raster row blocks to bound memory) instead of millions of ``locate_cell``
calls; it is adequate because this field is a misfit *target*, never a term in
the momentum residual.
"""

import os

import numpy as np
from firedrake import Function

_VARIABLES = ("dhdt_smith", "dhdt_cpom")

# Pixel-to-cell reach as a multiple of the local cell scale sqrt(area); see
# the distance guard in load_dhdt_obs.
_REACH = float(os.environ.get("ISMIP7_DHDT_REACH", "0.75"))


def _obs_kit_path(data_root=None):
    r"""Locate the ISMIP7 observations MIPkit netCDF."""
    if data_root is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data_root = os.environ.get("ISMIP7_DATA_ROOT",
                                   os.path.join(here, "ISMIP7", "AIS"))
    explicit = os.environ.get("ISMIP7_OBS_KIT")
    if explicit:
        if not os.path.exists(explicit):
            raise FileNotFoundError(f"ISMIP7_OBS_KIT={explicit} does not exist")
        return explicit
    obs_dir = os.path.join(data_root, "obs", "mipkit")
    if os.path.isdir(obs_dir):
        cands = sorted(f for f in os.listdir(obs_dir)
                       if f.startswith("AntarcticaObsISMIP7") and f.endswith(".nc"))
        if cands:
            # sorted() puts v1.10 before v1.2 lexically; take the newest by
            # numeric version instead.
            def _ver(name):
                import re
                m = re.search(r"v(\d+)\.(\d+)", name)
                return (int(m.group(1)), int(m.group(2))) if m else (0, 0)
            return os.path.join(obs_dir, max(cands, key=_ver))
    raise FileNotFoundError(
        f"ISMIP7 observations MIPkit not found under {obs_dir}. Set "
        f"ISMIP7_OBS_KIT to the AntarcticaObsISMIP7-v*.nc path, or pull it "
        f"with antarctica/scripts/download_forcing.py."
    )


def _cache_rasters(variable, data_root=None, cache_dir=None):
    r"""Write NaN-free value and 0/1 coverage GeoTIFFs for ``variable``.

    The MIPkit is ~11 GB and carries NaN, which ``icepack.interpolate`` cannot
    distinguish from a real value. Rewriting the single 1 km slice we need as
    two small EPSG:3031 rasters costs a one-time pass and makes the sampling
    path identical to BedMachine's.
    """
    import netCDF4 as nc
    import rasterio
    from rasterio.transform import from_origin

    src = _obs_kit_path(data_root)
    if cache_dir is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cache_dir = os.path.join(here, "antarctica", "data", "dhdt_cache")
    os.makedirs(cache_dir, exist_ok=True)
    tag = f"{variable}_{os.path.basename(src).replace('.nc', '')}"
    val_fn = os.path.join(cache_dir, f"{tag}_value.tif")
    cov_fn = os.path.join(cache_dir, f"{tag}_valid.tif")
    if os.path.exists(val_fn) and os.path.exists(cov_fn):
        return val_fn, cov_fn

    with nc.Dataset(src) as d:
        if variable not in d.variables:
            raise KeyError(
                f"{variable!r} not in {src}; available dH/dt fields: "
                f"{[v for v in _VARIABLES if v in d.variables]}"
            )
        x = np.asarray(d.variables["x1km"][:], dtype="f8")
        y = np.asarray(d.variables["y1km"][:], dtype="f8")
        arr = d.variables[variable][:]
        if arr.ndim == 3:          # dhdt_cpom: mean over epochs
            arr = np.ma.mean(arr, axis=0)
        arr = np.ma.filled(np.ma.masked_invalid(arr), np.nan).astype("f8")

    valid = np.isfinite(arr).astype("f8")
    value = np.where(np.isfinite(arr), arr, 0.0)

    # GeoTIFF rows run north -> south; flip if the netCDF y ascends.
    if y[1] > y[0]:
        value = value[::-1, :]
        valid = valid[::-1, :]
        y = y[::-1]
    dx = abs(float(x[1] - x[0]))
    dy = abs(float(y[0] - y[1]))
    transform = from_origin(float(x[0]) - dx / 2.0,
                            float(y[0]) + dy / 2.0, dx, dy)
    prof = dict(driver="GTiff", height=value.shape[0], width=value.shape[1],
                count=1, dtype="float64", crs="EPSG:3031", transform=transform,
                compress="deflate", tiled=True)
    for fn, a in ((val_fn, value), (cov_fn, valid)):
        with rasterio.open(fn, "w", **prof) as dst:
            dst.write(a, 1)
    return val_fn, cov_fn


def load_dhdt_obs(Q_g, variable=None, data_root=None, min_coverage=0.5):
    r"""Observed mean dH/dt (m/yr ice equivalent) on the geometry space.

    ``Q_g`` must be a **DG0** (cell-wise constant) space: the binning below
    treats one dof as one cell, and ``assemble(TestFunction(Q_g)*dx)`` is the
    per-cell area only at degree 0.

    Returns ``(dhdt, mask)``, both Functions on ``Q_g``. ``mask`` is 1 where the
    observed coverage of the cell is at least ``min_coverage`` and 0 elsewhere;
    ``dhdt`` is 0 wherever ``mask`` is 0, so the two multiply cleanly into a
    misfit form. Multiply ``mask`` by a grounded indicator to restrict the
    constraint to grounded ice (recommended -- see module docstring).
    """
    import rasterio
    from firedrake import SpatialCoordinate, VectorFunctionSpace, assemble, dx
    from firedrake import TestFunction
    from scipy.spatial import cKDTree

    if Q_g.ufl_element().degree() != 0:
        raise ValueError(
            "load_dhdt_obs requires a DG0 geometry space (one dof per cell); "
            f"got degree {Q_g.ufl_element().degree()}. The pixel binning and "
            "the cell-area weights are only meaningful cell-wise."
        )

    variable = variable or os.environ.get("ISMIP7_DHDT_VAR", "dhdt_smith")
    if variable not in _VARIABLES:
        raise ValueError(f"dH/dt variable must be one of {_VARIABLES}, "
                         f"got {variable!r}")
    val_fn, cov_fn = _cache_rasters(variable, data_root=data_root)

    with rasterio.open(val_fn) as src:
        value = src.read(1)
        tr = src.transform
    with rasterio.open(cov_fn) as src:
        valid = src.read(1) > 0.5

    ny, nx = value.shape
    px, py = abs(tr.a), abs(tr.e)
    xs = tr.c + px * (np.arange(nx) + 0.5)
    ys = tr.f - py * (np.arange(ny) + 0.5)

    # cell centroids of the geometry space
    W = VectorFunctionSpace(Q_g.mesh(), "DG", 0)
    cen = Function(W).interpolate(SpatialCoordinate(Q_g.mesh())).dat.data_ro
    cen = cen.reshape(-1, 2)
    ncell = cen.shape[0]
    tree = cKDTree(cen)

    # Nearest-centroid alone would snap pixels lying OUTSIDE the mesh (open
    # ocean beyond the domain) onto boundary cells, inflating the integral --
    # measured at +145.9 Gt/yr against a raster truth of +127.0 before this
    # guard. Reject a pixel beyond the local cell scale.
    cell_area = assemble(TestFunction(Q_g) * dx).dat.data_ro
    reach = _REACH * np.sqrt(np.maximum(cell_area, 1e-30))

    comm = Q_g.mesh().comm
    if comm.size > 1:
        from mpi4py import MPI as _MPI
    else:
        _MPI = None

    # The raster is ~13M valid pixels; materializing coordinates, distances and
    # the MPI reduction buffer for all of them at once costs ~0.7 GB on EVERY
    # rank. Walk it in row blocks instead and accumulate the bin counts, which
    # bounds the footprint without changing the result. The blocks are cut on
    # RASTER ROWS, which every rank reads identically, so the per-block reduce
    # below is collective and aligned regardless of the mesh partition.
    rows_per_block = max(1, int(1_000_000 // max(nx, 1)))
    npix = np.zeros(ncell, dtype="f8")
    ssum = np.zeros(ncell, dtype="f8")
    for r0 in range(0, ny, rows_per_block):
        r1 = min(r0 + rows_per_block, ny)
        sel = valid[r0:r1]
        ii, jj = np.nonzero(sel)
        blk_vals = value[r0:r1][sel]
        pts = np.column_stack((xs[jj], ys[r0 + ii]))
        dist, owner = tree.query(pts, k=1)

        # UNDER MPI each rank holds only its own cells, so `cen` is a PARTIAL
        # set of centroids and the nearest LOCAL centroid is not the nearest
        # global one: without this every rank would absorb pixels belonging to
        # other ranks' cells, corrupting partition-boundary cells and -- worse
        # -- making the field depend on the rank count. Keep a pixel only on
        # the rank whose local nearest is also the global nearest.
        if _MPI is not None:
            gmin = np.empty_like(dist)
            comm.Allreduce([dist, _MPI.DOUBLE], [gmin, _MPI.DOUBLE],
                           op=_MPI.MIN)
            owns = dist <= gmin + 1e-9
        else:
            owns = np.ones(dist.shape, dtype=bool)

        keep = owns & (dist <= reach[owner])
        owner, blk_vals = owner[keep], blk_vals[keep]
        npix += np.bincount(owner, minlength=ncell)
        ssum += np.bincount(owner, weights=blk_vals, minlength=ncell)

    # coverage = observed pixel area / cell area
    cov_frac = npix * px * py / np.maximum(cell_area, 1e-30)

    dhdt = Function(Q_g, name="dhdt_obs")
    mask = Function(Q_g, name="dhdt_obs_mask")
    ok = (npix > 0) & (cov_frac >= min_coverage)
    mask.dat.data[:] = np.where(ok, 1.0, 0.0)
    dhdt.dat.data[:] = np.where(ok, ssum / np.maximum(npix, 1.0), 0.0)
    return dhdt, mask
