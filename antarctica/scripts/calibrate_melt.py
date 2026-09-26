#!/usr/bin/env python3
r"""ISMIP7 ocean-melt calibration (Burgard quadratic_mixed_slope, local TF):
the per-basin K the runs used before issue 26 was decided. A run now melts
with one K and a thermal-forcing offset per basin, the tracked calibration
(runconfig.MELT_CALIBRATION_DEFAULT, from select_melt_parameters.py), and
reads a file written here only when ISMIP7_K_PER_BASIN_NPZ names it.

Mesh: the section 4 MAP for the configured ISMIP7_FRICTION, named by
`icepack2_tools/naming.py` (ISMIP7_LC; ISMIP7_INV_H5 names a different MAP).
Only the mesh is read from it.

Geometry: the same ISMIP7_GEOMETRY_SPACE the forward reads (default dg0).

* `dg0` melts the cells the forward melts, through the forward's own path:
  bed and thickness sampled onto the cells (ISMIP7_RASTER_SAMPLE), the
  surface from flotation, the slope of `forcing.compute_sin_alpha`,
  thermal forcing and salinity at each centroid and its own draft,
  the forward's melt set `forcing.melt_receiving` (floating, `h > 0`),
  cell areas. A K fitted here is the K the forward applies, by construction.
  Ice-free cells with `haf <= 0` are left out, as the forward leaves them
  out: the observations cover real shelves only.
* `cg1` is the nodal calibration the earlier K files came from: BedMachine
  interpolated onto CG1 nodes with its raster surface and `mask == 3` as
  the floating mask, grad(draft) projected onto CG1 and capped at 5e-3,
  lumped-mass areas. With the same cap the DG0 fit reproduces these within
  about 10 percent per basin, 22 percent in basin 7
  (GEOMETRY_DISCRETIZATION.md); the slope convention is what separates
  them.

The K file records the geometry and the slope convention it was fitted
under, and `load_K_per_basin` warns once when a run melts under another.

Forcing: OI climatology TF + so (8 km, 60 m vertical) from the ISMIP7
meltMIP folder, sampled at the local ice-shelf draft.

Slope sin(alpha): the same ISMIP7_MELT_SLOPE the forward reads. `ant` (the
default, the ISMIP7 reference) is the constant ISMIP7_SIN_ALPHA_ANT on every
shelf, the slope the toolbox's K05, K50 and K95 were sampled with. `local` is
this mesh's draft slope, capped at ISMIP7_SIN_ALPHA_CAP when one is named; the
default is no cap under dg0, the forward's convention, and 5e-3 under cg1 to
suppress unstructured-mesh noise on the nodal slope.

Aggregation: melt at K = 1 is integrated to IMBIE2 basins (8 km labels,
nearest-neighbour onto the mesh), then compared against the per-basin
observation table that `_obs_csv` resolves (ISMIP7_MELT_OBS_CSV names one).
The run prints the table it opened and its integrated target, with a `[!]`
line when the default search ended on the older table.

Since melt is linear in K, the Term-1 optimum is closed form:

    K* = sum_b M_obs(b) * M_1(b) / sigma(b)^2
       / sum_b   M_1(b)^2          / sigma(b)^2

Output: antarctica/results/calibrated_K_per_basin_<LC>.npz. No run searches
for it: a forward reads it only when ISMIP7_K_PER_BASIN_NPZ names it.
ISMIP7_K_OUT names another destination; a bare filename resolves under
antarctica/results/.

Usage:
    ISMIP7_LC=2500 python antarctica/scripts/calibrate_melt.py
    ISMIP7_K_OUT=/scratch/check/K_2000.npz ISMIP7_LC=2000 \
        python antarctica/scripts/calibrate_melt.py
"""

import os, sys, glob
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
_PROJECT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _PROJECT)

import firedrake as fd
from firedrake import (
    Function, FunctionSpace, VectorFunctionSpace, CheckpointFile,
    assemble, dx,
)
from firedrake.petsc import PETSc

import rasterio
import icepack

from icepack2_tools.forcing import (quadratic_mixed_slope, compute_sin_alpha,
                                    is_floating, melt_receiving,
                                    melt_slope, sin_alpha_ant,
                                    SIN_ALPHA_ANT_DEFAULT, _K_PERCENTILES,
                                    _RHO_I)
K05, K50, K95 = _K_PERCENTILES
from icepack2_tools.geometry import sample_to_geometry
from icepack2_tools.melt_selection import read_melt_table, sample_nearest
from icepack2_tools.mpi_stats import (global_count, global_mean, global_range,
                                      global_size)
from icepack2_tools.naming import map_basename
from icepack2_tools.runconfig import (friction as _friction, lc as _lc,
                                      geometry_space as _geometry_space,
                                      obs_data_root, raster_sample)

DATA_ROOT = os.environ.get(
    "ISMIP7_DATA_ROOT", os.path.join(_PROJECT, "ISMIP7", "AIS")
)
MESH_DIR = os.path.join(_PROJECT, "antarctica", "mesh")
BEDMACHINE_DIR = os.path.join(obs_data_root(), "bedmachine")

LC = _lc()
INV_H5 = os.environ.get(
    "ISMIP7_INV_H5", os.path.join(MESH_DIR, map_basename(_friction(), LC))
)

CLIM_TF = os.path.join(DATA_ROOT, "meltMIP", "OI_Climatology_ismip8km_60m_tf_extrap.nc")
CLIM_SO = os.path.join(DATA_ROOT, "meltMIP", "OI_Climatology_ismip8km_60m_so_extrap.nc")
IMBIE2_NC = os.path.join(
    DATA_ROOT, "parameterisations", "ocean", "imbie2",
    "basin_numbers_ismip8km_v2.nc",
)
# Observed basal melt per IMBIE2 basin. The melt-calibration product re-released
# in July 2026 (Source Cooperative, ismip7-ais-melt-calibration) combines Paolo
# (2023), Davison (2023) and Adusumilli (2020) and raises the integrated target
# from 865 to 1067 Gt/yr, so the total-match K calibrated against the older
# Paolo+Adusumilli table is 23% low. Prefer the new table, fall back to the old
# one so a tree that predates the re-release still runs, and let
# ISMIP7_MELT_OBS_CSV name either explicitly. The new table is searched for in
# two places: where `download_mirror.py` lands it, and beside the old table,
# where a site that stages files by hand tends to put it, as IU Quartz did.
_OBS_CSV_NEW = "Melt_Paolo_Davison_Adusumilli_imbie2.csv"
_OBS_CSV_CANDIDATES = (
    os.path.join(DATA_ROOT, "meltobs", _OBS_CSV_NEW),
    os.path.join(DATA_ROOT, "parameterisations", "ocean", "meltobs",
                 _OBS_CSV_NEW),
    os.path.join(DATA_ROOT, "parameterisations", "ocean", "meltobs",
                 "Melt_Paolo_Err_Adusumilli_imbie2_v3.csv"),
)


def _obs_csv():
    r"""Path to the per-basin melt observations, newest available first."""
    named = os.environ.get("ISMIP7_MELT_OBS_CSV")
    if named:
        return named
    for path in _OBS_CSV_CANDIDATES:
        if os.path.exists(path):
            return path
    return _OBS_CSV_CANDIDATES[-1]


OBS_CSV = _obs_csv()


def _announce_obs_table():
    r"""Name the table in use, loudly when the search fell through to the old one.

    Runs before the existence checks in `main`, so a tree carrying neither table
    is told every place that was searched ahead of the FileNotFoundError, which
    names the old table alone.
    """
    PETSc.Sys.Print(f"  Observation table: {OBS_CSV}")
    if os.environ.get("ISMIP7_MELT_OBS_CSV") or OBS_CSV != _OBS_CSV_CANDIDATES[-1]:
        return
    searched = " nor ".join(_OBS_CSV_CANDIDATES[:-1])
    PETSc.Sys.Print(
        f"  [!] {_OBS_CSV_NEW} was found at neither {searched}, so this "
        f"calibrates against the older Paolo and Adusumilli table, whose "
        f"integrated target is 865.0 Gt/yr where the July 2026 table has "
        f"1067.4. Stage the new table at one of those paths, or name a table "
        f"with ISMIP7_MELT_OBS_CSV."
    )


def _k_out():
    r"""Where the calibration is written.

    The default is ``antarctica/results/calibrated_K_per_basin_<LC>.npz``,
    which a forward reads only when ISMIP7_K_PER_BASIN_NPZ names it.
    ISMIP7_K_OUT names another destination. Like ISMIP7_MAP_OUT, a bare
    filename resolves under the default directory.
    """
    out_dir = os.path.join(_PROJECT, "antarctica", "results")
    named = os.environ.get("ISMIP7_K_OUT")
    if not named:
        return os.path.join(out_dir, f"calibrated_K_per_basin_{LC}.npz")
    path = named if os.path.dirname(named) else os.path.join(out_dir, named)
    # np.savez appends the suffix when it is absent; keep the printed path true.
    return path if path.endswith(".npz") else path + ".npz"


# The geometry the melt is evaluated on, the same knob the forward reads.
# ``dg0`` fits K through the forward's own cell by cell melt path (section
# ``forward_geometry``), so the forward applies the melt its K was fitted to;
# ``cg1`` is the nodal calibration the earlier K files came from.
GEOMETRY = _geometry_space()


def default_slope_cap(space):
    r"""The draft-slope cap a calibration applies when none is named.

    The nodal CG1 slope carries unstructured-mesh noise, so the CG1
    calibration caps it at 5e-3. The forward applies no cap to its cell slope
    (``forcing.compute_sin_alpha``), so a DG0 calibration fits the slope the
    forward melts with, uncapped, unless ``ISMIP7_SIN_ALPHA_CAP`` names one.
    """
    return float("inf") if space == "dg0" else 5e-3


SIN_ALPHA_CAP = float(os.environ.get("ISMIP7_SIN_ALPHA_CAP",
                                     default_slope_cap(GEOMETRY)))
# The slope convention, the same knob the forward reads: ``ant`` is one
# constant sin(alpha) (the protocol's reference, no cap applies); ``local``
# is this mesh's draft slope, capped as above.
MELT_SLOPE = melt_slope()
# The ice to seawater density ratio simulation.py builds the DG0 surface with.
RHO_RATIO = 917.0 / 1024.0


def _load_mesh():
    PETSc.Sys.Print(f"  Loading mesh from: {INV_H5}")
    with CheckpointFile(INV_H5, "r") as chk:
        mesh = chk.load_mesh()
    return mesh


def _bedmachine_path():
    m = glob.glob(os.path.join(BEDMACHINE_DIR, "*.nc"))
    if not m:
        raise FileNotFoundError(BEDMACHINE_DIR)
    return m[0]


def _interp_bedmachine(mesh, Q):
    bm = _bedmachine_path()
    PETSc.Sys.Print(f"  Interpolating BedMachine onto mesh: {bm}")
    bed = icepack.interpolate(rasterio.open(f"netcdf:{bm}:bed"), Q)
    thk = icepack.interpolate(rasterio.open(f"netcdf:{bm}:thickness"), Q)
    sur = icepack.interpolate(rasterio.open(f"netcdf:{bm}:surface"), Q)
    msk = icepack.interpolate(rasterio.open(f"netcdf:{bm}:mask"), Q)
    return bed, thk, sur, msk


def _grid_interp(nc_path, var, mx, my, draft=None):
    r"""Nearest-neighbour interpolation onto mesh nodes via scipy.

    Loads the field once into memory (avoids xarray.interp memory blowups
    we saw with the 8km mesh); NaN in coverage gaps and points off the grid
    read as 0, as the forward reads the climatology
    (``melt_selection.sample_nearest``, which can keep them NaN instead).
    """
    return sample_nearest(nc_path, var, mx, my, draft=draft, fill=0.0)


def _compute_sin_alpha(mesh, h, s):
    Q = h.function_space()
    V = VectorFunctionSpace(mesh, "CG", 1)
    draft = Function(Q).interpolate(s - h)
    grad_draft = fd.project(fd.grad(draft), V)
    g = grad_draft.dat.data_ro
    gmag = np.sqrt(g[:, 0] ** 2 + g[:, 1] ** 2)
    return gmag / np.sqrt(1.0 + gmag * gmag)


def calibration_geometry(mesh):
    r"""The melt inputs on CG1 nodes, as the earlier calibrations built them:
    BedMachine interpolated with its raster surface and mask, grad(draft)
    projected onto CG1, nodal lumped-mass areas."""
    Q = FunctionSpace(mesh, "CG", 1)
    bed, thk, sur, msk = _interp_bedmachine(mesh, Q)
    h_np, s_np, mask_np = thk.dat.data_ro, sur.dat.data_ro, msk.dat.data_ro
    PETSc.Sys.Print(f"  BedMachine: h min={h_np.min():.1f}  med={np.median(h_np):.1f}  "
                    f"max={h_np.max():.1f}; mask: floating(=3) "
                    f"{(np.round(mask_np)==3).sum()} nodes")
    return {
        "x": mesh.coordinates.dat.data_ro[:, 0],
        "y": mesh.coordinates.dat.data_ro[:, 1],
        "draft": np.minimum(s_np - h_np, 0.0),
        "sin_a": (np.full(h_np.shape[0], sin_alpha_ant()) if MELT_SLOPE == "ant"
                  else _compute_sin_alpha(mesh, thk, sur)),
        # Authoritative floating mask from BedMachine: mask == 3
        "floating": np.round(mask_np).astype(int) == 3,
        "area": assemble(fd.TestFunction(Q) * dx).dat.data_ro,
        "dofs": "nodes",
    }


def forward_cells(mesh):
    r"""The forward's cold-start geometry on the DG0 cells of ``mesh``, as
    ``simulation.setup_model`` builds it for the budd and regularized_coulomb
    laws (no thickness floor): bed and thickness sampled onto the cells
    (``ISMIP7_RASTER_SAMPLE``) and the surface from flotation. The one
    builder the calibrations and ``check_melt_bound.py`` share, so the cells
    a fit sums over and the cells a check melts are the same.

    Returns ``Q``, ``Q_g`` and ``V``, the Functions ``b``, ``h`` and ``s``,
    and the centroids ``x`` and ``y``."""
    Q = FunctionSpace(mesh, "CG", 1)
    Q_g = FunctionSpace(mesh, "DG", 0)
    bm = _bedmachine_path()
    PETSc.Sys.Print(f"  Sampling BedMachine onto DG0 cells ({raster_sample()}): {bm}")
    b_dg = sample_to_geometry(rasterio.open(f"netcdf:{bm}:bed"), Q_g, Q,
                              method=raster_sample())
    h_dg = sample_to_geometry(rasterio.open(f"netcdf:{bm}:thickness"), Q_g, Q,
                              method=raster_sample())
    s_dg = Function(Q_g).interpolate(
        fd.max_value(b_dg + h_dg, (1.0 - RHO_RATIO) * h_dg))
    xy = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    return {"Q": Q, "Q_g": Q_g, "V": VectorFunctionSpace(mesh, "CG", 1),
            "b": b_dg, "h": h_dg, "s": s_dg,
            "x": xy[:, 0].copy(), "y": xy[:, 1].copy()}


def forward_geometry(mesh):
    r"""The melt inputs on DG0 cells, as the forward melts them: the cells of
    `forward_cells`, the cell slope of ``forcing.compute_sin_alpha``, the
    forcing at each cell centroid and its own draft, the forward's melt set
    and cell areas.

    The melt set is ``forcing.melt_receiving``, ``haf <= 0`` on cells with
    ``h > 0``, the one the forward's callbacks melt. An ice-free cell also has
    ``haf <= 0``; the observations cover real shelves only, so it is left out
    of the fit, and the forward leaves it out of the melt."""
    c = forward_cells(mesh)
    b_np, h_np, s_np = c["b"].dat.data_ro, c["h"].dat.data_ro, c["s"].dat.data_ro
    afloat = is_floating(s_np, b_np)
    floating = melt_receiving(s_np, b_np, h_np)
    comm = mesh.comm
    h_lo, h_hi = global_range(h_np, comm)
    PETSc.Sys.Print(f"  BedMachine on cells: h min={h_lo:.1f}  "
                    f"mean={global_mean(h_np, comm):.1f}  max={h_hi:.1f}; "
                    f"floating (haf <= 0, h > 0) "
                    f"{global_count(floating, comm)} cells, "
                    f"{global_count(afloat & ~(h_np > 0), comm)} ice-free "
                    f"haf <= 0 cells left out")
    return {
        "x": c["x"],
        "y": c["y"],
        "draft": np.minimum(s_np - h_np, 0.0),
        "sin_a": compute_sin_alpha({"Q": c["Q"], "V": c["V"], "Q_g": c["Q_g"],
                                    "h": c["h"], "s": c["s"]}),
        "floating": floating,
        "area": assemble(fd.TestFunction(c["Q_g"]) * dx).dat.data_ro,
        "dofs": "cells",
    }


def fit_per_basin_K(basin, melt_1, area, bids_obs, M_obs, sigma_obs):
    r"""Per-basin K from melt at K = 1.

    ``melt_1`` is the melt at K = 1 (m/yr of ice) on dofs with areas ``area``
    (m^2), zero where not floating, and ``basin`` the IMBIE2 basin of each
    dof. Melt is linear in K, so each basin's K is the ratio of the observed
    to the modelled total, the Term-1 weighted scalar is closed form, and the
    total-match scalar is the ratio of the sums. Returns the modelled totals
    at K = 1 in Gt/yr and the three K, with NaN for a basin no dof melts in.
    """
    kgyr = melt_1 * area * float(_RHO_I)
    M_1 = np.array([float(kgyr[basin == bid].sum()) / 1e12 for bid in bids_obs])
    w = 1.0 / np.maximum(sigma_obs, 1e-6) ** 2
    K_star = float(np.sum(M_obs * M_1 * w)) / max(float(np.sum(M_1 * M_1 * w)), 1e-30)
    K_total = float(np.sum(M_obs)) / max(float(np.sum(M_1)), 1e-30)
    K_basin = np.where(M_1 > 1e-12, M_obs / np.where(M_1 > 1e-12, M_1, 1.0), np.nan)
    return {"M_1": M_1, "K_star": K_star, "K_total": K_total, "K_basin": K_basin}


def _load_obs():
    r"""Basin ids, observed melt and its uncertainty, both in Gt/yr, from
    the table ``_obs_csv`` resolved (``melt_selection.read_melt_table``
    locates the columns by header name, so either published layout reads).
    """
    return read_melt_table(OBS_CSV)


def main():
    PETSc.Sys.Print("=== ISMIP7 melt calibration (Burgard quadratic_mixed_slope) ===")
    _announce_obs_table()
    K_out = _k_out()
    PETSc.Sys.Print(f"  K output: {K_out}")
    for p in (INV_H5, CLIM_TF, CLIM_SO, IMBIE2_NC, OBS_CSV):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    # Read ahead of the mesh work: the integrated target identifies the table
    # at a glance, and a malformed one fails here instead of after BedMachine.
    bids_obs, M_obs, sigma_obs = _load_obs()
    PETSc.Sys.Print(f"  Integrated target: {float(np.sum(M_obs)):.1f} Gt/yr "
                    f"over {len(bids_obs)} basins")

    mesh = _load_mesh()
    # num_vertices()/num_cells() count this rank's plex, halo included; the
    # coordinate dofs and the owned cell set are reduced to global totals.
    PETSc.Sys.Print(f"  Mesh: {global_size(mesh.coordinates)} vertices, "
                    f"{mesh.comm.allreduce(mesh.cell_set.size)} cells")
    how = ("the forward's cell by cell melt path" if GEOMETRY == "dg0"
           else "CG1 nodes")
    PETSc.Sys.Print(f"  Geometry: {GEOMETRY} ({how})")

    g = forward_geometry(mesh) if GEOMETRY == "dg0" else calibration_geometry(mesh)
    dofs = g["dofs"]
    mesh_x, mesh_y, draft = g["x"], g["y"], g["draft"]
    PETSc.Sys.Print(f"  Draft range: {draft.min():.0f} .. {draft.max():.0f} m")

    PETSc.Sys.Print("  Interpolating TF climatology...")
    tf = _grid_interp(CLIM_TF, "tf", mesh_x, mesh_y, draft=draft)
    PETSc.Sys.Print("  Interpolating so climatology...")
    sal = _grid_interp(CLIM_SO, "so", mesh_x, mesh_y, draft=draft)
    PETSc.Sys.Print("  Interpolating IMBIE2 basins...")
    # basinNumber is in [0..19]; NaN in coverage gaps gets filled with 0 by
    # the nearest-neighbour interp. Distinguish by also requiring that the
    # dof sits over BedMachine ice (floating or grounded).
    basin_raw = _grid_interp(IMBIE2_NC, "basinNumber", mesh_x, mesh_y)
    basin = np.round(basin_raw).astype(int)

    PETSc.Sys.Print(f"  TF range: {tf.min():.2f} .. {tf.max():.2f} K  "
                    f"S range: {sal.min():.2f} .. {sal.max():.2f} PSU")

    sin_a = g["sin_a"]
    if MELT_SLOPE == "ant":
        PETSc.Sys.Print(f"  Slope: ant, sin(alpha) = {sin_alpha_ant():g} on every "
                        f"{dofs[:-1]} (the protocol's mean Antarctic slope)")
    else:
        PETSc.Sys.Print(f"  Slope: local; sin(alpha) uncapped p50={np.median(sin_a):.2e} "
                        f"p95={np.quantile(sin_a, 0.95):.2e}; cap={SIN_ALPHA_CAP:g}")
        sin_a = np.minimum(sin_a, SIN_ALPHA_CAP)

    floating = g["floating"]
    PETSc.Sys.Print(f"  Floating {dofs}: {int(floating.sum())} / {len(floating)}")

    # Melt at K = 1 (m/yr ice equivalent)
    melt_1 = quadratic_mixed_slope(tf, sal, sin_a, K=1.0)
    melt_1 = np.where(floating, melt_1, 0.0)

    # Areas: nodal lumped mass on CG1, the cell area on DG0.
    area = g["area"]
    PETSc.Sys.Print(f"  Sum of {dofs[:-1]} areas: {float(area.sum()):.3e} m^2 "
                    f"(AIS area ~1.4e13); floating "
                    f"{float(area[floating].sum()) / 1e6:.1f} km^2")

    fit = fit_per_basin_K(basin, melt_1, area, bids_obs, M_obs, sigma_obs)
    M_1, K_star, K_total, K_basin = (fit["M_1"], fit["K_star"],
                                     fit["K_total"], fit["K_basin"])

    PETSc.Sys.Print("\n  basin |   M_obs    sigma   |  M_1(K=1)     | M_obs/M_1")
    PETSc.Sys.Print("  ------+--------------------+---------------+-----------")
    for i, bid in enumerate(bids_obs):
        ratio = M_obs[i] / M_1[i] if M_1[i] > 1e-12 else float("nan")
        PETSc.Sys.Print(
            f"  {bid:5d} | {M_obs[i]:8.2f} {sigma_obs[i]:7.2f}  | "
            f"{M_1[i]:11.3e}  | {ratio:10.3e}"
        )

    # Per-basin K: K_b = M_obs(b) / M_1(b). Linear, so this exactly fits
    # the basin-integrated obs. Written out for the control and projections,
    # so each basin melts with its own calibrated K.
    PETSc.Sys.Print("")
    PETSc.Sys.Print("  Per-basin K (M_obs / M_1, dimensionless):")
    for bid, kb in zip(bids_obs, K_basin):
        flag = ""
        if np.isfinite(kb):
            if kb < K05:    flag = "<K05"
            elif kb < K50:  flag = "in K05-K50"
            elif kb < K95:  flag = "in K50-K95"
            else:           flag = ">K95"
        PETSc.Sys.Print(f"    basin {bid:2d}: K_b = {kb:.3e}  {flag}")

    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"  Total obs:           {float(np.sum(M_obs)):8.1f} Gt/yr")
    PETSc.Sys.Print(f"  Total model @K=1:    {float(np.sum(M_1)):.3e} Gt/yr")
    PETSc.Sys.Print(f"  K* (Term-1 weighted) {K_star:.3e}")
    PETSc.Sys.Print(f"  K  (total-match)     {K_total:.3e}")
    PETSc.Sys.Print(f"  ISMIP7 toolbox (July 2026, constant slope "
                    f"{SIN_ALPHA_ANT_DEFAULT:g}): K05={K05:g} K50={K50:g} K95={K95:g}"
                    + ("" if MELT_SLOPE == "ant" else
                       "; not like for like with a local-slope fit"))

    rho_i_si = float(_RHO_I)
    total_at_K = float((quadratic_mixed_slope(tf, sal, sin_a, K=K_star)
                        * (floating).astype(float)
                        * area * rho_i_si).sum()) / 1e12
    PETSc.Sys.Print(f"  Sanity: integrated melt at K* = {total_at_K:.1f} Gt/yr")

    # The K on the mesh (one float per dof of the calibration geometry).
    K_field = np.zeros_like(basin, dtype=float)
    for bid, kb in zip(bids_obs, K_basin):
        if np.isfinite(kb):
            K_field[basin == bid] = kb
    K_field[~floating] = 0.0

    os.makedirs(os.path.dirname(K_out), exist_ok=True)
    # Provenance travels with the numbers. Two published observation tables are
    # in circulation and their integrated targets differ by 23%, so a K file
    # that does not name its own source cannot be told apart from the other
    # calibration once it is on disk. The forward reads basin_ids and K_basin,
    # plus sin_alpha_cap and geometry_space, which load_K_per_basin compares
    # against its own slope and geometry and warns about; the other entries
    # cost nothing.
    np.savez(
        K_out,
        basin_ids=bids_obs, K_basin=K_basin,
        M_obs=M_obs, M_1=M_1, sigma_obs=sigma_obs,
        K_star=K_star, K_total=K_total,
        basin_on_mesh=basin, K_field=K_field,
        obs_csv=os.path.basename(OBS_CSV), obs_total_gtyr=float(M_obs.sum()),
        mesh_source=os.path.basename(INV_H5),
        sin_alpha_cap=(SIN_ALPHA_CAP if MELT_SLOPE == "local" else float("inf")),
        melt_slope=MELT_SLOPE,
        sin_alpha_ant=(sin_alpha_ant() if MELT_SLOPE == "ant" else float("nan")),
        geometry_space=GEOMETRY, raster_sample=raster_sample(),
        floating_area_km2=float(area[floating].sum()) / 1e6,
        integrated_at_K_star_gtyr=total_at_K,
    )
    PETSc.Sys.Print(f"  Saved: {K_out}")

if __name__ == "__main__":
    main()
