#!/usr/bin/env python3
r"""ISMIP7 ocean-melt calibration (Burgard quadratic_mixed_slope, local TF).

Geometry: the highest-resolution icepack simulation mesh
(`inversion_icepack2_<lc>.h5`, default lc=2500), with bed / thickness /
surface / floating-mask **reinterpolated directly from BedMachine v4.1**
so thicknesses are not h_clamp'd and the floating mask comes from
BedMachine's authoritative `mask == 3`.

Forcing: OI climatology TF + so (8 km, 60 m vertical) from the ISMIP7
meltMIP folder, sampled at the local ice-shelf draft.

Slope sin(alpha): grad(s - h) projected to CG1 on the fine mesh, capped
at ISMIP7_SIN_ALPHA_CAP (default 5e-3) to suppress unstructured-mesh noise.

Aggregation: per-node melt is integrated to IMBIE2 basins (8 km labels,
nearest-neighbour onto the mesh), then compared against the
Paolo/Adusumilli per-basin observations.

Since melt is linear in K, the Term-1 optimum is closed form:

    K* = sum_b M_obs(b) * M_1(b) / sigma(b)^2
       / sum_b   M_1(b)^2          / sigma(b)^2

Usage:
    ISMIP7_LC=2500 python antarctica/scripts/calibrate_melt.py
"""

import os, sys, csv, glob
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

from icepack2_tools.forcing import quadratic_mixed_slope, _RHO_I
from icepack2_tools.runconfig import lc as _lc

DATA_ROOT = os.environ.get(
    "ISMIP7_DATA_ROOT", os.path.join(_PROJECT, "ISMIP7", "AIS")
)
MESH_DIR = os.path.join(_PROJECT, "antarctica", "mesh")
BEDMACHINE_DIR = os.path.join(_PROJECT, "antarctica", "data", "bedmachine")

LC = _lc()
INV_H5 = os.environ.get(
    "ISMIP7_INV_H5", os.path.join(MESH_DIR, f"inversion_icepack2_{LC}.h5")
)

CLIM_TF = os.path.join(DATA_ROOT, "meltMIP", "OI_Climatology_ismip8km_60m_tf_extrap.nc")
CLIM_SO = os.path.join(DATA_ROOT, "meltMIP", "OI_Climatology_ismip8km_60m_so_extrap.nc")
IMBIE2_NC = os.path.join(
    DATA_ROOT, "parameterisations", "ocean", "imbie2",
    "basin_numbers_ismip8km_v2.nc",
)
OBS_CSV = os.path.join(
    DATA_ROOT, "parameterisations", "ocean", "meltobs",
    "Melt_Paolo_Err_Adusumilli_imbie2_v3.csv",
)

SIN_ALPHA_CAP = float(os.environ.get("ISMIP7_SIN_ALPHA_CAP", "5e-3"))


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
    we saw with the 8km mesh) and uses RegularGridInterpolator with
    method='nearest'.
    """
    import xarray as xr
    from scipy.interpolate import RegularGridInterpolator

    ds = xr.open_dataset(nc_path)
    da = ds[var]
    # Sort axes ascending (RegularGridInterpolator requirement)
    x_arr = np.asarray(ds["x"].values)
    y_arr = np.asarray(ds["y"].values)
    if x_arr[0] > x_arr[-1]:
        x_arr = x_arr[::-1]; flip_x = True
    else:
        flip_x = False
    if y_arr[0] > y_arr[-1]:
        y_arr = y_arr[::-1]; flip_y = True
    else:
        flip_y = False

    zdim = None
    for d in da.dims:
        if d.lower() in ("z", "depth", "lev"):
            zdim = d; break

    if zdim is not None:
        z_arr = np.asarray(ds[zdim].values)
        if z_arr[0] > z_arr[-1]:
            z_arr = z_arr[::-1]; flip_z = True
        else:
            flip_z = False
        # Order axes (z, y, x) and load
        data = da.transpose(zdim, "y", "x").values.astype(np.float32)
        if flip_z: data = data[::-1, :, :]
        if flip_y: data = data[:, ::-1, :]
        if flip_x: data = data[:, :, ::-1]
        # Replace NaNs with 0 so nearest-neighbour returns 0 in coverage gaps
        data = np.nan_to_num(data, nan=0.0)
        interp = RegularGridInterpolator(
            (z_arr, y_arr, x_arr), data,
            method="nearest", bounds_error=False, fill_value=0.0,
        )
        # Clip draft into z range
        z_clipped = np.clip(np.asarray(draft), z_arr[0], z_arr[-1])
        pts = np.column_stack([z_clipped, my, mx])
        out = interp(pts)
    else:
        data = da.transpose("y", "x").values.astype(np.float32)
        if flip_y: data = data[::-1, :]
        if flip_x: data = data[:, ::-1]
        data = np.nan_to_num(data, nan=0.0)
        interp = RegularGridInterpolator(
            (y_arr, x_arr), data,
            method="nearest", bounds_error=False, fill_value=0.0,
        )
        pts = np.column_stack([my, mx])
        out = interp(pts)

    ds.close()
    return out


def _compute_sin_alpha(mesh, h, s):
    Q = h.function_space()
    V = VectorFunctionSpace(mesh, "CG", 1)
    draft = Function(Q).interpolate(s - h)
    grad_draft = fd.project(fd.grad(draft), V)
    g = grad_draft.dat.data_ro
    gmag = np.sqrt(g[:, 0] ** 2 + g[:, 1] ** 2)
    return gmag / np.sqrt(1.0 + gmag * gmag)


def _load_obs():
    bids, mobs, sobs = [], [], []
    with open(OBS_CSV) as f:
        r = csv.reader(f); next(r)
        for row in r:
            if not row or not row[1]: continue
            bids.append(int(row[0]))
            mobs.append(float(row[1]))
            sobs.append(float(row[3]))
    return np.array(bids), np.array(mobs), np.array(sobs)


def main():
    PETSc.Sys.Print("=== ISMIP7 melt calibration (Burgard quadratic_mixed_slope) ===")
    for p in (INV_H5, CLIM_TF, CLIM_SO, IMBIE2_NC, OBS_CSV):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    mesh = _load_mesh()
    Q = FunctionSpace(mesh, "CG", 1)
    PETSc.Sys.Print(f"  Mesh: {mesh.num_vertices()} vertices, "
                    f"{mesh.num_cells()} cells")

    bed, thk, sur, msk = _interp_bedmachine(mesh, Q)
    h_np = thk.dat.data_ro
    b_np = bed.dat.data_ro
    s_np = sur.dat.data_ro
    mask_np = msk.dat.data_ro

    PETSc.Sys.Print(f"  BedMachine: h min={h_np.min():.1f}  med={np.median(h_np):.1f}  "
                    f"max={h_np.max():.1f}; mask: floating(=3) "
                    f"{(np.round(mask_np)==3).sum()} nodes")

    mesh_x = mesh.coordinates.dat.data_ro[:, 0]
    mesh_y = mesh.coordinates.dat.data_ro[:, 1]

    draft = np.minimum(s_np - h_np, 0.0)
    PETSc.Sys.Print(f"  Draft range: {draft.min():.0f} .. {draft.max():.0f} m")

    PETSc.Sys.Print("  Interpolating TF climatology...")
    tf = _grid_interp(CLIM_TF, "tf", mesh_x, mesh_y, draft=draft)
    PETSc.Sys.Print("  Interpolating so climatology...")
    sal = _grid_interp(CLIM_SO, "so", mesh_x, mesh_y, draft=draft)
    PETSc.Sys.Print("  Interpolating IMBIE2 basins...")
    # basinNumber is in [0..19]; NaN in coverage gaps gets filled with 0 by
    # the nearest-neighbour interp. Distinguish by also requiring that the
    # node sits over BedMachine ice (floating or grounded).
    basin_raw = _grid_interp(IMBIE2_NC, "basinNumber", mesh_x, mesh_y)
    basin = np.round(basin_raw).astype(int)

    PETSc.Sys.Print(f"  TF range: {tf.min():.2f} .. {tf.max():.2f} K  "
                    f"S range: {sal.min():.2f} .. {sal.max():.2f} PSU")

    PETSc.Sys.Print("  Computing sin(alpha)...")
    sin_a = _compute_sin_alpha(mesh, thk, sur)
    PETSc.Sys.Print(f"  sin(alpha) uncapped p50={np.median(sin_a):.2e} "
                    f"p95={np.quantile(sin_a, 0.95):.2e}; cap={SIN_ALPHA_CAP:.0e}")
    sin_a = np.minimum(sin_a, SIN_ALPHA_CAP)

    # Authoritative floating mask from BedMachine: mask == 3
    floating = (np.round(mask_np).astype(int) == 3)
    PETSc.Sys.Print(f"  Floating nodes (BedMachine mask==3): "
                    f"{int(floating.sum())} / {len(floating)}")

    # Melt at K = 1 (m/yr ice equivalent)
    melt_1 = quadratic_mixed_slope(tf, sal, sin_a, K=1.0)
    melt_1 = np.where(floating, melt_1, 0.0)

    # Lumped mass: nodal area weights, one assemble only.
    PETSc.Sys.Print("  Assembling lumped mass...")
    v = fd.TestFunction(Q)
    lumped = assemble(v * dx).dat.data_ro
    total_area = float(lumped.sum())
    PETSc.Sys.Print(f"  Sum of nodal weights: {total_area:.3e} m^2 "
                    f"(AIS area ~1.4e13)")

    bids_obs, M_obs, sigma_obs = _load_obs()
    rho_i_si = float(_RHO_I)
    node_kgyr = melt_1 * lumped * rho_i_si

    PETSc.Sys.Print("\n  basin |   M_obs    sigma   |  M_1(K=1)     | M_obs/M_1")
    PETSc.Sys.Print("  ------+--------------------+---------------+-----------")
    M_1 = np.zeros_like(M_obs, dtype=float)
    for i, bid in enumerate(bids_obs):
        mask = (basin == bid)
        M_1[i] = float(node_kgyr[mask].sum()) / 1e12
        ratio = M_obs[i] / M_1[i] if M_1[i] > 1e-12 else float("nan")
        PETSc.Sys.Print(
            f"  {bid:5d} | {M_obs[i]:8.2f} {sigma_obs[i]:7.2f}  | "
            f"{M_1[i]:11.3e}  | {ratio:10.3e}"
        )

    w = 1.0 / np.maximum(sigma_obs, 1e-6) ** 2
    K_star = float(np.sum(M_obs * M_1 * w)) / max(float(np.sum(M_1 * M_1 * w)), 1e-30)
    K_total = float(np.sum(M_obs)) / max(float(np.sum(M_1)), 1e-30)

    # Per-basin K: K_b = M_obs(b) / M_1(b). Linear, so this exactly fits
    # the basin-integrated obs. We'll write it out for use by CTRL2015 /
    # projections so each basin gets its own (calibrated) effective K.
    K_basin = np.where(M_1 > 1e-12, M_obs / M_1, np.nan)

    PETSc.Sys.Print("")
    PETSc.Sys.Print("  Per-basin K (M_obs / M_1, dimensionless):")
    for bid, kb in zip(bids_obs, K_basin):
        flag = ""
        if np.isfinite(kb):
            if kb < 8.5e-5:   flag = "<K5"
            elif kb < 1.15e-4: flag = "in K5-K50"
            elif kb < 1.70e-4: flag = "in K50-K95"
            else:              flag = ">K95"
        PETSc.Sys.Print(f"    basin {bid:2d}: K_b = {kb:.3e}  {flag}")

    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"  Total obs:           {float(np.sum(M_obs)):8.1f} Gt/yr")
    PETSc.Sys.Print(f"  Total model @K=1:    {float(np.sum(M_1)):.3e} Gt/yr")
    PETSc.Sys.Print(f"  K* (Term-1 weighted) {K_star:.3e}")
    PETSc.Sys.Print(f"  K  (total-match)     {K_total:.3e}")
    PETSc.Sys.Print("  Burgard reference:   K5=8.5e-5  K50=1.15e-4  K95=1.70e-4")

    total_at_K = float((quadratic_mixed_slope(tf, sal, sin_a, K=K_star)
                        * (floating).astype(float)
                        * lumped * rho_i_si).sum()) / 1e12
    PETSc.Sys.Print(f"  Sanity: integrated melt at K* = {total_at_K:.1f} Gt/yr")

    # Save per-basin K and the K-on-mesh field (one float per CG1 node)
    K_field = np.zeros_like(basin, dtype=float)
    for bid, kb in zip(bids_obs, K_basin):
        if np.isfinite(kb):
            K_field[basin == bid] = kb
    K_field[~floating] = 0.0

    out_dir = os.path.join(_PROJECT, "antarctica", "results")
    os.makedirs(out_dir, exist_ok=True)
    K_out = os.path.join(out_dir, f"calibrated_K_per_basin_{LC}.npz")
    np.savez(
        K_out,
        basin_ids=bids_obs, K_basin=K_basin,
        M_obs=M_obs, M_1=M_1, sigma_obs=sigma_obs,
        K_star=K_star, K_total=K_total,
        basin_on_mesh=basin, K_field=K_field,
    )
    PETSc.Sys.Print(f"  Saved: {K_out}")


if __name__ == "__main__":
    main()
