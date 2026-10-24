r"""A/B the ways of putting BedMachine onto DG0 geometry cells.

    mpiexec -n 2 python antarctica/scripts/probe_raster_sampling.py MESH.msh

``vertex``    : icepack bilinear at the 3 vertices, L2-projected (3 px / cell)
``cell_mean`` : mean of the raster over the cell (geometry.raster_cell_mean)
``patch``     : hat-weighted patch mean, the lumped L2 projection onto CG1,
                projected (patchmean.py); conserves the raster's mass

Reports, for bed b, thickness H and surface s: the rms facet jump over all
interior facets and over GROUNDED interior facets (the DG0 driving stress is
entirely the facet jump in s, so this is the roughness the momentum balance
sees), the mean ice thickness along the exterior boundary where ice exists
(front <h>), total volume, and the size of the change in b itself. No solve.
"""
import glob
import math
import os
import sys
import time

import numpy as np
import rasterio
from firedrake import (And, Function, FunctionSpace, Mesh, assemble, avg,
                       conditional, dS, ds, dx, jump, max_value)
from firedrake.petsc import PETSc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from icepack2_tools.geometry import sample_bed_thickness, sample_to_geometry  # noqa: E402, F401

mesh_fn = sys.argv[1]
bm = os.environ.get("ISMIP7_BEDMACHINE") or sorted(glob.glob(os.path.join(
    os.path.dirname(__file__), "..", "data", "bedmachine", "*.nc")))[0]
mesh = Mesh(mesh_fn)
Q = FunctionSpace(mesh, "DG", 0)
Qc = FunctionSpace(mesh, "CG", 1)
rho = 917.0 / 1024.0
one = Function(Q).assign(1.0)
L_all = assemble(avg(one) * dS)
# Cell size indicator: cells wider than 10 km are "interior" (the coarse
# region of a 20 km-interior mesh); their shared facets are where the user's
# interior problem lives.
import ufl  # noqa: E402
_area = Function(Q).interpolate(ufl.CellVolume(mesh))
big = Function(Q).interpolate(conditional(_area > (10e3 ** 2) * 0.433, 1.0, 0.0))
L_int = assemble(big("+") * big("-") * dS)
n_big = int(mesh.comm.allreduce(float(big.dat.data_ro.sum())))


def stats(b, H):
    s = Function(Q).interpolate(max_value(b + H, (1.0 - rho) * H))
    g = Function(Q).interpolate(conditional(And(rho * H + b > 0.0, H > 1.0), 1.0, 0.0))
    L_g = assemble(avg(g) * dS)
    rms = lambda f: math.sqrt(assemble(jump(f) ** 2 * dS) / L_all)              # noqa: E731
    rms_g = lambda f: math.sqrt(assemble(avg(g) * jump(f) ** 2 * dS) / L_g)     # noqa: E731
    front = assemble(H * ds) / max(assemble(conditional(H > 1.0, 1.0, 0.0) * ds), 1.0)
    rms_i = lambda f: (math.sqrt(assemble(big("+") * big("-") * jump(f) ** 2 * dS) / L_int)  # noqa: E731
                       if L_int > 0 else float("nan"))
    return dict(rms_jb=rms(b), rms_js=rms(s), rms_jh=rms(H),
                rms_jb_g=rms_g(b), rms_js_g=rms_g(s),
                rms_jb_i=rms_i(b), rms_js_i=rms_i(s), rms_jh_i=rms_i(H),
                front_h=front, vol_km3=assemble(H * dx) / 1e9, s=s)


# BedMachine's own per-pixel mask (0 ocean, 1 ice-free land, 2 grounded ice,
# 3 floating ice, 4 Lake Vostok): its cell mean is the grounded / ice FRACTION
# of each cell, the reference either sampling method is trying to represent.
from icepack2_tools.geometry import raster_cell_mean  # noqa: E402
_mask_ds = rasterio.open(f"netcdf:{bm}:mask")
_win_mask = _mask_ds.read(1)
# Cheap route: sample the mask raster with the same lattice, via a float view
# of the predicate. raster_cell_mean reads the dataset; wrap a MemoryFile.
from rasterio.io import MemoryFile  # noqa: E402
def frac_of(pred):
    arr = pred(_win_mask).astype("f4")
    # A clean untiled GTiff profile: the NetCDF driver's block size (3334) is
    # not a multiple of 16 and GTiff refuses it.
    prof = dict(driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
                dtype="float32", crs=_mask_ds.crs, transform=_mask_ds.transform, tiled=False)
    with MemoryFile() as mf:
        with mf.open(**prof) as d:
            d.write(arr, 1)
        with mf.open() as d:
            return raster_cell_mean(d, Q)
f_grounded = frac_of(lambda m: (m == 2) | (m == 4))
n_cells = int(mesh.comm.allreduce(float(one.dat.data_ro.size)))

# The raster's own mass over the mesh, by the patch rule's quadrature: what
# the patch sampling conserves and the others miss.
from icepack2_tools import patchmean as P  # noqa: E402
from icepack2_tools.runconfig import patch_qdegree  # noqa: E402
_deg = patch_qdegree()
_Qs, _xq, _yq = P.quadrature_points(mesh, _deg)
_fq = Function(_Qs)
_fq.dat.data[:] = np.nan_to_num(P.raster_at(rasterio.open(f"netcdf:{bm}:thickness"), _xq, _yq))
vol_raster = assemble(_fq * P._measure((_deg, "default"))) / 1e9
PETSc.Sys.Print(f"  BedMachine thickness over the mesh (quadrature degree {_deg}): {vol_raster:.4e} km3")

# Methods to compare: the base samplings and, with _front, the front rule on
# top of them (geometry.sample_bed_thickness). Override with a comma list as
# the second argument.
methods = (sys.argv[2].split(",") if len(sys.argv) > 2
           else ["vertex", "cell_mean", "patch", "vertex_front", "patch_front"])
res = {}
for method in methods:
    t0 = time.time()
    b, H, counts = sample_bed_thickness(bm, Q, Qc, floor=0.0, method=method)
    if counts is not None:
        PETSc.Sys.Print(f"  {method:12s} front cells: {counts}")
    r = stats(b, H)
    r["t"] = time.time() - t0
    # DG0 flotation as the momentum balance sees it, vs BedMachine's fraction
    gmask = Function(Q).interpolate(conditional(And(rho * H + b > 0.0, H > 1.0), 1.0, 0.0))
    fl = Function(Q).interpolate(conditional(And(rho * H + b <= 0.0, H > 1.0), 1.0, 0.0))
    ng = int(mesh.comm.allreduce(float(gmask.dat.data_ro.sum())))
    nf = int(mesh.comm.allreduce(float(fl.dat.data_ro.sum())))
    # grounded by the method but mostly floating/ocean by BedMachine, and the reverse
    g_but_lt_half = int(mesh.comm.allreduce(float(((gmask.dat.data_ro > 0.5) & (f_grounded.dat.data_ro < 0.5)).sum())))
    f_but_gt_half = int(mesh.comm.allreduce(float(((fl.dat.data_ro > 0.5) & (f_grounded.dat.data_ro > 0.5)).sum())))
    r.update(ng=ng, nf=nf, g_wrong=g_but_lt_half, f_wrong=f_but_gt_half, gmask=gmask)
    r["b"], r["H"] = b, H
    res[method] = r
    PETSc.Sys.Print(
        f"  {method:12s}  rms|jump| all: b {r['rms_jb']:7.2f}  s {r['rms_js']:7.2f}  h {r['rms_jh']:7.2f} m"
        f"  | grounded: b {r['rms_jb_g']:7.2f}  s {r['rms_js_g']:7.2f} m"
        f"  | front<h> {r['front_h']:6.1f} m  vol {r['vol_km3']:.4e} km3 "
        f"({100 * (r['vol_km3'] / vol_raster - 1):+.3f}% of the raster)  [{r['t']:.1f}s]")
    PETSc.Sys.Print(
        f"  {'':12s}  INTERIOR (facets between two >10 km cells, {n_big} such cells): "
        f"rms|jump| b {r['rms_jb_i']:7.2f}  s {r['rms_js_i']:7.2f}  h {r['rms_jh_i']:7.2f} m")
    PETSc.Sys.Print(
        f"  {'':12s}  flotation: grounded {r['ng']} floating {r['nf']} cells of {n_cells}; "
        f"grounded-but-BedMachine<50% {r['g_wrong']}, floating-but-BedMachine>50% {r['f_wrong']}")

if not {"vertex", "cell_mean", "patch"} <= set(res):
    sys.exit(0)
db = Function(Q).assign(res["cell_mean"]["b"] - res["vertex"]["b"])
area = assemble(one * dx)
rms_db = math.sqrt(assemble(db ** 2 * dx) / area)
A_big = assemble(big * dx)
rms_db_i = math.sqrt(assemble(big * db ** 2 * dx) / A_big) if A_big > 0 else float("nan")
rms_db_m = math.sqrt(assemble((1 - big) * db ** 2 * dx) / max(area - A_big, 1.0))
mx = mesh.comm.allreduce(float(np.abs(db.dat.data_ro).max()) if db.dat.data_ro.size else 0.0, op=max)
flip = int(mesh.comm.allreduce(float((res["vertex"]["gmask"].dat.data_ro != res["cell_mean"]["gmask"].dat.data_ro).sum())))
PETSc.Sys.Print(f"  cells whose grounded flag differs between the two samplings: {flip}")
for k, r in res.items():
    r.pop("b"); r.pop("H"); r.pop("s"); r.pop("gmask")
v, c = res["vertex"], res["cell_mean"]
PETSc.Sys.Print(f"  cell_mean - vertex bed: rms {rms_db:.2f} m (interior cells {rms_db_i:.2f}, margin cells {rms_db_m:.2f}), max |diff| {mx:.1f} m")
PETSc.Sys.Print(f"  grounded rms jump in s: {v['rms_js_g']:.2f} -> {c['rms_js_g']:.2f} m "
                f"({100*(c['rms_js_g']/v['rms_js_g']-1):+.1f}%)")
p = res["patch"]
PETSc.Sys.Print(f"  patch vs vertex: grounded rms jump in s {v['rms_js_g']:.2f} -> {p['rms_js_g']:.2f} m "
                f"({100*(p['rms_js_g']/v['rms_js_g']-1):+.1f}%), b {v['rms_jb_g']:.2f} -> {p['rms_jb_g']:.2f} m; "
                f"volume {v['vol_km3']:.4e} -> {p['vol_km3']:.4e} km3 (raster {vol_raster:.4e})")
