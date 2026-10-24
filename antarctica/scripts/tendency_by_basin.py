#!/usr/bin/env python3
r"""Where a historical's initial mass tendency comes from, per basin.

A forward under the apparent mass balance carries ``a_ref_mb``, the field
that cancelled its t = 0 thickness tendency: a_ref = div(uH) - (SMB - melt)
on the start geometry, with the run's own upwind flux divergence. This
reads a state, rebuilds that divergence from the thickness and the velocity
the state holds (the same facet-upwind rule the transport uses), and from
the two recovers the forcing SMB - melt the run saw. It then evaluates the
divergence under the observed velocity the state carries too, so each
Mouginot basin gets:

- SMB - melt (the forcing),
- div(u H) under the model velocity and under the observed velocity,
- the imbalance SMB - melt - div under each, in Gt/yr over grounded ice,

against the Smith ICESat mass rate of the same basin where a
``compare_hist_dhdt.py`` JSON is given. A basin whose imbalance is positive
under the observed velocity too has too much forcing (or a bed that
cannot drain it); one positive only under the model velocity has a model
that drains too little. The observed velocity is zero where unobserved, so
a basin with much unobserved ice carries a biased observed divergence; the
unobserved share of each basin's area is printed.

    mpiexec -n 4 python antarctica/scripts/tendency_by_basin.py --kit KIT.nc \
        --smith dhdt.json --out tendency.json results/hist_cesm2_waccm_lcr75_2000_final.h5
"""
import argparse
import json
import os
import sys

import numpy as np

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_SCRIPTS)))
from icepack2_tools.runconfig import obs_data_root  # noqa: E402

DEFAULT_KIT = os.path.join(obs_data_root(), "..", "..", "ISMIP7", "AIS", "obs", "mipkit",
                           "AntarcticaObsISMIP7-v1.2.nc")
RHO_I, RHO_W = 917.0, 1024.0
GT_PER_M3 = RHO_I / 1.0e12
H_MIN = 1.0
# Mouginot and Rignot (2017) basins as the kit numbers them, inferred from
# their geography (11 = G-H, the Amundsen; 10 = F-G, Getz; 9 = Ep-F, Siple)
BASIN_NAMES = {1: "A-Ap", 2: "Ap-B", 3: "B-C", 4: "C-Cp", 5: "Cp-D", 6: "D-Dp", 7: "Dp-E", 8: "E-Ep",
               9: "Ep-F", 10: "F-G", 11: "G-H", 12: "H-Hp", 13: "Hp-I", 14: "I-Ipp", 15: "Ipp-J",
               16: "J-Jpp", 17: "Jpp-K", 18: "K-A"}


def chk_thickness(path, mesh):
    import firedrake as fd
    with fd.CheckpointFile(path, "r") as chk:
        return chk.load_function(mesh, name="thickness")


def divergence(mesh, H, u):
    r"""Per-cell net outflow [m^3/yr] of the facet-upwind flux of ``H`` by ``u``."""
    from firedrake import (FacetNormal, TestFunction, assemble, conditional, dS, ds, dot, gt)
    v = TestFunction(H.function_space())
    n = FacetNormal(mesh)
    un = dot(u("+"), n("+"))
    q = conditional(gt(un, 0.0), H("+"), H("-")) * un
    unb = dot(u, n)
    form = q * (v("+") - v("-")) * dS + conditional(gt(unb, 0.0), H * unb, 0.0) * v * ds
    return assemble(form).dat.data_ro.copy()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("state")
    ap.add_argument("--kit", default=DEFAULT_KIT)
    ap.add_argument("--smith", default=None, help="compare_hist_dhdt.py JSON for the observed basin rates")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    import firedrake as fd
    from firedrake.petsc import PETSc
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    with fd.CheckpointFile(args.state, "r") as chk:
        mesh = chk.load_mesh()
        H0 = chk.load_function(mesh, name="H_init")
        b = chk.load_function(mesh, name="bed")
        u = chk.load_function(mesh, name="velocity")
        uo = chk.load_function(mesh, name="velocity_obs")
        try:
            aref = chk.load_function(mesh, name="a_ref_mb")
        except Exception:
            aref = None   # a run without the reference: the divergences alone
    Q0 = H0.function_space()
    area = fd.assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro
    # the reference was built on the start geometry, so the divergences are
    # evaluated there too
    # the start geometry for a state carrying the reference (it was built
    # there); the evolved geometry otherwise, so the divergence is the flux
    # the state's own velocity carries through its own thickness
    H_div = H0 if aref is not None else chk_thickness(args.state, mesh)
    div_m = divergence(mesh, H_div, u)
    div_o = divergence(mesh, H_div, uo)
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro
    speed_o = fd.Function(Q0).project(fd.sqrt(fd.dot(uo, uo))).dat.data_ro
    hv, bv = H0.dat.data_ro, b.dat.data_ro
    haf = hv + np.where(bv < 0.0, bv, 0.0) * RHO_W / RHO_I
    grounded = (hv > H_MIN) & (haf > 0.0)
    aref_v = aref.dat.data_ro * area if aref is not None else np.zeros_like(area)
    parts = comm.gather((div_m, div_o, aref_v, area.copy(), grounded.copy(),
                         np.asarray(xy).copy(), (speed_o < 0.5).copy()), root=0)
    if comm.rank != 0:
        return
    div_m, div_o, aref_v, ar, gr, pts, unobs = (np.concatenate([p[i] for p in parts]) for i in range(7))
    import netCDF4
    kit = netCDF4.Dataset(args.kit)
    x = np.asarray(kit["x"][:], float)
    y = np.asarray(kit["y"][:], float)
    basins = np.asarray(kit["mouginot_basins"][:])
    ix = np.clip(np.rint((pts[:, 0] - x[0]) / (x[1] - x[0])).astype(int), 0, len(x) - 1)
    iy = np.clip(np.rint((pts[:, 1] - y[0]) / (y[1] - y[0])).astype(int), 0, len(y) - 1)
    bs = basins[iy, ix]
    smb_melt = div_m - aref_v            # m^3/yr per cell, the forcing the run saw
    smith = {}
    if args.smith:
        with open(args.smith) as f:
            d = json.load(f)[0]
        smith = {k: v.get("obs_gtyr") for k, v in d["basins"].items()}
    rows = []
    total = {"smb_melt": 0.0, "div_model": 0.0, "div_obs": 0.0, "imb_model": 0.0, "imb_obs": 0.0, "aref": 0.0}
    for bsn in sorted(np.unique(bs[gr])):
        m = gr & (bs == bsn)
        r = {"basin": int(bsn), "name": BASIN_NAMES.get(int(bsn), "?"),
             "area_km2": float(ar[m].sum() / 1e6),
             "unobserved_share": float(ar[m & unobs].sum() / ar[m].sum()),
             "smb_melt": float(smb_melt[m].sum() * GT_PER_M3),
             "div_model": float(div_m[m].sum() * GT_PER_M3),
             "div_obs": float(div_o[m].sum() * GT_PER_M3),
             "aref": float(aref_v[m].sum() * GT_PER_M3)}
        r["imb_model"] = r["smb_melt"] - r["div_model"]
        r["imb_obs"] = r["smb_melt"] - r["div_obs"]
        r["smith"] = smith.get(str(int(bsn)))
        rows.append(r)
        for k in total:
            total[k] += r[k]
    PETSc.Sys.Print(f"grounded ice, Gt/yr; a_ref total {total['aref']:+,.0f} (the budget's amb column at t=0)")
    PETSc.Sys.Print("basin  name   area_Mkm2 unobs  SMB-melt  div_model  div_obs  imb_model  imb_obs  Smith")
    for r in rows:
        PETSc.Sys.Print(f"{r['basin']:5d} {r['name']:>6s} {r['area_km2'] / 1e6:9.3f} {r['unobserved_share']:5.2f} "
                        f"{r['smb_melt']:9.0f} {r['div_model']:10.0f} {r['div_obs']:8.0f} {r['imb_model']:10.0f} "
                        f"{r['imb_obs']:8.0f} {('%+.0f' % r['smith']) if r['smith'] is not None else '-':>6s}")
    PETSc.Sys.Print(f"total {'':>6s} {'':>9s} {'':>5s} {total['smb_melt']:9.0f} {total['div_model']:10.0f} "
                    f"{total['div_obs']:8.0f} {total['imb_model']:10.0f} {total['imb_obs']:8.0f}")
    if args.out:
        with open(args.out, "w") as f:
            json.dump({"state": os.path.basename(args.state), "basins": rows, "total": total}, f, indent=1)


if __name__ == "__main__":
    main()
