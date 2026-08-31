#!/usr/bin/env python3
r"""Physical (thermomechanical) fluidity prior mean A_prior(x) for the ISMIP7
Antarctic inversion, following the mismip_time-dependent-da approach: the
inversion control becomes phi = log(A / A_prior) - a small deviation from a
physically-motivated fluidity - instead of log(A / A0_const), so the prior
regularizes the DEVIATION, not the amplitude (Recinos/fenics_ice framing).

This standalone driver computes and PLOTS A_prior for inspection. The same
`compute_fluidity_prior` is called inside inversion_icepack2.py at setup so the
inversion's prior mean is identical to what this script visualizes.

A_prior comes from a fixed-velocity depth-averaged enthalpy (Stefan) solve
(icepack2_tools/thermo_model.py, ported from the mismip thermo_model) driven by
the OBSERVED geometry/velocity and a mean-annual surface-temperature field.

Usage:
    OMP_NUM_THREADS=1 ISMIP7_LC=32000 \
      ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_320000_32000.msh \
      python antarctica/scripts/thermo_prior.py [--shear_amp 30] [--q_geo 50]
"""
import os
import sys
import argparse

import firedrake as fd
from firedrake import Constant, FunctionSpace, VectorFunctionSpace, Function, max_value
from firedrake.petsc import PETSc
import rasterio
import icepack

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(_ROOT))          # repo root
sys.path.insert(0, os.path.join(_ROOT, "scripts"))  # for simulation helpers if needed
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")

from icepack2_tools.thermo_model import compute_fluidity_prior, DEFAULTS
from icepack2_tools.dual_friction import weertman_anchor
from icepack2_tools.forcing import (load_racmo_smb_climatology,
                                    load_mean_annual_surface_temperature)

lc = int(os.environ.get("ISMIP7_LC", "32000"))
lc_coarse = int(os.environ.get("ISMIP7_LC_COARSE", str(lc * 10)))


def find_file(d, p):
    import glob
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def main():
    ap = argparse.ArgumentParser()
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k}", type=float, default=v)
    ap.add_argument("--max_picard", type=int, default=30)
    ap.add_argument("--tvar", default="tas", choices=["tas", "ts"])
    ap.add_argument("--outfile",
                    default=os.path.join(MESH_DIR, f"fluidity_prior_{lc}.h5"))
    ap.add_argument("--figfile",
                    default=os.path.join(_ROOT, "figs", f"fluidity_prior_{lc}.png"))
    args = ap.parse_args()
    P = {k: getattr(args, k) for k in DEFAULTS}

    mesh_fn = os.environ.get(
        "ISMIP7_MESH", os.path.join(MESH_DIR, f"antarctica_{lc_coarse}_{lc}.msh"))
    PETSc.Sys.Print(f"Loading mesh: {mesh_fn}")
    mesh = fd.Mesh(mesh_fn)
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)

    bm = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm}:bed"), Q)
    H = Function(Q, name="thickness").interpolate(
        max_value(icepack.interpolate(rasterio.open(f"netcdf:{bm}:thickness"), Q),
                  Constant(0.0)))
    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q, name="surface").interpolate(
        max_value(b + H, (Constant(1.0) - rho_ratio) * H))

    vel = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel}:VX"), rasterio.open(f"netcdf:{vel}:VY")),
        V, fillvalue=0.0)

    m_slide = float(os.environ.get("ISMIP7_M_SLIDE", "3.0"))
    C = weertman_anchor(H, s, u_obs, m_slide, Q)                # balance friction
    acc = load_racmo_smb_climatology(Q)
    T_srf = load_mean_annual_surface_temperature(Q, var=args.tvar)
    PETSc.Sys.Print(f"  T_srf ({args.tvar}) [{T_srf.dat.data_ro.min():.1f}, "
                    f"{T_srf.dat.data_ro.max():.1f}] K; "
                    f"acc [{acc.dat.data_ro.min():.2f}, {acc.dat.data_ro.max():.2f}] m/yr")

    A_prior = compute_fluidity_prior(u_obs, H, s, b, C, acc, T_srf, p=P,
                                     max_picard=args.max_picard)

    os.makedirs(os.path.dirname(args.outfile), exist_ok=True)
    with fd.CheckpointFile(args.outfile, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(A_prior, name="fluidity_prior")
        chk.save_function(T_srf, name="T_srf")
    PETSc.Sys.Print(f"Saved fluidity prior: {args.outfile}")

    if mesh.comm.size == 1:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import LogNorm
            os.makedirs(os.path.dirname(args.figfile), exist_ok=True)
            fig, axes = plt.subplots(1, 2, figsize=(13, 6))
            Alo = max(1.0, float(A_prior.dat.data_ro.min()))
            Ahi = float(A_prior.dat.data_ro.max())
            for ax, (f, ttl, kw) in zip(axes, [
                    (A_prior, "thermo fluidity prior A_prior", dict(norm=LogNorm(Alo, Ahi), cmap="magma")),
                    (T_srf, "mean-annual surface T (K)", dict(cmap="inferno"))]):
                col = fd.tripcolor(f, axes=ax, **kw)
                fig.colorbar(col, ax=ax, shrink=0.8)
                ax.set_title(ttl)
                ax.set_aspect("equal")
                ax.axis("off")
            fig.tight_layout()
            fig.savefig(args.figfile, dpi=150, bbox_inches="tight")
            PETSc.Sys.Print(f"Saved figure: {args.figfile}")
        except Exception as e:
            PETSc.Sys.Print(f"  (figure skipped: {e})")


if __name__ == "__main__":
    main()
