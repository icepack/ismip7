#!/usr/bin/env python3
r"""Ice flux through the mesh boundary of a state on a mesh without a buffer.

On Rice's meshes the front is the mesh boundary, so the calving flux the
budget books as ``outflux_gtyr`` is the flux through the exterior facets
that hold ice. This reads a MAP or a forward state and assembles that flux
under its own velocity and under the observed velocity it carries, over all
ice, and split into floating and grounded ice by the height above flotation,
with the front's ice-weighted mean thickness and normal speed. Read-only,
no solve; every number is an assembled integral, so any rank count reads
any file. The observed front, BedMachine's 500 m front under MEaSUREs,
carries 1,138 to 1,297 Gt/yr (icepack/ismip7#167).

    mpiexec -n 4 python antarctica/scripts/boundary_front_flux.py \
        --out flux.json map=MAP.h5 t2010=RUN_2000_t2010.0.h5 ...
"""
import argparse
import json
import os
import sys

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_SCRIPTS)))

import firedrake as fd  # noqa: E402
from firedrake import Constant, FacetNormal, assemble, conditional, dot, ds, gt, lt, max_value  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402

RHO_I, RHO_W = 917.0, 1024.0
GT_PER_M3 = RHO_I / 1.0e12
H_MIN = 1.0


def load(path):
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        f = {}
        for name in ("thickness", "bed", "velocity", "velocity_obs"):
            try:
                f[name] = chk.load_function(mesh, name=name)
            except Exception:
                f[name] = None
        t = float(chk.get_attr("/", "t_yr")) if chk.has_attr("/", "t_yr") else None
    return mesh, f, t


def flux_terms(mesh, h, b, u):
    r"""(total, floating, grounded) flux [m^3/yr], front length [m] and
    ice-weighted thickness and normal speed integrals on the exterior
    facets holding ice."""
    n = FacetNormal(mesh)
    ice = conditional(gt(h, H_MIN), 1.0, 0.0)
    haf = h + min_haf(b)
    floating = conditional(lt(haf, 0.0), 1.0, 0.0)
    un = max_value(dot(u, n), Constant(0.0))
    q = ice * h * un
    out = {
        "flux": assemble(q * ds), "flux_floating": assemble(floating * q * ds),
        "flux_grounded": assemble((1.0 - floating) * q * ds),
        "length": assemble(ice * ds), "h_int": assemble(ice * h * ds), "un_int": assemble(ice * un * ds),
        "length_floating": assemble(ice * floating * ds),
        "un_floating_int": assemble(ice * floating * un * ds),
    }
    return {k: float(v) for k, v in out.items()}


def min_haf(b):
    # height above flotation h_af = h + min(b, 0) rho_w / rho_i: ice floats
    # where h < -b rho_w / rho_i
    return conditional(lt(b, 0.0), b, Constant(0.0)) * Constant(RHO_W / RHO_I)


def check(label, path):
    mesh, f, t = load(path)
    h, b = f["thickness"], f["bed"]
    res = {"label": label, "file": os.path.basename(path), "t_yr": t}
    for key, u in (("model", f["velocity"]), ("obs", f["velocity_obs"])):
        if u is None:
            res[key] = None
            continue
        terms = flux_terms(mesh, h, b, u)
        res[key] = {
            "flux_gtyr": terms["flux"] * GT_PER_M3,
            "floating_gtyr": terms["flux_floating"] * GT_PER_M3,
            "grounded_gtyr": terms["flux_grounded"] * GT_PER_M3,
            "front_km": terms["length"] / 1e3,
            "mean_h_m": terms["h_int"] / terms["length"] if terms["length"] else None,
            "mean_un_myr": terms["un_int"] / terms["length"] if terms["length"] else None,
            "floating_front_km": terms["length_floating"] / 1e3,
            "floating_un_myr": (terms["un_floating_int"] / terms["length_floating"]
                                if terms["length_floating"] else None),
        }
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("states", nargs="+", help="label=path")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    results = []
    for item in args.states:
        label, path = item.split("=", 1) if "=" in item else (os.path.basename(item), item)
        if not os.path.exists(path):
            sys.exit(f"no such state: {path}")
        r = check(label, path)
        results.append(r)
        for key in ("model", "obs"):
            v = r[key]
            if v is None:
                PETSc.Sys.Print(f"  {label}: no {key} velocity")
                continue
            PETSc.Sys.Print(
                f"  {label} ({key}): front flux {v['flux_gtyr']:,.0f} Gt/yr (floating {v['floating_gtyr']:,.0f}, "
                f"grounded {v['grounded_gtyr']:,.0f}); front {v['front_km']:,.0f} km, mean h {v['mean_h_m']:.0f} m, "
                f"mean u.n {v['mean_un_myr']:.0f} m/yr (floating {v['floating_un_myr'] or 0:.0f} m/yr)")
        if r["model"] and r["obs"] and r["obs"]["flux_gtyr"]:
            PETSc.Sys.Print(f"  {label}: model / observed-velocity flux = "
                            f"{r['model']['flux_gtyr'] / r['obs']['flux_gtyr']:.2f}")
    if args.out and fd.COMM_WORLD.rank == 0:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=1)


if __name__ == "__main__":
    main()
