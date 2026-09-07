#!/usr/bin/env python3
r"""Model dH/dt vs observed (dhdt_smith) for one or more MAP checkpoints.

The transient inversion's payoff diagnostic: for each MAP, rebuild the SAME
one-step implicit-Euler DG0 upwind prognostic used by the inversion's dH/dt
term (source = RACMO SMB - OI-climatology per-basin-K melt, both at the
reference geometry) and score the resulting thickness tendency against the
observed mean map. Because the t=0 *velocity* misfit cannot distinguish a
velocity-only MAP from a transient one (both fit u), this tendency comparison
is the observable that can.

Usage:
    python antarctica/scripts/compare_dhdt.py LABEL=PATH [LABEL=PATH ...]

Environment: ISMIP7_DHDT_* and ISMIP7_K_* knobs as in the inversion.
Writes a PNG (obs / model / difference per MAP) to antarctica/figs/ and
prints a scoreboard. Read-only with respect to the MAPs.
"""

import os
import sys

import numpy as np
import firedrake as fd
from firedrake import (
    Constant, Function, FunctionSpace, VectorFunctionSpace, TestFunction,
    FacetNormal, assemble, dx, dS, ds, dot, jump,
)
from firedrake.petsc import PETSc

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.obs_dhdt import load_dhdt_obs
from icepack2_tools.runconfig import lc as _lc
from icepack2_tools.forcing import (
    load_racmo_smb_climatology, load_K_per_basin,
    make_climatology_ocean_callback,
)

DT = float(os.environ.get("ISMIP7_DHDT_DT", "1.0"))
CLIM0 = int(os.environ.get("ISMIP7_DHDT_CLIM_START", "2003"))
CLIM1 = int(os.environ.get("ISMIP7_DHDT_CLIM_END", "2019"))
RHO_GT = 917.0 / 1e12


def one_step_dhdt(path):
    r"""(mesh, dhdt_model, dhdt_obs, mask, cell_area, haf) for one MAP."""
    with fd.CheckpointFile(path, "r") as c:
        mesh = c.load_mesh()
        u0 = c.load_function(mesh, "velocity")
        H = c.load_function(mesh, "thickness")
        b = c.load_function(mesh, "bed")
        try:
            s = c.load_function(mesh, "surface")
        except Exception:
            rr = 917.0 / 1024.0
            s = Function(H.function_space()).interpolate(
                fd.max_value(b + H, (Constant(1.0) - rr) * H))

    Q_g = H.function_space()
    if Q_g.ufl_element().degree() != 0:
        raise SystemExit(f"{path}: geometry is not DG0; this diagnostic "
                         f"replicates the DG0 transport and would not be "
                         f"comparing like with like.")
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)

    dhdt_obs, cov = load_dhdt_obs(Q_g)
    smb = load_racmo_smb_climatology(Q_g, clim_start=CLIM0, clim_end=CLIM1)

    melt = Function(Q_g, name="ocean_melt_ref")
    if os.environ.get("ISMIP7_DHDT_MELT", "1") != "0":
        # Resolution comes from icepack2_tools.runconfig, the same owner the
        # inversion reads: a mismatch here picks a per-basin K calibrated at
        # another resolution, which perturbs the melt source and hence the
        # model tendency being scored.
        _k_lc = os.path.join(_ROOT, "results",
                             f"calibrated_K_per_basin_{_lc()}.npz")
        _k_2500 = os.path.join(_ROOT, "results",
                               "calibrated_K_per_basin_2500.npz")
        k_npz = os.environ.get("ISMIP7_K_PER_BASIN_NPZ",
                               _k_lc if os.path.exists(_k_lc) else _k_2500)
        W2 = VectorFunctionSpace(mesh, "DG", 0)
        xy = Function(W2).interpolate(
            fd.SpatialCoordinate(mesh)).dat.data_ro.reshape(-1, 2)
        K = load_K_per_basin(k_npz, xy[:, 0].copy(), xy[:, 1].copy(), fill=0.0)
        K = K * float(os.environ.get("ISMIP7_K_SCALE", "1.0"))
        ctx = {"mesh": mesh, "Q": Q, "V": V, "Q_g": Q_g,
               "geom_xy": (xy[:, 0].copy(), xy[:, 1].copy()),
               "h": H, "b": b, "s": s, "ocean_melt": melt}
        make_climatology_ocean_callback(K)(ctx, 0.0)

    h_next = Function(Q_g)
    w = TestFunction(Q_g)
    nrm = FacetNormal(mesh)
    un = dot(u0, nrm)
    unp = (un + abs(un)) / 2
    dt_c = Constant(DT)
    F_h = ((h_next - H) / dt_c * w * dx
           + (unp("+") * h_next("+") - unp("-") * h_next("-")) * jump(w) * dS
           + unp * h_next * w * ds
           - (smb - melt) * w * dx)
    fd.solve(F_h == 0, h_next, solver_parameters={
        "snes_type": "ksponly", "ksp_type": "gmres",
        "pc_type": "bjacobi", "sub_pc_type": "ilu", "ksp_rtol": 1e-10})

    dhdt_model = Function(Q_g, name="dhdt_model")
    dhdt_model.assign((h_next - H) / dt_c)

    area = assemble(TestFunction(Q_g) * dx).dat.data_ro.copy()
    hv, bv = H.dat.data_ro, b.dat.data_ro
    haf = hv + np.minimum(bv, 0.0) * (1024.0 / 917.0)
    return mesh, dhdt_model, dhdt_obs, cov, area, haf, hv


def score(label, mesh, dhdt_model, dhdt_obs, cov, area, haf, hv):
    r"""Scoreboard for one MAP.

    Every reported number is a GLOBAL reduction. Each rank owns only its own
    slice of the DG0 data, so a rank-local sum would print N different partial
    ice sheets under ``mpiexec -n N`` -- and this is the diagnostic that
    validation decisions are read off, so a plausible-looking partial integral
    is worse than no number at all.
    """
    comm = mesh.comm

    def gsum(a):
        return float(comm.allreduce(float(np.sum(a))))

    m = dhdt_model.dat.data_ro
    o = dhdt_obs.dat.data_ro
    sel = (cov.dat.data_ro > 0.5) & (haf > 0.0) & (hv > 1.0)
    w = area[sel]
    mm, oo = m[sel], o[sel]
    n_sel = int(comm.allreduce(int(sel.sum())))
    w_tot = gsum(w)
    gt = lambda f: gsum(f * w) * RHO_GT
    rms = float(np.sqrt(gsum(w * (mm - oo) ** 2) / w_tot))
    # area-weighted correlation
    wm, wo = gsum(w * mm) / w_tot, gsum(w * oo) / w_tot
    r = float(gsum(w * (mm - wm) * (oo - wo))
              / np.sqrt(gsum(w * (mm - wm) ** 2)
                        * gsum(w * (oo - wo) ** 2)))
    PETSc.Sys.Print(f"\n===== {label} (grounded+observed, n={n_sel}) =====")
    PETSc.Sys.Print(
        f"  integrated dH/dt   model {gt(mm):+8.1f}  obs {gt(oo):+8.1f} Gt/yr")
    PETSc.Sys.Print(f"  RMS(model - obs)   {rms:8.3f} m/yr")
    PETSc.Sys.Print(f"  area-wtd corr      {r:8.3f}")
    for lo, hi, lab in [(0, 50, "GL band 0-50 m"), (50, 200, "50-200 m"),
                        (200, 1e9, "interior >200 m")]:
        s2 = sel & (haf > lo) & (haf <= hi)
        w2 = area[s2]
        w2_tot = gsum(w2)
        if w2_tot <= 0.0:
            continue
        rms2 = float(np.sqrt(gsum(w2 * (m[s2] - o[s2]) ** 2) / w2_tot))
        PETSc.Sys.Print(
            f"    {lab:<18} model {gsum(m[s2]*w2)*RHO_GT:+8.1f} "
            f"obs {gsum(o[s2]*w2)*RHO_GT:+8.1f} Gt/yr   "
            f"RMS {rms2:6.3f} m/yr")
    return dict(rms=rms, corr=r, model_gt=gt(mm), obs_gt=gt(oo))


def main():
    runs = []
    for arg in sys.argv[1:]:
        if "=" not in arg:
            raise SystemExit(f"expected LABEL=PATH, got {arg!r}")
        runs.append(tuple(arg.split("=", 1)))
    if not runs:
        raise SystemExit(__doc__)

    panels = []
    for label, path in runs:
        PETSc.Sys.Print(f"[{label}] {path}")
        mesh, dm, do, cov, area, haf, hv = one_step_dhdt(path)
        score(label, mesh, dm, do, cov, area, haf, hv)
        panels.append((label, mesh, dm, do, cov))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from firedrake.pyplot import tripcolor
        n = len(panels)
        fig, axes = plt.subplots(n, 3, figsize=(15, 4.6 * n), squeeze=False)
        for i, (label, mesh, dm, do, cov) in enumerate(panels):
            diff = Function(dm.function_space()).assign(dm - do)
            msk = cov.dat.data_ro <= 0.5
            for f in (do, dm, diff):
                f.dat.data[msk] = np.nan
            for j, (f, t, lim) in enumerate([(do, "observed (Smith)", 2),
                                             (dm, f"model ({label})", 2),
                                             (diff, "model - obs", 2)]):
                ax = axes[i][j]
                c = tripcolor(f, axes=ax, cmap="RdBu_r", vmin=-lim, vmax=lim)
                ax.set_title(f"{t} [m/yr]")
                ax.set_aspect("equal")
                ax.set_xticks([]); ax.set_yticks([])
                fig.colorbar(c, ax=ax, shrink=0.75)
        out = os.path.join(_ROOT, "figs", "dhdt_compare.png")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        fig.savefig(out, dpi=170, bbox_inches="tight")
        PETSc.Sys.Print(f"\nSaved figure: {out}")
    except Exception as e:
        PETSc.Sys.Print(f"(figure skipped: {e})")


if __name__ == "__main__":
    main()
