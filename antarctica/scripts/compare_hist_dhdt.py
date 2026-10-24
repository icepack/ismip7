#!/usr/bin/env python3
r"""A historical run's thickness change against the Smith et al. ICESat dH/dt.

Read-only. For each state named, the modelled rate is the thickness change
since the run's start (``thickness`` minus the ``H_init`` the checkpoint
carries) over the elapsed time, put on the observations kit's 1 km grid by
nearest cell centre, and compared with ``dhdt_smith`` (ICESat to ICESat-2,
2003 to 2019, m/yr ice equivalent) where the record has data. Statistics
are over grounded ice (height above flotation > 0 in the state) and over
all ice, pixel-weighted: bias, rms, correlation, and the mass rate each
implies; then per Mouginot basin the model and observed mean rates and
mass rates, so a gain the record does not show can be placed. ``--plot``
draws observed, modelled and their difference per state.

    mpiexec -n 4 python antarctica/scripts/compare_hist_dhdt.py --out dhdt.json --plot dhdt.png \
        amb_on=results/hist_cesm2_waccm_lcr75_2000_final.h5 amb_off=results/hist_cesm2_waccm_lcr75noamb_2000_final.h5
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


def load_state(path):
    import firedrake as fd
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        h = chk.load_function(mesh, name="thickness")
        h0 = chk.load_function(mesh, name="H_init")
        b = chk.load_function(mesh, name="bed")
        t = float(chk.get_attr("/", "t_yr"))
        t0 = float(chk.get_attr("/", "t_start_yr")) if chk.has_attr("/", "t_start_yr") else None
    Q0 = h.function_space()
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro
    area = fd.assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro
    return mesh, h, h0, b, t, t0, np.asarray(xy), np.asarray(area)


def block_mean(model, obs, grounded, basins, n):
    r"""Block means over n x n pixels of the two fields, over the pixels where
    both are finite; a block is grounded where most of its model pixels are,
    and takes the basin of its centre pixel."""
    ny, nx = model.shape
    ny2, nx2 = ny // n, nx // n
    m = model[:ny2 * n, :nx2 * n].reshape(ny2, n, nx2, n)
    o = obs[:ny2 * n, :nx2 * n].reshape(ny2, n, nx2, n)
    g = grounded[:ny2 * n, :nx2 * n].reshape(ny2, n, nx2, n)
    b = basins[:ny2 * n, :nx2 * n].reshape(ny2, n, nx2, n)
    both = np.isfinite(m) & np.isfinite(o)
    cnt = both.sum(axis=(1, 3))
    with np.errstate(invalid="ignore", divide="ignore"):
        mm = np.where(both, m, 0.0).sum(axis=(1, 3)) / cnt
        om = np.where(both, o, 0.0).sum(axis=(1, 3)) / cnt
    ok = cnt >= max(1, (n * n) // 4)
    mm = np.where(ok, mm, np.nan)
    om = np.where(ok, om, np.nan)
    gb = (g & both).sum(axis=(1, 3)) * 2 > np.maximum(cnt, 1)
    bb = b[:, n // 2, :, n // 2]
    return mm, om, gb & ok, bb


def stats(model, obs, w):
    m = np.isfinite(model) & np.isfinite(obs)
    if m.sum() < 10:
        return {"n": int(m.sum())}
    d = model[m] - obs[m]
    return {"n": int(m.sum()), "bias_myr": float(d.mean()), "rms_myr": float(np.sqrt((d ** 2).mean())),
            "corr": float(np.corrcoef(model[m], obs[m])[0, 1]),
            "model_gtyr": float(model[m].sum() * w * GT_PER_M3), "obs_gtyr": float(obs[m].sum() * w * GT_PER_M3),
            "model_mean_myr": float(model[m].mean()), "obs_mean_myr": float(obs[m].mean())}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("states", nargs="+", help="label=path")
    ap.add_argument("--kit", default=DEFAULT_KIT)
    ap.add_argument("--t0", type=float, default=None, help="start year if the state does not record it (2003 for the historicals)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--plot", default=None, help="PNG: observed, modelled, difference per state")
    ap.add_argument("--block", type=float, default=0.0,
                    help="average both fields over blocks of this size [km] before the statistics and "
                         "the plot, so a cell-scale flux-divergence pattern does not hide the basin signal")
    args = ap.parse_args()
    from mpi4py import MPI
    from firedrake.petsc import PETSc
    comm = MPI.COMM_WORLD
    kit = x = y = obs = basins = None
    if comm.rank == 0:
        import netCDF4
        kit = netCDF4.Dataset(args.kit)
        x = np.asarray(kit["x1km"][:], float)
        y = np.asarray(kit["y1km"][:], float)
        obs = np.asarray(kit["dhdt_smith"][:], float)
        obs = np.where(np.isfinite(obs) & (np.abs(obs) < 1e3), obs, np.nan)
        basins = np.asarray(kit["mouginot_basins"][::2, ::2])
    results, panels = [], []
    for item in args.states:
        label, path = item.split("=", 1) if "=" in item else (os.path.basename(item), item)
        if not os.path.exists(path):
            sys.exit(f"no such state: {path}")
        mesh, h, h0, b, t, t0, xy, area = load_state(path)
        if t0 is None:
            t0 = args.t0
        if t0 is None:
            sys.exit(f"{path} records no start year; pass --t0")
        rate = (h.dat.data_ro - h0.dat.data_ro) / (t - t0)
        hv, bv = h.dat.data_ro, b.dat.data_ro
        haf = hv + np.where(bv < 0.0, bv, 0.0) * RHO_W / RHO_I
        grounded = (hv > H_MIN) & (haf > 0.0)
        ice = hv > H_MIN
        h0v = h0.dat.data_ro
        grounded0 = (h0v > H_MIN) & (h0v + np.where(bv < 0.0, bv, 0.0) * RHO_W / RHO_I > 0.0)
        parts = comm.gather((rate.copy(), grounded.copy(), ice.copy(), xy.copy(), area.copy(), grounded0.copy()), root=0)
        if comm.rank != 0:
            continue
        rate, grounded, ice, pts, ar, grounded0 = (np.concatenate([p[i] for p in parts]) for i in range(6))
        gl_area0, gl_area1 = float(ar[grounded0].sum() / 1e6), float(ar[grounded].sum() / 1e6)
        gl_lost, gl_gained = float(ar[grounded0 & ~grounded].sum() / 1e6), float(ar[~grounded0 & grounded].sum() / 1e6)
        from scipy.spatial import cKDTree
        X, Y = np.meshgrid(x, y)
        d, idx = cKDTree(pts).query(np.column_stack([X.ravel(), Y.ravel()]), workers=-1)
        inside = (d <= 1.5 * np.sqrt(ar[idx])).reshape(X.shape)
        idx = idx.reshape(X.shape)
        model = np.where(inside & ice[idx], rate[idx], np.nan)
        gr = inside & grounded[idx]
        fl = inside & ice[idx] & ~grounded[idx]
        w = abs(float((x[1] - x[0]) * (y[1] - y[0])))
        obs_s, basins_s = obs, basins
        floating_stats = stats(np.where(fl, model, np.nan), obs, w)
        if args.block > 0:
            # block means of both fields over the pixels where both are
            # finite; a block is grounded where most of its model pixels are
            n = max(1, int(round(args.block * 1e3 / abs(x[1] - x[0]))))
            model, obs_s, gr, basins_s = block_mean(model, obs, gr, basins, n)
            w *= n * n
        res = {"label": label, "file": os.path.basename(path), "t_yr": t, "t0": t0, "block_km": args.block,
               "grounded": stats(np.where(gr, model, np.nan), obs_s, w),
               "all_ice": stats(model, obs_s, w), "floating": floating_stats,
               "grounded_area_km2": {"start": gl_area0, "end": gl_area1, "ungrounded": gl_lost, "grounded": gl_gained},
               "basins": {}}
        for bsn in np.unique(basins_s[gr]):
            sel = gr & (basins_s == bsn)
            res["basins"][str(int(bsn))] = stats(np.where(sel, model, np.nan), obs_s, w)
        results.append(res)
        panels.append((label, model, obs_s, gr))
        g, f = res["grounded"], res["floating"]
        PETSc.Sys.Print(f"  {label}: {t0:.0f}-{t:.1f}, grounded vs Smith: bias {g['bias_myr']:+.3f} m/yr, rms {g['rms_myr']:.3f}, "
                        f"corr {g['corr']:.2f}; mass rate model {g['model_gtyr']:+,.0f} vs obs {g['obs_gtyr']:+,.0f} Gt/yr "
                        f"over the {g['n']:,} km2 both cover")
        PETSc.Sys.Print(f"  {label}: floating ice: model {f.get('model_gtyr', float('nan')):+,.0f} vs obs {f.get('obs_gtyr', float('nan')):+,.0f} Gt/yr "
                        f"(mean {f.get('model_mean_myr', float('nan')):+.3f} vs {f.get('obs_mean_myr', float('nan')):+.3f} m/yr); "
                        f"grounded area {gl_area0:,.0f} -> {gl_area1:,.0f} km2 (ungrounded {gl_lost:,.0f}, grounded {gl_gained:,.0f})")
        worst = sorted(res["basins"].items(), key=lambda kv: -abs(kv[1].get("model_gtyr", 0) - kv[1].get("obs_gtyr", 0)))[:5]
        for bsn, s in worst:
            PETSc.Sys.Print(f"    basin {bsn}: model {s['model_gtyr']:+,.0f} vs obs {s['obs_gtyr']:+,.0f} Gt/yr "
                            f"(mean {s['model_mean_myr']:+.3f} vs {s['obs_mean_myr']:+.3f} m/yr)")
    if comm.rank == 0 and args.out:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=1)
    if comm.rank == 0 and args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(panels)
        fig, axes = plt.subplots(n, 3, figsize=(16, 5.2 * n), squeeze=False)
        ext = [x[0] / 1e3, x[-1] / 1e3, y[0] / 1e3, y[-1] / 1e3]
        origin = "lower" if y[1] > y[0] else "upper"
        for i, (label, model, obs_p, gr) in enumerate(panels):
            tag = f" ({args.block:g} km means)" if args.block > 0 else ""
            for j, (title, field) in enumerate((("Smith ICESat dH/dt 2003-2019" + tag, obs_p),
                                                (f"{label} modelled{tag}", model),
                                                (f"{label} minus Smith{tag}", model - obs_p))):
                ax = axes[i, j]
                im = ax.imshow(np.where(np.isfinite(model), field, np.nan), extent=ext, origin=origin,
                               cmap="RdBu", vmin=-1.0, vmax=1.0)
                ax.set_title(title)
                ax.set_aspect("equal")
                ax.set_xlim(-2800, 2800)
                ax.set_ylim(-2400, 2400)
            fig.colorbar(im, ax=axes[i, :].tolist(), label="m/yr", shrink=0.8)
        fig.savefig(args.plot, dpi=110)
        PETSc.Sys.Print(f"  wrote {args.plot}")


if __name__ == "__main__":
    main()
