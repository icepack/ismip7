#!/usr/bin/env python3
r"""Where the exact cliff push acts, and what its version changes (issue #166).

    python antarctica/scripts/check_cliff_facets.py STATE.h5 [--out PREFIX]
        [--h-ice H] [--compare OTHER.h5] [--imbie2 PATH]

STATE is a MAP or a forward checkpoint (``thickness``, ``bed``, ``surface``,
and ``velocity`` and ``basal_stress`` when it holds a solved state). For every
facet between a cell holding ice (``H >= h_ice``, by default
``ISMIP7_FRONT_HMIN``) and one holding none, it classes the facet by the
ice-free neighbour's bed ``B`` against the ice side's base and surface:

==================  =========================================
class               neighbour
==================  =========================================
``ocean_deep``      seabed at or below the ice base
``ocean_shoal``     seabed between the ice base and sea level
``land_low``        land at or below the ice base
``land_partial``    land between the ice base and surface
``land_above``      land at or above the ice surface
==================  =========================================

Per class it reports the facets and their length, the edge thickness, the
push of ``exact_front`` version 1 (the free-cliff push) and version 2 (the
face above ``B``), both from ``dual_friction.cliff_push``, and the change
between them, in total and per IMBIE basin. With a solved state it compares
each grounded edge cell's change in push, summed over its front facets, with
the cell's own basal drag ``|tau_b| A``.

``--compare OTHER.h5`` names the same MAP solved under the other version, on
the same mesh: ``score_map.py --save-state`` on a copy of the MAP whose
``exact_front`` attribute is rewritten, since a forward refuses an
``ISMIP7_EXACT_FRONT`` that differs from its MAP. It adds the speed change by
distance to the nearest facet whose push the version changes, and the
grounded discharge (``region_budget.discharge``) of each, in total and per
basin. STATE is the reference the change is taken from.

Runs under MPI. Every statistic is gathered to rank 0, so any rank count
prints the same numbers. Writes PREFIX.txt, PREFIX_facets.npz and, with
matplotlib, PREFIX.png.
"""

import argparse
import os
import sys

import numpy as np
import firedrake as fd
from firedrake import (Constant, Function, FunctionSpace,
                       SpatialCoordinate, TestFunction, VectorFunctionSpace,
                       assemble, conditional, dS, dx, ge, gt, inner, sqrt)

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_SCRIPTS)))
sys.path.insert(0, _SCRIPTS)

from icepack2_tools.dual_friction import cliff_push, g, rho_I  # noqa: E402
from icepack2_tools.runconfig import front_hmin  # noqa: E402

CLASSES = ("ocean_deep", "ocean_shoal", "land_low", "land_partial", "land_above")
#: the classes whose push version 2 changes
CHANGED = ("ocean_shoal", "land_partial", "land_above")
TOL = 1e-6                     # m: a bed this close to the base counts as at it
DIST_BINS = (0.0, 2e3, 6e3, 20e3, np.inf)
YEAR_S = 365.25 * 86400.0


def _load(path):
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        f = {k: chk.load_function(mesh, name=k) for k in ("thickness", "bed", "surface")}
        for k in ("velocity", "basal_stress"):
            try:
                f[k] = chk.load_function(mesh, name=k)
            except (KeyError, RuntimeError, ValueError):
                pass
    return mesh, f


def _gather(comm, arr):
    parts = comm.gather(np.asarray(arr), root=0)
    return np.concatenate(parts) if comm.rank == 0 else None


def _wq(v, w, ps):
    r"""Weighted quantiles of ``v``."""
    if v.size == 0:
        return [np.nan for _ in ps]
    o = np.argsort(v)
    c = np.cumsum(w[o]) / w.sum()
    return [v[o][min(np.searchsorted(c, p), v.size - 1)] for p in ps]


def classify(B, base, s):
    r"""The class of each front facet from the neighbour's bed ``B`` and the
    ice side's base and surface (numpy arrays)."""
    ocean = B < 0.0
    return {
        "ocean_deep": ocean & (B <= base + TOL),
        "ocean_shoal": ocean & (B > base + TOL),
        "land_low": ~ocean & (B <= base + TOL),
        "land_partial": ~ocean & (B > base + TOL) & (B < s - TOL),
        "land_above": ~ocean & (B >= s - TOL),
    }


def facet_census(mesh, f, h_ice):
    r"""Per-facet arrays on rank 0 (None elsewhere) for every front facet."""
    H, b, s = f["thickness"], f["bed"], f["surface"]
    T = FunctionSpace(mesh, "HDiv Trace", 0)
    w = TestFunction(T)("+")
    ice = conditional(ge(H, h_ice), 1.0, 0.0)
    x = SpatialCoordinate(mesh)
    sides = (("+", "-"), ("-", "+"))

    def oriented(fn):
        return sum(ice(p) * (1.0 - ice(q)) * fn(p, q) for p, q in sides)

    def per_facet(expr):
        return assemble(expr * w * dS).dat.data_ro.copy()

    exprs = {
        "length": 1.0,
        "front": oriented(lambda p, q: 1.0),
        "H": oriented(lambda p, q: H(p)),
        "s": oriented(lambda p, q: s(p)),
        "B": oriented(lambda p, q: b(q)),
        "push1": oriented(lambda p, q: cliff_push(H(p), s(p), version=1)),
        "push2": oriented(lambda p, q: cliff_push(H(p), s(p), version=2, b_other=b(q))),
        "x": x[0],
        "y": x[1],
    }
    raw = {k: per_facet(e) for k, e in exprs.items()}
    L = raw["length"]
    sel = (L > 0.0) & (raw["front"] > 0.5 * L)
    out = {"length": L[sel]}
    for k in ("H", "s", "B", "push1", "push2", "x", "y"):
        out[k] = raw[k][sel] / L[sel]
    return {k: _gather(mesh.comm, v) for k, v in out.items()}


def edge_cell_ratio(mesh, f, h_ice):
    r"""For each grounded edge cell: the change in push summed over its front
    facets, ``sum |push2 - push1| L``, and its basal drag ``|tau_b| A`` (both
    MPa m^2), gathered to rank 0."""
    H, b, s = f["thickness"], f["bed"], f["surface"]
    Q0 = FunctionSpace(mesh, "DG", 0)
    phi = TestFunction(Q0)
    ice = conditional(ge(H, h_ice), 1.0, 0.0)
    # only the CHANGED classes, a bed above the base: elsewhere the two
    # versions agree and their difference is roundoff
    change = sum(ice(p) * (1.0 - ice(q)) * phi(p)
                 * conditional(gt(b(q), s(p) - H(p) + TOL), 1.0, 0.0)
                 * abs(cliff_push(H(p), s(p), version=2, b_other=b(q))
                       - cliff_push(H(p), s(p), version=1))
                 for p, q in (("+", "-"), ("-", "+")))
    dpush = assemble(change * dS).dat.data_ro.copy()
    tau = f["basal_stress"]
    drag = assemble(sqrt(inner(tau, tau)) * phi * dx).dat.data_ro.copy()
    area = assemble(phi * dx).dat.data_ro.copy()
    grounded = Function(Q0).interpolate(
        conditional(ge(s - H, b + Constant(TOL)), 0.0, 1.0)).dat.data_ro.copy()
    centroid = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro.copy()
    thick = Function(Q0).interpolate(H).dat.data_ro.copy()
    keep = (dpush > 0.0) & (grounded > 0.5)
    out = {"dpush": dpush[keep], "drag": drag[keep], "area": area[keep],
           "H": thick[keep], "x": centroid[keep, 0], "y": centroid[keep, 1]}
    return {k: _gather(mesh.comm, v) for k, v in out.items()}


def cell_state(mesh, f, h_ice, basins_of):
    r"""Cell centroids, areas, velocity, grounded/floating flags and the
    grounded discharge per basin of one solved state."""
    from region_budget import discharge, regions

    Q0 = FunctionSpace(mesh, "DG", 0)
    V0 = VectorFunctionSpace(mesh, "DG", 0)
    centroid = Function(V0).interpolate(SpatialCoordinate(mesh))
    u = Function(V0).interpolate(f["velocity"])
    area = assemble(TestFunction(Q0) * dx).dat.data_ro.copy()
    state = {"mesh": mesh, "thickness": f["thickness"], "bed": f["bed"],
             "velocity": f["velocity"]}
    reg = regions(state)
    reg["open"] = Function(Q0)
    reg["open"].dat.data[:] = 1.0 - reg["ice"].dat.data_ro
    q_total = discharge(state, reg)
    xy = centroid.dat.data_ro
    q_basin = {}
    grounded_all = reg["grounded"].dat.data_ro.copy()
    if basins_of is not None:
        cell_basin = basins_of(xy[:, 0], xy[:, 1])
        # the same list on every rank: discharge assembles collectively
        ids = np.unique(np.concatenate(mesh.comm.allgather(np.unique(cell_basin))))
        for k in ids:
            reg["grounded"].dat.data[:] = grounded_all * (cell_basin == k)
            q_basin[int(k)] = discharge(state, reg)
        reg["grounded"].dat.data[:] = grounded_all
    out = {"x": xy[:, 0], "y": xy[:, 1], "area": area,
           "ux": u.dat.data_ro[:, 0], "uy": u.dat.data_ro[:, 1],
           "grounded": grounded_all, "floating": reg["floating"].dat.data_ro.copy()}
    return {k: _gather(mesh.comm, v) for k, v in out.items()}, q_total, q_basin


def _basin_lookup(path):
    from icepack2_tools.forcing import _basin_on_mesh, imbie2_basin_path

    path = path or imbie2_basin_path()
    if not os.path.exists(path):
        return None, f"no IMBIE2 basin grid at {path}: per-basin rows skipped"
    return (lambda xs, ys: _basin_on_mesh(xs, ys, imbie2=path)), f"basins from {path}"


def census_lines(fc, basins_of):
    L = fc["length"]
    classes = classify(fc["B"], fc["s"] - fc["H"], fc["s"])
    dpush = fc["push2"] - fc["push1"]
    lines = [f"front facets: {L.size} ({L.sum() / 1e3:,.0f} km); push in MN/m "
             f"(1 MN/m over a 2 km cell is a 0.5 kPa basal stress)"]
    head = (f"{'class':<13s} {'facets':>7s} {'km':>7s} {'H p50/p90/p99 (m)':>20s} "
            f"{'v1 push p50/p90/p99':>22s} {'v2 push p50/p90/p99':>22s} "
            f"{'sum v1 (GN)':>12s} {'sum v2-v1 (GN)':>15s}")
    lines.append(head)
    for name in CLASSES:
        m = classes[name]
        w = L[m]
        Hq = _wq(fc["H"][m], w, (0.5, 0.9, 0.99))
        p1 = _wq(fc["push1"][m], w, (0.5, 0.9, 0.99))
        p2 = _wq(fc["push2"][m], w, (0.5, 0.9, 0.99))
        lines.append(
            f"{name:<13s} {m.sum():7d} {w.sum() / 1e3:7,.0f} "
            f"{'/'.join(f'{q:.0f}' for q in Hq):>20s} "
            f"{'/'.join(f'{q:.1f}' for q in p1):>22s} "
            f"{'/'.join(f'{q:.1f}' for q in p2):>22s} "
            f"{(fc['push1'][m] * w).sum() / 1e3:12,.1f} {(dpush[m] * w).sum() / 1e3:15,.1f}")
    changed = np.zeros(L.size, bool)
    for name in CHANGED:
        changed |= classes[name]
    for cut in (100.0, 200.0, 500.0):
        m = changed & (fc["H"] > cut)
        lines.append(f"changed facets with edge H > {cut:.0f} m: {m.sum()} "
                     f"({L[m].sum() / 1e3:,.0f} km), v2 - v1 {(dpush[m] * L[m]).sum() / 1e3:,.1f} GN")
    # ties broken by position, so every rank count lists the same facets
    order = [i for i in np.lexsort((fc["y"], fc["x"], dpush * L)) if changed[i]][:12]
    lines.append("largest changes (x km, y km, class, edge H m, B - base m, v1 -> v2 MN/m, length km):")
    for i in order:
        name = next(c for c in CLASSES if classes[c][i])
        lines.append(f"  {fc['x'][i] / 1e3:7.0f} {fc['y'][i] / 1e3:7.0f}  {name:<13s} "
                     f"{fc['H'][i]:6.0f} {fc['B'][i] - (fc['s'][i] - fc['H'][i]):7.0f} "
                     f"{fc['push1'][i]:8.1f} -> {fc['push2'][i]:8.1f}  {L[i] / 1e3:5.2f}")
    if basins_of is not None:
        basin = basins_of(fc["x"], fc["y"])
        lines.append("v2 - v1 by IMBIE basin, changed facets (GN; km of facet):")
        for k in np.unique(basin[changed]):
            m = changed & (basin == k)
            lines.append(f"  basin {k:3d}: {(dpush[m] * L[m]).sum() / 1e3:10,.1f}  ({L[m].sum() / 1e3:,.0f} km)")
    return lines, changed


def ratio_lines(ec):
    r = ec["dpush"] / np.maximum(ec["drag"], 1e-30)
    lines = [f"grounded edge cells whose push version 2 changes: {r.size}; "
             f"change in push / basal drag |tau_b| A, area-weighted p50/p90/p99 "
             f"{'/'.join(f'{q:.2f}' for q in _wq(r, ec['area'], (0.5, 0.9, 0.99)))}"]
    for cut in (0.25, 1.0, 4.0):
        m = r > cut
        lines.append(f"  ratio > {cut:g}: {m.sum()} cells, {ec['area'][m].sum() / 1e6:,.0f} km^2, "
                     f"edge H p50 {np.median(ec['H'][m]) if m.any() else np.nan:.0f} m")
    return lines


def compare_lines(ref, oth, fc, changed, q_ref, q_oth, qb_ref, qb_oth):
    from scipy.spatial import cKDTree

    tree = cKDTree(np.column_stack([oth["x"], oth["y"]]))
    dist, idx = tree.query(np.column_stack([ref["x"], ref["y"]]))
    if dist.max() > 1.0:
        raise SystemExit(f"--compare: the states are on different meshes (a centroid "
                         f"is {dist.max():.1f} m from its match)")
    du = np.hypot(oth["ux"][idx] - ref["ux"], oth["uy"][idx] - ref["uy"])
    speed = np.hypot(ref["ux"], ref["uy"])
    lines = [f"grounded discharge: reference {q_ref:,.1f} Gt/yr, other {q_oth:,.1f}, "
             f"change {q_oth - q_ref:+,.2f} ({100 * (q_oth - q_ref) / max(q_ref, 1e-9):+.3f} %)"]
    if qb_ref:
        lines.append("  by IMBIE basin (Gt/yr): reference, change, %")
        for k in sorted(qb_ref):
            d = qb_oth.get(k, 0.0) - qb_ref[k]
            pct = 100 * d / qb_ref[k] if abs(qb_ref[k]) > 1e-9 else np.nan
            lines.append(f"    basin {k:3d}: {qb_ref[k]:9.1f} {d:+9.2f} {pct:+8.3f}")
    walls = np.column_stack([fc["x"][changed], fc["y"][changed]])
    to_wall = (cKDTree(walls).query(np.column_stack([ref["x"], ref["y"]]))[0]
               if walls.size else np.full(ref["x"].size, np.inf))
    lines.append("speed change |du| by distance to the nearest changed facet "
                 "(area-weighted p50/p99 m/yr; km^2 with |du| > 10 m/yr; with |du| > 10 % of the speed):")
    for kind in ("grounded", "floating"):
        on = ref[kind] > 0.5
        for lo, hi in zip(DIST_BINS[:-1], DIST_BINS[1:]):
            m = on & (to_wall >= lo) & (to_wall < hi)
            a = ref["area"][m]
            big = m & (du > 10.0)
            rel = m & (du > 0.1 * np.maximum(speed, 1e-9))
            lines.append(
                f"  {kind:<8s} {lo / 1e3:4.0f}-{hi / 1e3:<4.0f} km: {a.sum() / 1e6:11,.0f} km^2, "
                f"|du| {'/'.join(f'{q:.2f}' for q in _wq(du[m], a, (0.5, 0.99)))}, "
                f"{ref['area'][big].sum() / 1e6:9,.0f}, {ref['area'][rel].sum() / 1e6:9,.0f}")
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("state")
    ap.add_argument("--out", default=None, help="output prefix (default: STATE without .h5)")
    ap.add_argument("--h-ice", type=float, default=None,
                    help="thickness that holds ice [m] (default ISMIP7_FRONT_HMIN)")
    ap.add_argument("--compare", default=None,
                    help="the same MAP solved under the other version, same mesh")
    ap.add_argument("--imbie2", default=None, help="IMBIE2 basin grid (default: the data root's)")
    args = ap.parse_args()
    h_ice = args.h_ice if args.h_ice is not None else front_hmin()
    out = args.out or os.path.splitext(args.state)[0] + "_cliff"
    comm = fd.COMM_WORLD
    basins_of, basin_note = _basin_lookup(args.imbie2)

    mesh, f = _load(args.state)
    fc = facet_census(mesh, f, h_ice)
    ec = edge_cell_ratio(mesh, f, h_ice) if "basal_stress" in f else None
    cmp = None
    if args.compare:
        if "velocity" not in f:
            raise SystemExit(f"--compare needs a solved state; {args.state} has no velocity")
        ref, q_ref, qb_ref = cell_state(mesh, f, h_ice, basins_of)
        mesh2, f2 = _load(args.compare)
        oth, q_oth, qb_oth = cell_state(mesh2, f2, h_ice, basins_of)
        cmp = (ref, oth, q_ref, q_oth, qb_ref, qb_oth)

    if comm.rank != 0:
        return
    lines = [os.path.basename(args.state), f"h_ice {h_ice:g} m, rho_I g = {rho_I * g:.4g} MPa/m, "
             f"{comm.size} ranks; {basin_note}"]
    census, changed = census_lines(fc, basins_of)
    lines += census
    if ec is not None:
        lines += ratio_lines(ec)
    if cmp is not None:
        ref, oth, q_ref, q_oth, qb_ref, qb_oth = cmp
        lines.append(f"compared with {os.path.basename(args.compare)}")
        lines += compare_lines(ref, oth, fc, changed, q_ref, q_oth, qb_ref, qb_oth)
    text = "\n".join(lines)
    print(text, flush=True)
    with open(out + ".txt", "w") as fh:
        fh.write(text + "\n")
    np.savez(out + "_facets.npz", changed=changed, **fc)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        L = fc["length"]
        d = np.abs(fc["push2"] - fc["push1"]) * L / 1e3          # GN
        m = changed & (d > 0)
        fig, ax = plt.subplots(figsize=(9, 8))
        ax.scatter(fc["x"][~changed] / 1e3, fc["y"][~changed] / 1e3, s=0.2, c="0.8")
        sc = ax.scatter(fc["x"][m] / 1e3, fc["y"][m] / 1e3, c=d[m], s=3, cmap="magma_r",
                        norm=LogNorm(vmin=max(d[m].min(), 1e-3), vmax=d[m].max()))
        ax.set_aspect("equal")
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        ax.set_title("Facets whose cliff push exact_front version 2 changes")
        fig.colorbar(sc, ax=ax, label="|v2 - v1| push times facet length (GN)")
        fig.savefig(out + ".png", dpi=130, bbox_inches="tight")
    except Exception as exc:  # the numbers are the result; a figure is a bonus
        print("no figure:", exc)


if __name__ == "__main__":
    main()
