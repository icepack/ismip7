#!/usr/bin/env python3
r"""Score a forward run's calving front against the Greene et al. (2022) record.

Read-only. For every state checkpoint named, the run's ice extent (the sign
of its level set) is compared with the Greene mask of the epoch nearest the
checkpoint's year (the masks are in the ISMIP7 observations kit, 24 epochs
1997-2021 on a 500 m grid). Both are put on the kit's grid at 1 km, and the
front of each is the edge of its ice against the ocean it is connected to:
the level set is anchored on the ice extent, so it also bounds every rock
outcrop and every ice-free pocket inland (the Transantarctic, Ellsworth and
Prince Charles mountains carry more of that boundary than the coast does),
while the record's masks draw only the coast. Filling what the open ocean
cannot reach makes the two comparable.

The residual on the model's front is the signed distance to the observed
coast, r > 0 where the model front lies seaward of the observed one (the
model keeps ice the record says is gone), r < 0 where it lies inland. Its
rms is the distance between the two fronts, its bias their net offset, and
the seaward and inland shares say whether a small rms hides a front that is
right on average and wrong everywhere. The same is reported the other way,
from the observed front to the model's (``obs_front``, on the mesh only:
the islands beyond its buffer are no front the run can have, and on a
mesh without a buffer the coast is scored only where it lies within the
run's extent): a shelf the model has lost shows in both, ice the model
keeps past the record shows mainly in the first. Both are given per Mouginot basin. A member of a calving
campaign is judged on the mean front rms over the epochs its run covers,
against the pinned-front and no-calving baselines. The ice-free area
inland that the ocean does not reach is reported too: it starts as the
outcrops and lakes, and a law whose rate is not gated on grounded ice can
grow it, since the level set bounds those edges as well. Each state also
carries the ice area of the run and of the record on the kit grid (ocean
filled, so inland rock counts as ice in both), and the summary gives their
change over the scored epochs: the record's area is nearly steady over
2003-2015 (+0.5 thousand km2/yr on 13.6 million), so a law tuned to the
continent as a whole should hold the area as the record does.

    mpiexec -n 8 python antarctica/scripts/score_front.py --out score.json \
        results/hist_cesm2_waccm_<tag>_2000_t20*.h5
"""
import argparse
import json
import os
import sys

import numpy as np

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)
from icepack2_tools.runconfig import obs_data_root  # noqa: E402

DEFAULT_KIT = os.path.join(obs_data_root(), "..", "..", "ISMIP7", "AIS", "obs", "mipkit",
                           "AntarcticaObsISMIP7-v1.2.nc")
STRIDE = 2  # the kit's 500 m grid, scored at 1 km


def year_of_days(days):
    import datetime
    t0 = datetime.date(1900, 1, 1)
    out = []
    for d in np.asarray(days):
        day = t0 + datetime.timedelta(days=int(d))
        start = datetime.date(day.year, 1, 1)
        length = (datetime.date(day.year + 1, 1, 1) - start).days
        out.append(day.year + (day - start).days / length)
    return np.array(out)


def stats(r, w):
    total = float(w.sum())
    if not total > 0:
        return {"n": 0}
    seaward, inland = r > 0, r < 0
    return {
        "n": int(len(r)),
        "rms_m": float(np.sqrt((r ** 2 * w).sum() / total)),
        "bias_m": float((r * w).sum() / total),
        "mae_m": float((np.abs(r) * w).sum() / total),
        "seaward_share": float(w[seaward].sum() / total),
        "inland_share": float(w[inland].sum() / total),
        "within_1km": float(w[np.abs(r) <= 1e3].sum() / total),
        "within_2km": float(w[np.abs(r) <= 2e3].sum() / total),
    }


def ocean_connected(ice):
    r"""The ice-free pixels the open ocean reaches, 4-connected from the
    grid's corner (far offshore on the kit grid)."""
    from scipy.ndimage import label
    lab, _ = label(~ice)
    return lab == lab[0, 0]


def front_weights(ice, ocean, dx):
    r"""Front length [m] carried by each ice pixel: one ``dx`` per edge it
    shares with an ocean pixel."""
    w = np.zeros(ice.shape, np.float32)
    w[1:, :] += ice[1:, :] & ocean[:-1, :]
    w[:-1, :] += ice[:-1, :] & ocean[1:, :]
    w[:, 1:] += ice[:, 1:] & ocean[:, :-1]
    w[:, :-1] += ice[:, :-1] & ocean[:, 1:]
    return w * dx


def signed_distance(ice, dx):
    r"""Exact Euclidean signed distance [m] to the edge of ``ice``, negative inside."""
    from scipy.ndimage import distance_transform_edt
    return np.where(ice, -distance_transform_edt(ice) * dx,
                    distance_transform_edt(~ice) * dx).astype(np.float32)


class Record:
    r"""The Greene masks and Mouginot basins of the kit on the 1 km grid,
    held on rank 0; each epoch's ocean-connected ice, its front and the
    signed distance to it on demand."""

    def __init__(self, path, comm):
        self.comm = comm
        self.kit = None
        if comm.rank == 0:
            import netCDF4
            self.kit = netCDF4.Dataset(path)
            self.x = np.asarray(self.kit["x"][::STRIDE], float)
            self.y = np.asarray(self.kit["y"][::STRIDE], float)
            self.dx = float(self.x[1] - self.x[0])
            self.years = year_of_days(self.kit["greene_mask_time"][:])
            self.basins = np.asarray(self.kit["mouginot_basins"][::STRIDE, ::STRIDE])
        self.years = comm.bcast(self.years if comm.rank == 0 else None, root=0)
        self._cache = {}

    def epoch_near(self, year, tolerance):
        k = int(np.argmin(np.abs(self.years - year)))
        return k if abs(self.years[k] - year) <= tolerance else None

    def epoch(self, k):
        r"""Rank 0 only: (ice, ocean, front weights, signed distance) of epoch ``k``."""
        if k not in self._cache:
            ice = np.asarray(self.kit["icemask_greene"][k, ::STRIDE, ::STRIDE]) > 0
            ocean = ocean_connected(ice)
            ice = ~ocean
            self._cache[k] = (ice, ocean, front_weights(ice, ocean, self.dx), signed_distance(ice, self.dx))
        return self._cache[k]


class ModelRaster:
    r"""Rank 0 only: the map from the record's grid to the run's cells
    (nearest cell centre; pixels farther than a cell from any centre are off
    the mesh and ocean), built once per mesh."""

    def __init__(self, record):
        self.record = record
        self.key = None

    def extent(self, phi, xy, area):
        key = (len(phi), float(xy[0, 0]), float(xy[-1, 1]))
        if key != self.key:
            from scipy.spatial import cKDTree
            X, Y = np.meshgrid(self.record.x, self.record.y)
            d, idx = cKDTree(xy).query(np.column_stack([X.ravel(), Y.ravel()]), workers=-1)
            inside = d <= 1.5 * np.sqrt(area[idx])
            self.idx = np.where(inside, idx, -1).reshape(X.shape)
            self.key = key
        ice = (self.idx >= 0) & (phi[self.idx] < 0.0)
        ocean = ocean_connected(ice)
        # ice-free pixels on the mesh the ocean does not reach: rock
        # outcrops and lakes at the start, and any ice a law removed inland
        # (the level set bounds those too, so a rate law can act there)
        inland = (self.idx >= 0) & ~ice & ~ocean
        return ~ocean, ocean, float(inland.sum()) * self.record.dx ** 2


def score_state(path, record, raster, tolerance, comm):
    import firedrake as fd
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        t = float(chk.get_attr("/", "t_yr"))
        phi = chk.load_function(mesh, name="levelset")
    k = record.epoch_near(t, tolerance)
    if k is None:
        return {"file": os.path.basename(path), "t_yr": t, "epoch": None} if comm.rank == 0 else None
    Q0 = phi.function_space()
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro
    area = fd.assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro
    # every rank's cells to rank 0, where the record is
    parts = comm.gather((phi.dat.data_ro.copy(), np.asarray(xy).copy(), area.copy()), root=0)
    if comm.rank != 0:
        return None
    pm, pts, ar = (np.concatenate([p[i] for p in parts]) for i in range(3))
    ice_m, ocean_m, inland_free = raster.extent(pm, pts, ar)
    ice_o, ocean_o, w_o, dist_o = record.epoch(k)
    w_m = front_weights(ice_m, ocean_m, record.dx)
    dist_m = signed_distance(ice_m, record.dx)
    # the observed coast only where the run has cells: the islands beyond
    # the mesh's buffer are tens of km from any front the run can have
    front_m, front_o = w_m > 0, (w_o > 0) & (raster.idx >= 0)
    # model front to the observed coast, and the observed coast to the
    # model's: r > 0 where the model front lies seaward of the observed one
    r_m, r_o = dist_o[front_m], -dist_m[front_o]
    out = {"file": os.path.basename(path), "t_yr": t, "epoch": int(k), "epoch_year": float(record.years[k]),
           "front_km": float(w_m[front_m].sum() / 1e3), "obs_front_km": float(w_o[front_o].sum() / 1e3),
           "inland_ice_free_km2": inland_free / 1e6,
           "ice_area_km2": float(ice_m.sum()) * record.dx ** 2 / 1e6,
           "obs_ice_area_km2": float(ice_o.sum()) * record.dx ** 2 / 1e6,
           "front": stats(r_m, w_m[front_m]), "obs_front": stats(r_o, w_o[front_o]), "basins": {}}
    b_m, b_o = record.basins[front_m], record.basins[front_o]
    for b in np.unique(b_m):
        m, o = b_m == b, b_o == b
        if m.sum() >= 5:
            out["basins"][str(int(b))] = {"front_km": float(w_m[front_m][m].sum() / 1e3),
                                          "front": stats(r_m[m], w_m[front_m][m]),
                                          "obs_front": stats(r_o[o], w_o[front_o][o])}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("states", nargs="+", help="state checkpoints of one run (any years)")
    ap.add_argument("--kit", default=DEFAULT_KIT)
    ap.add_argument("--tolerance", type=float, default=0.6, help="largest |t_yr - epoch year| scored [yr]")
    ap.add_argument("--out", default=None, help="JSON output")
    args = ap.parse_args()
    missing = [p for p in args.states if not os.path.exists(p)]
    if missing:
        sys.exit("no such state: " + ", ".join(missing))
    from firedrake.petsc import PETSc
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    record = Record(args.kit, comm)
    raster = ModelRaster(record)
    results = []
    for path in args.states:
        res = score_state(path, record, raster, args.tolerance, comm)
        if comm.rank == 0:
            results.append(res)
            if res.get("epoch") is None:
                PETSc.Sys.Print(f"  {res['file']}: t={res['t_yr']:.2f}, no Greene epoch within {args.tolerance} yr")
            else:
                f, g = res["front"], res["obs_front"]
                PETSc.Sys.Print(f"  {res['file']}: t={res['t_yr']:.2f} vs epoch {res['epoch_year']:.2f}: front rms "
                                f"{f['rms_m'] / 1e3:.2f} km, bias {f['bias_m'] / 1e3:+.2f} km, seaward {f['seaward_share']:.2f}, "
                                f"within 2 km {f['within_2km']:.2f}, front {res['front_km']:.0f} km "
                                f"(obs {res['obs_front_km']:.0f} km, obs->model rms {g['rms_m'] / 1e3:.2f} km), "
                                f"inland ice-free {res['inland_ice_free_km2']:.0f} km2, ice area "
                                f"{res['ice_area_km2']:.0f} km2 (obs {res['obs_ice_area_km2']:.0f})")
    if comm.rank == 0:
        scored = [r for r in results if r.get("epoch") is not None]
        rms = [r["front"]["rms_m"] for r in scored if r["front"].get("n")]
        summary = {"mean_rms_m": float(np.mean(rms)) if rms else None,
                   "final_rms_m": float(rms[-1]) if rms else None,
                   "mean_bias_m": float(np.mean([r["front"]["bias_m"] for r in scored if r["front"].get("n")])) if rms else None,
                   "mean_obs_rms_m": float(np.mean([r["obs_front"]["rms_m"] for r in scored if r["obs_front"].get("n")])) if rms else None,
                   "epochs": [r["epoch_year"] for r in scored],
                   "area_change_km2": float(scored[-1]["ice_area_km2"] - scored[0]["ice_area_km2"]) if scored else None,
                   "obs_area_change_km2": float(scored[-1]["obs_ice_area_km2"] - scored[0]["obs_ice_area_km2"]) if scored else None}
        PETSc.Sys.Print(f"summary: mean front rms {summary['mean_rms_m'] and summary['mean_rms_m'] / 1e3:.2f} km over "
                        f"{len(scored)} epochs, final {summary['final_rms_m'] and summary['final_rms_m'] / 1e3:.2f} km; "
                        f"ice area change {summary['area_change_km2'] or 0:+.0f} km2 "
                        f"(record {summary['obs_area_change_km2'] or 0:+.0f} km2)")
        if args.out:
            with open(args.out, "w") as f:
                json.dump({"summary": summary, "states": results}, f, indent=1)


if __name__ == "__main__":
    main()
