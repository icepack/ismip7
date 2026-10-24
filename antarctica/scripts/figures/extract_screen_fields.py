#!/usr/bin/env python3
r"""Pull what calving_screen_figures.py needs out of a screen's checkpoints.

Serial (one rank), one pass per member over the state at ``--year``:
cell means of the thickness, its start value, the speed and the level set,
plus per-region volume above flotation change, into ``fields_<tag>_<member>.npz``;
and once per screen the triangulation (vertices, cells, cell areas, cell
centres), the Mouginot basin region of every cell, and the extent raster of
the members against the Greene et al. (2022) mask nearest ``--year`` on the
observations kit's grid, into ``mesh_<tag>.npz``. Regions follow the IMBIE
grouping of the 18 basins (West, East, Peninsula).

    python extract_screen_fields.py --kit <front kit .nc> --tag calvdg1h --year 2013 \
        --out screen_data results/
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

RHO_ICE, RHO_WATER = 917.0, 1028.0
STRIDE = 2   # the kit's 500 m grid, used at 1 km as score_front.py does
# Mouginot/Rignot basin ids of the kit (A-A' = 1 anticlockwise to K-A = 18) by IMBIE region
REGION_OF = {**{b: "East Antarctica" for b in (1, 2, 3, 4, 5, 6, 7, 8, 17, 18)},
             **{b: "West Antarctica" for b in (9, 10, 11, 12, 16)},
             **{b: "Antarctic Peninsula" for b in (13, 14, 15)}}
REGION_NAMES = ["Antarctic Peninsula", "West Antarctica", "East Antarctica"]


def cell_order(Q0):
    return Q0.cell_node_map().values[:, 0]


def load_member(path):
    import firedrake as fd
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        t = float(chk.get_attr("/", "t_yr"))
        law = str(chk.get_attr("/", "calving_law")) if chk.has_attr("/", "calving_law") else ""
        h = chk.load_function(mesh, name="thickness")
        h0 = chk.load_function(mesh, name="H_init")
        phi = chk.load_function(mesh, name="levelset")
        u = chk.load_function(mesh, name="velocity")
        b = chk.load_function(mesh, name="bed")
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    order = cell_order(Q0)
    cm = lambda expr: fd.Function(Q0).interpolate(expr).dat.data_ro[order].astype(np.float32)  # noqa: E731
    out = dict(t_yr=t, calving_law=law, h=cm(h), h0=cm(h0), phi=cm(phi), bed=cm(b), speed=cm(fd.sqrt(fd.inner(u, u))))
    geom = None
    if True:
        V = mesh.coordinates.function_space()
        xy = mesh.coordinates.dat.data_ro
        cells = V.cell_node_map().values
        area = fd.assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro[order]
        centres = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro[order]
        geom = dict(x=xy[:, 0].astype(np.float64), y=xy[:, 1].astype(np.float64), cells=cells.astype(np.int32),
                    area=area.astype(np.float64), cx=centres[:, 0].astype(np.float64), cy=centres[:, 1].astype(np.float64))
    return out, geom


def vaf_gt(h, bed, area):
    r"""Volume above flotation of each cell [Gt]."""
    hf = np.maximum(-bed, 0.0) * RHO_WATER / RHO_ICE
    return np.maximum(h - hf, 0.0) * area * RHO_ICE / 1e12


class Kit:
    def __init__(self, path, year):
        import netCDF4
        from score_front import ocean_connected, year_of_days
        self.ds = netCDF4.Dataset(path)
        self.x = np.asarray(self.ds["x"][::STRIDE], float)
        self.y = np.asarray(self.ds["y"][::STRIDE], float)
        self.dx = float(self.x[1] - self.x[0])
        self.basins = np.asarray(self.ds["mouginot_basins"][::STRIDE, ::STRIDE]).astype(np.int8)
        years = year_of_days(self.ds["greene_mask_time"][:])
        k = int(np.argmin(np.abs(years - year)))
        self.epoch_year = float(years[k])
        ice = np.asarray(self.ds["icemask_greene"][k, ::STRIDE, ::STRIDE]) > 0
        self.obs_ice = ~ocean_connected(ice)
        self.ocean_connected = ocean_connected

    def pixel(self, px, py):
        i = np.clip(np.rint((py - self.y[0]) / (self.y[1] - self.y[0])).astype(int), 0, len(self.y) - 1)
        j = np.clip(np.rint((px - self.x[0]) / self.dx).astype(int), 0, len(self.x) - 1)
        return i, j

    def raster_index(self, cx, cy, area):
        from scipy.spatial import cKDTree
        X, Y = np.meshgrid(self.x, self.y)
        d, idx = cKDTree(np.column_stack([cx, cy])).query(np.column_stack([X.ravel(), Y.ravel()]), workers=-1)
        inside = d <= 1.5 * np.sqrt(area[idx])
        return np.where(inside, idx, -1).reshape(X.shape)

    def extent(self, idx, phi):
        ice = (idx >= 0) & (phi[idx] < 0.0)
        return ~self.ocean_connected(ice)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results")
    ap.add_argument("--kit", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--year", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--members", default=None, help="comma-separated subset")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    pat = re.compile(rf"hist_\w+?_{re.escape(args.tag)}_(\w+?)_\d+_t{args.year:.1f}\.h5$")
    paths = {}
    for p in sorted(glob.glob(os.path.join(args.results, f"hist_*_{args.tag}_*_t{args.year:.1f}.h5"))):
        m = pat.search(os.path.basename(p))
        if m and (not args.members or m.group(1) in args.members.split(",")):
            paths[m.group(1)] = p
    if not paths:
        sys.exit(f"no t{args.year:.1f} states for {args.tag} under {args.results}")
    print("members:", ", ".join(paths), flush=True)
    kit = Kit(args.kit, args.year)
    geom = idx = region = None
    count = None
    for name, p in paths.items():
        f, g = load_member(p)
        if geom is None:
            geom = g
            i, j = kit.pixel(g["cx"], g["cy"])
            basin = kit.basins[i, j]
            region = np.array([REGION_NAMES.index(REGION_OF[int(b)]) if int(b) in REGION_OF else -1 for b in basin], np.int8)
            idx = kit.raster_index(g["cx"], g["cy"], g["area"])
            count = np.zeros(idx.shape, np.int16)
            print("triangulation and raster built", flush=True)
        ice = kit.extent(idx, f["phi"])
        count += ice
        dv = vaf_gt(f["h"], f["bed"], geom["area"]) - vaf_gt(f["h0"], f["bed"], geom["area"])
        region_dvaf = np.array([dv[region == r].sum() for r in range(len(REGION_NAMES))])
        np.savez_compressed(os.path.join(args.out, f"fields_{args.tag}_{name}.npz"), t_yr=f["t_yr"], calving_law=f["calving_law"], h=f["h"], h0=f["h0"],
                            phi=f["phi"], speed=f["speed"], region_names=np.array(REGION_NAMES), region_dvaf_gt=region_dvaf,
                            ice_area_km2=float(ice.sum()) * kit.dx ** 2 / 1e6)
        print(f"{name}: t={f['t_yr']:.2f} dVAF by region (Gt) " +
              ", ".join(f"{n} {v:+.0f}" for n, v in zip(REGION_NAMES, region_dvaf)), flush=True)
    s = slice(None, None, 2)   # the extent raster at 2 km is plenty for a continent map
    np.savez_compressed(os.path.join(args.out, f"mesh_{args.tag}.npz"), x=geom["x"], y=geom["y"], cells=geom["cells"],
                        area=geom["area"], cx=geom["cx"], cy=geom["cy"], region=region,
                        count_x=kit.x[s], count_y=kit.y[s], count=count[s, s], obs_ice=kit.obs_ice[s, s],
                        obs_epoch_year=kit.epoch_year, t_yr=args.year, n_members=len(paths))
    print("wrote", os.path.join(args.out, f"mesh_{args.tag}.npz"), flush=True)


if __name__ == "__main__":
    main()
