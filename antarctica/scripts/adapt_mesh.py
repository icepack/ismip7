r"""Adapt the mesh of a forward checkpoint after Gudmundsson et al. (2012).

    mpiexec -n 4 python antarctica/scripts/adapt_mesh.py CHECKPOINT.h5 \
        --out-checkpoint NEW.h5 [--out-mesh NEW.msh] [--rebuild-aref]

Reads the checkpoint's mesh and state, builds the desired-element-size field
from the ISMIP7_ADAPT_* configuration (icepack2_tools.adapt_mesh), remeshes
the domain globally with gmsh on rank 0, transfers the state onto the new mesh
(MapFbetweenMeshes semantics: interpolation, ThickMin outside, bed re-sampled
from BedMachine, surface by flotation) and writes a checkpoint the forward can
restart from. The new mesh and its boundary-id sidecar land next to the old
mesh in antarctica/mesh/.

--rebuild-aref drops the frozen apparent-mass-balance reference so the forward
rebuilds it on the new mesh; only legitimate for a t=0 state (the initial
adaptation), never mid-run.
"""
import argparse
import glob
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

import firedrake as fd  # noqa: E402
import rasterio  # noqa: E402
from firedrake import Function, FunctionSpace  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402

from icepack2_tools.adapt_mesh import (AdaptMeshConfig, desired_element_size,  # noqa: E402
                                       remesh_global, transfer_state)
from icepack2_tools.geometry import sample_to_geometry  # noqa: E402
from icepack2_tools.runconfig import obs_data_root  # noqa: E402
from mesh_naming import next_adapted_mesh_name, resolve_outline_buffer  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MESH_DIR = os.path.join(HERE, "..", "mesh")
DATA_DIR = obs_data_root()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("checkpoint", nargs="?", default=None,
                    help="forward/MAP checkpoint; omit with --source-mesh --from-obs")
    ap.add_argument("--source-mesh", default=None,
                    help="with --mesh-only --from-obs: the .msh to evaluate the size field on (a scaffold only)")
    ap.add_argument("--from-obs", action="store_true",
                    help="size field from OBSERVATIONS: MEaSUReS velocity (masked where unobserved) and "
                         "BedMachine thickness/bed, instead of the checkpoint's model fields")
    ap.add_argument("--out-checkpoint", default=None)
    ap.add_argument("--out-mesh", default=None, help="default: <old basename>_adaptN.msh in antarctica/mesh")
    ap.add_argument("--rebuild-aref", action="store_true")
    ap.add_argument("--mesh-only", action="store_true",
                    help="build the size field and the new mesh (+ sidecar) from this checkpoint, no transfer: "
                         "a mesh for a fresh inversion (--out-checkpoint is then ignored)")
    ap.add_argument("--no-remesh", action="store_true",
                    help="identity test: keep the old mesh (copied under the new name) and only run the transfer")
    args = ap.parse_args()
    t0 = time.time()
    cfg = AdaptMeshConfig.from_env()
    PETSc.Sys.Print(f"adapt: config {cfg}")

    weight = None
    if args.checkpoint is None and not args.source_mesh:
        raise SystemExit("give a checkpoint, or --source-mesh with --mesh-only --from-obs")
    if args.source_mesh:
        if not (args.mesh_only and args.from_obs):
            raise SystemExit("--source-mesh only makes sense with --mesh-only --from-obs")
        mesh = fd.Mesh(args.source_mesh, name="firedrake_default")
        attrs = {"mesh_basename": os.path.basename(args.source_mesh), "geometry_space": "dg0"}
        H = b = u = None
    else:
        with fd.CheckpointFile(args.checkpoint, "r") as chk:
            mesh = chk.load_mesh()
            attrs = {k: chk.get_attr("/", k) for k in
                     ("mesh_basename", "buffer_m", "geometry_space", "raster_sample", "t_yr")
                     if chk.has_attr("/", k)}
            H = chk.load_function(mesh, name="thickness")
            b = chk.load_function(mesh, name="bed")
            try:
                u = chk.load_function(mesh, name="velocity")
            except Exception:
                u = None
    if args.from_obs:
        # Observations on the scaffold mesh: the inversion's own loaders.
        import icepack
        from firedrake import VectorFunctionSpace, conditional
        Qc = FunctionSpace(mesh, "CG", 1)
        Q0 = FunctionSpace(mesh, "DG", 0)
        Vc = VectorFunctionSpace(mesh, "CG", 1)
        bm = sorted(glob.glob(os.path.join(DATA_DIR, "bedmachine", "*.nc")))[0]
        vel = sorted(glob.glob(os.path.join(DATA_DIR, "velocity", "*.nc")))[0]
        H = sample_to_geometry(rasterio.open(f"netcdf:{bm}:thickness"), Q0, Qc, floor=0.0, method="vertex")
        b = sample_to_geometry(rasterio.open(f"netcdf:{bm}:bed"), Q0, Qc, method="vertex")
        u = icepack.interpolate((rasterio.open(f"netcdf:{vel}:VX"), rasterio.open(f"netcdf:{vel}:VY")), Vc, fillvalue=0.0)
        err = icepack.interpolate((rasterio.open(f"netcdf:{vel}:ERRX"), rasterio.open(f"netcdf:{vel}:ERRY")), Vc, fillvalue=0.0)
        weight = Function(Qc).interpolate(conditional(err[0] > 0.0, 1.0, 0.0))
        n_obs = int(mesh.comm.allreduce(float(weight.dat.data_ro.sum())))
        PETSc.Sys.Print(f"adapt: size field from OBSERVATIONS (MEaSUReS velocity on {n_obs} observed nodes, BedMachine geometry)")
    basename = str(attrs.get("mesh_basename", "")).replace(".msh", "")
    if not basename:
        raise RuntimeError("checkpoint has no mesh_basename attribute; cannot find its .msh/sidecar")
    old_msh = os.path.join(MESH_DIR, basename + ".msh")
    old_sidecar = os.path.join(MESH_DIR, f"boundary_ids_{basename}.json")
    for f in (old_msh, old_sidecar):
        if not os.path.exists(f):
            raise FileNotFoundError(f)
    # extract_ice_outline() reads ISMIP7_BUFFER_M. A scaffold (--source-mesh)
    # only supplies points, so the new mesh takes the buffer the caller names,
    # else the scaffold's tag. An adapted checkpoint's mesh keeps the old
    # mesh's buffer: its record, else its name's tag, and the environment only
    # for a legacy mesh with neither. Nothing named is an error: 0 and 20000
    # are each wrong for some legacy mesh.
    env_buffer = os.environ.get("ISMIP7_BUFFER_M")
    if args.source_mesh:
        old_buffer = resolve_outline_buffer(env_buffer, basename)
    else:
        old_buffer = resolve_outline_buffer(attrs.get("buffer_m"), basename, env_buffer)
    if old_buffer is None:
        raise SystemExit(
            f"adapt: {basename} records no buffer_m and its name carries no "
            f"_buffered<N> tag; set ISMIP7_BUFFER_M to the buffer it was built with")
    os.environ["ISMIP7_BUFFER_M"] = str(old_buffer)
    PETSc.Sys.Print(f"adapt: outline buffer {old_buffer:g} m")
    out_msh = args.out_mesh or os.path.join(
        MESH_DIR,
        next_adapted_mesh_name(basename, os.environ.get("ISMIP7_EXPERIMENT_NAME")) + ".msh")
    if os.path.realpath(out_msh) == os.path.realpath(old_msh):
        raise SystemExit(
            f"adapt: --out-mesh names the reference mesh ({out_msh}). Writing "
            f"it would destroy the mesh {basename} was built on and disable "
            f"the physical-group check against it; give a fresh name."
        )
    new_basename = os.path.splitext(os.path.basename(out_msh))[0]
    PETSc.Sys.Print(f"adapt: {basename} (t={attrs.get('t_yr', '?')}) -> {new_basename}")

    if args.no_remesh:
        import shutil
        if mesh.comm.rank == 0:
            shutil.copy(old_msh, out_msh)
            shutil.copy(old_sidecar, os.path.join(MESH_DIR, f"boundary_ids_{new_basename}.json"))
        mesh.comm.barrier()
        n_ele = mesh.comm.allreduce(FunctionSpace(mesh, "DG", 0).dof_dset.size)
        PETSc.Sys.Print("adapt: --no-remesh, transferring onto a fresh load of the same mesh")
    else:
        h_des, diag = desired_element_size(mesh, cfg, H, b, u=u, weight=weight)
        n_ele = remesh_global(mesh, h_des, cfg, out_msh, old_msh, old_sidecar)
    if args.mesh_only:
        n_old = mesh.comm.allreduce(FunctionSpace(mesh, "DG", 0).dof_dset.size)
        PETSc.Sys.Print(f"adapt: mesh only: {n_old} -> {n_ele} cells, wrote {out_msh} and its sidecar "
                        f"in {time.time() - t0:.0f} s")
        return
    if not args.out_checkpoint:
        raise SystemExit("--out-checkpoint is required unless --mesh-only")
    mesh_new = fd.Mesh(out_msh, name="firedrake_default")

    bm_fn = sorted(glob.glob(os.path.join(DATA_DIR, "bedmachine", "*.nc")))[0]
    method = str(attrs.get("raster_sample", "vertex"))

    def bed_sampler(Q_g, Qc):
        return sample_to_geometry(rasterio.open(f"netcdf:{bm_fn}:bed"), Q_g, Qc, method=method)

    def thickness_sampler(Q_g, Qc):
        return sample_to_geometry(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q_g, Qc,
                                  floor=0.0, method=method)

    audit = transfer_state(args.checkpoint, mesh_new, cfg, args.out_checkpoint,
                           new_basename, bed_sampler, rebuild_aref=args.rebuild_aref,
                           thickness_sampler=thickness_sampler)
    n_old = mesh.comm.allreduce(FunctionSpace(mesh, "DG", 0).dof_dset.size)   # owned cells only
    PETSc.Sys.Print(f"adapt: {n_old} -> {n_ele} cells; bed re-sampled ({method}); "
                    f"volume change {audit['volume_change_pct']:+.4f}%; "
                    f"wrote {args.out_checkpoint} in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
