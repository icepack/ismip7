r"""ISMIP7's mesh adaptation: the shared scheme plus this model's wiring.

The scheme (Gudmundsson et al. (2012): criteria, ``Error2EleSize``, relaxation and
ratio limits, ``GLrange``/``CFrange`` bands, gmsh global remeshing with the reference's
element-count control, ``MapFbetweenMeshes`` transfer, and the three DG0
additions) lives in :mod:`icepack_tools.adapt_mesh`, shared with the other
icepack2 projects. Everything ISMIP7 adds is here:

* the Antarctic domain for the remesh (BedMachine outline, buffer and
  boundary classification from :mod:`icepack2_tools.mesh`, so the new mesh
  carries the same physical groups and the calving BC lands on the same
  facets);
* :func:`transfer_state`, which knows this model's checkpoint: which fields
  exist, which are frozen references (``a_ref_mb``, ``N_ref``, ``H_init``),
  the geometry route, and the attributes a restart needs.

The config comes from ``ISMIP7_ADAPT_*`` (``ISMIP7_ADAPT_PRESET=ua`` for the
Antarctic sizes); see ``ADAPTIVE_MESH.md``.
"""

import os

import numpy as np
from firedrake import Constant, Function, FunctionSpace, assemble, conditional, ds, dx, max_value
from firedrake.petsc import PETSc
from mpi4py import MPI

from icepack_tools.geometry import cg1_lift

from icepack_tools.adapt_mesh import (  # noqa: F401  (re-exported for the scripts)
    CRITERIA,
    AdaptMeshConfig,
    Criterion,
    RHO_RATIO,
    calving_front_points,
    cross_mesh_transfer,
    current_element_size_nodal,
    desired_element_size,
    dirac_delta,
    error_to_ele_size,
    facet_midpoints_and_jumps,
    gather_points,
    grounding_line_points,
    nodal_distance_to,
    physical_divergence,
    physical_groups,
    preserve_front,
    rebuild_reference_pressure,
    surface_route_thickness,
)
from icepack_tools.adapt_mesh import remesh_global as _remesh_global

__all__ = ["AdaptMeshConfig", "Criterion", "desired_element_size", "remesh_global",
           "transfer_state", "antarctica_geometry_builder"]


def antarctica_geometry_builder(cfg, data_dir=None):
    r"""A ``build_geometry()`` callable for :func:`icepack_tools.adapt_mesh.
    remesh_global`: the ISMIP7 Antarctic outline, buffered per
    ``ISMIP7_BUFFER_M`` and classified into calving / other segments exactly
    as ``mesh_antarctica.py`` does, so the physical groups match the old mesh."""
    from .mesh import build_gmsh_geometry, classify_boundaries, extract_ice_outline, load_bedmachine_mask
    if data_dir is None:
        # BedMachine sits under ISMIP7_OBS_DATA_ROOT, like every other reader.
        from .runconfig import obs_data_root
        data_dir = obs_data_root()
    mask, x, y = load_bedmachine_mask(data_dir)
    outline = extract_ice_outline(mask, x, y)
    boundaries, names = classify_boundaries(outline, mask, x, y)

    def build():
        build_gmsh_geometry(boundaries, names, cfg.mesh_size_min, cfg.mesh_size_max)
    return build


def remesh_global(mesh, h_des, cfg, out_msh, old_msh, old_sidecar, log=PETSc.Sys.Print):
    r"""Remesh the Antarctic domain onto ``h_des``; the new mesh must expose
    the old mesh's physical groups, and the old sidecar is copied under the
    new basename."""
    builder = antarctica_geometry_builder(cfg) if mesh.comm.rank == 0 else None
    sidecar_out = os.path.join(os.path.dirname(out_msh),
                               "boundary_ids_" + os.path.splitext(os.path.basename(out_msh))[0] + ".json")
    return _remesh_global(mesh, h_des, cfg, out_msh, builder, reference_msh=old_msh,
                          sidecar_in=old_sidecar, sidecar_out=sidecar_out, log=log)


def transfer_state(chk_in, mesh_new, cfg, chk_out, new_msh_basename, bed_sampler,
                   rebuild_aref=False, thickness_sampler=None, log=PETSc.Sys.Print):
    r"""``MapFbetweenMeshes`` for this model's checkpoint; writes ``chk_out``.

    ``bed_sampler(Q_g_new, Q_cg_new)`` returns the bed on the new mesh from
    data (the ``DefineGeometry`` route); ``thickness_sampler`` likewise, used
    only for the initial adaptation (``rebuild_aref``), where the reference takes ALL
    geometry from data. Otherwise ``cfg.geometry`` selects the route.
    Frozen references: ``a_ref_mb`` becomes the transferred physical
    divergence (the forward rebuilds it), ``N_ref`` is rebuilt from the
    transferred ratio, ``H_init`` moves like the thickness. Returns audits.
    """
    import firedrake as fd
    with fd.CheckpointFile(chk_in, "r") as chk:
        mesh_old = chk.load_mesh()
        # smb_elevation_feedback: the next segment's restart guard reads it
        # (forcing.smb_feedback_restart_error), and an adapted continuation
        # is the same chain.
        attrs = {k: chk.get_attr("/", k) for k in
                 ("t_yr", "friction", "geometry_space", "mesh_basename", "lc", "lc_coarse",
                  "buffer_m", "raster_sample", "smb_elevation_feedback") if chk.has_attr("/", k)}
        old = {}
        for name in ("log_friction", "log_fluidity", "fluidity_prior", "thickness", "bed", "surface",
                     "velocity", "membrane_stress", "basal_stress", "H_init", "phi_eff", "C_w0",
                     "N_ref", "a_ref_mb", "levelset", "thickness_dg"):
            try:
                old[name] = chk.load_function(mesh_old, name=name)
            except Exception:
                pass
    log(f"  adapt: transferring {sorted(old)} ({cfg.transfer})")

    geom = attrs.get("geometry_space", "dg0")
    Q_g = FunctionSpace(mesh_new, "DG" if geom == "dg0" else "CG", 0 if geom == "dg0" else 1)
    Qc = FunctionSpace(mesh_new, "CG", 1)
    rr = Constant(RHO_RATIO)
    H_old = old["thickness"]
    s_old = old.get("surface") or Function(H_old.function_space()).interpolate(
        max_value(old["bed"] + H_old, (Constant(1.0) - rr) * H_old))

    new = {}
    b = bed_sampler(Q_g, Qc)
    if rebuild_aref and thickness_sampler is not None:
        H, route = thickness_sampler(Q_g, Qc), "data"        # first run-step: all geometry from data
    elif cfg.geometry == "bh-FROM-sBS":
        H, route = surface_route_thickness(s_old, b, mesh_new, Q_g), "bh-FROM-sBS"
    else:
        H = cross_mesh_transfer(H_old, mesh_new, cfg.transfer if geom == "dg0" else "interpolate", default=cfg.thick_min)
        H.dat.data[:] = np.maximum(H.dat.data_ro, 0.0)
        route = "bs-FROM-hBS"
    if geom == "dg0" and route != "data" and cfg.front_preserve:
        n_fixed = preserve_front(mesh_old, H_old, mesh_new, H)
        log(f"  adapt: front-preserving transfer set {n_fixed} boundary cells")
    s = Function(Q_g).interpolate(max_value(b + H, (Constant(1.0) - rr) * H))
    new["thickness"], new["bed"], new["surface"] = H, b, s
    log(f"  adapt: geometry route {route}")

    for name, f in old.items():
        if name in ("thickness", "bed", "surface"):
            continue
        if name == "a_ref_mb" and rebuild_aref:
            continue
        if name == "H_init" and rebuild_aref:
            new[name] = Function(Q_g).assign(H)                # the t=0 extent anchor IS the initial thickness
            continue
        if name == "levelset":
            # The shared contract (icepack_tools.adapt_mesh.cross_mesh_transfer):
            # a level set is carried as its P1 lift and lands back in DG0, so
            # the front position keeps its sub-cell accuracy instead of being
            # quantised to the old cell size by DG0-to-DG0 point evaluation.
            # Outside the old mesh is water (positive distance); the forward's
            # extent anchor re-solves the eikonal problem from the transferred
            # extent, as the calving project reinitialises after a remesh. No
            # run loads this field back yet: a calving forward rebuilds its
            # level set from the transferred thickness and H_init, so carrying
            # it here future-proofs the checkpoint contents rather than
            # changing a result today.
            new[name] = cross_mesh_transfer(cg1_lift(f), mesh_new, "interpolate", default=1e6,
                                            element=f.function_space().ufl_element())
            continue
        how = cfg.transfer if (name in ("a_ref_mb", "thickness_dg", "H_init") and geom == "dg0") else "interpolate"
        new[name] = cross_mesh_transfer(f, mesh_new, how, default=cfg.thick_min if name in ("H_init", "thickness_dg") else 0.0)

    comm = mesh_new.comm
    if "a_ref_mb" in old and not rebuild_aref and "velocity" in old:
        P_old = physical_divergence(mesh_old, H_old, old["velocity"], old["a_ref_mb"])
        new["phys_div"] = cross_mesh_transfer(P_old, mesh_new, "interpolate")   # smooth: the spikes live in a_ref
        new.pop("a_ref_mb", None)
        lo = comm.allreduce(float(P_old.dat.data_ro.min()) if P_old.dat.data_ro.size else np.inf, op=MPI.MIN)
        hi = comm.allreduce(float(P_old.dat.data_ro.max()) if P_old.dat.data_ro.size else -np.inf, op=MPI.MAX)
        log(f"  adapt: a_ref_mb replaced by the transferred physical divergence P in [{lo:.2f}, {hi:.2f}] m/yr, "
            f"net {assemble(P_old * dx) / 1e9:+.1f} Gt/yr; the forward rebuilds a_ref on this mesh")
    if "N_ref" in old and "N_ref" in new:
        new["N_ref"] = rebuild_reference_pressure(H_old, s_old, old["N_ref"], H, s, mesh_new)
        log("  adapt: N_ref rebuilt from the transferred N_eff/N_ref ratio on the new geometry")

    def _front(Hf):
        m = Hf.function_space().mesh()
        return assemble(Hf * ds(domain=m)) / max(assemble(conditional(Hf > 1.0, 1.0, 0.0) * ds(domain=m)), 1.0)
    m_old, m_new = assemble(H_old * dx) / 1e9, assemble(H * dx) / 1e9
    audit = {"volume_old_km3": m_old, "volume_new_km3": m_new,
             "volume_change_pct": 100.0 * (m_new - m_old) / max(m_old, 1e-30),
             "front_h_old": _front(H_old), "front_h_new": _front(H)}
    log(f"  adapt: front <h> {audit['front_h_old']:.1f} -> {audit['front_h_new']:.1f} m")
    log(f"  adapt: ice volume {m_old:.6e} -> {m_new:.6e} km3 ({audit['volume_change_pct']:+.4f}%)")

    with fd.CheckpointFile(chk_out, "w") as chk:
        chk.save_mesh(mesh_new)
        for name, f in new.items():
            f.rename(name)
            chk.save_function(f, name=name)
        for k, v in attrs.items():
            chk.set_attr("/", k, v)
        chk.set_attr("/", "mesh_basename", new_msh_basename)
        chk.set_attr("/", "adapted_from", os.path.basename(chk_in))
        if rebuild_aref:
            chk.set_attr("/", "adapted_initial", 1)
    return audit
