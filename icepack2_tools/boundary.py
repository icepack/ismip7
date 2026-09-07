r"""Boundary-id sidecar resolution and validation, shared by every solver.

The gmsh physical-line count depends on resolution AND on the outline buffer,
so a ``boundary_ids.json`` built for one mesh is NOT valid for another. The
2500 m sidecar (35 groups) paired with the 32 km mesh (639 groups) tagged
1596 km of a 32 000 km boundary, so 95% of the ice front ran stress-free
instead of carrying 1/2(rho_I g h^2 - rho_W g d^2) - legal and silent for an
entire campaign.

Every reader therefore goes through :func:`load_boundary_ids`, which prefers
the per-mesh sidecar ``boundary_ids_<mesh stem>.json`` written next to the mesh
by ``mesh_antarctica.py`` and HARD-ERRORS on any mismatch with the mesh it is
about to be used on.
"""

import json
import os

# Firedrake is imported inside load_boundary_ids so that the sidecar-name rule
# stays importable from the (deliberately Firedrake-free, fast) preflight.


def sidecar_path(mesh_dir, mesh_hint=None):
    r"""Path of the boundary-id sidecar to use.

    ``ISMIP7_BNDIDS`` wins; otherwise the per-mesh
    ``boundary_ids_<mesh stem>.json`` if it exists, where the stem comes from
    ``mesh_hint`` (a .msh path or bare basename) or ``ISMIP7_MESH``; otherwise
    the shared ``boundary_ids.json``.
    """
    explicit = os.environ.get("ISMIP7_BNDIDS")
    if explicit:
        return explicit
    hint = mesh_hint or os.environ.get("ISMIP7_MESH") or ""
    stem = os.path.splitext(os.path.basename(hint))[0]
    if stem:
        per_mesh = os.path.join(mesh_dir, f"boundary_ids_{stem}.json")
        if os.path.exists(per_mesh):
            return per_mesh
    return os.path.join(mesh_dir, "boundary_ids.json")


def load_boundary_ids(mesh, mesh_dir, mesh_hint=None, print_coverage=True):
    r"""Resolve, load and validate the sidecar for ``mesh``.

    Returns ``(bnd_ids, calving_ids, sidecar_fn)``. Raises ``ValueError`` if
    any exterior marker on the mesh is unclassified by the sidecar, or if any
    sidecar id is absent from the mesh: an unclassified marker gets NO
    calving-terminus back-pressure, which is exactly the failure that used to
    be silent.
    """
    from firedrake import Constant, assemble, ds
    from firedrake.petsc import PETSc

    sidecar_fn = sidecar_path(mesh_dir, mesh_hint)
    with open(sidecar_fn) as f:
        bnd_ids = json.load(f)
    calving_ids = tuple(bnd_ids["calving"])

    mesh_markers = set(int(i) for i in mesh.exterior_facets.unique_markers)
    tagged = set(int(i) for i in calving_ids) | set(
        int(i) for i in bnd_ids.get("other", [])
    )
    unclassified = mesh_markers - tagged
    absent = tagged - mesh_markers
    if unclassified or absent:
        raise ValueError(
            f"boundary_ids sidecar does not match this mesh.\n"
            f"  sidecar : {sidecar_fn}\n"
            f"  mesh has {len(mesh_markers)} exterior markers; sidecar names "
            f"{len(tagged)}\n"
            f"  {len(unclassified)} mesh marker(s) unclassified"
            f"{' (e.g. ' + str(sorted(unclassified)[:5]) + ')' if unclassified else ''}\n"
            f"  {len(absent)} sidecar id(s) absent from the mesh"
            f"{' (e.g. ' + str(sorted(absent)[:5]) + ')' if absent else ''}\n"
            f"  Unclassified markers get NO calving-terminus back-pressure.\n"
            f"  Regenerate with:  ISMIP7_MESH=<the .msh> python "
            f"antarctica/scripts/make_boundary_ids.py"
        )

    if print_coverage and calving_ids:
        len_all = float(assemble(Constant(1.0) * ds(domain=mesh)))
        len_cal = float(assemble(Constant(1.0) * ds(calving_ids, domain=mesh)))
        PETSc.Sys.Print(
            f"  Calving terminus BC on {len_cal/1e3:.0f} km of "
            f"{len_all/1e3:.0f} km exterior boundary "
            f"({100*len_cal/max(len_all, 1e-9):.0f}%), "
            f"{len(calving_ids)}/{len(mesh_markers)} markers "
            f"({os.path.basename(sidecar_fn)})"
        )
    return bnd_ids, calving_ids, sidecar_fn
