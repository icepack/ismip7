r"""Geometry representation helpers shared by the inversion and the forward.

The prognostic geometry (``h``, ``s``, ``b``) lives in the space selected by
``ISMIP7_GEOMETRY_SPACE`` - DG0 by default, so that the momentum solve and the
mass transport use one thickness field and the calving-terminus traction cannot
disagree with the boundary flux. See ``GEOMETRY_DISCRETIZATION.md``.

Three operations need care under DG0 and are collected here so the inversion
and the forward cannot drift apart:

* :func:`sample_to_geometry` - get a raster onto the geometry space as a CELL
  AVERAGE rather than a centroid point sample.
* :func:`cg1_lift` - a bounded, volume-preserving CG1 reconstruction, for the
  few places that genuinely need a pointwise gradient of a cell-wise field.
* :func:`surface_slope` - ``grad(s)`` that works for CG1 or DG0 surfaces.

Under ``dg1`` the geometry is cell-wise linear on vertex nodes
(:func:`geometry_function_space`); the vertex sample is then the CG1
interpolant itself, held in DG1, and the cell rules (the front cells, the lake
fix) decide on the cell and write all three of its dofs.
"""

import weakref

import numpy as np
from mpi4py import MPI

from firedrake import (
    Constant,
    Function,
    FunctionSpace,
    TestFunction,
    assemble,
    dx,
    grad,
    max_value,
)

from .cellview import CellView, geometry_degree


_LUMPED_MASS = weakref.WeakKeyDictionary()


def _lumped_mass(mesh, Q_cg):
    r"""Lumped CG1 mass vector for this mesh, assembled once.

    It depends only on the mesh, and cg1_lift now runs once per forcing
    callback (~1 per timestep), so re-assembling it every call is half the
    per-call cost for nothing. Keyed on the MESH, which lives for the whole
    run: the CG1 space object cg1_lift builds is transient (nothing outside
    holds it), so keying on that would evict the entry on every call. The
    value is the plain array, not the assembled Cofunction, which would keep
    a strong reference back to its own key and pin the mesh forever.
    """
    cached = _LUMPED_MASS.get(mesh)
    if cached is None:
        cached = assemble(TestFunction(Q_cg) * dx).dat.data_ro.copy()
        _LUMPED_MASS[mesh] = cached
    return cached


def cg1_lift(f):
    r"""Lumped-mass CG1 reconstruction of a cell-wise (DG0) field.

    Each CG1 node takes the area-weighted mean of the adjacent cell values.
    Volume-preserving (``int f dx`` is exact) and a convex combination, so it
    cannot overshoot - unlike an L2 projection, which does (projecting the DG0
    thermomechanical fluidity prior to CG1 produced a NEGATIVE fluidity, min
    -9.88 against a DG0 range of [1.0, 446.7]).

    It is NOT unbiased on the domain boundary, where the stencil is one-sided
    and pulls boundary nodes toward interior values. On an unbuffered mesh the
    boundary is the calving front, and the terminus traction goes as ``h^2``,
    so using this on the thickness inflated the front (105 -> 200 m at 32 km)
    and multiplied the outflux by ~4.7. Keep it out of the momentum residual
    and out of every flux; it is for diagnostics, fixed reference scalings, and
    parameterizations only.
    """
    mesh = f.function_space().mesh()
    Q_cg = FunctionSpace(mesh, "CG", 1)
    lumped = _lumped_mass(mesh, Q_cg)
    rhs = assemble(TestFunction(Q_cg) * f * dx)
    out = Function(Q_cg)
    out.dat.data[:] = rhs.dat.data_ro / lumped
    return out


def surface_slope(s):
    r"""``grad(s)`` valid for a CG1, DG0 or DG1 surface.

    A DG0 surface has an identically zero cell gradient - UFL folds it away -
    because its slope lives entirely in the inter-cell jumps, which the
    momentum balance picks up weakly through its ``jump(s, nu)`` facet term.
    A DG1 surface carries part of its slope in those jumps too. Callers that
    need a POINTWISE slope (a friction anchor, a melt parameterization) get
    one from a CG1 reconstruction of any discontinuous surface instead.
    """
    if s.function_space().ufl_element().family() == "Discontinuous Lagrange":
        return grad(cg1_lift(s))
    return grad(s)


def geometry_function_space(mesh, space):
    r"""The function space of the geometry under ``ISMIP7_GEOMETRY_SPACE``
    ``space``: CG1 for ``cg1``, DG0 for ``dg0``, and for ``dg1`` the DG1
    space with its nodes on the vertices.

    Firedrake's default DG1 variant puts the three nodes at interior points of
    the cell, where a nodal floor does not bound the field at the vertices and
    a vertex sample is not the CG1 interpolant. On the vertices a non-negative
    nodal vector is a non-negative field, the cell mean is the mean of the
    three nodal values, and the CG1 interpolant is held exactly.
    """
    if space == "cg1":
        return FunctionSpace(mesh, "CG", 1)
    if space == "dg0":
        return FunctionSpace(mesh, "DG", 0)
    if space == "dg1":
        return FunctionSpace(mesh, "DG", 1, variant="equispaced")
    raise ValueError(f"no geometry space {space!r}")


def _lattice_bary(n):
    r"""Barycentric centroids of the ``n**2`` equal-area sub-triangles of a
    uniform n-fold split of a triangle, shape ``(n*n, 3)``."""
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    up = (i + j) <= n - 1                     # n(n+1)/2 upright sub-triangles
    dn = (i + j) <= n - 2                     # n(n-1)/2 inverted ones
    a = np.concatenate([(i[up] + 1.0 / 3.0) / n, (i[dn] + 2.0 / 3.0) / n])
    b = np.concatenate([(j[up] + 1.0 / 3.0) / n, (j[dn] + 2.0 / 3.0) / n])
    return np.column_stack([a, b, 1.0 - a - b])


def _lattice_windows(datasets, tri):
    r"""For each raster, the window covering the triangles ``tri`` (with a
    two-pixel margin) as ``(array, row0, col0, transform)``. Masked / nodata
    pixels of a float raster become NaN; an integer raster is read as is."""
    from rasterio.windows import from_bounds

    flat = tri.reshape(-1, 2)
    out = []
    for dataset in datasets:
        bnd = dataset.bounds
        rx, ry = abs(dataset.res[0]), abs(dataset.res[1])
        win = from_bounds(max(flat[:, 0].min() - 2 * rx, bnd.left),
                          max(flat[:, 1].min() - 2 * ry, bnd.bottom),
                          min(flat[:, 0].max() + 2 * rx, bnd.right),
                          min(flat[:, 1].max() + 2 * ry, bnd.top),
                          transform=dataset.transform)
        win = win.round_lengths(op="ceil").round_offsets(op="floor")
        if np.issubdtype(np.dtype(dataset.dtypes[0]), np.integer):
            arr = dataset.read(1, window=win)
        else:
            arr = dataset.read(1, window=win, masked=True)
            arr = np.ma.filled(arr.astype("f4"), np.nan)
        out.append((arr, int(win.row_off), int(win.col_off), dataset.transform))
    return out


def cell_samples(datasets, Q_dg, reduce, nmax=64, chunk=2048, density=1.0, subset=None):
    r"""Look rasters up on an equal-area lattice in every owned cell of the
    DG0 space ``Q_dg``.

    Each cell is split into ``n**2`` equal-area sub-triangles with
    ``n = ceil(density * sqrt(cell_area / pixel_area))`` (capped at
    ``nmax``; the pixel is the first raster's), and every raster is looked
    up at each sub-triangle centroid by pixel containment. ``reduce(idx,
    values)`` is called for groups of cells sharing ``n``: ``idx`` indexes
    the owned cells and ``values`` holds one ``(len(idx), n*n)`` array per
    raster. ``subset``, owned cell indices, restricts the cells looked up.
    Returns the owned cell -> dof map, the indexing of ``idx``. Rank-local,
    no collectives.
    """
    from rasterio.transform import rowcol

    mesh = Q_dg.mesh()
    X = mesh.coordinates.dat.data_ro_with_halos[:, :2]
    cells = mesh.coordinates.cell_node_map().values          # owned cells -> nodes
    dofs = Q_dg.cell_node_map().values[:, 0]                 # owned cells -> dof
    nc = cells.shape[0]
    cell_ids = np.arange(nc) if subset is None else np.asarray(subset, dtype=int)
    if cell_ids.size == 0:
        return dofs
    tri = X[cells]                                           # (nc, 3, 2)
    v0, v1, v2 = tri[:, 0], tri[:, 1], tri[:, 2]
    area = 0.5 * np.abs((v1[:, 0] - v0[:, 0]) * (v2[:, 1] - v0[:, 1])
                        - (v2[:, 0] - v0[:, 0]) * (v1[:, 1] - v0[:, 1]))
    apix = abs(datasets[0].res[0] * datasets[0].res[1])
    n = np.clip(np.ceil(density * np.sqrt(area / apix)).astype(int), 1, nmax)
    windows = _lattice_windows(datasets, tri[cell_ids])

    for nn in np.unique(n[cell_ids]):
        sel = cell_ids[n[cell_ids] == nn]
        lam = _lattice_bary(int(nn))                         # (k, 3)
        for s0 in range(0, sel.size, chunk):
            idx = sel[s0:s0 + chunk]
            P = np.einsum("kb,cbd->ckd", lam, tri[idx])      # (c, k, 2)
            values = []
            for arr, row0, col0, transform in windows:
                rows, cols = rowcol(transform, P[..., 0].ravel(), P[..., 1].ravel())
                r = np.clip(np.asarray(rows) - row0, 0, arr.shape[0] - 1)
                c = np.clip(np.asarray(cols) - col0, 0, arr.shape[1] - 1)
                values.append(arr[r, c].reshape(idx.size, -1))
            reduce(idx, values)
    return dofs


def raster_cell_mean(dataset, Q_dg, floor=None, nmax=64, chunk=2048):
    r"""Mean of a rasterio raster over each cell of a DG0 space.

    Each owned cell is split into ``n**2`` equal-area sub-triangles with
    ``n = ceil(sqrt(cell_area / pixel_area))`` (capped at ``nmax``), the raster
    is looked up at every sub-triangle centroid by *pixel containment*, and the
    cell value is the mean. Sampling density therefore tracks pixel density: a
    2 km cell on 500 m BedMachine gets 9 samples for its ~7 pixels, a 20 km
    interior cell gets 729 for its ~1600. Contrast :func:`sample_to_geometry`
    with ``method="vertex"``, which reads three pixels per cell whatever its
    size.

    Denser sampling was expected to give a smoother DG0 field. It does not:
    measured, the cell mean is ROUGHER across the facet jumps that ARE the
    driving stress. See :func:`sample_to_geometry` for the numbers and for why
    ``"vertex"`` remains the default.

    One raster window covering the rank's cells is read. With a locality-
    preserving partition that window is small; with PETSc's ``simple``
    partitioner every rank's cells are scattered continent-wide and the window
    is the whole raster (711 MB float32 for BedMachine), which is tolerated
    rather than fixed here.

    Masked / nodata pixels are excluded from the mean; a cell with no valid
    sample is left NaN for the caller to fill.
    """
    means = np.full(Q_dg.mesh().coordinates.cell_node_map().values.shape[0], np.nan)

    def reduce(idx, values):
        with np.errstate(all="ignore"):
            means[idx] = np.nanmean(values[0], axis=1, dtype="f8")

    dofs = cell_samples([dataset], Q_dg, reduce, nmax=nmax, chunk=chunk)
    out = Function(Q_dg)
    if means.size == 0:
        return out
    if floor is not None:
        means = np.maximum(means, floor)
    out.dat.data[dofs] = means
    return out


def raster_patch_mean(dataset, Q_dg, Q_cg):
    r"""Mean of a raster over the cells around each node of a CG1 space.

    The cell means of :func:`raster_cell_mean` (valid pixels only) are
    averaged over the cells sharing each node, weighted by area and counting
    only the cells that hold a valid pixel; a node with none around it gets
    0. Returns ``(field, observed)``, ``observed`` a boolean array over the
    node dofs. Assembled, so collective.
    """
    cm = raster_cell_mean(dataset, Q_dg)
    ok = np.isfinite(cm.dat.data_ro)
    valid = Function(Q_dg)
    valid.dat.data[:] = ok.astype(float)
    cm.dat.data[~ok] = 0.0
    phi = TestFunction(Q_cg)
    num = assemble(phi * valid * cm * dx).dat.data_ro
    den = assemble(phi * valid * dx).dat.data_ro
    observed = den > 0.0
    out = Function(Q_cg)
    out.dat.data[:] = np.where(observed, num / np.where(observed, den, 1.0), 0.0)
    return out, observed


def sample_to_geometry(raster, Q_g, Q_cg, floor=None, method="vertex"):
    r"""Sample a raster onto the geometry space ``Q_g`` as a cell average.

    ``raster`` is a rasterio dataset (a ``netcdf:file:variable`` handle).

    ``method`` selects the DG0 cell value (``runconfig.RASTER_SAMPLES``):
    ``"vertex"`` projects the CG1 vertex interpolant (three pixels per cell);
    ``"cell_mean"`` is :func:`raster_cell_mean`, the raster's mean over the
    cell. The front sampling (``"vertex_front"``) samples one raster as
    ``"vertex"``; its front cells are :func:`sample_bed_thickness`'s. Under
    CG1 geometry there is no cell, so ``method`` is ignored.

    For CG1 geometry this is just the nodal interpolant. For DG0 it is NOT the
    obvious ``icepack.interpolate(raster, Q_g)``: a DG0 dof sits at the cell
    centroid, so that would take a ONE-POINT sample of a 500 m BedMachine
    raster per (at 32 km) 32 km cell. Since the DG0 driving stress is entirely
    the facet jump in ``s``, that sampling noise is read as slope. Measured at
    32 km, centroid sampling against the cell average:

        rms |jump s|     366 m  vs  257 m       (42% rougher)
        peakedness       2.04   vs  1.49
        front <h>        209 m  vs  153 m       (BedMachine ice front: 145
                                                 all / 167 floating, median 152)
        |driving force|  7.7%   vs  1.0% from the CG1 value

    and the rough version failed to converge in 200 Newton iterations. The L2
    projection of the CG1 interpolant IS the cell average, which is what a DG0
    field means, so use that.

    MEASURED: ``"cell_mean"`` is rougher, and ``"vertex"`` stays the default.
    The cell mean is the more faithful average of the raster, but faithfulness
    is not what the DG0 driving stress wants. Neighbouring cells share two of
    their three vertex samples, so ``"vertex"`` damps the jump between them by
    construction, while two independent cell means do not. Against vertex
    sampling the cell mean raised interior surface jumps by 6% and bed and
    thickness jumps by 35%, and at 2 km the momentum solve did not converge
    within 60 minutes. It does classify flotation better (at 32 km the
    misclassified fraction falls from 9.1% to 3.2%), which is why the knob is
    kept rather than removed. ``antarctica/scripts/probe_raster_sampling.py``
    reproduces the comparison.

    ``floor`` optionally clamps the field from below (thickness >= h_clamp).
    The two paths apply it at opposite ends of the averaging and so disagree
    when ``floor > 0``: ``"vertex"`` clamps the CG1 interpolant and then
    projects (clamp, then average), while ``"cell_mean"`` averages the raster
    over the cell and clamps that average (average, then clamp). Only the
    second can return exactly ``floor``; the first returns it only where every
    vertex is at or below the floor. They coincide wherever the raster already
    exceeds the floor, which is why the inversion sees no difference
    (``ISMIP7_H_CLAMP`` defaults to 0 and BedMachine thickness is
    non-negative); the gap is reachable from the budd_legacy cold start, whose
    ``h_clamp_init`` is 10 m.
    """
    import icepack
    from .runconfig import raster_base_method
    method = raster_base_method(method)
    dataset = raster
    raster_fn = lambda sp: icepack.interpolate(dataset, sp)  # noqa: E731

    is_dg0 = Q_g.ufl_element() != Q_cg.ufl_element()
    linear = is_dg0 and geometry_degree(Q_g) == 1
    if linear and method == "cell_mean":
        raise ValueError(
            "raster sampling cell_mean has no DG1 rule; a DG1 geometry samples "
            "with vertex or vertex_front")
    if is_dg0 and method == "cell_mean":
        out = raster_cell_mean(dataset, Q_g, floor=floor)
        bad = np.isnan(out.dat.data_ro)
        # The fill is a collective L2 projection, so the branch must be taken
        # by every rank or none: a rank-local `bad.any()` hangs under MPI when
        # the nodata cells (or the owned cells) are not spread over all ranks.
        if Q_g.mesh().comm.allreduce(bool(bad.any()), op=MPI.LOR):
            # No valid pixel under the cell (nodata): take the vertex value.
            fill = Function(Q_g).project(raster_fn(Q_cg))
            out.dat.data[bad] = fill.dat.data_ro[bad]
            if floor is not None:
                out.dat.data[bad] = np.maximum(out.dat.data_ro[bad], floor)
        return out
    if method == "patch":
        # The lumped L2 projection of the raster onto CG1 (patchmean.py):
        # the DG0 projection below then carries the raster's mass exactly.
        # Nodes whose patch the raster does not cover take the vertex value.
        from .patchmean import quadrature_points, raster_patch_mean, fill_uncovered
        from .runconfig import patch_qdegree
        degree = patch_qdegree()
        Qs, x, y = quadrature_points(Q_cg.mesh(), degree)
        field, _cov, _ = raster_patch_mean(dataset, Qs, x, y, Q_cg, degree)
        # Collective: every rank asks before any reads the fallback.
        n_bad = Q_cg.mesh().comm.allreduce(
            int((~np.isfinite(field.dat.data_ro)).sum()), op=MPI.SUM)
        if n_bad:
            fill_uncovered(field, raster_fn(Q_cg))
        if floor is not None:
            field = Function(Q_cg).interpolate(max_value(field, Constant(floor)))
        if not is_dg0:
            return field
        if Q_g.ufl_element().degree() == 1:
            # A vertex-node DG1 geometry (PR 169) holds the CG1 field exactly,
            # and its cell means are the vertex means, so the mass is kept.
            return Function(Q_g).interpolate(field)
        return Function(Q_g).project(field)
    if method not in ("vertex", "cell_mean"):
        raise ValueError(f"unknown raster sampling method {method!r}")

    field = raster_fn(Q_cg)
    if floor is not None:
        field = Function(Q_cg).interpolate(max_value(field, Constant(floor)))
    if not is_dg0:
        return field if isinstance(field, Function) else Function(Q_cg).interpolate(field)
    if linear:
        # The CG1 interpolant lies in the vertex-node DG1 space, so it is held
        # exactly; its cell means are the DG0 cell averages below.
        return Function(Q_g).interpolate(field)
    return Function(Q_g).project(field)


# BedMachine Antarctica mask values: 0 ocean, 1 ice-free land, 2 grounded ice,
# 3 floating ice, 4 Lake Vostok.
LAKE_MASK = 4


def raise_bed_to_lake_ice_base(b, H, bm_fn, Q_g, Q_cg, method="vertex"):
    r"""Set the bed to the ice base ``s - H`` under BedMachine's subglacial lake.

    Under ``mask == 4`` (Lake Vostok) BedMachine's ``bed`` is the lake FLOOR
    and its ``thickness`` the ice alone, so ``b + H`` -- the surface every part
    of this model builds -- sits below BedMachine's surface by the lake's water
    column: a bowl a median 266 m and up to 916 m deep over 15,200 km2, whose
    walls carry driving stresses near 1 MPa into the friction anchor and the
    momentum balance. The ice base ``s - H`` taken from BedMachine's own
    surface removes it, and the ice over the lake keeps its thickness.

    A cell takes the ice base when any of its vertices lies on the lake, so the
    correction covers the lake's whole footprint on the mesh; in a cell only
    partly over the lake the cell-averaged ``s - H`` carries only that part of
    the water column. Under DG1 the whole cell takes it, at each of its nodes.
    ``b`` and ``H`` must come from :func:`sample_to_geometry` with the same
    ``method``. Modifies ``b`` in place and returns the global number of cells
    (or nodes, under CG1 geometry) it changed.
    """
    import icepack
    import rasterio
    from firedrake import conditional, gt, lt

    s_bm = sample_to_geometry(
        rasterio.open(f"netcdf:{bm_fn}:surface"), Q_g, Q_cg, method=method)
    mask = icepack.interpolate(
        rasterio.open(f"netcdf:{bm_fn}:mask"), Q_cg, method="nearest")
    on_lake = Function(Q_cg).interpolate(
        conditional(lt(abs(mask - LAKE_MASK), 0.5), 1.0, 0.0))
    # Into DG0 the vertex indicator arrives as its vertex mean, positive when
    # any vertex is on the lake; under CG1 geometry it is the indicator itself.
    # Under DG1 the cell decides, through the same DG0 mean.
    view = CellView(Q_g) if Q_g.ufl_element() != Q_cg.ufl_element() else None
    linear = view is not None and view.degree == 1
    touches = Function(view.Q0 if linear else Q_g).interpolate(on_lake)
    new_b = Function(Q_g).interpolate(conditional(gt(touches, 0.0), s_bm - H, b))
    if linear:
        moved = np.abs(view.nodes(new_b) - view.nodes(b)).max(axis=1)
    else:
        moved = np.abs(new_b.dat.data_ro - b.dat.data_ro)
    changed = int(np.count_nonzero(moved > 1e-6))
    b.assign(new_b)
    return Q_g.mesh().comm.allreduce(changed, op=MPI.SUM)


# ── Front cells (issue #167) ──────────────────────────────────────────────
#
# Vertex sampling gives a cell the BedMachine front crosses the mean of its
# vertex samples, and an ocean vertex samples zero, so the cell holds a
# fraction of the front's thickness: on IU's 2 km buffered mesh the front band
# held 40 m against BedMachine's 163 m floating front and carried 343 of about
# 1,200 Gt/yr. Even on a mesh whose edges follow the front (the _frontbm
# meshes) a vertex on the front samples a blend of ice and water. Under the
# front sampling (``vertex_front``) BedMachine's own mask decides which cells
# hold ice at the marine front, and the front cells take its thickness and
# bed over their ice.

# A cell holds ice when at least this share of its lattice samples are ice:
# the area-preserving choice on a mesh that does not follow the front; on one
# that does the share is about 0 or 1.
FRONT_ICE_FRACTION = 0.5
# Cells whose pixel-density samples are mixed are looked up again on a lattice
# this many times as fine in each direction: at pixel density a 1 km cell gets
# four samples, too few to tell a cell 64 % ice from one 25 % ice.
FRONT_LATTICE_DENSITY = 4.0


def front_cells(H, b, bm_fn, frac=FRONT_ICE_FRACTION):
    r"""Rebuild the marine front of the DG0 or DG1 geometry ``H``, ``b`` from
    BedMachine's own mask, in place.

    Every owned cell is looked up on the :func:`cell_samples` lattice in
    BedMachine's mask, bed and thickness, and its samples classified
    (``obs_icemask.classify`` of BedMachine's ice). A cell holds ice when at
    least ``frac`` of its samples are ice. Then a cell without ice that holds
    a marine sample is water, ``H = 0``; a cell with ice that shares a vertex
    with a water cell, or holds a marine sample itself, is rebuilt: its
    ``H`` and ``b`` are BedMachine's means over the ice samples it holds.
    Every other cell is untouched, so interior ice, land margins and the open
    buffer keep their vertex samples.

    Under DG1 the cells are classified the same way and a water or rebuilt
    cell takes its value at all three nodes, so it is constant and its mean
    is the DG0 cell's value; every other cell keeps its linear sample.

    Returns the global counts, collectively.
    """
    import rasterio
    from .front import vertex_neighbours
    from .obs_icemask import ICE, MARINE, _BM_ICE, classify

    view = CellView(H.function_space())
    Q_g = view.Q0
    ncell = Q_g.mesh().coordinates.cell_node_map().values.shape[0]
    f_ice = np.zeros(ncell)
    marine = np.zeros(ncell, bool)
    h_ice = np.full(ncell, np.nan)
    b_ice = np.full(ncell, np.nan)

    def reduce(idx, values):
        bm_mask, bed, thk = values
        cls = classify(np.isin(bm_mask, _BM_ICE), bm_mask, bed)
        ice = cls == ICE
        held = ice & np.isfinite(thk) & (thk > 0)
        n = held.sum(axis=1)
        f_ice[idx] = ice.mean(axis=1)
        marine[idx] = (cls == MARINE).any(axis=1)
        with np.errstate(all="ignore"):
            h_ice[idx] = np.where(held, thk, 0.0).sum(axis=1) / n
            b_ice[idx] = np.where(held, bed, 0.0).sum(axis=1) / n

    with rasterio.open(f"netcdf:{bm_fn}:mask") as m, \
            rasterio.open(f"netcdf:{bm_fn}:bed") as bd, \
            rasterio.open(f"netcdf:{bm_fn}:thickness") as th:
        dofs = cell_samples([m, bd, th], Q_g, reduce)
        mixed = np.flatnonzero(((f_ice > 0.0) & (f_ice < 1.0))
                               | (marine & (f_ice > 0.0)))
        cell_samples([m, bd, th], Q_g, reduce, density=FRONT_LATTICE_DENSITY,
                     subset=mixed)

    n_own = Q_g.dof_dset.size

    def per_dof(a, empty):
        out = np.full(n_own, empty, dtype=a.dtype)
        out[dofs] = a
        return out

    F, MAR = per_dof(f_ice, 0.0), per_dof(marine, False)
    HI, BI = per_dof(h_ice, np.nan), per_dof(b_ice, np.nan)
    ice = F >= frac
    held = ice & np.isfinite(HI)
    water = ~ice & MAR
    rebuilt = held & (vertex_neighbours(Q_g)(water) | MAR)

    h0 = np.array(view.means(H))
    view.set_cells(H, water, 0.0)
    view.set_cells(H, rebuilt, HI)
    view.set_cells(b, rebuilt, BI)
    counts = {
        "rebuilt": int(rebuilt.sum()),
        "rebuilt_thicker": int((rebuilt & (HI > h0)).sum()),
        "water": int(water.sum()),
        "emptied": int((water & (h0 > 0.0)).sum()),
    }
    comm = Q_g.mesh().comm
    return {k: comm.allreduce(v, op=MPI.SUM) for k, v in counts.items()}


def front_changed_nodes(H_new, H_old, Q_cg, tol=1.0):
    r"""The owned dofs of the continuous space ``Q_cg`` that touch a cell
    whose thickness differs between ``H_new`` and ``H_old`` by more than
    ``tol`` m, as a boolean array: the nodes of the cells a front sampling
    rebuilt or emptied (issue #167). Collective."""
    from firedrake import Function, TestFunction, assemble, conditional, dx, gt
    changed = Function(H_new.function_space()).interpolate(
        conditional(gt(abs(H_new - H_old), tol), 1.0, 0.0))
    w = TestFunction(Q_cg)
    return np.asarray(assemble(w * changed * dx).dat.data_ro > 0.0)


def sample_bed_thickness(bm_fn, Q_g, Q_cg, floor=None, method="vertex"):
    r"""BedMachine's bed and thickness on the geometry space, as
    ``(b, H, counts)``.

    Both are :func:`sample_to_geometry` under ``method``, the thickness
    floored at ``floor``. Under the front sampling (``vertex_front``) the
    marine front cells are then rebuilt (:func:`front_cells`) and the floor
    applied again; ``counts`` are its global counts, None otherwise.
    Collective.
    """
    import rasterio
    from .runconfig import raster_front

    b = sample_to_geometry(
        rasterio.open(f"netcdf:{bm_fn}:bed"), Q_g, Q_cg, method=method)
    H = sample_to_geometry(
        rasterio.open(f"netcdf:{bm_fn}:thickness"), Q_g, Q_cg, floor=floor,
        method=method)
    if not raster_front(method):
        return b, H, None
    if Q_g.ufl_element() == Q_cg.ufl_element():
        raise ValueError(
            f"raster sampling {method} rebuilds DG0 front cells; this "
            f"geometry is CG1")
    counts = front_cells(H, b, bm_fn)
    if floor is not None:
        H.interpolate(max_value(H, Constant(floor)))
    return b, H, counts


def remove_unconnected_ice(H, b, Q_g, h_min=1.0):
    r"""Remove ice that no chain of ice cells connects to grounded ice.

    Every ice cell (``H > h_min``) reachable from a grounded cell through
    ice cells sharing a facet is kept; the rest (icebergs and detached
    shelf fragments) is set to zero thickness. Modifies ``H`` in place and
    returns the global count of cells emptied and their area. DG0 only;
    collective.
    """
    import ufl
    from .front import facet_neighbours
    from .grounding import height_above_flotation

    if Q_g.ufl_element().degree() != 0:
        raise ValueError("remove_unconnected_ice needs a DG0 geometry")
    mesh = Q_g.mesh()
    ice = np.asarray(H.dat.data_ro > h_min)
    haf = Function(Q_g).interpolate(height_above_flotation(H, b)).dat.data_ro
    grounded = ice & (haf > 0.0)
    neighbours_of = facet_neighbours(Q_g)
    reached = grounded.copy()
    while True:
        new = ice & ~reached & neighbours_of(reached)
        if not mesh.comm.allreduce(bool(new.any()), op=MPI.LOR):
            break
        reached |= new
    berg = ice & ~reached
    area = Function(Q_g).interpolate(ufl.CellVolume(mesh)).dat.data_ro
    n = mesh.comm.allreduce(int(berg.sum()), op=MPI.SUM)
    a = mesh.comm.allreduce(float(area[berg].sum()), op=MPI.SUM)
    if n:
        H.dat.data[berg] = 0.0
    return n, a
