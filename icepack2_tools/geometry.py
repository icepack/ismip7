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
"""

import weakref

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
    r"""``grad(s)`` valid for a CG1 *or* DG0 surface.

    A DG0 surface has an identically zero cell gradient - UFL folds it away -
    because its slope lives entirely in the inter-cell jumps, which the
    momentum balance picks up weakly through its ``jump(s, nu)`` facet term.
    Callers that need a POINTWISE slope (a friction anchor, a melt
    parameterization) get one from a CG1 reconstruction instead.
    """
    if s.function_space().ufl_element().degree() == 0:
        return grad(cg1_lift(s))
    return grad(s)


def sample_to_geometry(raster_fn, Q_g, Q_cg, floor=None):
    r"""Sample a raster onto the geometry space ``Q_g`` as a cell average.

    ``raster_fn(space)`` must return the raster interpolated onto ``space``
    (i.e. a closure over ``icepack.interpolate`` and the file).

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

    ``floor`` optionally clamps the field from below (thickness >= h_clamp)
    before averaging.
    """
    field = raster_fn(Q_cg)
    if floor is not None:
        field = Function(Q_cg).interpolate(max_value(field, Constant(floor)))
    if Q_g.ufl_element() == Q_cg.ufl_element():
        return field if isinstance(field, Function) else Function(Q_cg).interpolate(field)
    return Function(Q_g).project(field)
