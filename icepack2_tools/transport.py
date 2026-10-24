r"""The forward's thickness transport: implicit upwind on the geometry cells.

One backward-Euler step of ``dh/dt + div(h u) = src`` with the velocity held
fixed, written for a DG0 or DG1 thickness (``ISMIP7_GEOMETRY_SPACE``). Each
facet carries the upwind flux ``max(u.n, 0) h`` from the cell it leaves, the
domain boundary lets ice out and none in, and under DG1 the volume term
``-h u.grad(phi)`` completes the discretization of ``div(h u)``. That term
vanishes identically for a DG0 test function, so under DG0 it is left out
and the form is the finite-volume update the forward has always solved.

Testing with the constant 1, which both spaces hold, leaves only the boundary
outflow: the interior fluxes cancel in pairs, and the volume term integrates
``u.grad(1) = 0``. That is why the forward's transport residual (mass change
minus source minus outflow) is zero to the linear solver's tolerance in
either space.
"""

from firedrake import FacetNormal, dS, ds, dot, dx, grad, inner, jump

__all__ = ["upwind_divergence", "transport_form", "front_transport_form"]


def upwind_divergence(h, phi, u):
    r"""The upwind discretization of ``div(h u)`` tested against ``phi``.

    ``h`` is the thickness (a Function or TrialFunction on the transport
    space), ``phi`` a test function on the same space and ``u`` the velocity.
    """
    n = FacetNormal(phi.function_space().mesh())
    un = dot(u, n)
    un_plus = (un + abs(un)) / 2
    form = ((un_plus("+") * h("+") - un_plus("-") * h("-")) * jump(phi) * dS
            + un_plus * h * phi * ds)
    if phi.ufl_element().degree() > 0:
        form += -h * inner(u, grad(phi)) * dx
    return form


def transport_form(h, h_old, phi, u, dt, src):
    r"""One backward-Euler step from ``h_old`` over ``dt`` with source
    ``src``; solve ``lhs == rhs`` for ``h`` (a TrialFunction)."""
    return ((h - h_old) / dt * phi * dx + upwind_divergence(h, phi, u)
            - src * phi * dx)


def front_transport_form(h, h_old, phi, dt, src, front):
    r"""The same backward-Euler step with the calving front as a moving
    boundary, for a front an advected level set carries
    (``ISMIP7_FRONT_ANCHOR=advect``). ``front`` is an
    :class:`icepack_tools.front_transport.FrontTransport`: only the cells
    inside the front export ice, a cell beyond it receives what follows the
    front's own advance, and the rest of the ice crossing the front is the
    calving flux ``front.calving(h)``, which leaves the cells it crosses from.
    Tested against the constant 1, the form leaves the outflow across the
    domain edge (``front.outflow(h)``) and that calving flux, so the
    transport residual is mass change minus source plus both."""
    return ((h - h_old) / dt * phi * dx + front.divergence(h, phi)
            - src * phi * dx)

