"""
Solve the eikonal equation |grad(d)| = 1 to compute distance from the
calving front (ice edge) into the ice sheet.

Uses the Latent Variable Proximal Point (LVPP) method:
  Dokken, Farrell, Keith, Papadopoulos, Surowiec (2025)
  "The latent variable proximal point algorithm for variational
   problems with inequality constraints"
  Computer Methods in Applied Mechanics and Engineering.
  https://doi.org/10.1016/j.cma.2025.118181
  Section 5.1: eikonal equation.

The calving front is identified as the boundary of the region where
BedMachine thickness is zero (the buffer zone between the ice edge
and the mesh boundary).  d = 0 is enforced in this buffer via a
penalty term, connecting seamlessly to the DirichletBC at the mesh
boundary.  The distance field then ramps up from the ice edge.

Usage:
    from eikonal_distance import identify_calving_front, solve_eikonal_distance
    front_mask = identify_calving_front(H)
    d = solve_eikonal_distance(mesh, front_mask)
"""

import firedrake as fd
from firedrake import (
    Function,
    FunctionSpace,
    MixedFunctionSpace,
    Constant,
    TestFunction,
    split,
    inner,
    dot,
    sqrt,
    div,
    dx,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    DirichletBC,
    norm,
    conditional,
)
from firedrake.petsc import PETSc

from icepack2_tools.mpi_stats import (
    global_absmax, global_count, global_max, global_size,
)


def identify_calving_front(H, h_threshold=10.0):
    """Identify the calving front from ice thickness.

    Returns a mask that is 1.0 where thickness is at or below the
    clamp value (buffer zone / ocean) and 0.0 in the ice interior.
    The boundary of this mask is the calving front.

    Parameters
    ----------
    H : firedrake.Function
        Ice thickness (may be clamped to a minimum value).
    h_threshold : float
        Thickness below which nodes are considered buffer/ocean.
        Should match the clamp value used when interpolating H
        (default 10 m).

    Returns
    -------
    front_mask : firedrake.Function
        CG1 function: 1.0 in the buffer zone, 0.0 in the ice.
    """
    Q = H.function_space()
    front_mask = Function(Q, name="calving_front_mask")
    front_mask.interpolate(
        conditional(H < Constant(h_threshold + 0.1), Constant(1.0), Constant(0.0))
    )
    n_front = global_count(front_mask.dat.data_ro > 0.5, Q.mesh().comm)
    n_total = global_size(front_mask)
    PETSc.Sys.Print(
        f"  Calving front mask: {n_front}/{n_total} nodes in buffer (H < {h_threshold} m)"
    )
    return front_mask


def identify_grounding_line(H, b, rho_I=917.0, rho_W=1024.0, kH=1.0):
    """Identify the grounding line from thickness and bed.

    Returns a mask that is 1.0 at the grounding line (where height
    above flotation ≈ 0) and 0.0 away from it. Uses the smooth
    Dirac delta of height-above-flotation (Ua approach).

    Parameters
    ----------
    H : firedrake.Function
        Ice thickness.
    b : firedrake.Function
        Bed elevation.
    rho_I, rho_W : float
        Ice and water densities.
    kH : float
        Sharpness parameter (m⁻¹). 1/kH is the vertical width
        of the grounding zone.

    Returns
    -------
    gl_mask : firedrake.Function
        CG1 function: ~1 at the GL, ~0 away from it.
    """
    Q = H.function_space()
    h_f = (
        Constant(rho_W / rho_I) * sqrt(Constant(0.0) - b) * sqrt(Constant(0.0) - b)
    )  # avoid max_value in UFL
    # Simpler: height above flotation
    haf = H - Constant(rho_W / rho_I) * Function(Q).interpolate(
        conditional(-b > Constant(0.0), -b, Constant(0.0))
    )
    # Dirac delta: peaks at haf=0, width ~1/kH
    # delta(haf) = 0.5 * kH * sech^2(kH * haf)
    gl_mask = Function(Q, name="grounding_line_mask")
    kH_c = Constant(kH)
    gl_mask.interpolate(
        Constant(0.5) * kH_c / (fd.cosh(kH_c * haf) * fd.cosh(kH_c * haf))
    )
    # Normalize to [0, 1] by the GLOBAL peak, so the mask is uniformly scaled
    # across ranks and a rank whose slice is all zero does not skip it.
    gl_max = global_max(gl_mask)
    if gl_max > 0:
        gl_mask.dat.data[:] /= gl_max

    n_gl = global_count(gl_mask.dat.data_ro > 0.1, Q.mesh().comm)
    PETSc.Sys.Print(f"  GL mask: {n_gl} nodes near grounding line (kH={kH})")
    return gl_mask


def solve_eikonal_distance(mesh, front_mask, max_its=40, tol=1e-4):
    """Solve |grad(d)| = 1 with d = 0 at the calving front.

    The front_mask marks the buffer zone (H ≈ 0) where d is forced
    to zero via a penalty term.  Since the buffer extends to the mesh
    boundary, this connects smoothly to the DirichletBC there.  The
    distance field ramps up from the ice edge at rate 1.

    Coordinates are temporarily non-dimensionalized (divided by a
    characteristic length L) so the LVPP iteration is well-conditioned.
    The result is returned in physical units (meters).

    Parameters
    ----------
    mesh : firedrake.Mesh
        The computational mesh.
    front_mask : firedrake.Function
        CG1 function with 1.0 in the buffer zone, 0.0 in the ice.
        Returned by ``identify_calving_front``.
    max_its : int
        Maximum LVPP iterations.
    tol : float
        Convergence tolerance on ||u^k - u^{k-1}||_L2.

    Returns
    -------
    firedrake.Function
        Distance field d(x) on CG1 space, in meters.
    """
    # Non-dimensionalize: scale coordinates to O(1). The extent is reduced
    # across ranks so every partition is scaled by the SAME factor; a
    # rank-local extent deforms the mesh inconsistently and an element
    # straddling a partition boundary is then assembled from two different
    # coordinate scalings.
    coords = mesh.coordinates
    coords_orig = coords.dat.data.copy()
    L = global_absmax(coords, mesh.comm)
    if not L > 0.0:
        raise ValueError(
            "degenerate mesh: all coordinates are zero, cannot "
            "non-dimensionalize"
        )
    PETSc.Sys.Print(f"  Scaling coordinates by L = {L/1e3:.0f} km")
    coords.dat.data[:] /= L

    try:
        d = _solve_lvpp(mesh, front_mask, max_its, tol)
        # Scale distance back to meters
        d.dat.data[:] *= L
    finally:
        # Always restore original coordinates
        coords.dat.data[:] = coords_orig

    PETSc.Sys.Print(f"  Distance field: max={global_max(d)/1e3:.1f} km")
    return d


def _solve_lvpp(mesh, front_mask, max_its, tol):
    """LVPP solve on the (already non-dimensionalized) mesh."""
    degree = 0
    Q_rt = FunctionSpace(mesh, "RT", degree + 1)
    V_dg = FunctionSpace(mesh, "DG", degree)

    Z = MixedFunctionSpace([V_dg, Q_rt])
    z = Function(Z)
    u, psi = split(z)

    u_old = Function(V_dg)
    psi_old = Function(Q_rt)

    z_test = TestFunction(Z)
    v, w = split(z_test)

    alpha = Constant(1.0)

    # Project front mask to DG0 and threshold
    front_ind = Function(V_dg)
    front_ind.project(front_mask)
    front_ind.interpolate(
        conditional(front_ind > Constant(0.1), Constant(1.0), Constant(0.0))
    )
    n_front = int((front_ind.dat.data_ro > 0.5).sum())
    PETSc.Sys.Print(f"  LVPP: {n_front} buffer cells for penalty")

    # LVPP variational form — equations (5.5a) and (5.5b)
    F = inner(psi / sqrt(1 + dot(psi, psi)), w) * dx(degree=4)
    F += inner(u, div(w)) * dx
    F += inner(div(psi), v) * dx
    F -= inner(div(psi_old), v) * dx
    F += alpha * inner(Constant(1.0), v) * dx

    # Penalty: d ≈ 0 in the buffer zone.
    # Only needs to dominate the source alpha (≤ 10).
    penalty = Constant(100.0)
    F += penalty * front_ind * inner(u, v) * dx

    bcs = DirichletBC(Z.sub(0), 0, "on_boundary")

    sp = {
        "snes_type": "newtonls",
        "snes_monitor": None,
        "snes_max_it": 100,
        "snes_linesearch_type": "l2",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 200,
        "snes_rtol": tol,
        "snes_atol": tol,
        "snes_divergence_tolerance": 1e10,
    }

    NLVP = NonlinearVariationalProblem(F, z, bcs=bcs)
    NLVS = NonlinearVariationalSolver(NLVP, solver_parameters=sp)

    u_fn, psi_fn = z.subfunctions

    for k in range(max_its):
        u_old.assign(u_fn)
        NLVS.solve()
        psi_old.assign(psi_fn)

        nrm = norm(u_fn - u_old, "L2")
        PETSc.Sys.Print(f"  LVPP it {k}: error = {nrm:.2e}, alpha = {float(alpha):.1f}")
        if nrm < tol:
            break
        if float(alpha) < 10:
            alpha.assign(sqrt(2) * alpha)

    # Project to CG1
    Q_cg = FunctionSpace(mesh, "CG", 1)
    d = Function(Q_cg, name="eikonal_distance")
    d.project(u_fn)
    d.interpolate(fd.max_value(d, Constant(0.0)))
    return d
