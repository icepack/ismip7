r"""DG1 geometry, the per-cell rules (``ISMIP7_GEOMETRY_SPACE=dg1``).

Under DG1 the forcing, the SMB-elevation feedback, the backdate and the
ISMIP7 output decide per cell on the cell means of the geometry, with the
surface rebuilt from them by the DG0 formula, so each is the DG0 rule of
those means and a DG0 melt calibration applies as fitted. The apparent mass
balance is the one rule that keeps the slopes: its reference is the DG1
field whose integral against every basis function is the t=0 divergence, so
it holds the t=0 state in its slopes as well as its means.

What these hold: the cell-local mass inverse; the DG1 reference holding a
transport step at the t=0 state, which a cell-constant one does not; the
forcing geometry, the draft slope and the feedback's surface change on a
DG1 run equal to the DG0 ones on its cell means; the melt calibration
record; the backdate on cell means; the inversion's dH/dt term, which
scores the cell means of a DG1 tendency on the grounded cells of the DG0
rule and steps with its forcing on the cells; the ISMIP7 output booking the
DG1 transport's own outflow.

Serial, small rectangle meshes.
"""

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
import icepack2_tools.dual_friction  # noqa: F401,E402

import numpy as np                                              # noqa: E402
import pytest                                                   # noqa: E402
from firedrake import (                                         # noqa: E402
    Constant,
    FacetNormal,
    Function,
    FunctionSpace,
    LinearVariationalProblem,
    LinearVariationalSolver,
    RectangleMesh,
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
    VectorFunctionSpace,
    as_vector,
    assemble,
    cos,
    dot,
    ds,
    dx,
    lhs,
    max_value,
    pi,
    rhs,
    sin,
    solve,
)

from icepack2_tools import geometry as G                        # noqa: E402
from icepack2_tools.cellview import CellView, local_mass_inverse  # noqa: E402
from icepack2_tools.forcing import (                            # noqa: E402
    SMBElevationFeedback,
    _check_melt_provenance,
    compute_sin_alpha,
    flotation_surface,
    forcing_geometry,
    melt_receiving,
)
from icepack2_tools.grounding import height_above_flotation   # noqa: E402
from icepack2_tools.obs_dhdt import (                           # noqa: E402
    CELL_MEAN_SOLVER,
    backdate_cells,
    backdate_thickness,
    cell_mean_form,
    grounded_cells,
)
from icepack2_tools.runconfig import dg1_refusals, melt_geometry_space  # noqa: E402
from icepack2_tools.transport import transport_form, upwind_divergence  # noqa: E402

RHO_RATIO = 917.0 / 1028.0
L = 60e3


def _mesh(n=8, jitter=0.0, seed=0):
    r"""A rectangle mesh, its interior vertices moved by up to ``jitter``
    of a cell so the cell areas differ."""
    mesh = RectangleMesh(n, n, L, L)
    if jitter:
        xy = mesh.coordinates.dat.data
        h = L / n
        inside = (xy[:, 0] > 1.0) & (xy[:, 0] < L - 1.0) & (xy[:, 1] > 1.0) & (xy[:, 1] < L - 1.0)
        rng = np.random.default_rng(seed)
        xy[inside] += rng.uniform(-jitter * h, jitter * h, size=(inside.sum(), 2))
    return mesh


def _ctx(mesh, space, h_expr, b_expr):
    r"""The forward's context as the forcing reads it, for one geometry
    space: ``h``, ``b``, ``s`` on the geometry space, and under DG1 the cell
    view with the forcing on the DG0 cells."""
    Q_g = G.geometry_function_space(mesh, space)
    h = Function(Q_g, name="thickness").interpolate(h_expr)
    b = Function(Q_g, name="bed").interpolate(b_expr)
    s = Function(Q_g).interpolate(
        max_value(b + h, (Constant(1.0) - RHO_RATIO) * h))
    view = CellView(Q_g) if space == "dg1" else None
    Q_f = view.Q0 if view is not None else Q_g
    return {"mesh": mesh, "Q": FunctionSpace(mesh, "CG", 1),
            "V": VectorFunctionSpace(mesh, "CG", 1), "Q_g": Q_g,
            "h": h, "b": b, "s": s, "H_init": h.copy(deepcopy=True),
            "rho_ratio": Constant(RHO_RATIO), "cell_view": view,
            "accum": Function(Q_f), "ocean_melt": Function(Q_f)}


def _dg0_of_means(ctx1):
    r"""The DG0 context whose cells hold the cell means of a DG1 one."""
    view = ctx1["cell_view"]
    mesh = ctx1["mesh"]
    Q0 = FunctionSpace(mesh, "DG", 0)
    h0, b0 = Function(Q0), Function(Q0)
    h0.dat.data[:] = view.means(ctx1["h"])
    b0.dat.data[:] = view.means(ctx1["b"])
    ctx0 = _ctx(mesh, "dg0", 0.0, 0.0)
    ctx0["h"].assign(h0)
    ctx0["b"].assign(b0)
    ctx0["s"].interpolate(max_value(b0 + h0, (Constant(1.0) - RHO_RATIO) * h0))
    return ctx0


def _straddling_geometry(mesh):
    r"""A sheet that thins and goes afloat across x = 30 km on a bed
    sloping below sea level, so some DG1 cells hold both grounded and
    floating vertices."""
    x, y = SpatialCoordinate(mesh)
    h = 1200.0 - 0.015 * x + 40.0 * sin(2 * pi * y / L)
    b = 300.0 - 0.03 * x + 30.0 * cos(2 * pi * y / L)
    return h, b


# ── what a dg1 run does not carry ────────────────────────────────────────

def test_dg1_names_each_setting_it_does_not_carry():
    assert dg1_refusals() == []
    asked = dict(legacy_transport=True, subelement=True,
                 geometry_transfer=True, timing_cache=True)
    lines = dg1_refusals(**asked)
    assert len(lines) == len(asked)
    for name in asked:
        assert len(dg1_refusals(**{name: asked[name]})) == 1


# ── the cell-local mass inverse ──────────────────────────────────────────

@pytest.mark.parametrize("space", ["dg0", "dg1"])
def test_the_mass_inverse_returns_the_field_its_integrals_came_from(space):
    rng = np.random.default_rng(4)
    mesh = _mesh(jitter=0.3)
    Q = G.geometry_function_space(mesh, space)
    view = CellView(Q)
    f = Function(Q)
    f.dat.data[:] = rng.normal(size=f.dat.data.size)
    area = assemble(TestFunction(view.Q0) * dx).dat.data_ro
    assert np.ptp(area) > 0.1 * area.mean()
    out = view.solve_mass(assemble(f * TestFunction(Q) * dx), Function(Q), area)
    np.testing.assert_allclose(out.dat.data_ro, f.dat.data_ro, rtol=1e-10, atol=1e-12)
    # the mean of the coefficients is the cell mean of the integrand
    r = rng.normal(size=(5, 3))
    a = rng.uniform(1.0, 3.0, size=5)
    np.testing.assert_allclose(local_mass_inverse(r, a).mean(axis=1), r.sum(axis=1) / a)


# ── the apparent mass balance ────────────────────────────────────────────

def _flow(mesh):
    x, y = SpatialCoordinate(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    return Function(V).interpolate(
        as_vector((300.0 * (1 + 0.5 * y / L), 60.0 * sin(pi * x / L))))


def _step(h0, u, src, dt=1.0):
    Q = h0.function_space()
    h1 = Function(Q)
    F = transport_form(TrialFunction(Q), h0, TestFunction(Q), u, Constant(dt), src)
    LinearVariationalSolver(LinearVariationalProblem(lhs(F), rhs(F), h1),
                            solver_parameters={"ksp_type": "preonly",
                                               "pc_type": "lu"}).solve()
    return h1


def test_a_dg1_reference_holds_the_t0_state_in_its_slopes():
    r"""The reference the forward builds under DG1 cancels the t=0
    divergence against every basis function, so a transport step from the
    t=0 state with it as the only source returns that state, slopes and
    all. A cell-constant reference, the DG0 rule on the cell means, cancels
    each cell's net flux alone and lets the slopes move."""
    mesh = _mesh(n=12)
    Q1 = G.geometry_function_space(mesh, "dg1")
    view = CellView(Q1)
    x, y = SpatialCoordinate(mesh)
    h0 = Function(Q1).project(800.0 + 300.0 * sin(2 * pi * x / L) * cos(pi * y / L))
    u = _flow(mesh)
    phi = TestFunction(Q1)
    area = assemble(TestFunction(view.Q0) * dx).dat.data_ro
    a_ref = view.solve_mass(assemble(upwind_divergence(h0, phi, u)), Function(Q1), area)
    h1 = _step(h0, u, a_ref)
    np.testing.assert_allclose(h1.dat.data_ro, h0.dat.data_ro, rtol=0, atol=1e-8)

    flat = Function(Q1)
    view.set_cells(flat, np.ones(area.size, dtype=bool), view.means(a_ref))
    h1_flat = _step(h0, u, flat)
    moved = np.abs(view.nodes(h1_flat) - view.nodes(h0)).max()
    assert moved > 1.0


# ── the forcing reads the cell means ─────────────────────────────────────

def test_dg1_forcing_reads_the_dg0_forcing_of_its_cell_means():
    mesh = _mesh(n=10)
    ctx1 = _ctx(mesh, "dg1", *_straddling_geometry(mesh))
    ctx0 = _dg0_of_means(ctx1)
    got, want = forcing_geometry(ctx1), forcing_geometry(ctx0)
    for g, w in zip(got, want):
        np.testing.assert_allclose(g, w, rtol=1e-13, atol=1e-9)
    h, b, s = got
    np.testing.assert_array_equal(melt_receiving(s, b, h),
                                  melt_receiving(want[2], want[1], want[0]))
    # in a cell that straddles flotation the surface is the formula on the
    # means, which the mean of the nodal surface is not
    nodal = ctx1["cell_view"].means(ctx1["s"])
    straddles = np.abs(nodal - s) > 1.0
    assert straddles.any()
    np.testing.assert_allclose(s, flotation_surface(b, h, RHO_RATIO), rtol=1e-14)
    # and the forcing lives on the cells, one dof each
    assert ctx1["accum"].dat.data_ro.size == h.size


@pytest.mark.parametrize("slope", ["ant", "local"])
def test_the_dg1_draft_slope_is_the_dg0_slope_of_its_cell_means(slope, monkeypatch):
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", slope)
    mesh = _mesh(n=10)
    ctx1 = _ctx(mesh, "dg1", *_straddling_geometry(mesh))
    ctx0 = _dg0_of_means(ctx1)
    got, want = compute_sin_alpha(ctx1), compute_sin_alpha(ctx0)
    assert got.shape == want.shape == ctx1["accum"].dat.data_ro.shape
    np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-15)


def test_the_dg1_feedback_measures_the_surface_of_the_cell_means():
    mesh = _mesh(n=10)
    h_expr, b_expr = _straddling_geometry(mesh)
    ctx1 = _ctx(mesh, "dg1", h_expr, b_expr)
    # zero at a cold start, where the thickness is the initial one
    np.testing.assert_array_equal(SMBElevationFeedback.surface_change(ctx1), 0.0)
    x, _ = SpatialCoordinate(mesh)
    ctx1["h"].interpolate(h_expr - 0.002 * x)
    ctx0 = _dg0_of_means(ctx1)
    ctx0["H_init"].dat.data[:] = ctx1["cell_view"].means(ctx1["H_init"])
    got = SMBElevationFeedback.surface_change(ctx1)
    assert np.abs(got).max() > 10.0
    np.testing.assert_allclose(got, SMBElevationFeedback.surface_change(ctx0),
                               rtol=1e-12, atol=1e-9)


def test_a_dg1_run_melts_with_a_dg0_calibration(monkeypatch, tmp_path):
    r"""A calibration records the geometry space it melted; a dg1 run melts
    the dg0 melt of its cell means, so a dg0 file applies and a cg1 file is
    refused, as it is under dg0."""
    npz = str(tmp_path / "deltaT.npz")
    for run, fitted, ok in (("dg1", "dg0", True), ("dg1", "cg1", False),
                            ("dg0", "dg0", True), ("dg0", "dg1", False)):
        monkeypatch.setenv("ISMIP7_GEOMETRY_SPACE", run)
        data = {"melt_slope": "ant", "geometry_space": fitted}
        if ok:
            _check_melt_provenance(npz, data, "calibrate_deltaT.py",
                                   "ISMIP7_DELTAT_PER_BASIN_NPZ", strict=True)
        else:
            with pytest.raises(ValueError, match="geometry"):
                _check_melt_provenance(npz, data, "calibrate_deltaT.py",
                                       "ISMIP7_DELTAT_PER_BASIN_NPZ", strict=True)
    assert [melt_geometry_space(s) for s in ("dg0", "cg1", "dg1")] == ["dg0", "cg1", "dg0"]


# ── the backdate ─────────────────────────────────────────────────────────

def test_a_dg1_backdate_moves_each_cell_by_its_means_change():
    rng = np.random.default_rng(5)
    mesh = _mesh(n=10)
    Q1 = G.geometry_function_space(mesh, "dg1")
    view = CellView(Q1)
    h_expr, b_expr = _straddling_geometry(mesh)
    h = Function(Q1).interpolate(h_expr)
    b = Function(Q1).interpolate(b_expr)
    n = view.Q0.dof_dset.size
    dhdt = rng.normal(0.0, 2.0, size=n)
    observed = (rng.random(n) < 0.8).astype(float)
    means0 = np.array(view.means(h))
    # one grounded cell thickening by more than it holds: undone to zero
    haf = means0 - np.maximum(-np.asarray(view.means(b)), 0.0) / RHO_RATIO
    grounded = np.flatnonzero((haf > 0.0) & (observed > 0.5))
    dhdt[grounded[0]] = 2.0 * means0[grounded[0]] / 12.0
    nodes0 = view.nodes(h)

    new, changed = backdate_cells(view, h, b, dhdt, observed, 12.0, RHO_RATIO)

    want, want_changed = backdate_thickness(means0, view.means(b), dhdt, observed,
                                            12.0, RHO_RATIO)
    np.testing.assert_array_equal(changed, want_changed)
    np.testing.assert_allclose(view.means(h), want, rtol=1e-12, atol=1e-9)
    nodes = view.nodes(h)
    assert nodes.min() >= 0.0
    np.testing.assert_array_equal(nodes[~changed], nodes0[~changed])
    np.testing.assert_array_equal(nodes[grounded[0]], 0.0)
    # elsewhere the slope is kept: each vertex moved by the cell's change
    kept = changed & (nodes0 + (want - means0)[:, None] >= 0.0).all(axis=1)
    np.testing.assert_allclose(nodes[kept] - nodes0[kept],
                               np.broadcast_to((want - means0)[kept, None], nodes[kept].shape),
                               atol=1e-9)
    assert (~changed).any() and kept.any()


# ── the inversion's dH/dt term ───────────────────────────────────────────

def _cell_means_by_solve(f, Q0):
    mean = Function(Q0)
    solve(cell_mean_form(mean, f) == 0, mean, solver_parameters=CELL_MEAN_SOLVER)
    return mean


def test_the_dhdt_term_scores_the_cell_means_of_a_dg1_tendency():
    r"""The DG0 solve the term runs is the cell mean to roundoff, on cells of
    unequal area, and a slope that keeps every mean leaves it unchanged."""
    rng = np.random.default_rng(6)
    mesh = _mesh(n=8, jitter=0.3)
    Q1 = G.geometry_function_space(mesh, "dg1")
    view = CellView(Q1)
    f = Function(Q1)
    f.dat.data[:] = rng.normal(0.0, 3.0, size=f.dat.data.size)
    mean = _cell_means_by_solve(f, view.Q0)
    np.testing.assert_allclose(mean.dat.data_ro, view.means(f), rtol=1e-12, atol=1e-13)
    sloped = f.copy(deepcopy=True)
    tilt = rng.normal(0.0, 5.0, size=(view.Q0.dof_dset.size, 3))
    view.write_nodes(sloped, view.nodes(f) + tilt - tilt.mean(axis=1, keepdims=True))
    assert np.abs(view.nodes(sloped) - view.nodes(f)).max() > 1.0
    np.testing.assert_allclose(_cell_means_by_solve(sloped, view.Q0).dat.data_ro,
                               mean.dat.data_ro, rtol=1e-12, atol=1e-12)


def test_the_dhdt_grounded_cells_are_the_dg0_rule_on_the_cell_means():
    mesh = _mesh(n=10)
    h_expr, b_expr = _straddling_geometry(mesh)
    Q1 = G.geometry_function_space(mesh, "dg1")
    view = CellView(Q1)
    H1 = Function(Q1).interpolate(h_expr)
    b1 = Function(Q1).interpolate(b_expr)
    Q0 = view.Q0
    H0, b0 = Function(Q0), Function(Q0)
    H0.dat.data[:] = view.means(H1)
    b0.dat.data[:] = view.means(b1)
    # under DG0 the rule is the one the term always applied, on the dofs
    haf0 = Function(Q0).interpolate(height_above_flotation(H0, b0)).dat.data_ro
    want = np.where((haf0 > 0.0) & (H0.dat.data_ro > 1.0), 1.0, 0.0)
    np.testing.assert_array_equal(grounded_cells(CellView(Q0), H0, b0).dat.data_ro, want)
    np.testing.assert_array_equal(grounded_cells(view, H1, b1).dat.data_ro, want)
    haf1 = view.nodes(Function(Q1).interpolate(height_above_flotation(H1, b1)))
    assert ((haf1 > 0.0).any(axis=1) & ~(haf1 > 0.0).all(axis=1)).any()
    assert 0 < want.sum() < want.size


def test_a_dg1_dhdt_step_forced_on_the_cells_conserves_mass():
    mesh = _mesh(n=10)
    Q1 = G.geometry_function_space(mesh, "dg1")
    view = CellView(Q1)
    x, y = SpatialCoordinate(mesh)
    h0 = Function(Q1).project(800.0 + 300.0 * sin(2 * pi * x / L) * cos(pi * y / L))
    u = _flow(mesh)
    src = Function(view.Q0).interpolate(0.3 + 0.2 * cos(pi * x / L))
    dt = 2.0
    h1 = _step(h0, u, src, dt=dt)
    un = dot(u, FacetNormal(mesh))
    outflow = assemble((un + abs(un)) / 2 * h1 * ds)
    change = assemble(h1 * dx) - assemble(h0 * dx)
    assert change == pytest.approx(dt * (assemble(src * dx) - outflow), rel=1e-10)


# ── the ISMIP7 output ────────────────────────────────────────────────────

def test_the_output_books_the_dg1_transports_own_outflow(tmp_path):
    from icepack2_tools.ismip7_output import AnnualOutput
    mesh = _mesh(n=8)
    Q1 = G.geometry_function_space(mesh, "dg1")
    view = CellView(Q1)
    x, y = SpatialCoordinate(mesh)
    h = Function(Q1).interpolate(600.0 + 0.004 * x + 100.0 * sin(pi * y / L))
    u = _flow(mesh)
    annual = AnnualOutput(mesh, view.Q0, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)
    annual.start_year(view.mean_function(h))
    annual.begin_step()
    zero = Constant(0.0)
    annual.book_advance(0.5, zero, zero, None, h, u,
                        np.ones(view.Q0.dof_dset.size, dtype=bool))
    annual.commit_step()
    un = dot(u, FacetNormal(mesh))
    outflow = assemble((un + abs(un)) / 2 * h * ds)
    booked = -float((annual.year_acc["licalvf"] * annual.cell_area).sum())
    assert outflow > 0.0
    assert booked == pytest.approx(0.5 * outflow, rel=1e-12)
    annual.close()
