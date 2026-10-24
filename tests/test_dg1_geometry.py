r"""DG1 geometry (``ISMIP7_GEOMETRY_SPACE=dg1``): cell-wise linear thickness,
surface and bed on vertex nodes, discontinuous across cells.

What these hold: the space and its provenance (a ``_dg1`` MAP, refused by
every other space); the cell view the per-cell rules read, which is the
identity under DG0; the positivity floor, the only limiter, which keeps cell
means; BedMachine sampled as the CG1 interpolant with the DG0 cell values as
its means, and front cells rebuilt whole; the upwind transport at second
order and conserving mass; the driving stress at second order with both its
cell and facet terms live; the exact-front push on piecewise-constant DG1
equal to DG0's; and the Budd shelf gate exact at every point of a floating
DG1 cell.

Serial, small rectangle meshes, a synthetic BedMachine-shaped NetCDF.
"""

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import (                      # noqa: E402
    budd_nhat,
    budd_nhat_ungated,
    budd_reference_pressure,
    build_rc_residual,
    effective_pressure,
)

import os                                                       # noqa: E402
import subprocess                                               # noqa: E402
from pathlib import Path                                        # noqa: E402

import netCDF4                                                  # noqa: E402
import numpy as np                                              # noqa: E402
import pytest                                                   # noqa: E402
import rasterio                                                 # noqa: E402
from finat.ufl import FiniteElement                             # noqa: E402
from firedrake import (                                         # noqa: E402
    Constant,
    Function,
    FunctionSpace,
    LinearVariationalProblem,
    LinearVariationalSolver,
    RectangleMesh,
    SpatialCoordinate,
    TensorFunctionSpace,
    TestFunction,
    TrialFunction,
    VectorFunctionSpace,
    as_vector,
    assemble,
    conditional,
    cos,
    dot,
    dx,
    exp,
    grad,
    inner,
    lhs,
    lt,
    max_value,
    pi,
    rhs,
    sin,
)
from icepack2.constants import (                                # noqa: E402
    gravity as g,
    ice_density as RHO_I,
    water_density as RHO_W,
)

from icepack2_tools import geometry as G                        # noqa: E402
from icepack2_tools.cellview import CellView, clamp_cells, floor_fill  # noqa: E402
from icepack2_tools.front import clamp_thickness                # noqa: E402
from icepack2_tools.naming import map_geom_tag                  # noqa: E402
from icepack2_tools.runconfig import (                          # noqa: E402
    geometry_space,
    geometry_space_mismatch,
    warm_start_across_spaces,
)
from icepack2_tools.transport import transport_form, upwind_divergence  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
R = RHO_I / RHO_W


@pytest.fixture
def space_env(monkeypatch):
    def set_space(value):
        monkeypatch.setenv("ISMIP7_GEOMETRY_SPACE", value)
    return set_space


# ── the space and its provenance ─────────────────────────────────────────

def test_dg1_is_a_geometry_space_with_its_own_map_tag(space_env):
    tags = {}
    for space in ("dg0", "cg1", "dg1"):
        space_env(space)
        assert geometry_space() == space
        tags[space] = map_geom_tag()
    assert tags == {"dg0": "_dg0", "cg1": "", "dg1": "_dg1"}
    space_env("dg2")
    with pytest.raises(ValueError, match="'dg1'"):
        geometry_space()


@pytest.mark.parametrize("space, tag", [("dg0", "_dg0"), ("cg1", ""), ("dg1", "_dg1")])
def test_the_runners_name_a_map_as_the_forward_finds_it(space, tag):
    script = REPO / "antarctica" / "scripts" / "ismip7_names.sh"
    out = subprocess.run(
        ["bash", "-c", f"source '{script}'; ismip7_geom_tag"],
        env={**os.environ, "ISMIP7_GEOMETRY_SPACE": space},
        capture_output=True, text=True, check=True).stdout
    assert out == tag


def test_a_dg1_state_seeds_only_a_dg1_run():
    assert geometry_space_mismatch("dg1", "dg1", "map.h5") is None
    for run, recorded in (("dg1", "dg0"), ("dg1", "cg1"), ("dg1", None),
                          ("dg0", "dg1"), ("cg1", "dg1")):
        refusal = geometry_space_mismatch(run, recorded, "map.h5")
        assert refusal and "map.h5" in refusal and run in refusal
    # dg0 and cg1 keep the rules they had
    for run, recorded in (("dg0", "cg1"), ("cg1", "dg0"), ("dg0", None)):
        assert geometry_space_mismatch(run, recorded, "map.h5") is None


def test_an_inversion_takes_another_spaces_controls_alone(monkeypatch):
    r"""A DG0 MAP's controls may start a DG1 inversion and the reverse; its
    geometry may not, and asking for it is refused. Within a space, and
    between DG0 and CG1, nothing changes."""
    monkeypatch.delenv("ISMIP7_WARM_START_GEOMETRY", raising=False)
    assert warm_start_across_spaces("dg1", "dg0", "vf.h5") is True
    assert warm_start_across_spaces("dg0", "dg1", "map.h5") is True
    for run, recorded in (("dg1", "dg1"), ("dg0", "dg0"), ("dg0", "cg1"), ("cg1", None)):
        assert warm_start_across_spaces(run, recorded, "map.h5") is False
    monkeypatch.setenv("ISMIP7_WARM_START_GEOMETRY", "0")
    assert warm_start_across_spaces("dg1", "dg0", "vf.h5") is True
    monkeypatch.setenv("ISMIP7_WARM_START_GEOMETRY", "1")
    with pytest.raises(ValueError, match="ISMIP7_WARM_START_GEOMETRY"):
        warm_start_across_spaces("dg1", "dg0", "vf.h5")
    assert warm_start_across_spaces("dg1", "dg1", "map.h5") is False


def _redistribute(monkeypatch, *argv):
    scripts = str(REPO / "antarctica" / "scripts")
    monkeypatch.syspath_prepend(scripts)
    import redistribute_checkpoint
    monkeypatch.setattr("sys.argv", ["redistribute_checkpoint.py", *argv])
    redistribute_checkpoint.main()


def _dg1_full_state(path):
    from firedrake import CheckpointFile
    mesh = RectangleMesh(4, 4, 4e3, 4e3)
    Q_g = G.geometry_function_space(mesh, "dg1")
    Q = FunctionSpace(mesh, "CG", 1)
    x, y = SpatialCoordinate(mesh)
    H = Function(Q_g).interpolate(500.0 + 0.01 * x)
    s = Function(Q_g).interpolate(400.0 + 0.002 * y)
    with CheckpointFile(str(path), "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(Function(Q, name="log_friction"))
        chk.save_function(Function(Q, name="log_fluidity"))
        chk.save_function(H, name="N_ref_thickness")
        chk.save_function(s, name="N_ref_surface")
        chk.set_attr("/", "friction", "budd")
        chk.set_attr("/", "geometry_space", "dg1")
    return H, s


def test_a_dg1_rewrite_keeps_the_budd_reference_geometry(tmp_path, monkeypatch):
    from firedrake import CheckpointFile
    path = tmp_path / "map_dg1.h5"
    H, s = _dg1_full_state(path)
    _redistribute(monkeypatch, "--input", str(path))
    with CheckpointFile(str(path), "r") as chk:
        mesh = chk.load_mesh()
        for name, want in (("N_ref_thickness", H), ("N_ref_surface", s)):
            got = chk.load_function(mesh, name)
            assert np.allclose(np.sort(got.dat.data_ro),
                               np.sort(want.dat.data_ro))


def test_dg1_refuses_to_publish_a_timing_cache(tmp_path, monkeypatch):
    path = tmp_path / "map_dg1.h5"
    _dg1_full_state(path)
    manifest = tmp_path / "cache.json"
    with pytest.raises(ValueError, match="dg1 does not carry timing caches"):
        _redistribute(monkeypatch, "--input", str(path), "--output",
                      str(tmp_path / "cache.h5"), "--manifest", str(manifest))
    assert not manifest.exists() and not (tmp_path / "cache.h5").exists()


# ── the cell view and the positivity floor ───────────────────────────────

def _mesh(n=6, length=6e3):
    return RectangleMesh(n, n, length, length)


def test_dg1_nodes_sit_on_the_vertices():
    mesh = _mesh()
    Q1 = G.geometry_function_space(mesh, "dg1")
    nodes = Function(VectorFunctionSpace(mesh, Q1.ufl_element())).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro
    verts = mesh.coordinates.dat.data_ro
    gap = np.linalg.norm(nodes[:, None, :] - verts[None, :, :], axis=2).min(axis=1)
    assert gap.max() < 1e-9
    with pytest.raises(ValueError, match="vertices"):
        CellView(FunctionSpace(mesh, "DG", 1))
    with pytest.raises(ValueError, match="continuous"):
        CellView(FunctionSpace(mesh, "CG", 1))


def test_cell_means_are_the_cell_averages():
    mesh = _mesh()
    x, y = SpatialCoordinate(mesh)
    Q1 = G.geometry_function_space(mesh, "dg1")
    f = Function(Q1).project(exp(x / 3e3) * cos(y / 2e3))
    view = CellView(Q1)
    average = Function(view.Q0).project(f).dat.data_ro
    np.testing.assert_allclose(view.means(f), average, rtol=1e-12)


def test_the_cell_view_under_dg0_is_the_array_itself():
    rng = np.random.default_rng(1)
    Q0 = FunctionSpace(_mesh(), "DG", 0)
    view = CellView(Q0)
    f = Function(Q0)
    f.dat.data[:] = rng.normal(size=f.dat.data.size)
    assert view.means(f) is f.dat.data_ro or np.shares_memory(view.means(f), f.dat.data_ro)
    mask = rng.random(f.dat.data.size) < 0.3
    expected = f.dat.data_ro.copy()
    expected[mask] = 7.0
    view.set_cells(f, mask, 7.0)
    np.testing.assert_array_equal(f.dat.data_ro, expected)
    # the positivity floor with one node per cell is the DG0 floor, bitwise
    h = rng.normal(50.0, 80.0, size=200)
    free = rng.random(200) < 0.2
    np.testing.assert_array_equal(
        clamp_cells(h.copy()[:, None], 10.0, free, None)[:, 0],
        clamp_thickness(h.copy(), 10.0, free, None))


def test_the_positivity_floor_keeps_means_and_lifts_vertices():
    rng = np.random.default_rng(2)
    v = rng.normal(20.0, 40.0, size=(500, 3))
    exempt = rng.random(500) < 0.1
    before = v.copy()
    mean0 = before.mean(axis=1)
    floor = 5.0
    floors = np.where(exempt, 0.0, floor)
    out = clamp_cells(v, floor, exempt)
    assert np.all(out >= floors[:, None])
    low = mean0 < floors
    np.testing.assert_allclose(out[low], np.broadcast_to(floors[low, None], out[low].shape))
    np.testing.assert_allclose(out.mean(axis=1)[~low], mean0[~low], rtol=0, atol=1e-12)
    untouched = before.min(axis=1) >= floors
    np.testing.assert_array_equal(out[untouched], before[untouched])
    # a slope scaled exactly to the floor: the lowest vertex sits on it
    scaled = ~low & ~untouched
    assert scaled.any()
    np.testing.assert_allclose(out[scaled].min(axis=1), floors[scaled], atol=1e-12)


def test_the_floor_adds_exactly_the_fill_its_means_call_for():
    rng = np.random.default_rng(3)
    mesh = _mesh(n=12)
    Q1 = G.geometry_function_space(mesh, "dg1")
    view = CellView(Q1)
    h = Function(Q1)
    h.dat.data[:] = rng.normal(20.0, 40.0, size=h.dat.data.size)
    area = assemble(TestFunction(view.Q0) * dx).dat.data_ro
    exempt = rng.random(area.size) < 0.1
    floor = 5.0
    fill = (floor_fill(view.means(h), floor, exempt, None) * area).sum()
    before = assemble(h * dx)
    nodes = view.nodes(h)
    view.write_nodes(h, clamp_cells(nodes.copy(), floor, exempt, None))
    assert fill > 0.0
    np.testing.assert_allclose(assemble(h * dx) - before, fill, rtol=1e-12)
    # a slope scaling that let a scaled cell's mean drop 1 m shows as the gap
    out = clamp_cells(nodes.copy(), floor, exempt, None)
    scaled = (nodes.min(axis=1) < np.where(exempt, 0.0, floor)) \
        & (nodes.mean(axis=1) >= np.where(exempt, 0.0, floor))
    assert scaled.any()
    out[scaled] -= 1.0
    view.write_nodes(h, out)
    gap = assemble(h * dx) - before - fill
    np.testing.assert_allclose(gap, -area[scaled].sum(), rtol=1e-9)


# ── BedMachine onto DG1 ──────────────────────────────────────────────────

PX = 500.0
XS = np.arange(-2000.0, 22000.0 + 1.0, PX)
YS = XS[::-1].copy()
SHELF_H, GROUNDED_H = 200.0, 800.0


def _bedmachine(path):
    r"""Grounded ice on a sloping bed for x < 6 km, a 200 m shelf to the
    front at x = 12 km, ocean beyond (the geometry of test_front_cells.py
    with a bed slope, so the DG1 sample has slopes to keep)."""
    XX, YY = np.meshgrid(XS, YS)
    grounded = XX < 6000.0
    ice = XX < 12000.0
    mask = np.where(grounded, 2, np.where(ice, 3, 0))
    thk = np.where(grounded, GROUNDED_H - 0.02 * XX, np.where(ice, SHELF_H, 0.0))
    bed = np.where(grounded, 100.0 + 0.01 * YY, -500.0)
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", len(XS))
        d.createDimension("y", len(YS))
        for name, vals in (("x", XS), ("y", YS)):
            var = d.createVariable(name, "f8", (name,))
            var[:] = vals
            var.standard_name = f"projection_{name}_coordinate"
            var.units = "m"
        for name, typ, vals in (("mask", "i1", mask), ("thickness", "f4", thk),
                                ("bed", "f4", bed),
                                ("surface", "f4", np.where(grounded, bed + thk, 0.1 * thk))):
            var = d.createVariable(name, typ, ("y", "x"))
            var[:] = vals.astype(typ)
    return str(path)


def _bm_spaces():
    mesh = RectangleMesh(20, 20, 20e3, 20e3)
    return (FunctionSpace(mesh, "CG", 1), FunctionSpace(mesh, "DG", 0),
            G.geometry_function_space(mesh, "dg1"))


def test_a_dg1_vertex_sample_is_the_cg1_interpolant_with_the_dg0_means(tmp_path):
    bm = _bedmachine(tmp_path / "bm.nc")
    Q, Q0, Q1 = _bm_spaces()
    b1, H1, _ = G.sample_bed_thickness(bm, Q1, Q, floor=0.0, method="vertex")
    b0, H0, _ = G.sample_bed_thickness(bm, Q0, Q, floor=0.0, method="vertex")
    H_cg = G.sample_to_geometry(rasterio.open(f"netcdf:{bm}:thickness"), Q, Q, floor=0.0)
    np.testing.assert_allclose(H1.dat.data_ro,
                               Function(Q1).interpolate(H_cg).dat.data_ro, rtol=1e-12)
    view = CellView(Q1)
    np.testing.assert_allclose(view.means(H1), H0.dat.data_ro, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(view.means(b1), b0.dat.data_ro, rtol=1e-9, atol=1e-9)
    with pytest.raises(ValueError, match="cell_mean"):
        G.sample_to_geometry(rasterio.open(f"netcdf:{bm}:bed"), Q1, Q, method="cell_mean")


def test_dg1_front_cells_are_rebuilt_whole_with_the_dg0_values(tmp_path):
    bm = _bedmachine(tmp_path / "bm.nc")
    Q, Q0, Q1 = _bm_spaces()
    b1, H1, counts1 = G.sample_bed_thickness(bm, Q1, Q, floor=0.0, method="vertex_front")
    b0, H0, counts0 = G.sample_bed_thickness(bm, Q0, Q, floor=0.0, method="vertex_front")
    assert counts1 == counts0 and counts1["rebuilt"] > 0
    view = CellView(Q1)
    np.testing.assert_allclose(view.means(H1), H0.dat.data_ro, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(view.means(b1), b0.dat.data_ro, rtol=1e-9, atol=1e-9)
    xc = Function(VectorFunctionSpace(Q0.mesh(), "DG", 0)).interpolate(
        SpatialCoordinate(Q0.mesh())).dat.data_ro[:, 0]
    nodes = view.nodes(H1)
    front = (xc > 11e3) & (xc < 12e3)
    np.testing.assert_allclose(nodes[front], SHELF_H)       # constant, the shelf
    assert np.all(nodes[xc > 12e3] == 0.0)
    grounded = xc < 5e3
    assert np.all(np.ptp(nodes[grounded], axis=1) > 1.0)    # slopes kept inland


# ── transport ────────────────────────────────────────────────────────────

LT = 100e3
U0 = 300.0


def _flow(mesh):
    x, y = SpatialCoordinate(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = Function(V).interpolate(as_vector((U0 * (1 + 0.5 * y / LT), 60.0 * sin(pi * x / LT))))
    h_exact = 500.0 + 300.0 * sin(2 * pi * x / LT) * cos(pi * y / LT)
    return u, h_exact


def _steady_error(space, n):
    r"""L2 error of the steady upwind solution of ``div(h u) = f``, with
    ``f`` manufactured from a smooth thickness that vanishes on the inflow
    sides (x = 0 and y = 0), where the boundary lets no ice in."""
    mesh = RectangleMesh(n, n, LT, LT)
    x, y = SpatialCoordinate(mesh)
    u, _ = _flow(mesh)
    h_exact = (500.0 * sin(pi * x / (2 * LT)) * sin(pi * y / (2 * LT))
               * (1.0 + 0.3 * cos(2 * pi * x / LT)))
    f = h_exact * (u[0].dx(0) + u[1].dx(1)) + dot(u, grad(h_exact))
    Q = G.geometry_function_space(mesh, space)
    h = Function(Q)
    phi = TestFunction(Q)
    F = upwind_divergence(TrialFunction(Q), phi, u) - f * phi * dx(degree=8)
    LinearVariationalSolver(
        LinearVariationalProblem(lhs(F), rhs(F), h),
        solver_parameters={"ksp_type": "preonly", "pc_type": "lu"}).solve()
    return np.sqrt(assemble((h - h_exact) ** 2 * dx(degree=8)))


def test_dg1_transport_is_second_order_where_dg0_is_first():
    rates = {}
    for space in ("dg0", "dg1"):
        e = [_steady_error(space, n) for n in (8, 16, 32)]
        rates[space] = np.log2(e[1] / e[2])
    assert 0.8 < rates["dg0"] < 1.3
    assert rates["dg1"] > 1.7


@pytest.mark.parametrize("space", ["dg0", "dg1"])
def test_a_transport_step_conserves_mass_through_inflow_and_outflow(space):
    mesh = RectangleMesh(12, 12, LT, LT)
    u, h_exact = _flow(mesh)
    Q = G.geometry_function_space(mesh, space)
    h_old = Function(Q).project(h_exact)
    src = Function(Q).project(2.0 + sin(pi * SpatialCoordinate(mesh)[1] / LT))
    h = Function(Q)
    dt = Constant(5.0)
    F = transport_form(TrialFunction(Q), h_old, TestFunction(Q), u, dt, src)
    LinearVariationalSolver(LinearVariationalProblem(
        lhs(F), rhs(F), h),
        solver_parameters={"ksp_type": "preonly", "pc_type": "lu"}).solve()
    n = mesh.coordinates.function_space().mesh()
    from firedrake import FacetNormal, ds
    un = dot(u, FacetNormal(n))
    outflow = assemble((un + abs(un)) / 2 * h * ds)
    change = assemble(h * dx) - assemble(h_old * dx)
    budget = float(dt) * (assemble(src * dx) - outflow)
    assert change == pytest.approx(budget, rel=1e-10)


def test_the_dg0_transport_form_is_the_finite_volume_update():
    r"""Under DG0 the shared form assembles exactly the finite-volume update
    the forward solved before DG1 existed."""
    mesh = RectangleMesh(8, 8, LT, LT)
    u, h_exact = _flow(mesh)
    Q0 = FunctionSpace(mesh, "DG", 0)
    h_old = Function(Q0).project(h_exact)
    src = Function(Q0).assign(1.5)
    dt = Constant(2.0)
    phi, h = TestFunction(Q0), TrialFunction(Q0)
    from firedrake import FacetNormal, dS, ds, jump
    un = dot(u, FacetNormal(mesh))
    up = (un + abs(un)) / 2
    old = ((h - h_old) / dt * phi * dx
           + (up("+") * h("+") - up("-") * h("-")) * jump(phi) * dS
           + up * h * phi * ds - src * phi * dx)
    new = transport_form(h, h_old, phi, u, dt, src)
    A_old, A_new = assemble(lhs(old)), assemble(lhs(new))
    np.testing.assert_array_equal(A_old.M.values, A_new.M.values)
    np.testing.assert_array_equal(assemble(rhs(old)).dat.data_ro,
                                  assemble(rhs(new)).dat.data_ro)


# ── driving stress ───────────────────────────────────────────────────────

def _residual(mesh, H, s, b, exact_front=0, H_cell=None, calving_ids=None):
    r"""The velocity block of the production residual at the zero state:
    the driving stress (and the front push) alone."""
    dg0 = FiniteElement("DG", "triangle", 0)
    Z = (VectorFunctionSpace(mesh, "CG", 1)
         * TensorFunctionSpace(mesh, dg0, symmetry=True)
         * VectorFunctionSpace(mesh, dg0))
    Q = FunctionSpace(mesh, "CG", 1)
    one = Function(Q).assign(1.0)
    F = build_rc_residual(
        Function(Z), Function(Q), Function(Q), H=H, s=s, b=b, C_w0=one, A4_base=one,
        n_flow=Constant(3.0), n_flow_val=3.0, m_slide=3.0, tau_c=0.1, alpha=1e-4,
        H_ref=100.0, k_lim=0.0, calving_ids=calving_ids, exact_front=exact_front,
        front_hmin=1.0, H_cell=H_cell)
    return assemble(F).subfunctions[0].dat.data_ro.copy()


def _driving_stress_error(space, n):
    r"""Relative error of the driving force on a smooth test field ``w``:
    the production residual's velocity block applied to the CG1 interpolant
    of ``w``, against ``-rho_I g int H grad(s) . w`` for the exact fields.
    H and s are L2-projected, so a DG1 surface jumps across every facet."""
    mesh = RectangleMesh(n, n, LT, LT)
    x, y = SpatialCoordinate(mesh)
    H_ex = 1000.0 + 300.0 * cos(pi * x / LT) * cos(pi * y / LT)
    s_ex = 1500.0 + 400.0 * sin(2 * pi * x / LT) * sin(pi * y / LT)
    Q = G.geometry_function_space(mesh, space)
    H = Function(Q).project(H_ex)
    s = Function(Q).project(s_ex)
    b = Function(Q).assign(-2000.0)
    r = _residual(mesh, H, s, b)
    w_ex = as_vector((sin(pi * x / LT) * sin(pi * y / LT), x * y / LT ** 2))
    w = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(w_ex)
    exact = float(assemble(-RHO_I * g * H_ex * inner(grad(s_ex), w_ex) * dx(degree=10)))
    return abs((r * w.dat.data_ro).sum() - exact) / abs(exact)


def test_the_dg1_driving_stress_is_second_order():
    r"""Both terms of the driving stress live and second order together, as
    CG1's cell term alone is; DG0's facet term alone is first order
    (GEOMETRY_DISCRETIZATION.md). Measured: 64 x 64 gives rates 1.99, 2.00
    and 0.88 for DG1, CG1 and DG0, with DG0 still approaching 1."""
    rates = {}
    for space in ("dg0", "dg1", "cg1"):
        e = [_driving_stress_error(space, n) for n in (16, 32)]
        rates[space] = np.log2(e[0] / e[1])
    assert rates["dg1"] > 1.8 and rates["cg1"] > 1.8
    assert rates["dg0"] < 1.2


@pytest.mark.parametrize("version", [1, 2])
@pytest.mark.parametrize("b_water", [-1000.0, 300.0])
def test_exact_front_on_piecewise_constant_dg1_is_dg0s(version, b_water):
    r"""A grounded slab beside water and beside rock, with a film on the water
    side: piecewise-constant DG1 fields push exactly as the DG0 cells do,
    with the cell means deciding which cells hold ice."""
    L, B, W, NX = 20e3, 10e3, 4e3, 10
    mesh = RectangleMesh(int(NX * (L + B) / L), 2, L + B, W)
    x, _ = SpatialCoordinate(mesh)
    Q0 = FunctionSpace(mesh, "DG", 0)
    H0 = Function(Q0).interpolate(conditional(lt(x, L), 1600.0, 0.4))
    b0 = Function(Q0).interpolate(conditional(lt(x, L), -300.0, b_water))
    s0 = Function(Q0).interpolate(max_value(b0 + H0, (1 - R) * H0))
    Q1 = G.geometry_function_space(mesh, "dg1")
    H1, b1, s1 = (Function(Q1).interpolate(f) for f in (H0, b0, s0))
    r0 = _residual(mesh, H0, s0, b0, exact_front=version)
    view = CellView(Q1)
    r1 = _residual(mesh, H1, s1, b1, exact_front=version, H_cell=view.mean_function(H1))
    np.testing.assert_allclose(r1, r0, rtol=1e-10, atol=1e-10 * np.abs(r0).max())
    with pytest.raises(ValueError, match="H_cell"):
        _residual(mesh, H1, s1, b1, exact_front=version)


# ── the Budd shelf gate on DG1 ───────────────────────────────────────────

def test_the_budd_gate_zeroes_every_point_of_a_floating_dg1_cell():
    r"""A shelf thinning across its cells, with the flotation surface at the
    nodes and a surface offset far above roundoff: ``N`` inside each cell is
    the roundoff residue the old sign test passed, so the ungated law reaches
    the cap there, and the HAF gate is exactly zero at every point. Grounded
    cells keep their friction, and a cell the grounding line crosses is gated
    point by point."""
    mesh = RectangleMesh(6, 2, 3e3, 1e3)
    x, _ = SpatialCoordinate(mesh)
    Q1 = G.geometry_function_space(mesh, "dg1")
    H = Function(Q1).interpolate(900.0 - 0.2 * x)
    # floating for x < 2 km (bed at -2000 m), grounded beyond (bed at -300 m),
    # constant on each cell, so the cells on x = 2 km hold one bed each
    b0 = Function(FunctionSpace(mesh, "DG", 0)).interpolate(
        conditional(lt(x, 2e3), -2000.0, -300.0))
    b = Function(Q1).interpolate(b0)
    s = Function(Q1).interpolate(max_value(b + H, (1 - R) * H + 1e-9))
    probe = FunctionSpace(mesh, "DG", 3)
    xs = Function(probe).interpolate(x).dat.data_ro
    N = effective_pressure(H, s)
    gated = Function(probe).interpolate(budd_nhat(N, None, H, b, nhat_floor=0.02, nhat_cap=3.0))
    ungated = Function(probe).interpolate(
        budd_nhat_ungated(N, None, H, nhat_floor=0.02, nhat_cap=3.0))
    shelf = xs < 2e3 - 1.0
    assert np.all(gated.dat.data_ro[shelf] == 0.0)
    assert np.allclose(ungated.dat.data_ro[shelf], 3.0)
    # on the -300 m bed the ice floats from 900 - 0.2 x = (rho_W/rho_I) 300,
    # inside the last column of cells
    x_gl = (900.0 - 300.0 / R) / 0.2
    assert 2.5e3 < x_gl < 3e3
    assert np.all(gated.dat.data_ro[(xs > 2e3 + 1.0) & (xs < x_gl - 1.0)] > 0.0)
    assert np.all(gated.dat.data_ro[xs > x_gl + 1.0] == 0.0)


def test_the_dg1_budd_reference_starts_n_hat_where_the_inversion_does():
    r"""Grounded DG1 cells whose bed crosses sea level, beside a grounding
    line: ``N`` has kinks inside these cells, so a nodal interpolant of it as
    ``N_ref`` moves ``N_hat`` off the inversion's ``N / N`` inside them, and
    the reference evaluated point by point on the same fields does not."""
    mesh = RectangleMesh(8, 2, 8e3, 2e3)
    x, _ = SpatialCoordinate(mesh)
    Q1 = G.geometry_function_space(mesh, "dg1")
    b = Function(Q1).interpolate(150.0 - 0.07 * x)          # sea level at 2.1 km
    H = Function(Q1).interpolate(700.0 - 0.06 * x)
    s = Function(Q1).interpolate(max_value(b + H, (1 - R) * H))
    N = effective_pressure(H, s)
    probe = FunctionSpace(mesh, "DG", 3)

    def nhat(N_ref):
        return Function(probe).interpolate(
            budd_nhat(N, N_ref, H, b, nhat_floor=0.02, nhat_cap=3.0)).dat.data_ro

    inversion = nhat(None)
    pointwise = nhat(budd_reference_pressure(H, s))
    nodal = nhat(Function(Q1).interpolate(budd_reference_pressure(H, s)))
    # the two N expressions are built separately, so the form compiler may
    # order their arithmetic differently: equal to roundoff, not bitwise
    np.testing.assert_allclose(pointwise, inversion, rtol=1e-12, atol=0.0)
    assert np.abs(nodal - inversion).max() > 0.01
