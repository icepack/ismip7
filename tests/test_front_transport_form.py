r"""The forward's transport under an advected front (ISMIP7_FRONT_ANCHOR=advect).

``icepack2_tools.transport.front_transport_form`` is the forward's
backward-Euler step with the shared
:class:`icepack_tools.front_transport.FrontTransport` as its divergence: only
the cells inside the front export ice, and the ice that crosses the front
without following it is the calving flux. The forward's mass check reads the
step's mass change as the source minus the outflow across the domain edge
minus that calving flux, on a DG0 and on a DG1 thickness, so that identity is
what is tested here, for a free, a held and a retreating front. The ISMIP7
output books the outflow of the exporting cells only.

The shared class's own cases (a free advance keeps its ice, a barrier holds,
enclosed water joins the ice, icebergs go, serially and on three ranks) are
``icepack_tools/test/front_transport_test.py``.
"""
import numpy as np
import pytest

fd = pytest.importorskip("firedrake")
pytest.importorskip("icepack_tools.front_transport")

from firedrake import (Constant, Function, FunctionSpace, TestFunction,  # noqa: E402
                       TrialFunction, VectorFunctionSpace, assemble, dx)

from icepack_tools.front_transport import FrontTransport  # noqa: E402
from icepack_tools.levelset import LevelSet  # noqa: E402

from icepack2_tools.transport import front_transport_form  # noqa: E402

LENGTH, WIDTH, NX, NY = 100e3, 20e3, 20, 4
FRONT0, H0, U, SMB, DT = 50e3, 300.0, 500.0, 0.3, 0.5


def _strip(degree, **level_set):
    mesh = fd.RectangleMesh(NX, NY, LENGTH, WIDTH)
    x = fd.SpatialCoordinate(mesh)
    Q0 = FunctionSpace(mesh, "DG", 0)
    Q = Q0 if degree == 0 else FunctionSpace(
        mesh, fd.FiniteElement("DG", mesh.ufl_cell(), 1, variant="equispaced"))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        fd.as_vector((Constant(U), Constant(0.0))))
    h0 = Function(Q0).interpolate(
        fd.conditional(x[0] < FRONT0, Constant(H0), Constant(0.0)))
    h = Function(Q).interpolate(h0)
    ls = LevelSet(mesh, h0, law="prescribed", h_min=1.0, anchor="advect", **level_set)
    ls.phi.interpolate(x[0] - Constant(FRONT0))
    ls.update_cell_fields()
    return mesh, Q, u, h, ls


@pytest.mark.parametrize("degree", [0, 1])
@pytest.mark.parametrize("rate, calves", [(0.0, False), (U, True), (U + 300.0, True)])
def test_one_step_closes_the_mass_budget(degree, rate, calves):
    mesh, Q, u, h, ls = _strip(degree)
    front = FrontTransport(ls, u)
    h_old = Function(Q).assign(h)
    src = Function(Q).assign(SMB)
    dt = Constant(DT)
    F = front_transport_form(TrialFunction(Q), h_old, TestFunction(Q), dt, src, front)
    front.advance(DT, Constant(rate))
    fd.solve(fd.lhs(F) == fd.rhs(F), h,
             solver_parameters={"ksp_type": "preonly", "pc_type": "lu"})
    m0, m1 = float(assemble(h_old * dx)), float(assemble(h * dx))
    source = float(assemble(src * dx)) * DT
    out = float(assemble(front.outflow(h))) * DT
    calving = float(assemble(front.calving(h))) * DT
    assert abs((m1 - m0) - (source - out - calving)) < 1e-10 * m0
    reach = H0 * U * WIDTH * DT
    if calves:
        # the ice reaching the front this step, within the first-order
        # scheme's spread
        assert calving == pytest.approx(reach, rel=0.1)
    else:
        assert abs(calving) < 1e-12 * reach
    # each exporting cell's share sums to the total
    cells = assemble(front.calving(h, TestFunction(FunctionSpace(mesh, "DG", 0))))
    assert float(cells.dat.data_ro.sum()) * DT == pytest.approx(calving, rel=1e-12, abs=1e-6)


@pytest.mark.parametrize("degree", [0, 1])
def test_a_retreating_front_is_not_a_substep_error(degree):
    """The substep controller reads the transport's own increment and leaves
    out the cells the front crossed. Read after the front's rules, the cut of
    a cell the front retreats off (it keeps the part still behind the front)
    follows a gain from the accumulation: a reversal that grows, by more than
    the controller accepts at its ceiling (20 m), at any step. That is how
    the 2 km buffered calving screen of 16 October stalled within 0.15 yr."""
    from icepack2_tools.cellview import CellView
    from icepack2_tools.substep import SubstepController
    # reinitialised on model time, as the forward does (ISMIP7_FRONT_REINIT_YR)
    mesh, Q, u, h, ls = _strip(degree, reinit_interval=2.5)
    front = FrontTransport(ls, u)
    view = CellView(Q)
    h_old = Function(Q).assign(h)
    src = Function(Q)
    tau = 0.05
    dt = Constant(tau)
    F = front_transport_form(TrialFunction(Q), h_old, TestFunction(Q), dt, src, front)
    solver = fd.LinearVariationalSolver(
        fd.LinearVariationalProblem(fd.lhs(F), fd.rhs(F), h),
        solver_parameters={"ksp_type": "preonly", "pc_type": "lu"})
    rate = Constant(U + 2500.0)                      # retreat at 2.5 km/yr
    after_rules = SubstepController(tol=1.0)
    transport = SubstepController(tol=1.0)
    worst = {"after_rules": 0.0, "transport": 0.0}
    cut = 0
    for _ in range(60):                              # 3 yr: the front crosses a cell
        before = np.array(h.dat.data_ro)
        front.advance(tau, rate)
        crossed = front.inside_old != (front.chi_in.dat.data_ro > 0.5)
        # accumulation on the ice only, as the forward forces it
        src.interpolate(SMB * front.chi_in)
        h_old.assign(h)
        solver.solve()
        h_transport = np.array(h.dat.data_ro)
        # the forward's rules (simulation._front_rules), on the whole cell
        keep, band, touch = front.rules()
        cut += int(np.count_nonzero((keep > 0.0) & (keep < 1.0)))
        view.write_nodes(h, view.nodes(h) * keep[:, None])
        give, back, _ = front.excess(np.array(view.means(h)), band, touch)
        view.write_nodes(h, view.nodes(h) * (1.0 - give)[:, None])
        view.add_cells(h, back)
        worst["after_rules"] = max(worst["after_rules"], after_rules.observe(
            before, h.dat.data_ro, tau))
        worst["transport"] = max(worst["transport"], transport.observe(
            before, h_transport, tau, skip=view.node_mask(crossed)))
    assert cut > 0, "the front never retreated off a cell"
    assert worst["after_rules"] > 20.0
    assert worst["transport"] < 0.1


def test_the_output_books_the_outflow_of_the_exporting_cells_only(tmp_path):
    from icepack2_tools.ismip7_output import AnnualOutput
    mesh = fd.RectangleMesh(4, 1, 4.0, 1.0)
    Q = FunctionSpace(mesh, "DG", 0)
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        fd.as_vector((Constant(1.0), Constant(0.0))))
    h = Function(Q).assign(1.0)
    xc = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    inside = Function(Q)
    inside.dat.data[:] = np.where(xc[:, 0] < 2.0, 1.0, 0.0)   # the outflow edge is beyond the front
    grounded = np.ones(Q.dof_dset.size, dtype=bool)
    zero = Constant(0.0)
    booked = {}
    for label, exporting in (("all", None), ("inside", inside)):
        annual = AnnualOutput(mesh, Q, str(tmp_path / label / "annual.h5"),
                              str(tmp_path / label / "scalars.csv"),
                              first_year=2015, rho_ratio=917.0 / 1024.0)
        annual.start_year(h)
        annual.begin_step()
        annual.book_advance(1.0, zero, zero, None, h, u, grounded, exporting=exporting)
        annual.commit_step()
        booked[label] = -float((annual.year_acc["licalvf"] * annual.cell_area).sum())
        annual.close()
    # every cell exports: the unit thickness leaves across the unit edge at x = 4
    assert booked["all"] == pytest.approx(1.0, rel=1e-12)
    # the cells at the edge are beyond the front and export nothing
    assert booked["inside"] == 0.0
