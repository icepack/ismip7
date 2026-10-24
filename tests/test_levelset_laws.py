r"""What each calving law does to the front, and what one step removes.

The level-set unit tests were lost (issue #35). These rebuild the part that
belongs here: the laws as ISMIP7 configures them (``ISMIP7_CALVING`` naming a
law in ``icepack_tools.calving``, ``ISMIP7_CALVING_PARAMS`` its parameters,
``ISMIP7_CALVING_MODULE`` a law of our own) driving ``icepack2_tools.levelset``
(the ``extent`` anchoring default), the masks the transport applies
afterwards, the drag gate the momentum residual reads, and the step-size
behaviour of the retreat. The laws themselves and the shared class's
operators are tested in ``icepack_tools``; what is checked here is the
behaviour this repository depends on.

Serial, one small unit mesh, no data files and no momentum solve: a velocity
and a thickness are imposed, so every number below follows from the law and
the geometry alone.
"""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from firedrake import (Constant, Function, FunctionSpace, SpatialCoordinate,
                       TensorFunctionSpace, UnitSquareMesh, VectorFunctionSpace,
                       as_vector)
from icepack2.constants import ice_density, water_density
from icepack_tools.calving import FrontState

from icepack2_tools import runconfig
from icepack2_tools.levelset import LevelSet, initial_distance

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Ice on the left half of a unit domain, open water beyond, so the extent has
#: an interior front. Cells are 1/8 wide; on the default diagonal every front
#: cell is a triangle with one 1/8 facet on x = 0.5 and area 1/128.
N = 8
H_ICE = 300.0
HMIN = 1.0
FRONT_LEN = 1.0 / N
CELL_AREA = 1.0 / (2 * N * N)
#: An extending flow u = (STRAIN x, 0), so e_xx = STRAIN and the von Mises
#: stress is set by the strain rate rather than by a regulariser.
STRAIN = 0.1
#: A bed deep enough that every ice cell floats, and one shallow enough that
#: every ice cell is grounded.
BED_FLOATING = -2000.0
BED_GROUNDED = -100.0
SIGMA_MAX_GROUNDED = 1.0
SIGMA_MAX_FLOATING = 0.15


def _state(law, bed=BED_FLOATING, **kwargs):
    mesh = UnitSquareMesh(N, N)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    h = Function(Q0).interpolate(H_ICE * (x < 0.5))
    b = Function(Q0).interpolate(Constant(bed))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((STRAIN * x, Constant(0.0))))
    ls = LevelSet(mesh, h, law=law, h_min=HMIN, **kwargs)
    return mesh, Q0, h, b, u, ls


@pytest.fixture
def calving_env(monkeypatch):
    r"""A clean calving configuration, set through the environment the way a
    run is configured."""
    for knob in ("ISMIP7_CALVING", "ISMIP7_CALVING_PARAMS", "ISMIP7_CALVING_MODULE",
                 *runconfig.RETIRED_CALVING_KNOBS):
        monkeypatch.delenv(knob, raising=False)
    return monkeypatch


def _configured_law(env, name, params=""):
    env.setenv("ISMIP7_CALVING", name)
    if params:
        env.setenv("ISMIP7_CALVING_PARAMS", params)
    return runconfig.calving_law_object()


def _front_state(mesh, h, b, u, ls, A=None, n=None):
    M = Function(TensorFunctionSpace(mesh, "DG", 0, symmetry=True))
    return FrontState(u, M, None, h, b, ls, A=A, n=n,
                      rho_i=ice_density, rho_w=water_density)


def _xc(mesh):
    return Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro[:, 0]


def test_the_extent_anchor_is_the_default():
    r"""ISMIP7 anchors on the ice extent rather than on a fixed contour,
    because on a sheet with melt and a buffered ocean the thickness is the
    authority on where ice is."""
    _, _, _, _, _, ls = _state("none")
    assert ls.anchor == "extent"


def test_a_front_cell_is_one_with_an_ice_free_neighbour():
    mesh, _, h, _, _, ls = _state("none")
    front = ls.front_len.dat.data_ro > 0.0
    xc = _xc(mesh)
    # the front runs at x = 0.5, so every front cell sits in the column
    # immediately inside it
    assert front.any()
    assert np.all(xc[front] < 0.5)
    assert np.all(xc[front] > 0.5 - 2.0 / N)


def test_fixed_holds_the_front_and_removes_nothing_from_it():
    r"""`fixed` is the control configuration: the front cannot move, so the
    rate is zero and only the cells beyond the t=0 extent are flagged."""
    mesh, Q0, h, b, u, ls = _state("fixed")
    rate = ls.advance(0.1, u)
    assert rate == 0.0
    beyond, frac = ls.calving_masks()
    assert frac is None, "a held front sheds no sub-cell mass"
    xc = _xc(mesh)
    assert np.all(beyond == (xc > 0.5))


@pytest.mark.parametrize("bed, sigma_max", [
    (BED_FLOATING, SIGMA_MAX_FLOATING),
    (BED_GROUNDED, SIGMA_MAX_GROUNDED),
])
def test_the_configured_law_drives_the_front(calving_env, bed, sigma_max):
    r"""``ISMIP7_CALVING=vonmises_strain`` with its thresholds in
    ``ISMIP7_CALVING_PARAMS``, evaluated on the forward's kind of front state
    and handed to the level set as its rate:
    ``c = |u| sqrt(3) A^(-1/n) e~^(1/n) / sigma_max``. Under uniaxial
    extension ``e~ = e_xx / sqrt(2)``, and the threshold is the floating one
    over the deep bed and the grounded one over the shallow bed."""
    law = _configured_law(calving_env, "vonmises_strain",
                          f"sigma_max_gr={SIGMA_MAX_GROUNDED}, sigma_max_fl={SIGMA_MAX_FLOATING}")
    assert law.front_mode == "prescribed"
    mesh, Q0, h, b, u, ls = _state(law.front_mode, bed=bed)
    A, n, dt = 20.0, 3.0, 0.1
    state = _front_state(mesh, h, b, u, ls, A=Constant(A), n=n)
    rate = ls.advance(dt, u, rate=law.rate(state, 0.0))
    front = ls.front_len.dat.data_ro > 0.0
    xc = _xc(mesh)
    sigma = np.sqrt(3.0) * A ** (-1.0 / n) * (STRAIN / np.sqrt(2.0)) ** (1.0 / n)
    expected = STRAIN * xc[front] * sigma / sigma_max
    assert np.allclose(ls.c_cell.dat.data_ro[front], expected, rtol=1e-6)
    assert rate == pytest.approx(expected.mean(), rel=1e-6)
    beyond, frac = ls.calving_masks()
    assert np.allclose(frac[front], expected * dt * FRONT_LEN / CELL_AREA,
                       rtol=1e-6)
    assert np.all(frac[~front] == 0.0)
    assert not beyond[h.dat.data_ro > HMIN].any(), \
        "one small step must not carry the front past a whole ice cell"


def test_the_restart_check_rebuilds_the_gate_a_level_set_front_writes():
    r"""A level-set front writes its own ocean-drag gate (drag only more than
    a cell diameter beyond its front) and a checkpoint records its residual
    under that gate. ``front.levelset_drag_gate`` is what a restart rebuilds
    from the saved level set before it measures the loaded state: it has to
    equal the level set's own gate, at the start and after the front has
    retreated inside the t=0 extent, where the cold start's t=0-extent rule
    (``ocean_drag_cells``) leaves the retreated cells drag-free and the two
    gates disagree (every production restart of October 2026 read ||F||
    about 1e8 there against a recorded 1e-3)."""
    from firedrake import CellDiameter
    from icepack2_tools.front import (facet_neighbours, levelset_drag_gate,
                                      ocean_drag_cells)
    mesh, Q0, h, b, u, ls = _state("prescribed")
    diam = Function(Q0).interpolate(CellDiameter(mesh)).dat.data_ro
    extent0 = h.dat.data_ro >= HMIN
    gate = levelset_drag_gate(ls.phi.dat.data_ro, diam)
    assert np.array_equal(gate, ls.drag_mask.dat.data_ro)
    assert gate.sum() > 0, "the open half of the square carries drag"
    # the ice retreats by two cell columns: the t=0 extent is unchanged, the
    # level set follows the thickness on its next advance
    xc = _xc(mesh)
    h.dat.data[(xc > 0.25) & (xc < 0.5)] = 0.0
    ls.advance(0.1, u, rate=Function(Q0).interpolate(Constant(0.0)))
    gate = levelset_drag_gate(ls.phi.dat.data_ro, diam)
    assert np.array_equal(gate, ls.drag_mask.dat.data_ro)
    ice = h.dat.data_ro >= HMIN
    cold = ocean_drag_cells(ice, facet_neighbours(Q0), extent0)
    retreated = extent0 & ~ice
    assert not cold[retreated].any(), "the t=0-extent rule never drags inside the t=0 extent"
    assert gate[retreated].any(), "the level set drags the retreated cells beyond its strip"
    assert not np.array_equal(gate, cold.astype(float))


def test_the_shed_fraction_is_the_swept_area_over_the_cell():
    r"""``min(1, c dt L / A)``: the fraction a front cell loses is the area the
    front sweeps through it, which is what makes the tally a volume. Every
    front cell here has L = 1/8 and A = 1/128, so ``0.05 * 0.5 * 16 = 0.4``."""
    mesh, Q0, h, b, u, ls = _state("prescribed")
    c = Function(Q0).interpolate(Constant(0.05))
    ls.advance(0.5, u, rate=c)
    beyond, frac = ls.calving_masks()
    front = ls.front_len.dat.data_ro > 0.0
    assert front.sum() == N
    assert np.allclose(ls.front_len.dat.data_ro[front], FRONT_LEN)
    assert np.allclose(ls.cell_area[front], CELL_AREA)
    assert np.allclose(ls.c_cell.dat.data_ro[front], 0.05)
    assert np.allclose(frac[front], 0.4)
    assert np.all(frac[~front] == 0.0)
    assert not beyond.any()


def test_retreat_scales_with_the_step():
    r"""Two half steps remove what one full step does, to the order of the
    scheme: a run that halves its timestep must not calve a different amount.

    The masks are applied the way the transport applies them: cells beyond
    the front are emptied, then each front cell loses ``frac`` of what it
    still holds. One step of fraction ``f`` sheds ``f H``; two of ``f / 2``
    shed ``(f / 2)(2 - f / 2) H``, which differs from it at ``O(f^2)``.
    """
    c0, dt = 0.02, 0.4

    def shed(steps):
        mesh, Q0, h, b, u, ls = _state("prescribed")
        c = Function(Q0).interpolate(Constant(c0))
        total = 0.0
        for _ in range(steps):
            ls.advance(dt / steps, u, rate=c)
            beyond, frac = ls.calving_masks()
            data = h.dat.data
            total += float((data[beyond] * ls.cell_area[beyond]).sum())
            data[beyond] = 0.0
            removed = data * frac
            total += float((removed * ls.cell_area).sum())
            data -= removed
        return total

    f = c0 * dt * FRONT_LEN / CELL_AREA
    one, two = shed(1), shed(2)
    assert one == pytest.approx(N * H_ICE * CELL_AREA * f, rel=1e-9)
    assert two == pytest.approx(
        N * H_ICE * CELL_AREA * (f / 2) * (2 - f / 2), rel=1e-9)
    assert two == pytest.approx(one, rel=0.05)


def test_a_zero_rate_law_leaves_the_front_where_it_was():
    mesh, Q0, h, b, u, ls = _state("prescribed")
    before = ls.phi.dat.data_ro.copy()
    ls.advance(0.1, u, rate=Function(Q0).interpolate(Constant(0.0)))
    assert np.allclose(ls.phi.dat.data_ro, before, atol=1e-9)


def test_the_drag_gate_is_off_at_the_front_and_on_beyond_it():
    r"""The floor-cell ocean drag must not reach the cells the front runs
    through, or it damps the ice rather than the water. The gate is what the
    momentum residual holds, so it is checked against the distance itself."""
    mesh, Q0, h, b, u, ls = _state("prescribed")
    ls.advance(0.1, u, rate=Constant(0.02))
    gate = ls.drag_mask.dat.data_ro
    phi = ls.phi.dat.data_ro
    diam = ls.cell_diam.dat.data_ro
    assert np.all(gate[phi <= diam] == 0.0)
    assert np.all(gate[phi > diam] == 1.0)


def test_an_unknown_law_is_refused(calving_env):
    with pytest.raises(ValueError, match="front law"):
        _state("nonesuch")
    with pytest.raises(ValueError, match="ISMIP7_CALVING=nonesuch: unknown calving law"):
        _configured_law(calving_env, "nonesuch")


def test_a_law_without_its_input_is_refused(calving_env):
    r"""Both refusals matter: a missing rate or a missing rheology would
    otherwise be read as a zero rate, and the front would silently hold."""
    mesh, _, h, b, u, ls = _state("prescribed")
    with pytest.raises(ValueError, match="rate"):
        ls.advance(0.1, u)
    law = _configured_law(calving_env, "vonmises_strain")
    with pytest.raises(ValueError, match="fluidity"):
        law.rate(_front_state(mesh, h, b, u, ls), 0.0)


def test_floating_only_confines_the_law_to_floating_ice(calving_env):
    r"""``ISMIP7_CALVING_FLOATING_ONLY=1`` gates the configured law on
    grounded ice (a marine grounded cliff no longer calves at the rate cap:
    Kemp Land under ssp370, 22 Oct 2026); the recorded description says so,
    so a resume under the other rule is refused; the knob takes only 1 or 0."""
    law = _configured_law(calving_env, "hfb", "sigma_max=0.15")
    assert not law.grounded_gate and "grounded ice does not calve" not in law.describe()
    calving_env.setenv("ISMIP7_CALVING_FLOATING_ONLY", "1")
    gated = runconfig.calving_law_object()
    assert gated.grounded_gate and "grounded ice does not calve" in gated.describe()
    assert gated.describe() != law.describe()
    # a grounded cliff under a tensile membrane stress: the law's own rule
    # calves it, the gated law does not
    mesh, Q0, h, b, u, ls = _state(gated.front_mode, bed=BED_GROUNDED)
    M = Function(TensorFunctionSpace(mesh, "DG", 0, symmetry=True)).interpolate(
        Constant(((1.0, 0.0), (0.0, 0.0))))
    grounded = FrontState(u, M, None, h, b, ls, rho_i=ice_density, rho_w=water_density)
    ice = _xc(mesh) < 0.5
    assert np.all(Function(Q0).interpolate(law.rate(grounded, 0.0)).dat.data_ro[ice] > 0.0)
    assert np.all(Function(Q0).interpolate(gated.rate(grounded, 0.0)).dat.data_ro[ice] == 0.0)
    calving_env.setenv("ISMIP7_CALVING_FLOATING_ONLY", "0")
    assert not runconfig.calving_law_object().grounded_gate
    calving_env.setenv("ISMIP7_CALVING_FLOATING_ONLY", "yes")
    with pytest.raises(ValueError, match="ISMIP7_CALVING_FLOATING_ONLY"):
        runconfig.calving_law_object()


def test_the_parameters_are_the_laws_own(calving_env):
    r"""``ISMIP7_CALVING_PARAMS`` is checked against the law it configures when
    the law is made, which the forward does before its MAP load: a misspelt
    or impossible parameter fails in seconds, not after the initial solve."""
    law = _configured_law(calving_env, "hfb", "sigma_max=0.15,mode=zero_stress")
    assert (law.p["sigma_max"], law.p["mode"]) == (0.15, "zero_stress")
    assert "sigma_max=0.15" in law.describe()
    calving_env.setenv("ISMIP7_CALVING_PARAMS", "sigma_maks=0.15")
    with pytest.raises(ValueError, match="ISMIP7_CALVING=hfb: .*unknown parameter"):
        runconfig.calving_law_object()
    calving_env.setenv("ISMIP7_CALVING_PARAMS", "mode=nye")
    with pytest.raises(ValueError, match="mode"):
        runconfig.calving_law_object()


def test_no_law_means_no_parameters_and_no_level_set(calving_env):
    assert runconfig.calving_law() == "none"
    assert runconfig.calving_law_object() is None
    calving_env.setenv("ISMIP7_CALVING_PARAMS", "sigma_max_fl=0.2")
    with pytest.raises(ValueError, match="configure no law"):
        runconfig.calving_law()


def test_the_retired_threshold_knobs_are_refused_under_a_law(calving_env):
    r"""The von Mises thresholds used to be two knobs of their own. Set with a
    law configured they would now be ignored, so they are refused with the
    spelling that replaces them; under ``none`` they never meant anything."""
    calving_env.setenv("ISMIP7_CALVING_SIGMA_MAX_FLOATING", "0.2")
    assert runconfig.calving_law() == "none"
    calving_env.setenv("ISMIP7_CALVING", "vonmises")
    with pytest.raises(ValueError, match="ISMIP7_CALVING_PARAMS=sigma_max_fl="):
        runconfig.calving_law()


def test_a_law_of_our_own_joins_the_registry_once(calving_env, tmp_path):
    r"""``ISMIP7_CALVING_MODULE`` registers a law from a file before the name
    is looked up, and a second lookup in the same process does not register
    it twice (the forward resolves the law at startup and again in the run)."""
    path = tmp_path / "undercut_law.py"
    path.write_text(textwrap.dedent('''
        from firedrake import Constant
        from icepack_tools import calving


        @calving.register
        class Undercut(calving.Law):
            name = "undercut_test"
            defaults = {"c": 25.0}

            def rate(self, model, t):
                return Constant(float(self.p["c"]))
    '''))
    from icepack_tools import calving
    try:
        calving_env.setenv("ISMIP7_CALVING_MODULE", str(path))
        law = _configured_law(calving_env, "undercut_test", "c=40")
        assert law.p["c"] == 40
        assert runconfig.calving_law_object().p["c"] == 40
        calving_env.setenv("ISMIP7_CALVING_MODULE", str(tmp_path / "missing.py"))
        with pytest.raises(ValueError, match="not a file"):
            runconfig.calving_law_object()
    finally:
        calving.LAWS.pop("undercut_test", None)


def test_resolving_no_law_stays_free_of_firedrake():
    r"""``runconfig`` is imported by the preflight gate and the core report,
    which must stay fast; under the default ``none`` resolving the law must
    not reach the registry, and so not Firedrake."""
    code = ("import os, sys\n"
            "os.environ.pop('ISMIP7_CALVING', None)\n"
            "from icepack2_tools import runconfig\n"
            "assert runconfig.calving_law_object() is None\n"
            "assert 'firedrake' not in sys.modules, 'resolving none imported firedrake'\n")
    r = subprocess.run([sys.executable, "-c", code], cwd=REPO, capture_output=True, text=True,
                       env=dict(os.environ, PYTHONPATH=REPO + os.pathsep
                                + os.environ.get("PYTHONPATH", "")))
    assert r.returncode == 0, r.stderr


def test_fixed_anchors_on_the_t0_extent_not_the_restarted_one():
    r"""A resumed run must not re-freeze the front where it had already
    retreated to: that would permanently bar cells a continuous run of the
    same length would keep. ``phi_init`` from the t=0 thickness is what the
    forward passes for exactly this reason.
    """
    mesh = UnitSquareMesh(N, N)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    h0 = Function(Q0).interpolate(H_ICE * (x < 0.5))       # the t=0 extent
    h_now = Function(Q0).interpolate(H_ICE * (x < 0.25))   # already retreated
    phi0 = initial_distance(mesh, h0, h_min=HMIN)
    ls = LevelSet(mesh, h_now, law="fixed", h_min=HMIN, phi_init=phi0)
    # the forward advances before it asks for the masks, and the advance
    # re-solves phi from the current extent, so the bar must survive it
    u = Function(VectorFunctionSpace(mesh, "CG", 1))
    assert ls.advance(0.1, u) == 0.0
    beyond, _ = ls.calving_masks()
    xc = _xc(mesh)
    # the bar is the ORIGINAL front, so the band between them is still open
    assert np.all(beyond == (xc > 0.5))
    assert not beyond[(xc > 0.25) & (xc < 0.5)].any()
