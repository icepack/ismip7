r"""The SMB-elevation feedback: ``SMB(t) = SMB_base(t) + dacabfdz(t) (s - s_ref)``.

Decided on icepack/ismip7 issue 116 and on by default. These pin the
arithmetic, the reference surface, the gates that refuse a run on a missing
or shifted gradient, and the rule that a chain carries the feedback from its
cold start or not at all. The 32 km comparison of core 7 with and without it
is in the run records (``antarctica/runlog/*-i116on.json`` and
``*-i116off.json``).

The reader tests build trees of empty files: ``available_years`` reads names
only, which is all the gates ask about.
"""
import os

import numpy as np
import pytest

from icepack2_tools.forcing import (
    OCX, OCX_ATMOSPHERE_SOURCE, SMB_FEEDBACK_MARKER, ISMIP7Atmosphere,
    SMBElevationFeedback, build_smb_feedback, feedback_mode,
    flotation_surface, make_forcing_callback, smb_feedback_banner,
    smb_feedback_restart_error, smb_kgm2s_to_myr,
)
from icepack2_tools.runconfig import smb_elevation_feedback

RHO_RATIO = 917.0 / 1024.0
G = 1e-8                          # dacabfdz, kg m-2 s-1 per m (the files' p99 is 4.5e-8)
G_MYR = smb_kgm2s_to_myr(G)       # m/yr of ice per m


class _Field:
    r"""Stands in for a Firedrake Function: the feedback touches ``.dat`` only."""

    class _Dat:
        def __init__(self, values):
            self.data = np.array(values, dtype=float)

        @property
        def data_ro(self):
            return self.data

    def __init__(self, values):
        self.dat = self._Dat(values)


class _Gradient:
    r"""Stands in for ISMIP7Atmosphere: ``dacabfdz`` per year, and the years
    it was asked for."""
    esm, scenario = "CESM2-WACCM", "ssp585"

    def __init__(self, per_year):
        self.per_year = per_year
        self.reads = []

    def get_smb_gradient(self, year, mesh_x, mesh_y):
        self.reads.append(year)
        return np.full(len(mesh_x), self.per_year[year])

    def get_smb(self, year, mesh_x, mesh_y, anomaly=True):
        return np.full(len(mesh_x), 0.25)


def _ctx(b, H_init, h):
    n = len(b)
    return {"geom_xy": (np.zeros(n), np.zeros(n)), "b": _Field(b),
            "H_init": _Field(H_init), "h": _Field(h), "rho_ratio": RHO_RATIO,
            "accum": _Field(np.zeros(n))}


# Three cells: grounded (bed 100 m), afloat (bed -800 m, 400 m of ice) and
# grounded at t=0 but lifting off once 100 m thinner (bed -300 m, 400 m).
B = [100.0, -800.0, -300.0]
H0 = [1000.0, 400.0, 400.0]


# --- the knob ---------------------------------------------------------------

@pytest.mark.parametrize("value, on", [(None, True), ("1", True), (" 1 ", True),
                                       ("0", False), ("", False)])
def test_the_feedback_is_on_unless_turned_off(monkeypatch, value, on):
    if value is None:
        monkeypatch.delenv("ISMIP7_SMB_ELEVATION_FEEDBACK", raising=False)
    else:
        monkeypatch.setenv("ISMIP7_SMB_ELEVATION_FEEDBACK", value)
    assert smb_elevation_feedback() is on


@pytest.mark.parametrize("value", ["yes", "true", "on", "2", "dacabfdz"])
def test_the_knob_has_one_spelling_each_way(monkeypatch, value):
    monkeypatch.setenv("ISMIP7_SMB_ELEVATION_FEEDBACK", value)
    with pytest.raises(ValueError, match="ISMIP7_SMB_ELEVATION_FEEDBACK"):
        smb_elevation_feedback()


# --- the correction ---------------------------------------------------------

def test_no_surface_change_means_no_feedback():
    r"""At a cold start h is H_init, so the correction is exactly zero on
    every kind of cell: the ``balance`` apparent-MB reference folds the first
    step's SMB into a_ref and relies on it."""
    fb = SMBElevationFeedback(_Gradient({2015: 4.5e-8}))
    corr = fb.correction(_ctx(B, H0, H0), 2015)
    assert np.array_equal(corr, np.zeros(3))


def test_the_reference_is_rebuilt_from_H_init_and_ignores_the_stored_surface():
    r"""A cold start can load its surface from the MAP file, which need not
    match the flotation formula; the feedback never reads it."""
    ctx = _ctx(B, H0, H0)
    ctx["s"] = _Field([5000.0, -5000.0, 1.0])
    assert np.array_equal(SMBElevationFeedback(_Gradient({2015: G})).correction(ctx, 2015),
                          np.zeros(3))


def test_a_lowering_surface_takes_the_gradients_sign():
    r"""100 m of thinning: the grounded surface drops 100 m, the floating one
    drops by its freeboard share, and the cell that lifts off drops to its
    flotation surface. A positive dacabfdz lowers the SMB of all three."""
    fb = SMBElevationFeedback(_Gradient({2015: G}))
    h = [h0 - 100.0 for h0 in H0]
    corr = fb.correction(_ctx(B, H0, h), 2015)
    ds_grounded = -100.0
    ds_floating = -(1.0 - RHO_RATIO) * 100.0
    ds_lifted = (1.0 - RHO_RATIO) * 300.0 - (-300.0 + 400.0)
    assert np.allclose(corr, G_MYR * np.array([ds_grounded, ds_floating, ds_lifted]),
                       rtol=1e-12, atol=0.0)
    assert (corr < 0).all()


def test_the_units_are_metres_of_ice_per_year():
    r"""The file's kg m-2 s-1 per m, times 100 m of lowering, is
    1e-8 * 100 * 31556926 / 917 m/yr, the conversion acabf takes (issue #29)."""
    fb = SMBElevationFeedback(_Gradient({2015: G}))
    corr = fb.correction(_ctx([100.0], [1000.0], [900.0]), 2015)
    assert corr[0] == pytest.approx(-G * 100.0 * 31556926.0 / 917.0, rel=1e-12)


def test_the_gradient_is_read_once_per_forcing_year():
    grad = _Gradient({2015: G, 2016: 2 * G})
    fb = SMBElevationFeedback(grad)
    ctx = _ctx([100.0], [1000.0], [900.0])
    first = fb.correction(ctx, 2015)
    fb.correction(ctx, 2015)
    second = fb.correction(ctx, 2016)
    assert grad.reads == [2015, 2016]
    assert second[0] == pytest.approx(2 * first[0], rel=1e-12)


def test_the_flotation_surface_is_the_forwards():
    b, h = np.array([100.0, -800.0]), np.array([1000.0, 400.0])
    assert np.allclose(flotation_surface(b, h, RHO_RATIO),
                       [1100.0, (1.0 - RHO_RATIO) * 400.0])


# --- the forcing callback ---------------------------------------------------

def test_the_callback_adds_the_feedback_to_the_forced_smb():
    base = np.array([0.1, 0.2, 0.3])
    atm = _Gradient({2015: G})
    callback = make_forcing_callback(atm=atm, smb_baseline=base,
                                     smb_feedback=SMBElevationFeedback(atm))
    h = [h0 - 100.0 for h0 in H0]
    ctx = _ctx(B, H0, h)
    callback(ctx, 2015.1)
    expect = 0.25 + base + SMBElevationFeedback(_Gradient({2015: G})).correction(ctx, 2015)
    assert np.allclose(ctx["accum"].dat.data, expect, rtol=1e-14, atol=0.0)


def test_the_controls_fixed_smb_is_rewritten_so_the_feedback_never_compounds():
    r"""The control assigns its climatology once and builds its callback with
    no atmosphere. Adding the feedback onto ``accum`` itself would stack every
    step's correction on the last; the callback writes baseline + feedback
    from the baseline every step instead."""
    base = np.array([0.1, 0.2, 0.3])
    fb = SMBElevationFeedback(_Gradient({2015: G}))
    callback = make_forcing_callback(smb_baseline=base, smb_feedback=fb)
    ctx = _ctx(B, H0, [h0 - 10.0 for h0 in H0])
    for k in range(1, 6):
        callback(ctx, 2015.0 + 0.1 * k)
    ctx["h"] = _Field([h0 - 100.0 for h0 in H0])
    callback(ctx, 2015.6)
    expect = base + SMBElevationFeedback(_Gradient({2015: G})).correction(ctx, 2015)
    assert np.allclose(ctx["accum"].dat.data, expect, rtol=1e-14, atol=0.0)


def test_without_the_feedback_the_control_writes_its_climatology():
    base = np.array([0.1, 0.2, 0.3])
    ctx = _ctx(B, H0, H0)
    make_forcing_callback(smb_baseline=base)(ctx, 2015.1)
    assert np.array_equal(ctx["accum"].dat.data, base)


def test_a_feedback_with_no_smb_to_add_to_is_refused():
    with pytest.raises(ValueError, match="needs an SMB to add to"):
        make_forcing_callback(smb_feedback=SMBElevationFeedback(_Gradient({})))


# --- the gates: a missing or shifted gradient refuses the run ---------------

def _gradient_tree(root, years, esm="CESM2-WACCM", scenario="ssp585",
                   product="SDBN1-8000m", version="v2"):
    if scenario == OCX:
        d = os.path.join(root, OCX, esm, product, "dacabfdz", version)
    else:
        d = os.path.join(root, esm, scenario, product, "dacabfdz", version)
    os.makedirs(d, exist_ok=True)
    for y in years:
        open(os.path.join(
            d, f"dacabfdz_AIS_{esm}_{scenario}_{product}_{version}_{y}.nc"), "wb").close()


def _atm(root, esm="CESM2-WACCM", scenario="ssp585"):
    return ISMIP7Atmosphere(data_root=str(root), esm=esm, scenario=scenario)


def test_a_full_series_and_the_one_year_bridge_are_covered(tmp_path):
    _gradient_tree(tmp_path, range(2015, 2300))        # 2015..2299; 2300 is bridged
    atm = _atm(tmp_path)
    assert atm.coverage_problem(2015, 2300, "dacabfdz") is None
    assert atm.require_years(2015, 2300, ("dacabfdz",)) == (2015, 2299)


@pytest.mark.parametrize("years, says", [
    ([], "no files"),
    ([y for y in range(2015, 2301) if y != 2100], "(2100..2100)"),
    (range(2015, 2298), "missing (2299..2300)"),
    (range(2016, 2301), "missing (2015..2015)"),
])
def test_an_absent_short_or_holed_gradient_refuses_the_run(tmp_path, years, says):
    r"""The reader turns an absent variable into zeros, which would be a run
    with no feedback that says it has one; any other gap raises mid-run."""
    if years:
        _gradient_tree(tmp_path, years)
    atm = _atm(tmp_path)
    assert says in atm.coverage_problem(2015, 2300, "dacabfdz")
    with pytest.raises(FileNotFoundError) as e:
        SMBElevationFeedback(atm).check(2015, 2300)
    assert says in str(e.value)
    assert "ISMIP7_SMB_ELEVATION_FEEDBACK=0" in str(e.value)
    assert "data/CESM2-WACCM/ssp585/SDBN1-8000m/dacabfdz/" in str(e.value)


def test_with_the_feedback_off_nothing_is_checked(monkeypatch, tmp_path):
    monkeypatch.setenv("ISMIP7_SMB_ELEVATION_FEEDBACK", "0")
    fb = build_smb_feedback(_atm(tmp_path), 2015, 2300)      # an empty tree
    assert fb is None and feedback_mode(fb) == "off"
    monkeypatch.delenv("ISMIP7_SMB_ELEVATION_FEEDBACK")
    with pytest.raises(FileNotFoundError):
        build_smb_feedback(_atm(tmp_path), 2015, 2300)


def test_the_shifted_ocx_gradient_is_refused(tmp_path):
    r"""The OCX ``dacabfdz`` v1 was spatially shifted (discussion #45). With
    no v2 beside it the reader falls back to it, so the feedback refuses."""
    _gradient_tree(tmp_path, range(1979, 2026), esm=OCX_ATMOSPHERE_SOURCE,
                   scenario=OCX, version="v1")
    atm = _atm(tmp_path, esm=OCX_ATMOSPHERE_SOURCE, scenario=OCX)
    assert "is v1 on disk" in atm.version_problem("dacabfdz")
    with pytest.raises(FileNotFoundError, match="v2 or newer"):
        SMBElevationFeedback(atm).check(1979, 2025)
    _gradient_tree(tmp_path, range(1979, 2026), esm=OCX_ATMOSPHERE_SOURCE,
                   scenario=OCX, version="v2")
    assert atm.version_problem("dacabfdz") is None
    SMBElevationFeedback(atm).check(1979, 2025)


def test_the_version_floor_is_for_the_ocx_gradient_only(tmp_path):
    r"""The ESM gradients at v1 (MRI-ESM2-0 keeps one beside its v2) are
    usable, per the focus group; only the OCX v1 is wrong."""
    _gradient_tree(tmp_path, range(2015, 2301), esm="MRI-ESM2-0",
                   product="GEMB-SDBN1-8000m", version="v1")
    assert _atm(tmp_path, esm="MRI-ESM2-0").version_problem("dacabfdz") is None


def test_the_banner_names_the_gradient_the_run_reads(tmp_path):
    _gradient_tree(tmp_path, range(2015, 2301))
    fb = SMBElevationFeedback(_atm(tmp_path))
    banner = smb_feedback_banner(fb)
    assert banner.startswith(SMB_FEEDBACK_MARKER)
    assert "dacabfdz from CESM2-WACCM ssp585 SDBN1-8000m v2" in banner
    assert feedback_mode(fb) == "dacabfdz"
    assert smb_feedback_banner(None) == (
        f"{SMB_FEEDBACK_MARKER} off (ISMIP7_SMB_ELEVATION_FEEDBACK=0)")
    (line,) = fb.provenance()
    assert "atmosphere dacabfdz CESM2-WACCM ssp585 SDBN1-8000m v2" in line


# --- a chain carries the feedback from its cold start or not at all ---------

@pytest.mark.parametrize("recorded, requested, adapted, refused", [
    (None, "dacabfdz", False, True),         # a checkpoint from before the feedback
    (None, "off", False, False),
    ("dacabfdz", "dacabfdz", False, False),
    ("dacabfdz", "off", False, True),
    ("off", "dacabfdz", False, True),
    ("off", "off", False, False),
    (None, "dacabfdz", True, False),         # an adapted t=0 state is a new start
    ("dacabfdz", "off", True, False),
    ("dacabfdz", None, False, False),        # a caller that applies no forcing
    (None, None, False, False),
])
def test_the_restart_rule(recorded, requested, adapted, refused):
    problem = smb_feedback_restart_error(recorded, requested, adapted_initial=adapted)
    assert (problem is not None) is refused
    if refused:
        assert "ISMIP7_SMB_ELEVATION_FEEDBACK=" in problem
        assert "from its cold start" in problem


# --- on a real DG0 geometry --------------------------------------------------

def test_on_a_dg0_mesh_the_reference_is_the_forwards_surface_and_the_log_adds_up():
    r"""The arrays the feedback reads are the geometry dofs, and ``s_ref``
    from H_init is the forward's own surface: the UFL expression simulation.py
    interpolates, cell by cell. The yearly line's net is the integral of the
    correction."""
    fd = pytest.importorskip("firedrake")
    # 800 km across, so the net is Gt/yr and not a rounded zero; the second
    # half is a marine bed on which 60 m of thinning floats some cells
    mesh = fd.RectangleMesh(4, 4, 800e3, 800e3)
    Q = fd.FunctionSpace(mesh, "DG", 0)
    x, _ = fd.SpatialCoordinate(mesh)
    b = fd.Function(Q).interpolate(fd.conditional(x < 400e3, 200.0, -600.0))
    H_init = fd.Function(Q).interpolate(fd.Constant(500.0) + x / 2000.0)
    h = fd.Function(Q).interpolate(H_init - 60.0)
    rho_ratio = fd.Constant(RHO_RATIO)
    s_ufl = fd.Function(Q).interpolate(
        fd.max_value(b + H_init, (fd.Constant(1.0) - rho_ratio) * H_init))
    assert np.allclose(flotation_surface(b.dat.data_ro, H_init.dat.data_ro, RHO_RATIO),
                       s_ufl.dat.data_ro, rtol=1e-14, atol=0.0)

    xy = fd.Function(fd.VectorFunctionSpace(mesh, Q.ufl_element())).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    accum = fd.Function(Q)
    ctx = {"mesh": mesh, "geom_xy": (xy[:, 0].copy(), xy[:, 1].copy()), "b": b,
           "H_init": H_init, "h": h, "rho_ratio": rho_ratio, "accum": accum}
    lines = []
    fb = SMBElevationFeedback(_Gradient({2015: G}), log=lines.append)
    corr = fb.correction(ctx, 2015)
    ds = (flotation_surface(b.dat.data_ro, h.dat.data_ro, RHO_RATIO)
          - s_ufl.dat.data_ro)
    assert np.allclose(corr, G_MYR * ds, rtol=1e-12, atol=0.0)
    field = fd.Function(Q)
    field.dat.data[:] = corr
    net = float(fd.assemble(field * fd.dx)) * 917.0 / 1e12
    (line,) = lines
    assert line.strip().startswith("dacabfdz feedback 2015: net")
    assert abs(net) > 1.0 and f"net {net:+.2f} Gt/yr" in line
    assert SMB_FEEDBACK_MARKER not in line
    fb.correction(ctx, 2015)
    assert len(lines) == 1                         # once per forcing year
