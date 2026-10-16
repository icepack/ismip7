r"""The relaxed initial state (icepack2_tools.relaxation): the year the
relaxation runs, what it refuses, what its end state records, what a
re-inversion takes from it, and what a forward records it started from."""
import os
import sys

import numpy as np
import pytest

from icepack2_tools.handoff import OBJECTIVE_KEYS, OBJECTIVE_RECORD_KEYS
from icepack2_tools.relaxation import (
    END_STATE_ATTR, GEOMETRY_METHOD_RELAXED, RELAX_KEYS, anchor_ratio_counts,
    end_state_problems, inherited_geometry, init_state, is_relaxed,
    relax_backdate_years, relax_dt, relax_environment_problems,
    describe_relaxation, relax_experiment_name, relax_forcing,
    relax_fssa_attrs, relax_window, relaxation_attrs, relaxed_map_name,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MAP = {
    "misfit_norm": "logvel", "log_vel_weight": 85380.4, "log_vel_eps": 1.0,
    "gamma_theta": 24.93, "friction": "regularized_coulomb", "lake_ice_base": 1,
    "drag_gate": "vertex", "h_visc_floor": 2.5,
    "objective_total": 2592.0, "objective_iteration": 1036.0,
    "geometry_source_method": "target-native-bedmachine-cell-average-v1",
}


def _attrs(**over):
    return {**relaxation_attrs(
        source_map="/maps/inversion_icepack2_rc_n3_dg0_2000.h5", source_sha256="ab12",
        source_attrs={**MAP, **over}, t_start=2014.0, t_end=2015.0, dt=0.0125,
        backdate_years=1.0, forcing="ocx protocol"),
        **relax_fssa_attrs(1.0, "step")}


def test_the_year_is_2014_to_2015_at_half_the_production_step(monkeypatch):
    for k in ("ISMIP7_RELAX_START", "ISMIP7_RELAX_DT", "ISMIP7_GEOMETRY_BACKDATE",
              "ISMIP7_RELAX_FORCING"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("ISMIP7_DT", "0.025")
    assert relax_window() == (2014.0, 2015.0)
    assert relax_dt(2014.0, 2015.0) == pytest.approx(0.0125)
    assert relax_backdate_years(2014.0) == 1.0
    assert relax_forcing() == "ocx"
    monkeypatch.setenv("ISMIP7_DT", "0.05")         # follows the production step
    assert relax_dt(2014.0, 2015.0) == pytest.approx(0.025)


def test_the_window_and_step_are_refused_outside_their_rules(monkeypatch):
    monkeypatch.setenv("ISMIP7_RELAX_START", "2015")
    with pytest.raises(ValueError, match="before the 2015 geometry"):
        relax_window()
    monkeypatch.setenv("ISMIP7_RELAX_START", "2001")
    with pytest.raises(ValueError, match="2003"):
        relax_window()
    monkeypatch.setenv("ISMIP7_RELAX_DT", "0.3")
    with pytest.raises(ValueError, match="does not divide"):
        relax_dt(2014.0, 2015.0)
    monkeypatch.setenv("ISMIP7_GEOMETRY_BACKDATE", "12")
    with pytest.raises(ValueError, match="undoes 1 years"):
        relax_backdate_years(2014.0)
    monkeypatch.setenv("ISMIP7_GEOMETRY_BACKDATE", "1")
    assert relax_backdate_years(2014.0) == 1.0
    monkeypatch.setenv("ISMIP7_RELAX_FORCING", "ctrl")
    with pytest.raises(ValueError, match="ISMIP7_RELAX_FORCING"):
        relax_forcing()


def test_the_year_runs_free_of_the_apparent_mb_with_a_pinned_front():
    ok = dict(apparent_mb=None, output=False, calving="none", fixed_front=True)
    assert relax_environment_problems(**ok) == []
    for change, word in ((dict(apparent_mb="balance"), "ISMIP7_APPARENT_MB=0"),
                         (dict(output=True), "ISMIP7_OUTPUT=0"),
                         (dict(calving="vonmises"), "ISMIP7_CALVING=none"),
                         (dict(fixed_front=False), "ISMIP7_FIXED_FRONT=1")):
        problems = relax_environment_problems(**{**ok, **change})
        assert len(problems) == 1 and word in problems[0]


def test_the_end_state_carries_the_map_s_objective_and_leaves_its_value_behind():
    out = _attrs()
    assert out[END_STATE_ATTR] == 1
    for key in ("misfit_norm", "log_vel_weight", "gamma_theta", "lake_ice_base",
                "drag_gate", "h_visc_floor"):
        assert out[key] == MAP[key]
    assert set(out) - set(OBJECTIVE_KEYS) - {END_STATE_ATTR} == set(RELAX_KEYS)
    assert not set(out) & set(OBJECTIVE_RECORD_KEYS)
    assert out["relax_source_objective_total"] == 2592.0
    assert out["relax_source_map"] == "/maps/inversion_icepack2_rc_n3_dg0_2000.h5"
    assert (out["relax_t_start"], out["relax_t_end"], out["relax_dt"]) == (2014.0, 2015.0, 0.0125)


def test_the_end_state_records_the_stabilization_its_year_ran_under():
    assert relax_fssa_attrs(1, "step") == {
        "relax_fssa_theta": 1.0, "relax_fssa_reference": "step"}
    # off: the forward resolves no reference
    assert relax_fssa_attrs(0.0, None) == {
        "relax_fssa_theta": 0.0, "relax_fssa_reference": "none"}
    assert "FSSA theta 1, step reference" in describe_relaxation(_attrs())
    off = {**_attrs(), **relax_fssa_attrs(0.0, None)}
    assert "FSSA off" in describe_relaxation(off)
    # a relaxation made before the record held it
    older = {k: v for k, v in _attrs().items() if not k.startswith("relax_fssa_")}
    assert "FSSA unrecorded" in describe_relaxation(older)


def test_bytes_and_numpy_attributes_come_through_as_values():
    out = _attrs(misfit_norm=b"logvel", log_vel_weight=np.float64(5.0), lake_ice_base=np.int64(1))
    assert out["misfit_norm"] == "logvel" and type(out["log_vel_weight"]) is float
    assert type(out["lake_ice_base"]) is int


def test_a_relaxed_map_is_not_relaxed_again():
    with pytest.raises(ValueError, match="relaxed already"):
        _attrs(geometry_source_method=GEOMETRY_METHOD_RELAXED)
    with pytest.raises(ValueError, match="relaxed already"):
        _attrs(**{END_STATE_ATTR: 1})


RELAXED_OUT = "/maps/final/inversion_icepack2_rc_n3_dg0_2000_relax2014.h5"


def test_an_end_state_seeds_a_re_inversion_only_finished_and_on_its_mesh():
    state = {**_attrs(), "t_yr": 2015.0, "stalled": 0}
    ok = dict(same_mesh=True, geometry_taken=True, map_out=RELAXED_OUT)
    assert end_state_problems(state, **ok) == []
    short = {**state, "t_yr": 2014.6}
    assert "finish the relaxation" in end_state_problems(short, **ok)[0]
    assert "stalled" in end_state_problems({**state, "stalled": 1}, **ok)[0]
    assert "ISMIP7_MESH=checkpoint" in end_state_problems(
        state, **{**ok, "same_mesh": False})[0]
    assert "not taken" in end_state_problems(state, **{**ok, "geometry_taken": False})[0]


def test_a_re_inversion_writes_only_the_relaxed_map_name():
    state = {**_attrs(), "t_yr": 2015.0, "stalled": 0}
    # the production MAP the forwards load by default, which it would replace
    (why,) = end_state_problems(state, same_mesh=True, geometry_taken=True,
                                map_out="inversion_icepack2_rc_n3_dg0_2000.h5")
    assert "ISMIP7_MAP_OUT" in why and "_relax2014.h5" in why
    # another start year names another MAP
    assert end_state_problems(
        state, same_mesh=True, geometry_taken=True,
        map_out=RELAXED_OUT.replace("relax2014", "relax2013"))


def test_a_re_inversion_records_the_end_state_it_took_its_geometry_from():
    state = {**_attrs(), "t_yr": 2015.0, "geometry_source": "/maps/x.h5",
             "geometry_source_method": "checkpoint-native-v1"}
    assert is_relaxed(state)
    rec = inherited_geometry(state, geometry_taken=True,
                             warm_basename="relax_x_2000_final.h5", warm_sha256="cd34")
    assert rec["geometry_source_method"] == GEOMETRY_METHOD_RELAXED
    assert rec["geometry_source"] == rec["relax_state"] == "relax_x_2000_final.h5"
    assert rec["relax_state_sha256"] == "cd34"
    assert rec["relax_source_sha256"] == "ab12"
    assert (rec["relax_fssa_theta"], rec["relax_fssa_reference"]) == (1.0, "step")
    assert not set(rec) & set(OBJECTIVE_KEYS)
    # the next link of the chain resumes the re-inversion's own checkpoint
    link2 = inherited_geometry({**rec, "log_vel_weight": 1.0}, geometry_taken=True,
                               warm_basename="map_relax2014.h5")
    assert link2 == rec
    with pytest.raises(RuntimeError, match="does not take that geometry"):
        inherited_geometry(rec, geometry_taken=False, warm_basename="map_relax2014.h5")


def test_an_observed_warm_start_keeps_this_run_s_own_record():
    assert inherited_geometry(MAP, geometry_taken=True, warm_basename="m.h5") is None
    assert inherited_geometry({}, geometry_taken=False, warm_basename="m.h5") is None


def test_a_forward_records_which_initial_state_it_started_from():
    relaxed = {"geometry_source_method": GEOMETRY_METHOD_RELAXED}
    assert init_state(MAP) == "observed"
    # a relaxed MAP gives its controls alone, on its own mesh as on any other
    assert init_state(relaxed) == "relaxed-controls"


def test_names():
    src = "/maps/inversion_icepack2_rc_n3_dg0_2000_int5000.h5"
    assert relaxed_map_name(src, 2014.0) == "inversion_icepack2_rc_n3_dg0_2000_int5000_relax2014.h5"
    assert relax_experiment_name(src) == "relax_inversion_icepack2_rc_n3_dg0_2000_int5000"
    assert relax_experiment_name(src, "smoke").endswith("_int5000_smoke")


def test_the_anchor_ratio_counts_grounded_cells_only():
    c_old = np.array([1.0, 1.0, 1.0, 1.0, 2.0])
    c_new = np.array([1.0, 1.2, 1.5, 0.5, 9.0])
    grounded = np.array([True, True, True, True, False])
    n, n1, n3, s, lo, hi = anchor_ratio_counts(c_new, c_old, grounded)
    assert (n, n1, n3) == (4, 3, 2)
    assert s == pytest.approx(abs(np.log(1.2)) + np.log(1.5) + np.log(2.0))
    assert (lo, hi) == (pytest.approx(np.log(0.5)), pytest.approx(np.log(1.5)))
    assert anchor_ratio_counts(c_new, c_old, np.zeros(5, bool))[:4] == (0, 0, 0, 0.0)


# --- the driver refuses before the setup -----------------------------------

def _driver_env(monkeypatch, tmp_path):
    for k in ("ISMIP7_RELAX_START", "ISMIP7_RELAX_DT", "ISMIP7_GEOMETRY_BACKDATE",
              "ISMIP7_CALVING", "ISMIP7_CALVING_PARAMS", "ISMIP7_RESTART"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("ISMIP7_APPARENT_MB", "0")
    monkeypatch.setenv("ISMIP7_OUTPUT", "0")
    monkeypatch.setenv("ISMIP7_FIXED_FRONT", "1")
    monkeypatch.setenv("ISMIP7_RELAX_FORCING", "none")
    monkeypatch.setenv("ISMIP7_AUTO_RESUME", "0")


def _driver():
    pytest.importorskip("firedrake")
    import importlib
    sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))
    return importlib.import_module("relaxation.run")


def test_the_driver_refuses_the_apparent_mb_before_reading_the_map(monkeypatch, tmp_path):
    run = _driver()
    _driver_env(monkeypatch, tmp_path)
    fake = tmp_path / "map.h5"
    fake.write_bytes(b"not read")
    monkeypatch.setenv("ISMIP7_INVERSION", str(fake))
    monkeypatch.setenv("ISMIP7_APPARENT_MB", "1")
    with pytest.raises(RuntimeError, match="ISMIP7_APPARENT_MB=0"):
        run.main()
    monkeypatch.delenv("ISMIP7_INVERSION")
    with pytest.raises(RuntimeError, match="ISMIP7_INVERSION must name"):
        run.main()


def test_the_driver_refuses_a_relaxed_map_before_the_setup(monkeypatch, tmp_path):
    h5py = pytest.importorskip("h5py")
    run = _driver()
    _driver_env(monkeypatch, tmp_path)
    relaxed = tmp_path / "map_relax2014.h5"
    with h5py.File(relaxed, "w") as h:
        h["/"].attrs["geometry_source_method"] = GEOMETRY_METHOD_RELAXED
    monkeypatch.setenv("ISMIP7_INVERSION", str(relaxed))
    with pytest.raises(ValueError, match="relaxed already"):
        run.main()


def test_the_driver_resumes_only_its_own_end_state(monkeypatch, tmp_path):
    h5py = pytest.importorskip("h5py")
    from icepack2_tools.runconfig import file_sha256
    run = _driver()
    _driver_env(monkeypatch, tmp_path)
    source = tmp_path / "map.h5"
    with h5py.File(source, "w") as h:
        h["/"].attrs["geometry_source_method"] = MAP["geometry_source_method"]
    monkeypatch.setenv("ISMIP7_INVERSION", str(source))
    # a MAP re-inverted from a relaxation of this MAP carries its sha256 too
    relaxed = tmp_path / "map_relax2014.h5"
    with h5py.File(relaxed, "w") as h:
        h["/"].attrs["relax_source_sha256"] = file_sha256(str(source))
        h["/"].attrs["geometry_source_method"] = GEOMETRY_METHOD_RELAXED
    monkeypatch.setenv("ISMIP7_RESTART", str(relaxed))
    with pytest.raises(RuntimeError, match="not a relaxation"):
        run.main()


def test_a_warm_start_s_geometry_is_taken_only_under_its_own_sampling(monkeypatch):
    r"""A taken geometry keeps the sampling it was built with and the MAP
    records the run's, so forcing one sampled another way is refused."""
    from icepack2_tools.runconfig import warm_start_geometry
    monkeypatch.delenv("ISMIP7_WARM_START_GEOMETRY", raising=False)
    same = dict(same_mesh=True, same_lake=True)
    assert warm_start_geometry(**same, warm_sampling="vertex_front",
                               run_sampling="vertex_front")
    assert not warm_start_geometry(**same, warm_sampling="vertex",
                                   run_sampling="vertex_front")
    assert not warm_start_geometry(same_mesh=False, same_lake=True,
                                   warm_sampling="vertex", run_sampling="vertex")
    monkeypatch.setenv("ISMIP7_WARM_START_GEOMETRY", "1")
    assert warm_start_geometry(same_mesh=True, same_lake=False,
                               warm_sampling="vertex", run_sampling="vertex")
    with pytest.raises(ValueError, match="ISMIP7_RASTER_SAMPLE=vertex"):
        warm_start_geometry(**same, warm_sampling="vertex",
                            run_sampling="vertex_front")
    monkeypatch.setenv("ISMIP7_WARM_START_GEOMETRY", "0")
    assert not warm_start_geometry(**same, warm_sampling="vertex_front",
                                   run_sampling="vertex_front")
