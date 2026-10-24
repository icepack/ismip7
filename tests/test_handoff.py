r"""One objective across restarts (icepack2_tools.handoff): the settings
comparison, the accepted-iterate lookup, and the objective handoff gap.
Pure Python, no mesh."""
import numpy as np
import pytest

from icepack2_tools.handoff import (
    OBJECTIVE_KEYS,
    SUBELEMENT_SCHEME_VERSIONS,
    accepted_evaluation,
    check_subelement_record,
    frozen_in_control,
    handoff_gap,
    objective_mismatches,
)


def _settings(**over):
    base = {
        "misfit_norm": "sigma", "misfit_scale": 918512.0,
        "log_vel_weight": 85380.44865839917, "log_vel_eps": 1.0,
        "dhdt_weight": 1.0, "dhdt_net_sigma": 10.0, "dhdt_filter_l": 0.0,
        "dhdt_speed_ratio": 30.0, "dhdt_speed_slow": 10.0, "dhdt_speed_fast": 100.0,
        "prior_form": "bilaplacian", "gamma_theta": 1e5, "gamma_phi": 1e5,
        "prior_sigma_theta": 0.3, "prior_sigma_phi": 0.3, "prior_rho": 7500.0,
        "prior_sigma_alpha": 0.1828, "prior_rho_theta": 7500.0, "friction_c_ref": 0.0334,
        "friction_control": "sqrt", "friction": "budd", "n_flow": 3.0,
        "geometry_space": "dg0", "friction_anchor_length": 20000.0,
        "lake_ice_base": 1, "fluidity_prior_origin": "thermomechanical",
        "grad_precond": "mass_consistent",
        "subelement_friction": 1, "subelement_scheme": "sep1",
        "subelement_scheme_version": 2, "exact_front": 0,
        "fluidity_control": "all", "drag_gate": "vertex", "h_visc_floor": 1.0,
        "phi_grounded": "inverted", "raster_sample": "vertex",
        "vel_sample": "point", "u_obs_max": 0.0, "dhdt_reweight": 0.0,
        "dhdt_reweight_scale": 1.0, "kill_bergs": 0,
        "velocity_sample": "point",
    }
    base.update(over)
    return base


def test_identical_settings_have_no_mismatch():
    assert objective_mismatches(_settings(), _settings()) == []


def test_every_objective_setting_is_compared():
    for key in OBJECTIVE_KEYS:
        current = _settings()
        current[key] = "changed" if isinstance(current[key], str) else current[key] * 2 + 1
        out = objective_mismatches(_settings(), current)
        assert len(out) == 1 and out[0].startswith(f"{key}:"), key


def test_numbers_are_compared_as_numbers_and_missing_keys_are_not_mismatches():
    recorded = _settings(prior_rho="7500.0")          # an attribute read back as text
    assert objective_mismatches(recorded, _settings()) == []
    old = _settings()
    del old["misfit_scale"], old["prior_sigma_alpha"]  # an older checkpoint
    assert objective_mismatches(old, _settings()) == []
    # the auto weight re-derived on a restart is exactly the drift to catch
    out = objective_mismatches(_settings(), _settings(log_vel_weight=56646.35))
    assert out == ["log_vel_weight: 85380.44865839917 -> 56646.35"]


def test_a_sub_element_record_without_a_scheme_was_sep2():
    legacy = _settings(subelement_friction=1)
    del legacy["subelement_scheme"], legacy["subelement_scheme_version"]
    out = objective_mismatches(legacy, _settings(subelement_friction=1, subelement_scheme="sep1"))
    assert out == ["subelement_scheme: 'sep2' -> 'sep1'",
                   "subelement_scheme_version: 1 -> 2"]
    assert objective_mismatches(legacy, _settings(subelement_friction=1,
                                                  subelement_scheme="sep2",
                                                  subelement_scheme_version=1)) == []
    # without sub-element friction the scheme selects nothing
    off = dict(subelement_friction=0)
    plain = _settings(**off)
    del plain["subelement_scheme"]
    assert objective_mismatches(plain, _settings(**off)) == []
    assert objective_mismatches(_settings(subelement_scheme="sep2", **off),
                                _settings(subelement_scheme="sep1", **off)) == []
    no_record = _settings(subelement_scheme="sep2")
    del no_record["subelement_friction"]
    assert objective_mismatches(no_record, _settings(subelement_scheme="sep1", **off)) == []


def test_the_versions_are_the_forms_this_code_builds():
    # bump a scheme's version with any change to its residual at the same
    # control fields (icepack2_tools.handoff, the comment over the table)
    assert SUBELEMENT_SCHEME_VERSIONS == {"sep2": 1, "sep1": 2}


@pytest.mark.parametrize("scheme, control, stands_for", [
    ("sep2", "log", 1), ("sep2", "sqrt", 1),
    # both SEP1 forms agree where theta is zero
    ("sep1", "sqrt", 2), ("sep1", "exp", 2),
    # a log-control SEP1 record cannot say which form it was
    ("sep1", "log", "unrecorded"),
])
def test_a_sub_element_record_without_a_version(scheme, control, stands_for):
    old = _settings(subelement_scheme=scheme, friction_control=control)
    del old["subelement_scheme_version"]
    current = _settings(subelement_scheme=scheme, friction_control=control,
                        subelement_scheme_version=SUBELEMENT_SCHEME_VERSIONS[scheme])
    out = objective_mismatches(old, current)
    if stands_for == SUBELEMENT_SCHEME_VERSIONS[scheme]:
        assert out == []
    else:
        assert out == [f"subelement_scheme_version: {stands_for!r} -> "
                       f"{SUBELEMENT_SCHEME_VERSIONS[scheme]}"]
    # without sub-element friction on both sides the version selects nothing
    old["subelement_friction"] = 0
    assert objective_mismatches(old, dict(current, subelement_friction=0)) == []


def test_another_version_of_the_scheme_is_another_objective():
    out = objective_mismatches(_settings(subelement_scheme_version=1), _settings())
    assert out == ["subelement_scheme_version: 1 -> 2"]


def test_a_forward_rebuilds_only_the_version_this_code_builds():
    assert check_subelement_record("sep2", None, "log") == 1
    assert check_subelement_record("sep1", None, "sqrt") == 2
    assert check_subelement_record("sep1", np.int64(2), "log") == 2
    with pytest.raises(RuntimeError, match="before the scheme's version was recorded"):
        check_subelement_record("sep1", None, "log", source="old.h5")
    with pytest.raises(RuntimeError, match="version 1; this code builds version 2"):
        check_subelement_record("sep1", 1, "sqrt")
    with pytest.raises(RuntimeError, match="sub-element scheme 'sep3'"):
        check_subelement_record("sep3", 1, "log")


def test_the_accepted_evaluation_is_the_one_matching_the_objective():
    ring = [(1.0e9, "trial a"), (7.5e8, "accepted"), (7.9e8, "trial b")]
    assert accepted_evaluation(ring, 7.5e8)[1] == "accepted"
    assert accepted_evaluation(ring, 7.5e8 * (1 + 1e-12))[1] == "accepted"
    assert accepted_evaluation(ring, 7.0e8) is None
    assert accepted_evaluation([], 1.0) is None


def test_the_handoff_gap_is_relative():
    assert handoff_gap(1.0e9, 1.0e9) == 0.0
    assert handoff_gap(1.0e9, 1.001e9) == pytest.approx(1e-3)
    assert np.isfinite(handoff_gap(0.0, 1.0))


def test_auto_weight_frozen_only_within_one_control():
    sqrt_record = _settings()
    assert frozen_in_control(sqrt_record, "prior_sigma_alpha", "sqrt") == 0.1828
    assert frozen_in_control(sqrt_record, "prior_sigma_alpha", "exp") is None
    exp_record = _settings(friction_control="exp", prior_sigma_alpha=1.0)
    assert frozen_in_control(exp_record, "prior_sigma_alpha", "sqrt") is None
    assert frozen_in_control(exp_record, "friction_c_ref", "exp") == 0.0334
    legacy = {"prior_sigma_alpha": 0.2}
    assert frozen_in_control(legacy, "prior_sigma_alpha", "sqrt") is None
    assert frozen_in_control(_settings(friction_c_ref=0.0), "friction_c_ref", "sqrt") is None


def test_a_new_cliff_push_version_is_an_objective_change():
    r"""A MAP inverted under exact_front 1 warm-starts a version 2 inversion
    only non-strictly: the push at rock walls changed (issue #166)."""
    assert objective_mismatches(_settings(exact_front=1), _settings(exact_front=2)) == \
        ["exact_front: 1 -> 2"]
