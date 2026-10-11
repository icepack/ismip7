r"""``ISMIP7_GEOMETRY_RAMP``: the cold start's geometry ramp from a relaxed
MAP's own geometry (``runconfig.geometry_ramp_steps``)."""
import pytest

from icepack2_tools.runconfig import geometry_ramp_steps


def test_the_ramp_is_off_by_default(monkeypatch):
    monkeypatch.delenv("ISMIP7_GEOMETRY_RAMP", raising=False)
    assert geometry_ramp_steps() == 0
    monkeypatch.setenv("ISMIP7_GEOMETRY_RAMP", "")
    assert geometry_ramp_steps() == 0


def test_the_ramp_takes_a_step_count(monkeypatch):
    monkeypatch.setenv("ISMIP7_GEOMETRY_RAMP", " 10 ")
    assert geometry_ramp_steps() == 10


@pytest.mark.parametrize("raw", ["-1", "two", "2.5"])
def test_a_bad_count_stops_the_run(monkeypatch, raw):
    monkeypatch.setenv("ISMIP7_GEOMETRY_RAMP", raw)
    with pytest.raises(ValueError):
        geometry_ramp_steps()


def test_the_seed_is_unset_by_default(monkeypatch):
    from icepack2_tools.runconfig import seed_state
    monkeypatch.delenv("ISMIP7_SEED_STATE", raising=False)
    assert seed_state() is None
    monkeypatch.setenv("ISMIP7_SEED_STATE", "  ")
    assert seed_state() is None
    monkeypatch.setenv("ISMIP7_SEED_STATE", " /x/run_1000_t2004.0.h5 ")
    assert seed_state() == "/x/run_1000_t2004.0.h5"
