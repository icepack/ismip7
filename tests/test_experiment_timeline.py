r"""The protocol timeline the core drivers hand to the forward.

``t = Y.0`` is 1 January of year Y (see ``AnnualOutput``: "the year being
accumulated, its Jan 1 has passed"), so a run covering years A through B runs
from ``A.0`` to ``B+1.0``. Discussions #16 and #20: a projection covers
2015-2300 and a historical covers its start through 2014, which puts the
historical's final checkpoint at 2015.0, exactly where the projection begins.

Get the end wrong by one and two things break at once: the handoff checkpoint
lands at 2014.0, so the projection's first forcing call asks the SSP tree for
year 2014, which precedes that scenario and is an error; and the banked ISMIP7
series ends a year early, so the submission carries 2015-2299 under filenames
the readiness document says read 2015-2300.

Each driver is executed for real with ``run_core_experiment`` captured, so
this asserts the period the driver actually requests.
"""
import importlib
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")
sys.path.insert(0, SCRIPTS)

pytest.importorskip("firedrake")
experiment = importlib.import_module("experiment")     # noqa: E402

#: driver module -> the calendar years the protocol says it covers. The
#: historicals start in 2003 from the 2015 geometry with the Smith et al.
#: (2020) thinning undone (issue #117, the 25 September meeting).
HISTORICAL = {
    "historical.cesm_waccm": (2003, 2014),
    "historical.mri_esm2": (2003, 2014),
}
PROJECTIONS = {
    "projections.ssp126_cesm_waccm": (2015, 2300),
    "projections.ssp126_mri_esm2": (2015, 2300),
    "projections.ssp370_cesm_waccm": (2015, 2100),
    "projections.ssp370_mri_esm2": (2015, 2100),
    "projections.ssp585_cesm_waccm": (2015, 2300),
    "projections.ssp585_mri_esm2": (2015, 2300),
}


def _period(monkeypatch, module_name):
    r"""The (t_start, t_end) the driver hands run_core_experiment."""
    captured = {}

    def capture(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(experiment, "run_core_experiment", capture)
    mod = importlib.import_module(module_name)
    monkeypatch.setattr(mod, "run_core_experiment", capture, raising=False)
    # the drivers call it under `if __name__ == "__main__"`, so replay that
    # entry point rather than trusting the module-level constants
    code = open(mod.__file__).read()
    exec(compile(code, mod.__file__, "exec"),
         {"__name__": "__main__", "__file__": mod.__file__,
          "run_core_experiment": capture, "os": os, "sys": sys})
    return captured["t_start_default"], captured["t_end_default"]


@pytest.mark.parametrize("module_name,years", sorted(PROJECTIONS.items()))
def test_a_projection_runs_through_the_last_submitted_year(monkeypatch, module_name, years):
    first, last = years
    t_start, t_end = _period(monkeypatch, module_name)
    assert t_start == float(first)
    # t_end is 1 January of the year AFTER the last one covered, so the final
    # year_end fires at t_end and banks `last`
    assert t_end == float(last + 1)


@pytest.mark.parametrize("module_name,years", sorted(HISTORICAL.items()))
def test_a_historical_runs_through_the_last_submitted_year(monkeypatch, module_name, years):
    first, last = years
    t_start, t_end = _period(monkeypatch, module_name)
    assert t_start == float(first)
    assert t_end == float(last + 1)


def test_the_ocx_driver_uses_the_same_end_year_convention():
    r"""OCX does not go through run_core_experiment, so it carries its own
    T_START/T_END. It has to follow the same rule or core 11 banks a year
    fewer than its title claims."""
    ocx = importlib.import_module("projections.ocx")
    first, last = 2003, 2025          # from the backdated 2003 geometry (issue #117) to the OCX forcing's end
    assert ocx.T_START == float(first)
    assert ocx.T_END == float(last + 1)


def test_the_historical_hands_off_where_the_projection_starts(monkeypatch):
    r"""run_simulation takes t_start from the restart checkpoint, so the
    historical's end IS the projection's first step. If they disagree the
    projection simulates a year its own scenario tree does not cover."""
    hist_ends = {_period(monkeypatch, m)[1] for m in HISTORICAL}
    proj_starts = {_period(monkeypatch, m)[0] for m in PROJECTIONS}
    assert hist_ends == proj_starts == {2015.0}
