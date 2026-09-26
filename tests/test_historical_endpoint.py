r"""A control or projection branches only from a historical that finished.

A historical chain rewrites ``hist_<esm>_<lc>_final.h5`` at the end of every
job, so a chain that stopped early leaves one holding the year it reached.
``simulation.historical_endpoint`` is what both drivers ask before branching.
"""
import importlib.util
import os

import pytest
from firedrake import CheckpointFile, UnitSquareMesh

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")


@pytest.fixture
def simulation(tmp_path, monkeypatch):
    import sys
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    spec = importlib.util.spec_from_file_location(
        "simulation", os.path.join(SCRIPTS, "simulation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS_DIR", str(tmp_path))
    return mod


def _endpoint(sim, t_yr, name="hist_cesm2_waccm_prod_1000_final.h5"):
    path = os.path.join(sim.RESULTS_DIR, name)
    with CheckpointFile(path, "w") as chk:
        chk.save_mesh(UnitSquareMesh(2, 2))
        if t_yr is not None:
            chk.set_attr("/", "t_yr", t_yr)
    return path


def test_no_historical_means_no_endpoint(simulation):
    assert simulation.historical_endpoint("cesm2_waccm", "_prod", 2015.0, 1000) is None


def test_a_finished_historical_is_the_endpoint(simulation):
    path = _endpoint(simulation, 2015.0)
    assert simulation.historical_endpoint("cesm2_waccm", "_prod", 2015.0, 1000) == path


def test_a_historical_that_stopped_early_is_refused(simulation):
    _endpoint(simulation, 1950.3)
    with pytest.raises(RuntimeError, match="t_yr=1950.3, short of the 2015 handoff"):
        simulation.historical_endpoint("cesm2_waccm", "_prod", 2015.0, 1000)


def test_an_endpoint_without_a_year_is_refused(simulation):
    _endpoint(simulation, None)
    with pytest.raises(RuntimeError, match="no t_yr"):
        simulation.historical_endpoint("cesm2_waccm", "_prod", 2015.0, 1000)
