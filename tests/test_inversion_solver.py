"""MUMPS's analysis under the inversion and the forward: the PT-Scotch fallback, the
ISMIP7_MUMPS_ANALYSIS knob, and the cache fingerprint it leaves alone."""
import os
import sys
from pathlib import Path

import pytest

from icepack2_tools import solverconfig as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "antarctica" / "scripts"))


def test_mumps_ordering_follows_the_petsc_build(monkeypatch):
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: False)
    opts = sc._mumps_options("x_")
    assert "x_mat_mumps_icntl_29" not in opts and opts["x_pc_factor_mat_solver_type"] == "mumps"
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: True)
    monkeypatch.delenv("ISMIP7_MUMPS_ANALYSIS", raising=False)
    opts = sc._mumps_options("x_")
    assert opts["x_mat_mumps_icntl_28"] == 2 and opts["x_mat_mumps_icntl_29"] == 1
    monkeypatch.setenv("ISMIP7_MUMPS_ANALYSIS", "sequential")
    assert "x_mat_mumps_icntl_28" not in sc._mumps_options("x_")


# solver_configuration_fingerprint(solver_provenance(mode)) at 0e267b2, in a
# clean environment: the identity every published cache was prepared under.
PREPARED_CACHE_FINGERPRINTS = {
    "scpc_mumps": "ae74e2957ed02c3def0a033c0947fe169e324aa6b56a41688f6b5e34ed56b4e8",
    "schur_mumps": "df55bc6c0bd794fb078b40cbcd4398fef0211486501d83ddd172cfb3d896958c",
}


@pytest.mark.parametrize("mode", sorted(PREPARED_CACHE_FINGERPRINTS))
@pytest.mark.parametrize("ptscotch,analysis", [
    (True, "parallel"), (False, "parallel"), (True, "sequential")])
def test_mumps_analysis_leaves_the_cache_fingerprint_alone(
        monkeypatch, mode, ptscotch, analysis):
    """The analysis follows the build and the knob, and the record says which ran;
    the fingerprint stays the one existing caches were prepared under."""
    import timing_campaign as tc
    for name in list(os.environ):
        if name.startswith("ISMIP7_"):
            monkeypatch.delenv(name)
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: ptscotch)
    monkeypatch.setenv("ISMIP7_MUMPS_ANALYSIS", analysis)
    provenance = sc.solver_provenance(mode)
    recorded = {key: value for key, value in provenance["diagnostic_petsc_options"].items()
                if key.endswith("mat_mumps_icntl_28")}
    assert bool(recorded) == (ptscotch and analysis == "parallel")
    assert tc.solver_configuration_fingerprint(provenance) == PREPARED_CACHE_FINGERPRINTS[mode]


@pytest.mark.parametrize("ptscotch", [True, False])
def test_mumps_analysis_knob_is_checked(monkeypatch, ptscotch):
    """A bad value is refused on every build, including one without PT-Scotch."""
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: ptscotch)
    monkeypatch.setenv("ISMIP7_MUMPS_ANALYSIS", "distributed")
    with pytest.raises(ValueError):
        sc.mumps_analysis()
    with pytest.raises(ValueError):
        sc._mumps_options("x_")
    with pytest.raises(ValueError):
        sc.diagnostic_solver_parameters("scpc_mumps")
