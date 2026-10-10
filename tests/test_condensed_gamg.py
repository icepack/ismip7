r"""scpc_gamg iterates on the condensed system, and its tuning never reaches MUMPS.

The first scpc_gamg lane (quartz job 10517922, 2500/25000 x 16) ran one GAMG
V-cycle per outer FGMRES iteration on the matrix-free mixed system: 0.257 s per
iteration, 4894 of them, 126 s/step against condensed MUMPS's 33. The mode now
puts the Krylov method on the assembled condensed system and hands GAMG the
rigid-body modes. Every knob of that lands in the record's PETSc options, and
none of it may move the scpc_mumps options: those fingerprint every prepared
cache, which a GAMG lane starts from.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "antarctica" / "scripts"))

import timing_campaign as tc  # noqa: E402
from build_timing_matrix import SNES_COUNT_MARK, _solver_work  # noqa: E402
from icepack2_tools.solverconfig import (  # noqa: E402
    diagnostic_solver_parameters,
    solver_provenance,
)

KNOBS = (
    "ISMIP7_CONDENSED_KSP_TYPE",
    "ISMIP7_CONDENSED_KSP_ATOL_FACTOR",
    "ISMIP7_CONDENSED_KSP_RTOL",
    "ISMIP7_CONDENSED_KSP_RESTART",
    "ISMIP7_CONDENSED_NEAR_NULLSPACE",
    "ISMIP7_CONDENSED_PETSC_OPTIONS",
)


@pytest.fixture(autouse=True)
def no_knobs(monkeypatch):
    for name in KNOBS + ("ISMIP7_KSP_RTOL", "ISMIP7_KSP_MAXIT"):
        monkeypatch.delenv(name, raising=False)


def condensed(mode="scpc_gamg"):
    prefix = "condensed_field_"
    return {
        key[len(prefix):]: value
        for key, value in diagnostic_solver_parameters(mode).items()
        if key.startswith(prefix)
    }


def test_the_krylov_method_is_on_the_condensed_system():
    options = condensed()
    assert options["ksp_type"] == "fgmres"  # the coarse solve is GMRES
    assert diagnostic_solver_parameters("scpc_gamg")["ksp_type"] == "fgmres"
    assert options["ksp_gmres_restart"] == 100
    assert options["pc_type"] == "gamg"
    # Antarctica is held by drag: the rotation bought 2 % of the V-cycles and
    # made each 14 % dearer.
    assert options["near_nullspace"] == "none"


def test_the_inner_solve_stops_where_the_outer_tolerance_is_met(monkeypatch):
    r"""FGMRES hands the preconditioner unit vectors and the elimination is
    exact, so the outer relative residual after one iteration is the inner
    absolute residual. An absolute inner tolerance just under the outer
    relative one is the loosest that leaves the outer solve one iteration
    (19 V-cycles a solve on quartz where a relative 1e-7 spent 36)."""
    outer = diagnostic_solver_parameters("scpc_gamg")["ksp_rtol"]
    options = condensed()
    assert options["ksp_atol"] == pytest.approx(0.5 * outer)
    assert options["ksp_atol"] < outer
    assert options["ksp_rtol"] <= 1e-10  # the relative test never decides
    # It follows the outer tolerance rather than standing beside it.
    monkeypatch.setenv("ISMIP7_KSP_RTOL", "1e-4")
    assert condensed()["ksp_atol"] == pytest.approx(5e-5)
    monkeypatch.setenv("ISMIP7_CONDENSED_KSP_ATOL_FACTOR", "0.1")
    assert condensed()["ksp_atol"] == pytest.approx(1e-5)


def test_the_single_v_cycle_configuration_is_one_knob_away(monkeypatch):
    monkeypatch.setenv("ISMIP7_CONDENSED_KSP_TYPE", "preonly")
    options = condensed()
    assert options["ksp_type"] == "preonly"
    assert not any(key.startswith("ksp_") and key != "ksp_type" for key in options)
    monkeypatch.setenv("ISMIP7_CONDENSED_NEAR_NULLSPACE", "rigid_body")
    assert condensed()["near_nullspace"] == "rigid_body"


def test_a_rung_s_extra_options_are_applied_last_and_recorded(monkeypatch):
    monkeypatch.setenv(
        "ISMIP7_CONDENSED_PETSC_OPTIONS",
        "pc_gamg_threshold=0.02 -mg_levels_ksp_max_it=4 pc_gamg_asm_use_agg ksp_rtol=1e-8",
    )
    options = condensed()
    assert options["pc_gamg_threshold"] == "0.02"
    assert options["mg_levels_ksp_max_it"] == "4"
    assert options["pc_gamg_asm_use_agg"] is None
    assert options["ksp_rtol"] == "1e-8"
    recorded = solver_provenance("scpc_gamg")["diagnostic_petsc_options"]
    assert recorded["condensed_field_pc_gamg_threshold"] == "0.02"


@pytest.mark.parametrize("name, value", [
    ("ISMIP7_CONDENSED_NEAR_NULLSPACE", "rotations"),
    ("ISMIP7_CONDENSED_PETSC_OPTIONS", "condensed_field_pc_gamg_threshold=0.02"),
    ("ISMIP7_CONDENSED_PETSC_OPTIONS", "=0.02"),
])
def test_a_misspelt_knob_is_refused(monkeypatch, name, value):
    monkeypatch.setenv(name, value)
    with pytest.raises(ValueError):
        diagnostic_solver_parameters("scpc_gamg")


def test_no_gamg_knob_moves_the_cache_solver_s_fingerprint(monkeypatch):
    before = tc.solver_configuration_fingerprint(solver_provenance("scpc_mumps"))
    from icepack2_tools.solverconfig import _have_ptscotch
    expected = {
        "mat_type": "aij", "ksp_type": "preonly", "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    if _have_ptscotch():
        expected.update({"mat_mumps_icntl_28": 2, "mat_mumps_icntl_29": 1})
    assert condensed("scpc_mumps") == expected
    monkeypatch.setenv("ISMIP7_CONDENSED_KSP_TYPE", "gmres")
    monkeypatch.setenv("ISMIP7_CONDENSED_KSP_RTOL", "1e-9")
    monkeypatch.setenv("ISMIP7_CONDENSED_NEAR_NULLSPACE", "none")
    monkeypatch.setenv("ISMIP7_CONDENSED_PETSC_OPTIONS", "pc_gamg_threshold=0.02")
    after = tc.solver_configuration_fingerprint(solver_provenance("scpc_mumps"))
    assert after == before
    # ...while the lane's own record does tell the rungs apart.
    assert tc.solver_configuration_fingerprint(
        solver_provenance("scpc_gamg")
    ) != before


def test_solver_work_counts_the_line_search_s_solves_when_the_record_has_them():
    r"""SNES's linear-iteration count leaves out the solve NLEQ-ERR makes for its
    simplified Newton step. Job 10517922 recorded 1872 iterations and made 4894."""
    summary = {"count": 10, "snes_iterations_total": 112, "linear_iterations_total": 1872}
    assert _solver_work({"diagnostic_solve_summary": summary}) == f"11.2 × 16.7{SNES_COUNT_MARK}"
    summary["condensed_iterations_total"] = 4894
    assert _solver_work({"diagnostic_solve_summary": summary}) == "11.2 × 43.7"


def test_a_matrix_free_jacobian_is_frozen_at_the_newton_iterate_by_default(monkeypatch):
    r"""scpc_* apply the Jacobian matrix-free, at whatever the state Function
    holds; the line search had moved it to its trial point before solving for
    the simplified Newton step (1134 of a MUMPS lane's 1246 outer iterations).
    The smoke script (`make solver-smoke`) holds the solver itself to this."""
    from icepack2_tools.solverconfig import linearization_state

    monkeypatch.delenv("ISMIP7_FREEZE_LINEARIZATION", raising=False)
    assert linearization_state("scpc_mumps") == "frozen"
    assert linearization_state("scpc_gamg") == "frozen"
    # An assembled Jacobian never followed the state.
    assert linearization_state("full_mumps") == "assembled"
    assert linearization_state("schur_gamg") == "assembled"
    assert solver_provenance("scpc_gamg")["linearization_state"] == "frozen"
    before = tc.solver_configuration_fingerprint(solver_provenance("scpc_mumps"))

    monkeypatch.setenv("ISMIP7_FREEZE_LINEARIZATION", "0")
    assert linearization_state("scpc_mumps") == "live"
    assert linearization_state("full_mumps") == "assembled"
    assert solver_provenance("scpc_mumps")["linearization_state"] == "live"
    # The record says which; the caches, which are converged states, do not care.
    assert tc.solver_configuration_fingerprint(solver_provenance("scpc_mumps")) == before
