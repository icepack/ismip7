r"""Single source of truth for transient PETSc solver configuration.

This module is deliberately pure Python.  The simulation imports the PETSc
option dictionaries from here, while timing/provenance tools import the same
functions without paying the Firedrake import cost or duplicating defaults.

The named modes distinguish two materially different meanings of "MUMPS":

``schur_mumps``
    MUMPS factors PETSc's approximate velocity Schur preconditioner only.
``full_mumps``
    MUMPS factors the complete mixed ``(u, M, tau)`` Jacobian.  This is the
    memory-heavy reference configuration used before the scalable-solver work.

``scpc_*`` uses Firedrake Slate to eliminate the two DG0 stress fields exactly
cell by cell.  The remaining assembled CG1 velocity system is preconditioned
by GAMG or MUMPS.  These modes avoid relying on the block size inferred for an
AIJ submatrix, which is the unresolved weakness in PETSc ``selfp`` here.
"""

import math
import os


# A forward driver invoked outside the managed launchers must fall back to the
# established full-Jacobian reference, never to an unqualified development PC.
# The launchers name their own: the timing Makefile exports its campaign's
# solver, and the cluster forward runner (batch_runners/projection.sbatch)
# exports scpc_gamg for production. This default does not follow them, because
# it is not only the forward's: the inversion, whose taped solves follow
# ISMIP7_INVERSION_LINEAR_SOLVER whatever is set here, stamps the mode it
# resolves on its MAP, and redistribute_checkpoint.py fingerprints a
# published cache with it.
DIAGNOSTIC_SOLVER_DEFAULT = "full_mumps"
DIAGNOSTIC_SOLVER_MODES = (
    "schur_gamg",
    "schur_mumps",
    "scpc_gamg",
    "scpc_mumps",
    "full_mumps",
)
DIAGNOSTIC_SOLVER_ALIASES = {
    # Backward-compatible names used by timing runs before the modes were
    # disambiguated.  Reports record both the requested and canonical names.
    "iterative": "schur_gamg",
    "mumps": "schur_mumps",
}
# The inversion's taped solves: the annotated forward, the adjoint tlm_adjoint
# solves against it, and the solve that publishes the MAP's state. A knob of
# its own, so the scpc_mumps the timing campaign exports for the lane's
# ISMIP7_DIAGNOSTIC_LINEAR_SOLVER leaves the taped solve where it was. The
# approximate-Schur modes are out: their selfp preconditioner was never
# qualified for the mixed system. scpc_gamg is the default since the final
# Quartz round of issues 156 and 157 (README, "Inversion solver"): under the
# direct forward it took full_mumps's iterates 1.35 to 1.5 times faster an
# evaluation at 2 km and 4.6 times at 1 km, with half the memory.
INVERSION_SOLVER_DEFAULT = "scpc_gamg"
INVERSION_SOLVER_MODES = ("full_mumps", "scpc_mumps", "scpc_gamg")
# The inversion's taped forward under scpc_* solves its Newton corrections to
# a relative 1e-8, not the transient's 1e-6. Each evaluation starts from the
# last one's state at new controls, as far as ||F|| = 9e10 from the solution
# at 32 km, and there NLEQ-ERR, which judges a step by the norm of the Newton
# correction, could not use corrections solved to 1e-6: at a trial point the
# exact condensation took in 6 Newton iterations it stalled at
# ||F|| = 7.9e10 for 200. At 1e-8 it took the same 6. Backtracking on ||F||
# instead (bt) also took 6 there and costs less where it works (64 s an
# evaluation against 95 s at 2 km), but it damps nearly every step on the
# synthetic slabs: 113 Newton iterations against 4 under the cell-wise law
# and 20 against 4 under the sub-element law, with an exact LU as with GAMG
# (README, "Inversion solver"). ISMIP7_INVERSION_SNES_LINESEARCH and
# ISMIP7_INVERSION_KSP_RTOL override both. The transient is unchanged: its
# steps start next to the solution.
INVERSION_KSP_RTOL_DEFAULT = "1e-8"

SNES_TYPE_DEFAULT = "newtonls"
SNES_LINESEARCH_DEFAULT = "nleqerr"
SNES_RTOL_DEFAULT = "1e-8"
SNES_ATOL_DEFAULT = "1e-50"
SNES_STOL_DEFAULT = "0"
# PETSc 3.25 uses -1 for PETSC_DETERMINE, which restores the default 1e4
# residual-growth cutoff.  PETSC_UNLIMITED is -3.  The timing scouts must let
# the nonlinear method reach its real convergence/max-iteration outcome rather
# than report DIVERGED_DTOL after the first large nleqerr trial step.
SNES_DIVERGENCE_TOL_DEFAULT = "-3"
SNES_MAXIT_DEFAULT = "200"
SNES_ATOL_SCALE_DEFAULT = "100"
SNES_RESTART_FAILURE_ATOL_SCALE_DEFAULT = "1e-6"
SNES_KSP_EW_DEFAULT = "0"
SNES_MONITOR_DEFAULT = "0"
SNES_LOG_DEFAULT = "stdout"
SOLVER_VIEW_DEFAULT = "0"
# scpc_* apply a matrix-free Jacobian, which is the form's action at whatever
# the state Function holds. The NLEQ-ERR line search evaluates the residual at
# its trial point, which Firedrake writes into that Function, and then solves
# for its simplified Newton step: against J(x_trial), with a condensed system
# assembled at x_k. Frozen (the default), the Jacobian is built on a copy of
# the state that is refreshed only when SNES re-forms it
# (preconditioners.frozen_linearization). ``0`` restores the live state every
# lane before 2026-09-19 ran with.
FREEZE_LINEARIZATION_DEFAULT = "1"
# The inversion's publishing solve (final_solve_parameters): a solve that
# starts at an already-converged state must exit at iteration 0, and one that
# does not must stay short and loud rather than grind to the shared 200.
FINAL_SNES_STOL_DEFAULT = "1e-8"
FINAL_SNES_MAXIT_DEFAULT = "50"
FINAL_KSP_MAXIT_DEFAULT = "50"
# The inversion's objective evaluations: first one untaped Newton solve at the
# full exponents from the last converged state, then a taped solve that starts
# converged; one that fails goes to the failed-trial rescue. The in-tape
# 5-stage n=1->3 ladder this replaces restarted at n=1 from the n=3 state:
# every evaluation of the 2 km sub-element chains began at ||F|| 1e11-1e13
# and took ~200 Newton iterations, and its n=2.5->3 stage is where the
# failed trials diverged
# (snes_recinos_sqrtC_pattyn_noprior_subelement_nodhdt.log, 30 Sep 2026).
DIRECT_FORWARD_DEFAULT = "1"
# Successful direct solves took 1-17 Newton iterations on the 2 km and 32 km
# meshes; one that has not converged by 30 is a lost trial (30 Sep 2026).
DIRECT_FORWARD_MAXIT_DEFAULT = "30"
# A direct trial whose residual has grown this many times over its start is
# lost: fail it and let the line search backtrack rather than run the cap out
# (2 km SEP2, NOTS 1692389: 50 iterations to ||F|| 5e33, 289 s).
DIRECT_FORWARD_DTOL_DEFAULT = "1e6"
# Rungs of the continuation ladder the failed-trial rescue may climb: each
# failed rung cost 1-3 h on the 2 km mesh and ~80% of rescues failed anyway
# (2 km chains, 27-30 Sep 2026). 0 disables the rescue.
TRIAL_RESCUE_RUNGS_DEFAULT = "1"

KSP_RTOL_DEFAULT = "1e-6"
KSP_MAXIT_DEFAULT = "1000"
# scpc_gamg iterates on the assembled condensed velocity system, not on the
# matrix-free mixed one. The elimination is exact, so with one V-cycle per
# outer FGMRES iteration (the first configuration) every Krylov iteration pays
# a mixed-Jacobian action and three Slate sweeps over the mesh around its
# V-cycle: 0.19 s against 0.03 s for an iteration on the AIJ system. Flexible
# because the coarse solve is itself GMRES. ``preonly`` restores the single
# V-cycle. Quartz, 2500/25000 x 16, frozen linearization, s/step: single
# V-cycle 79.0; inner solve to rtol 1e-7 27.6; to the absolute tolerance below
# 19.3; condensed MUMPS 13.2.
CONDENSED_KSP_TYPE_DEFAULT = "fgmres"
# The inner solve stops on an ABSOLUTE residual, a fraction of the outer
# relative tolerance. FGMRES hands the preconditioner unit vectors and the
# elimination is exact, so the outer relative residual after one iteration IS
# the inner absolute residual (both 0.01527 after the first V-cycle of job
# 10517922): this is the loosest inner solve that still leaves the outer
# FGMRES one iteration. A relative 1e-7 spent 36 V-cycles a solve where the
# outer tolerance is met after 19. The relative test is parked out of reach.
CONDENSED_KSP_ATOL_FACTOR_DEFAULT = "0.5"
CONDENSED_KSP_RTOL_DEFAULT = "1e-12"
# Solves ran to 85 iterations at the first configuration's convergence rate;
# a restart inside that range stalls them.
CONDENSED_KSP_RESTART_DEFAULT = "100"
# Near-nullspace handed to GAMG for the condensed operator: ``none`` (GAMG's
# default, the two translations) or ``rigid_body`` (adds the in-plane
# rotation, which the membrane operator does not see on a floating shelf).
# The rotation halved the V-cycles of a synthetic shelf problem and did
# nothing for Antarctica, most of which is held by drag: 2 % fewer V-cycles,
# each 14 % dearer on the denser coarse grids (jobs 10520141 / 10520247).
CONDENSED_NEAR_NULLSPACES = ("rigid_body", "none")
CONDENSED_NEAR_NULLSPACE_DEFAULT = "none"
TRANSPORT_KSP_RTOL_DEFAULT = "1e-10"
TRANSPORT_KSP_MAXIT_DEFAULT = "500"
MASS_RESIDUAL_TOL_GT_DEFAULT = "5e-5"
CONTINUATION_STEPS_DEFAULT = "8"
RESCUE_MAXIT_DEFAULT = "600"
SUBCYCLES_DEFAULT = "1,4,16"
# Adaptive substepping (icepack2_tools.substep): off by default. The
# tolerance is metres of backward-Euler local error in the DG0 thickness.
SUBSTEP_ADAPT_DEFAULT = "0"
SUBSTEP_TOL_DEFAULT = "1.0"
SUBSTEP_INIT_DEFAULT = "1"
SUBSTEP_MAX_DEFAULT = "64"
SUBSTEP_QUIET_DEFAULT = "20"
SUBSTEP_HMIN_DEFAULT = "10"
# Free-surface stabilization (icepack2_tools.fssa): theta, 0 = off, and the
# velocity the surface change is measured from (fssa.REFERENCES). On by
# default since 6 October 2026.
FSSA_THETA_DEFAULT = "1"
FSSA_REFERENCE_DEFAULT = "auto"
RESCUE_ENABLED_DEFAULT = "1"


def _env(name, default):
    return os.environ.get(name, default)


def _enabled(name, default="0"):
    return _env(name, default).strip().lower() not in {
        "", "0", "false", "no", "off",
    }


def requested_diagnostic_solver():
    return _env(
        "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER", DIAGNOSTIC_SOLVER_DEFAULT
    ).strip().lower()


def snes_monitor_enabled():
    return _enabled("ISMIP7_SNES_MONITOR", SNES_MONITOR_DEFAULT)


def solver_view_enabled():
    return _enabled("ISMIP7_SOLVER_VIEW", SOLVER_VIEW_DEFAULT)


def linearization_state(mode=None):
    r"""Where a mode's Jacobian is linearized while the line search runs:
    ``assembled`` (an AIJ matrix, fixed at the Newton iterate by construction),
    or for the matrix-free scpc_* modes ``frozen`` at the iterate or ``live``
    (following the state Function to the line search's trial point)."""
    mode = diagnostic_solver_mode(mode)
    if not mode.startswith("scpc_"):
        return "assembled"
    frozen = _enabled("ISMIP7_FREEZE_LINEARIZATION", FREEZE_LINEARIZATION_DEFAULT)
    return "frozen" if frozen else "live"


def diagnostic_solver_mode(requested=None):
    r"""Canonical mode for ``requested`` (default: the environment's)."""
    if requested is None:
        requested = requested_diagnostic_solver()
    requested = str(requested).strip().lower()
    mode = DIAGNOSTIC_SOLVER_ALIASES.get(requested, requested)
    if mode not in DIAGNOSTIC_SOLVER_MODES:
        choices = ", ".join(DIAGNOSTIC_SOLVER_MODES)
        aliases = ", ".join(DIAGNOSTIC_SOLVER_ALIASES)
        raise ValueError(
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER must be one of "
            f"{choices} (legacy aliases: {aliases}), not {requested!r}"
        )
    return mode


def diagnostic_solver_label(mode=None):
    mode = diagnostic_solver_mode() if mode is None else mode
    return {
        "schur_gamg": "approximate-schur-gamg",
        "schur_mumps": "approximate-schur-mumps-ptscotch",
        "scpc_gamg": "exact-local-condensation-gamg",
        "scpc_mumps": "exact-local-condensation-mumps-ptscotch",
        "full_mumps": "full-jacobian-mumps",
    }[mode]


def _nonlinear_options():
    options = {
        "snes_type": _env("ISMIP7_SNES_TYPE", SNES_TYPE_DEFAULT),
        "snes_rtol": float(_env("ISMIP7_SNES_RTOL", SNES_RTOL_DEFAULT)),
        "snes_atol": float(_env("ISMIP7_SNES_ATOL", SNES_ATOL_DEFAULT)),
        "snes_max_it": int(_env("ISMIP7_SNES_MAXIT", SNES_MAXIT_DEFAULT)),
        "snes_linesearch_type": _env(
            "ISMIP7_SNES_LINESEARCH", SNES_LINESEARCH_DEFAULT
        ),
        "snes_divergence_tolerance": float(_env(
            "ISMIP7_SNES_DIVERGENCE_TOL", SNES_DIVERGENCE_TOL_DEFAULT
        )),
        "snes_stol": float(_env("ISMIP7_SNES_STOL", SNES_STOL_DEFAULT)),
    }
    if _enabled("ISMIP7_SNES_KSP_EW", SNES_KSP_EW_DEFAULT):
        options["snes_ksp_ew"] = None
    return options


def _outer_ksp_options():
    return {
        "ksp_type": "fgmres",
        "ksp_rtol": float(_env("ISMIP7_KSP_RTOL", KSP_RTOL_DEFAULT)),
        "ksp_max_it": int(_env("ISMIP7_KSP_MAXIT", KSP_MAXIT_DEFAULT)),
    }


def _gamg_options(prefix=""):
    return {
        f"{prefix}pc_type": "gamg",
        # Do not collapse the coarse problem onto one MPI rank and factor it.
        f"{prefix}pc_gamg_parallel_coarse_grid_solver": None,
        f"{prefix}mg_levels_ksp_type": "chebyshev",
        f"{prefix}mg_levels_pc_type": "jacobi",
        f"{prefix}mg_coarse_ksp_type": "gmres",
        f"{prefix}mg_coarse_ksp_rtol": 1e-2,
        f"{prefix}mg_coarse_ksp_max_it": 50,
        f"{prefix}mg_coarse_pc_type": "jacobi",
    }


def condensed_near_nullspace():
    name = _env(
        "ISMIP7_CONDENSED_NEAR_NULLSPACE", CONDENSED_NEAR_NULLSPACE_DEFAULT
    ).strip().lower()
    if name not in CONDENSED_NEAR_NULLSPACES:
        raise ValueError(
            "ISMIP7_CONDENSED_NEAR_NULLSPACE must be one of "
            f"{CONDENSED_NEAR_NULLSPACES}, not {name!r}"
        )
    return name


def _extra_condensed_options(prefix):
    r"""``ISMIP7_CONDENSED_PETSC_OPTIONS="pc_gamg_threshold=0.02 mg_levels_ksp_max_it=4"``:
    further options of the condensed GAMG solve, for a tuning rung. They are
    applied last and, like every other entry, land in the record's
    ``diagnostic_petsc_options``."""
    options = {}
    for token in _env("ISMIP7_CONDENSED_PETSC_OPTIONS", "").split():
        name, sep, value = token.lstrip("-").partition("=")
        if not name or name.startswith(prefix):
            raise ValueError(
                "ISMIP7_CONDENSED_PETSC_OPTIONS takes unprefixed name=value "
                f"entries (or a bare flag), not {token!r}"
            )
        options[f"{prefix}{name}"] = value if sep else None
    return options


def _condensed_gamg_options(prefix):
    r"""GAMG on SCPC's condensed velocity system: the Krylov method that
    iterates on it (see CONDENSED_KSP_TYPE_DEFAULT), the near-nullspace
    ``ISMIP7SCPC`` attaches to it, and a rung's extra options."""
    params = _gamg_options(prefix)
    ksp_type = _env(
        "ISMIP7_CONDENSED_KSP_TYPE", CONDENSED_KSP_TYPE_DEFAULT
    ).strip().lower()
    params[f"{prefix}ksp_type"] = ksp_type
    if ksp_type != "preonly":
        outer_rtol = float(_env("ISMIP7_KSP_RTOL", KSP_RTOL_DEFAULT))
        params.update({
            f"{prefix}ksp_atol": outer_rtol * float(_env(
                "ISMIP7_CONDENSED_KSP_ATOL_FACTOR",
                CONDENSED_KSP_ATOL_FACTOR_DEFAULT,
            )),
            f"{prefix}ksp_rtol": float(_env(
                "ISMIP7_CONDENSED_KSP_RTOL", CONDENSED_KSP_RTOL_DEFAULT
            )),
            f"{prefix}ksp_max_it": int(_env(
                "ISMIP7_KSP_MAXIT", KSP_MAXIT_DEFAULT
            )),
        })
        if ksp_type.endswith("gmres"):
            params[f"{prefix}ksp_gmres_restart"] = int(_env(
                "ISMIP7_CONDENSED_KSP_RESTART", CONDENSED_KSP_RESTART_DEFAULT
            ))
    # Not a PETSc option: ISMIP7SCPC reads it from the options database.
    params[f"{prefix}near_nullspace"] = condensed_near_nullspace()
    params.update(_extra_condensed_options(prefix))
    return params


def _have_ptscotch():
    from petsc4py import PETSc
    return bool(PETSc.Sys.hasExternalPackage("ptscotch"))


def mumps_analysis():
    r"""``ISMIP7_MUMPS_ANALYSIS``: ``parallel`` (PT-Scotch distributed
    analysis, where PETSc has it) or ``sequential`` (MUMPS's own). The
    inversion selects ``sequential`` for its condensed factorizations: under
    the distributed analysis the per-solve adjoint factorizations aborted the
    2 km chains in MUMPS_LOAD_RECV_MSGS (7 of 17 links on 6 Oct, none of the
    links run sequentially)."""
    value = _env("ISMIP7_MUMPS_ANALYSIS", "parallel").strip().lower()
    if value not in ("parallel", "sequential"):
        raise ValueError("ISMIP7_MUMPS_ANALYSIS must be parallel or sequential")
    return value


def _mumps_options(prefix=""):
    opts = {
        f"{prefix}pc_type": "lu",
        f"{prefix}pc_factor_mat_solver_type": "mumps",
    }
    # Read the knob first, so a bad value fails on every build.
    analysis = mumps_analysis()
    if analysis == "parallel" and _have_ptscotch():
        # Distributed analysis with PT-Scotch nested dissection.
        opts[f"{prefix}mat_mumps_icntl_28"] = 2
        opts[f"{prefix}mat_mumps_icntl_29"] = 1
    # Without PT-Scotch the distributed analysis cannot run: MUMPS fails the
    # factorization and the condensed preconditioner returns NaN at its first
    # application (DIVERGED_NANORINF on a PETSc built without it). MUMPS's own
    # sequential analysis and ordering is then the right default.
    return opts


def diagnostic_solver_parameters(mode=None):
    r"""Return the exact PETSc options used by the mixed diagnostic solve.

    ``mode`` names a solver other than the environment's; every other knob is
    still read from the environment."""
    mode = diagnostic_solver_mode(mode)
    params = _nonlinear_options()

    if mode == "full_mumps":
        # Complete mixed-Jacobian reference.  The three numerical controls are
        # retained from the pre-scalable-solver production configuration.
        params.update({
            "mat_type": "aij",
            # Match the pre-scalable-solver production reference exactly:
            # GMRES around a full MUMPS LU (normally one Krylov iteration).
            "ksp_type": "gmres",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "mat_mumps_icntl_14": 400,
            "mat_mumps_icntl_24": 1,
            "mat_mumps_cntl_3": 1e-12,
        })
        return params

    params.update(_outer_ksp_options())

    if mode.startswith("schur_"):
        params.update({
            "mat_type": "aij",
            "pc_type": "fieldsplit",
            "pc_fieldsplit_type": "schur",
            "pc_fieldsplit_schur_fact_type": "full",
            "pc_fieldsplit_0_fields": "1,2",
            "pc_fieldsplit_1_fields": "0",
            "fieldsplit_0_ksp_type": "preonly",
            "fieldsplit_0_pc_type": "bjacobi",
            "fieldsplit_0_sub_ksp_type": "preonly",
            "fieldsplit_0_sub_pc_type": "ilu",
            "pc_fieldsplit_schur_precondition": "selfp",
            "fieldsplit_1_mat_schur_complement_ainv_type": "blockdiag",
            "fieldsplit_1_ksp_type": "preonly",
        })
        if mode == "schur_gamg":
            params.update(_gamg_options("fieldsplit_1_"))
        else:
            params.update(_mumps_options("fieldsplit_1_"))
        return params

    # SCPC sees a matrix-free mixed operator, uses Slate to invert the two DG0
    # fields exactly on each cell, and assembles only the condensed CG1 system.
    params.update({
        "mat_type": "matfree",
        "pmat_type": "matfree",
        "pc_type": "python",
        "pc_python_type": "icepack2_tools.preconditioners.ISMIP7SCPC",
        "pc_sc_eliminate_fields": "1,2",
        "condensed_field_mat_type": "aij",
        "condensed_field_ksp_type": "preonly",
    })
    if mode == "scpc_gamg":
        params.update(_condensed_gamg_options("condensed_field_"))
    else:
        params.update(_mumps_options("condensed_field_"))
    return params


def inversion_solver_mode(requested=None):
    r"""``ISMIP7_INVERSION_LINEAR_SOLVER``, canonical: the solver of the
    inversion's taped forward, its adjoint and its publishing solve."""
    if requested is None:
        requested = _env("ISMIP7_INVERSION_LINEAR_SOLVER", INVERSION_SOLVER_DEFAULT)
    mode = str(requested).strip().lower()
    if mode not in INVERSION_SOLVER_MODES:
        raise ValueError(
            "ISMIP7_INVERSION_LINEAR_SOLVER must be one of "
            f"{', '.join(INVERSION_SOLVER_MODES)}, not {requested!r}"
        )
    return mode


def inversion_state_parameters(mode=None):
    r"""PETSc options of the inversion's taped forward under ``mode``.

    ``full_mumps`` is the inversion's own reference, kept as it was before the
    knob existed: the shared SNES options around a GMRES-wrapped MUMPS LU of
    the whole mixed Jacobian, with MUMPS printing its error return (INFOG(1),
    the workspace or pivot code) so a factorisation that fails is named and
    does not reach SNES only as DIVERGED_LINEAR_SOLVE (job 1612624).  The
    scpc_* modes are the forward's own options for that mode with the outer
    Krylov tolerance ``ISMIP7_INVERSION_KSP_RTOL`` (default 1e-8,
    INVERSION_KSP_RTOL_DEFAULT; the condensed solve's absolute tolerance
    follows it as the forward's follows ISMIP7_KSP_RTOL) and the line search
    ``ISMIP7_INVERSION_SNES_LINESEARCH``, by default the shared one."""
    mode = inversion_solver_mode(mode)
    if mode != "full_mumps":
        params = diagnostic_solver_parameters(mode)
        rtol = float(_env("ISMIP7_INVERSION_KSP_RTOL", INVERSION_KSP_RTOL_DEFAULT))
        params["ksp_rtol"] = rtol
        if "condensed_field_ksp_atol" in params:
            params["condensed_field_ksp_atol"] = rtol * float(_env(
                "ISMIP7_CONDENSED_KSP_ATOL_FACTOR",
                CONDENSED_KSP_ATOL_FACTOR_DEFAULT,
            ))
        if mode == "scpc_gamg":
            # A rung's extra condensed options still come last, as in the
            # forward (_extra_condensed_options).
            params.update(_extra_condensed_options("condensed_field_"))
        params["snes_linesearch_type"] = _env(
            "ISMIP7_INVERSION_SNES_LINESEARCH", params["snes_linesearch_type"])
        return params
    params = _nonlinear_options()
    params.update({
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 400,  # working memory increase
        "mat_mumps_icntl_24": 1,  # detect null pivots
        "mat_mumps_cntl_3": 1e-12,  # null pivot threshold
        "mat_mumps_icntl_4": 1,
    })
    return params


def inversion_adjoint_parameters(params):
    r"""The adjoint's options, from the taped forward's ``params``.

    The adjoint is one linear solve, so every ``snes_*`` entry goes. That
    matters twice. An absolute tolerance sized for the forward residual lets
    the adjoint exit at iteration 0 whenever ||dJ/du|| is small, returning a
    zero adjoint and a gradient that is the prior's alone (job 10432790). And
    under ``mat_type: matfree`` tlm_adjoint hands these options to a
    LinearVariationalSolver, whose ``ksponly`` default the forward's
    ``newtonls`` would replace with a line-searched Newton solve of a linear
    problem; the assembled path ignores SNES options either way. The outer
    Krylov tolerance is the forward's relative one, and the condensed solve's
    absolute tolerance is scale-free here too: FGMRES hands the preconditioner
    unit vectors."""
    adjoint = {k: v for k, v in params.items() if not k.startswith("snes_")}
    if adjoint.get("mat_type") == "matfree":
        adjoint["snes_type"] = "ksponly"
    return adjoint


def transport_solver_parameters():
    return {
        "ksp_type": "gmres",
        "ksp_rtol": float(_env(
            "ISMIP7_TRANSPORT_KSP_RTOL", TRANSPORT_KSP_RTOL_DEFAULT
        )),
        "ksp_max_it": int(_env(
            "ISMIP7_TRANSPORT_KSP_MAXIT", TRANSPORT_KSP_MAXIT_DEFAULT
        )),
        # Firedrake's LinearVariationalSolver also checks the wrapping SNES,
        # but this makes a failed inner KSP an error at its point of origin.
        "ksp_error_if_not_converged": None,
        "pc_type": "bjacobi",
        "sub_ksp_type": "preonly",
        "sub_pc_type": "ilu",
    }


def mass_residual_tol_gt():
    r"""Absolute fail-loud tolerance for a mass-budget residual [Gt]."""
    value = float(_env(
        "ISMIP7_MASS_RESIDUAL_TOL_GT", MASS_RESIDUAL_TOL_GT_DEFAULT
    ))
    if value < 0:
        raise ValueError("ISMIP7_MASS_RESIDUAL_TOL_GT must be nonnegative")
    return value


def continuation_steps():
    return int(_env("ISMIP7_CONTINUATION_STEPS", CONTINUATION_STEPS_DEFAULT))


def rescue_max_it():
    return int(_env("ISMIP7_RESCUE_MAXIT", RESCUE_MAXIT_DEFAULT))


def rescue_enabled():
    r"""Whether a failed direct transient solve may enter the rescue ladder."""
    return _enabled("ISMIP7_RESCUE_ENABLED", RESCUE_ENABLED_DEFAULT)


def fixed_front_enabled():
    r"""Whether ice advected beyond the t=0 extent is removed as calving
    (``ISMIP7_FIXED_FRONT``). Boolean-parsed like every other switch here:
    ``0``/``false``/``off`` disable it (it used to be presence-based)."""
    return _enabled("ISMIP7_FIXED_FRONT", "0")


def subcycles():
    values = tuple(
        int(value) for value in _env(
            "ISMIP7_SUBCYCLES", SUBCYCLES_DEFAULT
        ).split(",")
    )
    if not values or any(value < 1 for value in values):
        raise ValueError("ISMIP7_SUBCYCLES must be a comma-separated list of positive integers")
    return values


def substep_settings():
    r"""Adaptive substepping of each macro step (``ISMIP7_SUBSTEP_ADAPT``), or
    None when off. With it on, the fixed ``ISMIP7_SUBCYCLES`` retry list is
    replaced by a substep count chosen from the thickness error estimate:
    ``ISMIP7_SUBSTEP_TOL`` (m), starting at ``ISMIP7_SUBSTEP_INIT`` substeps,
    at most ``ISMIP7_SUBSTEP_MAX``, halving after ``ISMIP7_SUBSTEP_QUIET``
    quiet macro steps, over cells at least ``ISMIP7_SUBSTEP_HMIN`` m thick."""
    if not _enabled("ISMIP7_SUBSTEP_ADAPT", SUBSTEP_ADAPT_DEFAULT):
        return None
    return {
        "tol": float(_env("ISMIP7_SUBSTEP_TOL", SUBSTEP_TOL_DEFAULT)),
        "m_init": int(_env("ISMIP7_SUBSTEP_INIT", SUBSTEP_INIT_DEFAULT)),
        "m_max": int(_env("ISMIP7_SUBSTEP_MAX", SUBSTEP_MAX_DEFAULT)),
        "quiet_steps": int(_env("ISMIP7_SUBSTEP_QUIET", SUBSTEP_QUIET_DEFAULT)),
        "hmin": float(_env("ISMIP7_SUBSTEP_HMIN", SUBSTEP_HMIN_DEFAULT)),
    }


def fssa_theta():
    r"""``ISMIP7_FSSA_THETA``: weight of the free-surface stabilization of
    the forward's lagged thickness-velocity coupling (icepack2_tools.fssa);
    1 (the default) makes the lagged step stable at any size, 0 leaves the
    term out of the momentum balance."""
    value = float(_env("ISMIP7_FSSA_THETA", FSSA_THETA_DEFAULT))
    if value < 0:
        raise ValueError("ISMIP7_FSSA_THETA must be nonnegative")
    return value


def forward_fssa_theta(restart_metadata=None):
    r"""The stabilization weight a forward steps with: :func:`fssa_theta`,
    except on a restart from a checkpoint stepped without the stabilization
    (``restart_metadata`` with no ``fssa_tau`` record), which keeps it off
    unless ``ISMIP7_FSSA_THETA`` is set, so a chain keeps the momentum balance
    it began with. A prepared timing or map-check state (one that records a
    ``timing_cache_role``) was never stepped and starts like a cold start."""
    if (restart_metadata is not None
            and restart_metadata.get("timing_cache_role") is None
            and restart_metadata.get("fssa_tau") is None
            and "ISMIP7_FSSA_THETA" not in os.environ):
        return 0.0
    return fssa_theta()


def fssa_reference():
    r"""``ISMIP7_FSSA_REFERENCE``: the velocity the stabilization measures the
    surface change from, ``auto`` (the default), ``start`` or ``step``
    (icepack2_tools.fssa.resolve_reference, which also checks the value)."""
    return (_env("ISMIP7_FSSA_REFERENCE", FSSA_REFERENCE_DEFAULT)
            or FSSA_REFERENCE_DEFAULT).strip().lower()


def snes_atol_scale():
    return float(_env("ISMIP7_SNES_ATOL_SCALE", SNES_ATOL_SCALE_DEFAULT))


def snes_restart_failure_atol_scale():
    return float(_env(
        "ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE",
        SNES_RESTART_FAILURE_ATOL_SCALE_DEFAULT,
    ))


def nonlinear_solver_options():
    r"""Shared SNES options (type, tolerances, line search) for any momentum solve.

    Diagnostic and inversion solver builders layer their linear options on
    top. The condensed inversion modes may override the line search through
    ``ISMIP7_INVERSION_SNES_LINESEARCH``.
    """
    return _nonlinear_options()


def direct_forward_enabled():
    r"""Whether an inversion forward first tries one Newton solve at the full
    exponents from the last converged state (``ISMIP7_DIRECT_FORWARD``)."""
    return _enabled("ISMIP7_DIRECT_FORWARD", DIRECT_FORWARD_DEFAULT)


def direct_forward_max_it():
    return int(_env("ISMIP7_DIRECT_FORWARD_MAXIT", DIRECT_FORWARD_MAXIT_DEFAULT))


def direct_forward_parameters(params):
    r"""The direct solve's options, from the taped solve's ``params`` under
    any inversion solver mode: the relative test against the trial's own
    initial residual (no ``snes_atol``: an absolute floor at the accepted
    residual stopped small control steps before the state responded), the
    live step-size exit of ``final_solve_bounds``, at most
    ``ISMIP7_DIRECT_FORWARD_MAXIT`` Newton iterations, and a lost trial once
    the residual grows ``ISMIP7_DIRECT_FORWARD_DTOL`` times."""
    out = {k: v for k, v in params.items() if k != "snes_atol"}
    out.update(final_solve_bounds())
    out["snes_max_it"] = direct_forward_max_it()
    out["snes_divergence_tolerance"] = float(
        _env("ISMIP7_DIRECT_FORWARD_DTOL", DIRECT_FORWARD_DTOL_DEFAULT))
    return out


def trial_rescue_rungs():
    r"""Rungs of the continuation ladder a failed line-search trial may climb
    (``ISMIP7_TRIAL_RESCUE_RUNGS``); 0 sends it straight to backtracking."""
    rungs = int(_env("ISMIP7_TRIAL_RESCUE_RUNGS", TRIAL_RESCUE_RUNGS_DEFAULT))
    if rungs < 0:
        raise ValueError("ISMIP7_TRIAL_RESCUE_RUNGS must be >= 0")
    return rungs


def final_solve_bounds():
    r"""``snes_stol``/``snes_max_it`` for a solve that may start at a converged
    state: the step-size exit must be live (a floor-level residual cannot pass
    the relative test) and the iteration count bounded."""
    return {
        "snes_stol": float(_env("ISMIP7_FINAL_SNES_STOL", FINAL_SNES_STOL_DEFAULT)),
        "snes_max_it": int(_env("ISMIP7_FINAL_SNES_MAXIT", FINAL_SNES_MAXIT_DEFAULT)),
    }


def final_solve_parameters(base, fnorm_ref, *, viewer=None):
    r"""Options for the solve that publishes an already-converged mixed state.

    ``fnorm_ref`` is the residual norm the preceding converged forward solve
    reached.  Restarting Newton from that state leaves ``||F||`` at the
    rounding floor, where the relative test can never pass, ``snes_stol=0``
    disables the step-size exit, and every iteration is another full
    factorisation.  ``snes_atol = snes_atol_scale() * fnorm_ref`` makes such a
    solve report ``CONVERGED_FNORM_ABS`` at iteration 0 (the transient restart
    fast path uses the same rule), while a solve from a neighbouring state
    still converges through the normal relative path.  Iteration counts are
    bounded and the converged reason is always printed: this solve decides
    what gets published and must never be silent.  ``base`` is not modified.
    """
    params = dict(base)
    params.update(final_solve_bounds())
    params.update({
        "ksp_max_it": int(_env("ISMIP7_FINAL_KSP_MAXIT", FINAL_KSP_MAXIT_DEFAULT)),
        # An LU with perturbed null pivots is not an exact inverse; GMRES must
        # error out rather than iterate to PETSc's default 10000.
        "ksp_error_if_not_converged": None,
        "snes_monitor": viewer,
        "snes_converged_reason": viewer,
    })
    try:
        fnorm_ref = float(fnorm_ref)
    except (TypeError, ValueError):
        fnorm_ref = float("nan")
    if math.isfinite(fnorm_ref) and fnorm_ref > 0.0:
        params["snes_atol"] = snes_atol_scale() * fnorm_ref
    else:
        params.pop("snes_atol", None)
    return params


def effective_solver_env():
    r"""Effective solver knobs for the committed core-run environment block."""
    return {
        "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": DIAGNOSTIC_SOLVER_DEFAULT,
        "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER_CANONICAL": diagnostic_solver_mode(),
        "ISMIP7_SNES_TYPE": SNES_TYPE_DEFAULT,
        "ISMIP7_SNES_LINESEARCH": SNES_LINESEARCH_DEFAULT,
        "ISMIP7_SNES_RTOL": SNES_RTOL_DEFAULT,
        "ISMIP7_SNES_ATOL": SNES_ATOL_DEFAULT,
        "ISMIP7_SNES_STOL": SNES_STOL_DEFAULT,
        "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
        "ISMIP7_SNES_MAXIT": SNES_MAXIT_DEFAULT,
        "ISMIP7_SNES_ATOL_SCALE": SNES_ATOL_SCALE_DEFAULT,
        "ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE": (
            SNES_RESTART_FAILURE_ATOL_SCALE_DEFAULT
        ),
        "ISMIP7_SNES_KSP_EW": SNES_KSP_EW_DEFAULT,
        "ISMIP7_SNES_MONITOR": SNES_MONITOR_DEFAULT,
        "ISMIP7_SNES_LOG": SNES_LOG_DEFAULT,
        "ISMIP7_SOLVER_VIEW": SOLVER_VIEW_DEFAULT,
        "ISMIP7_KSP_RTOL": KSP_RTOL_DEFAULT,
        "ISMIP7_KSP_MAXIT": KSP_MAXIT_DEFAULT,
        "ISMIP7_TRANSPORT_KSP_RTOL": TRANSPORT_KSP_RTOL_DEFAULT,
        "ISMIP7_TRANSPORT_KSP_MAXIT": TRANSPORT_KSP_MAXIT_DEFAULT,
        "ISMIP7_MASS_RESIDUAL_TOL_GT": MASS_RESIDUAL_TOL_GT_DEFAULT,
        "ISMIP7_CONTINUATION_STEPS": CONTINUATION_STEPS_DEFAULT,
        "ISMIP7_RESCUE_MAXIT": RESCUE_MAXIT_DEFAULT,
        "ISMIP7_RESCUE_ENABLED": RESCUE_ENABLED_DEFAULT,
        "ISMIP7_SUBCYCLES": SUBCYCLES_DEFAULT,
        "ISMIP7_SUBSTEP_ADAPT": SUBSTEP_ADAPT_DEFAULT,
        "ISMIP7_SUBSTEP_TOL": SUBSTEP_TOL_DEFAULT,
        "ISMIP7_SUBSTEP_INIT": SUBSTEP_INIT_DEFAULT,
        "ISMIP7_SUBSTEP_MAX": SUBSTEP_MAX_DEFAULT,
        "ISMIP7_SUBSTEP_QUIET": SUBSTEP_QUIET_DEFAULT,
        "ISMIP7_SUBSTEP_HMIN": SUBSTEP_HMIN_DEFAULT,
        "ISMIP7_FSSA_THETA": FSSA_THETA_DEFAULT,
        "ISMIP7_FSSA_REFERENCE": FSSA_REFERENCE_DEFAULT,
    }


def solver_provenance(mode=None):
    r"""JSON-serializable complete effective solver configuration.

    ``mode`` gives the configuration this environment would have under another
    diagnostic solver: a timing lane under one solver uses it to fingerprint
    the solver its initial-state cache was prepared with."""
    requested = requested_diagnostic_solver() if mode is None else mode
    mode = diagnostic_solver_mode(requested)
    return {
        "diagnostic_mode_requested": requested,
        "diagnostic_mode": mode,
        "diagnostic_label": diagnostic_solver_label(mode),
        "diagnostic_petsc_options": diagnostic_solver_parameters(mode),
        # Not in the cache fingerprint: it changes how a solve gets to F = 0,
        # not the state it converges to.
        "linearization_state": linearization_state(mode),
        "transport_petsc_options": transport_solver_parameters(),
        "mass_residual_tolerance_gt": mass_residual_tol_gt(),
        "continuation_steps": continuation_steps(),
        "rescue_max_it": rescue_max_it(),
        "rescue_enabled": rescue_enabled(),
        "subcycles": list(subcycles()),
        "substep_adapt": substep_settings(),
        "fssa_theta": fssa_theta(),
        "fssa_reference": fssa_reference(),
        "snes_atol_policy": {
            "initial": float(_env("ISMIP7_SNES_ATOL", SNES_ATOL_DEFAULT)),
            "post_convergence_scale": snes_atol_scale(),
            "restart_failure_residual_scale": snes_restart_failure_atol_scale(),
        },
        "monitoring": {
            "enabled": snes_monitor_enabled(),
            "view": solver_view_enabled(),
            "destination": _env("ISMIP7_SNES_LOG", SNES_LOG_DEFAULT),
        },
        "options_prefixes": {
            "diagnostic": "ismip7_diagnostic_",
            "transport": "ismip7_transport_",
        },
    }
