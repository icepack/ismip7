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

import os


# A forward driver invoked outside the managed launchers must fall back to the
# established full-Jacobian reference, never to an unqualified development PC.
# The timing Makefile explicitly exports scpc_gamg for scalable-solver work.
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

SNES_TYPE_DEFAULT = "newtonls"
SNES_LINESEARCH_DEFAULT = "nleqerr"
SNES_RTOL_DEFAULT = "1e-8"
SNES_ATOL_DEFAULT = "1e-50"
SNES_STOL_DEFAULT = "0"
SNES_DIVERGENCE_TOL_DEFAULT = "-1"
SNES_MAXIT_DEFAULT = "200"
SNES_ATOL_SCALE_DEFAULT = "100"
SNES_RESTART_FAILURE_ATOL_SCALE_DEFAULT = "1e-6"
SNES_KSP_EW_DEFAULT = "0"
SNES_MONITOR_DEFAULT = "0"
SNES_LOG_DEFAULT = "stdout"
SOLVER_VIEW_DEFAULT = "0"

KSP_RTOL_DEFAULT = "1e-6"
KSP_MAXIT_DEFAULT = "1000"
TRANSPORT_KSP_RTOL_DEFAULT = "1e-10"
TRANSPORT_KSP_MAXIT_DEFAULT = "500"
MASS_RESIDUAL_TOL_GT_DEFAULT = "5e-5"
CONTINUATION_STEPS_DEFAULT = "8"
RESCUE_MAXIT_DEFAULT = "600"
SUBCYCLES_DEFAULT = "1,4,16"
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


def diagnostic_solver_mode():
    requested = requested_diagnostic_solver()
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


def _mumps_options(prefix=""):
    return {
        f"{prefix}pc_type": "lu",
        f"{prefix}pc_factor_mat_solver_type": "mumps",
        # Distributed analysis with PT-Scotch nested dissection.
        f"{prefix}mat_mumps_icntl_28": 2,
        f"{prefix}mat_mumps_icntl_29": 1,
    }


def diagnostic_solver_parameters():
    r"""Return the exact PETSc options used by the mixed diagnostic solve."""
    mode = diagnostic_solver_mode()
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
        params.update(_gamg_options("condensed_field_"))
    else:
        params.update(_mumps_options("condensed_field_"))
    return params


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


def subcycles():
    values = tuple(
        int(value) for value in _env(
            "ISMIP7_SUBCYCLES", SUBCYCLES_DEFAULT
        ).split(",")
    )
    if not values or any(value < 1 for value in values):
        raise ValueError("ISMIP7_SUBCYCLES must be a comma-separated list of positive integers")
    return values


def snes_atol_scale():
    return float(_env("ISMIP7_SNES_ATOL_SCALE", SNES_ATOL_SCALE_DEFAULT))


def snes_restart_failure_atol_scale():
    return float(_env(
        "ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE",
        SNES_RESTART_FAILURE_ATOL_SCALE_DEFAULT,
    ))


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
    }


def solver_provenance():
    r"""JSON-serializable complete effective solver configuration."""
    return {
        "diagnostic_mode_requested": requested_diagnostic_solver(),
        "diagnostic_mode": diagnostic_solver_mode(),
        "diagnostic_label": diagnostic_solver_label(),
        "diagnostic_petsc_options": diagnostic_solver_parameters(),
        "transport_petsc_options": transport_solver_parameters(),
        "mass_residual_tolerance_gt": mass_residual_tol_gt(),
        "continuation_steps": continuation_steps(),
        "rescue_max_it": rescue_max_it(),
        "rescue_enabled": rescue_enabled(),
        "subcycles": list(subcycles()),
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
