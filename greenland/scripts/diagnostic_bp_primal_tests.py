#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2026 David Lilien <dlilien@iu.edu>
#
"""
Perform diagnostic solves to get a basic handle on parameters and solve speed.

Usage:
    diagnostic_solves.py --mesh-dir <mesh_dir> --mesh-kind <mesh_kind>

    diagnostic_solves.py --help
"""

import argparse
import datetime
import sys
from pathlib import Path
from typing import Any

import numpy as np


DIMENSION = 3
DEFAULT_MESH_DIR = Path(__file__).resolve().parent.parent / "meshes"
DEFAULT_FIG_DIR = Path(__file__).resolve().parent.parent / "figs"
DEFAULT_LOG_DIR = Path("test_logs")
MESH_KINDS = ("detailed", "promice", "simple", "buffered")

firedrake: Any | None = None
icepack: Any | None = None
PETSc: Any | None = None
g: Any | None = None
ρ_I: Any | None = None
ρ_W: Any | None = None
extract_surface: Any | None = None


def import_solver_stack() -> None:
    """Import Firedrake and Icepack after command-line arguments are parsed."""
    global PETSc, extract_surface, firedrake, g, icepack, ρ_I, ρ_W

    import firedrake as firedrake_module
    import icepack as icepack_module
    from firedrake.petsc import PETSc as petsc_module
    from icepack.constants import gravity
    from icepack.constants import ice_density
    from icepack.constants import water_density
    from icepackaccs import extract_surface as extract_surface_function

    firedrake = firedrake_module
    icepack = icepack_module
    PETSc = petsc_module
    g = gravity
    ρ_I = ice_density
    ρ_W = water_density
    extract_surface = extract_surface_function


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Perform diagnostic solves to get a basic handle on parameters and "
            "solve speed."
        )
    )
    parser.add_argument(
        "--mesh-dir",
        type=Path,
        default=DEFAULT_MESH_DIR,
        help="Directory containing initialized h5 files.",
    )
    parser.add_argument(
        "--mesh-kind",
        choices=MESH_KINDS,
        default="promice",
        help="Mesh type to use: detailed, promice, simple, or buffered.",
    )
    parser.add_argument("--twostep", action="store_true", help="Use two-step solver.")
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=DEFAULT_FIG_DIR,
        help="Directory for diagnostic figures.",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Directory for timing logs.",
    )
    return parser.parse_args()


def dirichlet_ids(mesh_kind: str) -> list[int]:
    """Return Dirichlet boundary IDs for a mesh kind."""
    if mesh_kind == "promice":
        return np.arange(2, 3153).tolist()
    if mesh_kind == "detailed":
        return np.arange(1, 2289).tolist()
    return np.arange(1, 64).tolist()


def solver_options(dids: list[int]) -> tuple[dict, dict]:
    """Build line-search and trust-region solver option dictionaries."""
    opts0 = {
        "dirichlet_ids": dids,
        "diagnostic_solver_type": "petsc",
        "diagnostic_solver_parameters": {
            "snes_type": "newtonls",
            "snes_line_search_type": "bt",
            "snes_linesearch_order": 2,
            "snes_linesearch_max_it": 2500,
            "snes_linesearch_damping": 0.05,
            "snes_max_it": 5000,
            "snes_stol": 1.0e-6,
            "snes_rtol": 1.0e-5,
            "ksp_type": "bcgs",
            "ksp_max_it": 2500,
            "ksp_rtol": 1.0e-8,
            "ksp_atol": 1.0e-4,
            "ksp_converged_maxits": True,
            "pc_type": "bjacobi",
            "pc_factor_mat_solver_type": "mumps",
            "pc_factor_shift_amount": 1.0e-10,
            "snes_monitor": None,
            "snes_linesearch_monitor": None,
        },
    }
    opts1 = {
        "dirichlet_ids": dids,
        "diagnostic_solver_type": "petsc",
        "diagnostic_solver_parameters": {
            "snes_type": "newtontr",
            "snes_tr_delta0": 1.0e5,
            "snes_tr_fallback_type": "dogleg",
            "snes_max_it": 5000,
            "snes_stol": 1.0e-8,
            "snes_rtol": 1.0e-8,
            "ksp_type": "bcgs",
            "ksp_max_it": 100000,
            "ksp_rtol": 1.0e-16,
            "ksp_atol": 1.0e-16,
            "pc_type": "bjacobi",
            "pc_hypre_type": "boomeramg",
            "pc_factor_mat_solver_type": "mumps",
            "pc_factor_shift_amount": 1.0e-10,
            "snes_monitor": None,
            "snes_linesearch_monitor": None,
        },
    }
    return opts0, opts1


def resolution_levels(comm_size: int) -> np.ndarray:
    """Choose mesh resolutions based on MPI size."""
    if comm_size == 1:
        return np.array([250 * 2**i for i in range(4, 8)])
    if comm_size < 3:
        return np.array([250 * 2**i for i in range(3, 8)])
    if comm_size < 13:
        return np.array([250 * 2**i for i in range(2, 6)])
    return np.array([250 * 2**i for i in range(0, 4)])


def checkpoint_path(mesh_dir: Path, mesh_kind: str, lc: int) -> Path:
    """Build the checkpoint path matching mesh_greenland.py names."""
    return mesh_dir / f"greenland_{mesh_kind}_{10 * lc}_{lc}.h5"


def save_input_figure(
    fig_dir, stem, H, u_obs_mag, C, C_var, firedrake, extract_surface
):
    """Save diagnostic input fields."""
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(14, 14))
    colors = firedrake.tripcolor(extract_surface(H), vmin=0, vmax=4000, axes=ax[0][0])
    plt.colorbar(colors, extend="max", label="Ice thickness [m]")
    colors = firedrake.tripcolor(
        extract_surface(u_obs_mag),
        cmap="turbo",
        norm=mcolors.LogNorm(vmin=5, vmax=2000),
        axes=ax[0][1],
    )
    plt.colorbar(colors, extend="max", label="Velocity (m/yr)")
    colors = firedrake.tripcolor(
        extract_surface(C), vmin=0, vmax=0.5, cmap="plasma", axes=ax[1][0]
    )
    plt.colorbar(colors, extend="max", label="Constant C")
    colors = firedrake.tripcolor(
        extract_surface(C_var), vmin=0, vmax=0.5, cmap="plasma", axes=ax[1][1]
    )
    plt.colorbar(colors, extend="max", label="Variable C")

    for axes_row in ax:
        for axes in axes_row:
            axes.axis("equal")
    fig.savefig(fig_dir / f"{stem}_inputs.png", dpi=300)
    plt.close(fig)


def save_velocity_figure(fig_dir, stem, u_obs_mag, u0, firedrake, extract_surface):
    """Save observed and modeled velocity comparison."""
    import matplotlib.pyplot as plt

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    colors = firedrake.tripcolor(
        extract_surface(u_obs_mag), vmin=0, vmax=1000, axes=ax[0]
    )
    plt.colorbar(colors, extend="max", label="Observed flow speed [m/yr]")
    colors = firedrake.tripcolor(
        extract_surface(firedrake.sqrt(u0[0] ** 2.0 + u0[1] ** 2.0)),
        vmin=0,
        vmax=1000,
        axes=ax[1],
    )
    plt.colorbar(colors, extend="max", label="Modeled flow speed [m/yr]")
    ax[0].axis("equal")
    ax[1].axis("equal")
    fig.savefig(fig_dir / f"{stem}_veltest.png", dpi=300)
    plt.close(fig)


def single_test(
    checkpoint_fn: Path,
    lc: int,
    twostep: bool,
    opts0: dict,
    opts1: dict,
    try_damping: np.ndarray,
    try_d0: np.ndarray,
    fig_dir: Path,
) -> tuple[float, float, float]:
    """Perform a single diagnostic solve."""

    with firedrake.CheckpointFile(str(checkpoint_fn), "r") as chk:
        mesh = chk.load_mesh(name=f"greenland_{DIMENSION}d")
        u_obs = chk.load_function(mesh, name=f"u_{DIMENSION}d")
        H_in = chk.load_function(mesh, name=f"H_{DIMENSION}d")
        b_in = chk.load_function(mesh, name=f"b_{DIMENSION}d")

    area = firedrake.Constant(
        firedrake.assemble(firedrake.Constant(1.0) * firedrake.ds_t(mesh))
    )
    Q = firedrake.FunctionSpace(mesh, "CG", 2, vfamily="R", vdegree=0)
    V = firedrake.VectorFunctionSpace(mesh, "CG", 1, dim=2, vfamily="GL", vdegree=4)

    h0 = firedrake.Function(Q).interpolate(b_in + H_in)
    H = firedrake.Function(Q).interpolate(H_in)
    h = h0.copy(deepcopy=True)

    # Smooth the surface elevation to reduce bumps etc.
    α = firedrake.Constant(2e3)
    J = (
        0.5
        * ((h - h0) ** 2 + α**2 * firedrake.inner(firedrake.grad(h), firedrake.grad(h)))
        * firedrake.dx
    )
    F = firedrake.derivative(J, h)
    firedrake.solve(F == 0, h)

    τ_D = firedrake.Function(Q).interpolate(
        ρ_I
        * g
        * H
        * firedrake.sqrt(firedrake.grad(H)[0] ** 2.0 + firedrake.grad(H)[1] ** 2.0)
    )
    u_obs_mag = firedrake.Function(Q).interpolate(
        firedrake.max_value(
            firedrake.sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), firedrake.Constant(1.0e-3)
        )
    )

    C_var = firedrake.Function(Q).interpolate(firedrake.sqrt(τ_D / u_obs_mag / 2.0))
    C_m = firedrake.assemble(C_var * firedrake.dx) / area
    C = firedrake.Function(Q).interpolate(firedrake.Constant(C_m))

    if firedrake.COMM_WORLD.size == 1:
        save_input_figure(
            fig_dir,
            checkpoint_fn.stem,
            H,
            u_obs_mag,
            C,
            C_var,
            firedrake,
            extract_surface,
        )

    T = firedrake.Constant(260)
    A0 = icepack.rate_factor(T)

    def linear_pos_friction(**kwargs):
        """Evaluate flotation-adjusted linear sliding friction."""
        p_W = ρ_W * g * firedrake.max_value(0, kwargs["thickness"] - kwargs["surface"])
        p_I = ρ_I * g * kwargs["thickness"]
        ϕ = 1 - p_W / p_I
        return (
            ϕ
            * kwargs["C"] ** 2.0
            * firedrake.inner(kwargs["velocity"], kwargs["velocity"])
        )

    model = icepack.models.HybridModel(friction=linear_pos_friction)
    PETSc.Sys.Print(f"Beginning initial velocity solve {lc:d}")

    elapsed = 0.0
    damping_used = np.nan
    delta0_used = np.nan
    u0 = firedrake.Function(V).interpolate(firedrake.Constant(1.0e-1) * u_obs)

    if twostep:
        for damping in try_damping:
            opts0["diagnostic_solver_parameters"]["snes_linesearch_damping"] = damping
            solver0 = icepack.solvers.FlowSolver(model, **opts0)
            PETSc.Sys.Print(f"Trying LS solve with damping = {damping:e}")
            try:
                start_time = datetime.datetime.now()
                u0 = solver0.diagnostic_solve(
                    velocity=u0, fluidity=A0, C=C, surface=h, thickness=H
                )
                end_time = datetime.datetime.now()
                elapsed += (end_time - start_time).total_seconds()
                damping_used = damping
                break
            except firedrake.exceptions.ConvergenceError:
                pass
        else:
            PETSc.Sys.Print("LS Failed")

    PETSc.Sys.Print(f"Starting TR solve on {firedrake.COMM_WORLD.size:d} processes")
    for d0 in try_d0:
        opts1["diagnostic_solver_parameters"]["snes_tr_delta0"] = d0
        solver1 = icepack.solvers.FlowSolver(model, **opts1)
        PETSc.Sys.Print(f"Trying TR solve with d0 = {d0:e}")
        try:
            start_time = datetime.datetime.now()
            u0 = solver1.diagnostic_solve(
                velocity=u0, fluidity=A0, C=C, surface=h, thickness=H
            )
            end_time = datetime.datetime.now()
            elapsed += (end_time - start_time).total_seconds()
            delta0_used = d0
            break
        except firedrake.exceptions.ConvergenceError:
            pass
    else:
        PETSc.Sys.Print("No solution found!")

    PETSc.Sys.Print(
        f"Solve on {firedrake.COMM_WORLD.size:d} processes took {int(elapsed):d} s"
    )

    if firedrake.COMM_WORLD.size == 1:
        save_velocity_figure(
            fig_dir, checkpoint_fn.stem, u_obs_mag, u0, firedrake, extract_surface
        )

    return elapsed, damping_used, delta0_used


def run_tests(
    mesh_dir: Path, mesh_kind: str, twostep: bool, fig_dir: Path, log_dir: Path
) -> None:
    """Run diagnostic solves over the selected resolution levels."""
    opts0, opts1 = solver_options(dirichlet_ids(mesh_kind))
    lcs = resolution_levels(firedrake.COMM_WORLD.size)
    log_dir.mkdir(parents=True, exist_ok=True)
    timing_fn = log_dir / f"timings_{mesh_kind}_n{firedrake.COMM_WORLD.size:d}.txt"
    try_d0 = np.round(np.array([10 ** -(i / 2) * 1e7 for i in range(11)]), decimals=5)
    try_damping = np.round(
        np.array([0.5 * 10 ** (-i / 2) for i in range(11)]), decimals=5
    )

    with timing_fn.open("w") as fout:
        fout.write("nproc, lc (m), time (s), LS damping, TR delta0\n")

    for lc in lcs:
        checkpoint_fn = checkpoint_path(mesh_dir, mesh_kind, int(lc))
        if not checkpoint_fn.exists():
            print(f"Skipping missing checkpoint: {checkpoint_fn}")
            continue

        elapsed, damping, delta0 = single_test(
            checkpoint_fn,
            int(lc),
            twostep,
            opts0,
            opts1,
            try_damping,
            try_d0,
            fig_dir,
        )
        with timing_fn.open("a") as fout:
            fout.write(
                f"{firedrake.COMM_WORLD.size:d}, {int(lc):d}, {elapsed:f}, "
                f"{damping:e}, {delta0:e}\n"
            )


def main() -> None:
    """Run diagnostic solves from command-line arguments."""
    args = parse_args()
    sys.argv = [sys.argv[0]]
    import_solver_stack()
    run_tests(
        args.mesh_dir.expanduser(),
        args.mesh_kind,
        args.twostep,
        args.fig_dir.expanduser(),
        args.log_dir.expanduser(),
    )


if __name__ == "__main__":
    main()
