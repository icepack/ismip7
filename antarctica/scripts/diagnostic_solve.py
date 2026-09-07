#!/usr/bin/env python3
"""
Diagnostic velocity solve for Antarctica using icepack2.

Loads observational data onto the mesh, sets up an initial friction
estimate, and solves the SSA momentum balance using icepack2's dual
formulation with continuation from linear to Glen's flow law.

The grounding zone is handled seamlessly via height-above-flotation:
friction vanishes where ice is floating, and the surface elevation
transitions smoothly between grounded (s = b + h) and floating
(s = (1 - rho_I/rho_W) * h).

Usage:
    python scripts/diagnostic_solve.py
"""

import os
import sys
import glob
import numpy as np

import firedrake
from firedrake import (
    Constant,
    Function,
    max_value,
    sqrt,
    inner,
    grad,
    derivative,
    dx,
)
from firedrake.petsc import PETSc

import rasterio
import icepack
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
    glen_flow_law,
)
import colorcet as cc
import matplotlib.pyplot as plt

from mesh_naming import get_buffer_m, mesh_filename

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")

sys.path.insert(0, os.path.dirname(_ROOT))
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.runconfig import lc as _lc, lc_coarse as _lc_coarse

lc = _lc()
lc_coarse = _lc_coarse()
buffer_m = get_buffer_m()


def find_file(directory, pattern):
    matches = glob.glob(os.path.join(directory, pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching '{pattern}' in {directory}")
    return matches[0]


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    # ── Load Mesh ──────────────────────────────────────────────────────
    mesh_fn = os.environ.get("ISMIP7_MESH", mesh_filename(lc_coarse, lc, buffer_m))
    PETSc.Sys.Print(f"Loading mesh: {mesh_fn}")
    mesh = firedrake.Mesh(mesh_fn)
    PETSc.Sys.Print(f"  {mesh.num_vertices()} vertices, {mesh.num_cells()} cells")

    # Boundary classification, hard-checked against this mesh (a sidecar built
    # for another mesh leaves most of the front without terminus back-pressure).
    bnd_ids, calving_ids, _ = load_boundary_ids(mesh, MESH_DIR, mesh_hint=mesh_fn)
    other_ids = tuple(bnd_ids["other"])
    PETSc.Sys.Print(f"  {len(calving_ids)} calving + {len(other_ids)} other boundaries")

    # ── Function Spaces ────────────────────────────────────────────────
    Q = firedrake.FunctionSpace(mesh, "CG", 1)
    V = firedrake.VectorFunctionSpace(mesh, "CG", 1)

    dg0 = firedrake.FiniteElement("DG", "triangle", 0)
    Sigma = firedrake.TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = firedrake.VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    # ── Load Observational Data ────────────────────────────────────────
    PETSc.Sys.Print("Loading BedMachine data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    bed_src = rasterio.open(f"netcdf:{bm_fn}:bed")
    thick_src = rasterio.open(f"netcdf:{bm_fn}:thickness")

    b = icepack.interpolate(bed_src, Q)
    H_raw = icepack.interpolate(thick_src, Q)
    H = Function(Q).interpolate(max_value(H_raw, Constant(10.0)))

    # Surface elevation: seamless grounded/floating transition
    # s = max(b + H, (1 - rho_I/rho_W) * H)
    rho_ratio = Constant(917.0 / 1024.0)
    s_expr = max_value(b + H, (Constant(1.0) - rho_ratio) * H)
    s = Function(Q).interpolate(s_expr)

    # Smooth surface to reduce numerical noise (Tikhonov regularization)
    PETSc.Sys.Print("Smoothing surface elevation...")
    s0 = s.copy(deepcopy=True)
    alpha = Constant(2e3)
    J_smooth = 0.5 * ((s - s0) ** 2 + alpha**2 * inner(grad(s), grad(s))) * dx
    F_smooth = derivative(J_smooth, s)
    firedrake.solve(F_smooth == 0, s)

    PETSc.Sys.Print("Loading velocity data...")
    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    vx_src = rasterio.open(f"netcdf:{vel_fn}:VX")
    vy_src = rasterio.open(f"netcdf:{vel_fn}:VY")
    u_obs = icepack.interpolate((vx_src, vy_src), V, fillvalue=0.0)

    # ── Friction Estimate ──────────────────────────────────────────────
    PETSc.Sys.Print("Computing initial friction estimate...")

    u_speed = Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    )

    # Effective pressure ratio (height above flotation)
    # phi = 1 where fully grounded, -> 0 at grounding line, 0 where floating
    p_W = rho_W * g * max_value(Constant(0.0), -b)
    p_I = rho_I * g * H
    phi = Function(Q).interpolate(max_value(Constant(1.0) - p_W / p_I, Constant(0.01)))

    # ── icepack2 Momentum Balance ──────────────────────────────────────
    PETSc.Sys.Print("Setting up icepack2 momentum balance...")

    z = Function(Z)
    # Initial velocity guess: fraction of observed
    z.sub(0).interpolate(Constant(0.1) * u_obs)

    # Rheology scales
    n = Constant(glen_flow_law)
    m = Constant(glen_flow_law)
    tau_c = Constant(0.1)  # yield stress scale (MPa)
    eps_c = Constant(0.01)  # strain rate scale (1/yr)
    u_c = Constant(100.0)  # velocity scale (m/yr)

    u, M, tau = firedrake.split(z)

    # Sliding coefficient: K = u_c / (phi * tau_c)^m
    # Where floating (phi -> 0): K -> large (free sliding)
    # Where grounded: K is finite (friction active)
    K = u_c / (phi * tau_c) ** n

    fields = {
        "velocity": u,
        "membrane_stress": M,
        "basal_stress": tau,
        "thickness": H,
        "surface": s,
    }
    rheology = {
        "flow_law_exponent": n,
        "flow_law_coefficient": eps_c / tau_c**n,
        "sliding_exponent": m,
        "sliding_coefficient": K,
    }

    # Build the action functional
    L = (
        model.minimization.viscous_power(**fields, **rheology)
        + model.minimization.friction_power(**fields, **rheology)
        + model.minimization.momentum_balance(**fields)
        + model.minimization.calving_terminus(**fields, outflow_ids=calving_ids)
    )
    F = derivative(L, z)

    # No Dirichlet BCs needed: icepack2 is well-posed at h=0, so velocity
    # naturally vanishes where thickness goes to zero at land boundaries.
    # Calving fronts are handled by calving_terminus (ocean back-pressure).

    # Solver parameters
    sparams = {
        "snes_type": "newtonls",
        "snes_max_it": 200,
        "snes_linesearch_type": "nleqerr",
        "snes_monitor": None,
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    pparams = {"form_compiler_parameters": {"quadrature_degree": 4}}

    problem = firedrake.NonlinearVariationalProblem(F, z, **pparams)
    solver = firedrake.NonlinearVariationalSolver(problem, solver_parameters=sparams)

    # Continuation: step exponent from linear (1) to Glen's law (3)
    PETSc.Sys.Print("Solving with continuation...")
    for exponent in np.linspace(1.0, glen_flow_law, 5):
        n.assign(exponent)
        m.assign(exponent)
        PETSc.Sys.Print(f"  n = {exponent:.2f}")
        solver.solve()

    PETSc.Sys.Print("Solve complete!")

    # ── Plot Results ───────────────────────────────────────────────────
    u_sol = z.subfunctions[0]
    u_sol_mag = Function(Q).interpolate(sqrt(u_sol[0] ** 2 + u_sol[1] ** 2))

    # Manual positioning: [left, bottom, width, height] in figure coords
    # Three panels with generous spacing, colorbars half-height
    pw = 0.25  # panel width
    ph = 0.75  # panel height
    bot = 0.10  # bottom margin
    gap1 = 0.03  # gap between obs and mod
    gap2 = 0.10  # bigger gap before diff panel
    cb_w = 0.012  # colorbar width
    cb_gap = 0.015  # gap between panel and its colorbar
    cb_h = ph * 0.5  # colorbar is half panel height
    cb_bot = bot + (ph - cb_h) / 2  # vertically centered

    x0 = 0.04
    x1 = x0 + pw + gap1
    x_cb1 = x1 + pw + cb_gap
    x2 = x_cb1 + cb_w + gap2
    x_cb2 = x2 + pw + cb_gap

    fig = plt.figure(figsize=(22, 7))
    ax0 = fig.add_axes([x0, bot, pw, ph])
    ax1 = fig.add_axes([x1, bot, pw, ph])
    cax_speed = fig.add_axes([x_cb1, cb_bot, cb_w, cb_h])
    ax2 = fig.add_axes([x2, bot, pw, ph])
    cax_diff = fig.add_axes([x_cb2, cb_bot, cb_w, cb_h])
    axes = [ax0, ax1, ax2]

    # Consistent spatial extent
    coords = mesh.coordinates.dat.data_ro
    xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
    ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
    pad_xy = 0.02 * max(xmax - xmin, ymax - ymin)
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlim(xmin - pad_xy, xmax + pad_xy)
        ax.set_ylim(ymin - pad_xy, ymax + pad_xy)

    # Turn off y-axis labels on panels 2 & 3 (shared y range)
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])

    # Tripcolor (underneath)
    firedrake.tripcolor(u_speed, vmin=0, vmax=1000, axes=ax0, cmap=cc.cm.CET_L19)
    ax0.set_title("Observed speed")

    colors_speed = firedrake.tripcolor(
        u_sol_mag, vmin=0, vmax=1000, axes=ax1, cmap=cc.cm.CET_L19
    )
    ax1.set_title("Modeled speed")

    diff = Function(Q).interpolate(u_speed - u_sol_mag)
    colors_diff = firedrake.tripcolor(
        diff, vmin=-500, vmax=500, axes=ax2, cmap=cc.cm.CET_CBTD1
    )
    ax2.set_title("Observed - Modeled")

    # Mesh overlay on top
    for ax in axes:
        firedrake.triplot(
            mesh,
            axes=ax,
            interior_kw={"linewidth": 0.1, "alpha": 0.3, "color": "k"},
            boundary_kw={"linewidth": 1.0, "color": "k"},
        )

    # Colorbars — half height, vertically centered
    fig.colorbar(colors_speed, cax=cax_speed, label="m/yr")
    fig.colorbar(colors_diff, cax=cax_diff, label="m/yr")
    out_fn = os.path.join(FIG_DIR, f"diagnostic_solve_{lc}.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    PETSc.Sys.Print(f"Saved: {out_fn}")

    # Save checkpoint
    chk_fn = os.path.join(MESH_DIR, f"diagnostic_{lc}.h5")
    with firedrake.CheckpointFile(chk_fn, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(u_sol, name="u")
        chk.save_function(H, name="H")
        chk.save_function(b, name="b")
        chk.save_function(s, name="s")
    PETSc.Sys.Print(f"Saved: {chk_fn}")


if __name__ == "__main__":
    main()
