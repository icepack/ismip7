#!/usr/bin/env python3
"""
Grounding line retreat sensitivity using adjoint + eikonal distance.

Computes the adjoint sensitivity dJ/dH (J = total GL flux) and projects
onto a smooth GL retreat perturbation defined by the eikonal distance
from the grounding line.

The perturbation represents uniform thinning near the GL:
  ΔH(x) = -ε · φ(d_GL(x))
where d_GL is the eikonal distance from the GL and φ is a smooth
profile concentrating the thinning near the GL.

This avoids the H→0 singularity of calving front experiments: the
perturbation is smooth, the GL stays well-resolved, and the adjoint
linearization is valid.

The shape derivative gives dJ/dε — the GL flux response to unit
thinning concentrated at the grounding line.

Usage:
    python scripts/gl_sensitivity.py
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import firedrake as fd
from firedrake import (
    Constant,
    Function,
    max_value,
    sqrt,
    inner,
    grad,
    derivative,
    dx,
    split,
    assemble,
    Mesh,
    FunctionSpace,
    VectorFunctionSpace,
    TensorFunctionSpace,
    FiniteElement,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
)
from tlm_adjoint.firedrake import (
    reset_manager,
    start_manager,
    stop_manager,
    clear_caches,
    compute_gradient,
    Functional,
    EquationSolver,
)
from firedrake.petsc import PETSc

import os, sys, glob

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PROJECT = os.path.dirname(_ROOT)
sys.path.insert(0, _PROJECT)
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.mpi_stats import global_max, global_range, global_sum
from icepack2_tools.eikonal import identify_grounding_line, solve_eikonal_distance
import rasterio, icepack
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as grav,
    glen_flow_law,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")
RESULTS_DIR = os.path.join(_ROOT, "results")
SHAPEFILE = os.path.expanduser("~/data/shapefiles/IceShelf_Antarctica_v02.shp")

from mesh_naming import get_buffer_m, mesh_filename
from icepack2_tools.runconfig import lc as _lc, lc_coarse as _lc_coarse

lc = _lc()
lc_coarse = _lc_coarse()
buffer_m = get_buffer_m()

# GL thinning profile widths (km) — how far inland the thinning extends
THINNING_WIDTHS_KM = [5.0, 10.0, 20.0]
# Shelves to validate with forward FD
VALIDATE_SHELVES = ["Ross", "Pine_Island", "Getz", "Thwaites"]
MIN_SHELF_AREA = 5000e6


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def build_shelf_masks(shapefile, coords, min_area):
    """Build per-shelf boolean masks from the MEaSUREs shapefile."""
    gdf = gpd.read_file(shapefile)
    gdf = gdf[gdf.geometry.area > min_area].copy()
    gdf = gdf.sort_values("NAME").reset_index(drop=True)
    merge_groups = {
        "Ross": ["Ross_West", "Ross_East"],
        "Getz": ["Getz", "Getz_1", "Getz_2"],
        "Abbot": [
            "Abbot",
            "Abbot_1",
            "Abbot_2",
            "Abbot_3",
            "Abbot_4",
            "Abbot_5",
            "Abbot_6",
        ],
    }
    for merged_name, parts in merge_groups.items():
        part_mask = gdf["NAME"].isin(parts)
        if part_mask.sum() > 1:
            merged_geom = gdf.loc[part_mask, "geometry"].unary_union
            merged_row = gpd.GeoDataFrame(
                {"NAME": [merged_name], "geometry": [merged_geom]}, crs=gdf.crs
            )
            gdf = gpd.GeoDataFrame(
                pd.concat([gdf[~part_mask], merged_row], ignore_index=True), crs=gdf.crs
            )
    points = gpd.GeoSeries([Point(x, y) for x, y in coords], crs=gdf.crs)
    masks = {}
    for _, row in gdf.iterrows():
        name = row["NAME"]
        buffered = row.geometry.buffer(500)
        mask = np.array([buffered.contains(p) for p in points])
        if mask.sum() > 0:
            masks[name] = mask
    return masks


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    mesh_fn = os.environ.get("ISMIP7_MESH", mesh_filename(lc_coarse, lc, buffer_m))
    PETSc.Sys.Print(f"Loading mesh: {mesh_fn}")
    mesh = Mesh(mesh_fn)
    PETSc.Sys.Print(f"  {mesh.num_vertices()} vertices, {mesh.num_cells()} cells")

    # Sidecar resolved (per-mesh preferred, parametric fallback) and
    # HARD-CHECKED against this mesh: an id absent from the mesh makes
    # ds(id) integrate to zero, i.e. silently wrong physics with no crash.
    use_calving_terminus = os.environ.get("ISMIP7_NO_CALVING_TERMINUS") is None
    bnd_ids, calving_ids, bndids_fn = load_boundary_ids(
        mesh, MESH_DIR, mesh_hint=mesh_fn,
        print_coverage=use_calving_terminus,
    )

    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    # ── Data ──
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:bed"), Q)
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "10.0"))
    H = Function(Q).interpolate(
        max_value(
            icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q),
            Constant(h_clamp),
        )
    )
    PETSc.Sys.Print(f"  H clamp: {h_clamp} m")

    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q).interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))

    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V,
        fillvalue=0.0,
    )

    A0 = Function(Q).interpolate(Constant(icepack.rate_factor(Constant(260.0))))
    n_glen = Constant(glen_flow_law)
    tau_c = Constant(0.1)
    u_c = Constant(100.0)
    phi_eff = Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rho_W
            * grav
            * max_value(Constant(0.0), -b)
            / (rho_I * grav * max_value(H, Constant(1.0))),
            Constant(0.01),
        )
    )
    K_base = u_c / (phi_eff * tau_c) ** n_glen

    sparams = {
        "snes_type": "newtonls",
        "snes_max_it": 200,
        "snes_linesearch_type": "nleqerr",
        "snes_divergence_tolerance": -1,
        "snes_stol": 0.0,
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 200,
        "mat_mumps_icntl_24": 1,
        "mat_mumps_cntl_3": 1e-12,
    }
    fc_params = {"quadrature_degree": 4}

    # ── Load MAP ──
    inv_fn = os.path.join(MESH_DIR, f"inversion_icepack2_{lc}.h5")
    PETSc.Sys.Print(f"Loading MAP: {inv_fn}")
    with fd.CheckpointFile(inv_fn, "r") as chk:
        chk_mesh = chk.load_mesh()
        theta_chk = chk.load_function(chk_mesh, name="log_friction")
        phi_chk = chk.load_function(chk_mesh, name="log_fluidity")
    theta_map = Function(Q, name="theta_map")
    theta_map.dat.data[:] = theta_chk.dat.data_ro
    phi_map = Function(Q, name="phi_map")
    phi_map.dat.data[:] = phi_chk.dat.data_ro

    K_map_expr = K_base * fd.exp(-n_glen * theta_map)
    A_map_expr = A0 * fd.exp(phi_map)

    # ── Shelf masks ──
    PETSc.Sys.Print("Building shelf masks...")
    coords = mesh.coordinates.dat.data_ro
    shelf_mask_arrays = build_shelf_masks(SHAPEFILE, coords, MIN_SHELF_AREA)
    shelf_mask_fns = {}
    for name, arr in shelf_mask_arrays.items():
        fn = Function(Q)
        fn.dat.data[:] = arr.astype(float)
        shelf_mask_fns[name] = fn

    # ── Grounding line mask and eikonal distance ──
    PETSc.Sys.Print("\nComputing GL eikonal distance...")
    gl_mask = identify_grounding_line(H, b)
    d_gl = solve_eikonal_distance(mesh, gl_mask)
    d_gl_arr = d_gl.dat.data_ro
    H_arr = H.dat.data_ro

    # Grounded indicator (smooth Heaviside)
    # Grounding indicator via smoothed height-above-flotation
    from firedrake import conditional

    s_float = Function(Q).interpolate(
        b + (rho_W / rho_I) * max_value(-b, Constant(0.0))
    )
    haf = Function(Q).interpolate(s - s_float)
    haf0 = haf.copy(deepcopy=True)
    fd.solve(
        derivative(
            0.5
            * ((haf - haf0) ** 2 + Constant(100.0) ** 2 * inner(grad(haf), grad(haf)))
            * dx,
            haf,
        )
        == 0,
        haf,
    )
    grounded_indicator = Function(Q).interpolate(
        conditional(haf > 0, Constant(1.0), Constant(0.0))
    )
    grounded_arr = haf.dat.data_ro > 0

    # ══════════════════════════════════════════════════════════════════════
    # Forward reference solve
    # ══════════════════════════════════════════════════════════════════════
    PETSc.Sys.Print("\n=== Forward reference solve ===")
    h_ctrl = H.copy(deepcopy=True)
    z_adj = Function(Z)
    z_adj.sub(0).interpolate(Constant(0.1) * u_obs)

    s_ctrl = Function(Q).interpolate(
        max_value(b + h_ctrl, (Constant(1.0) - rho_ratio) * h_ctrl)
    )

    u_s, M_s, tau_s = split(z_adj)
    flds = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": h_ctrl,
        "surface": s_ctrl,
    }
    rh = {
        "flow_law_exponent": n_glen,
        "flow_law_coefficient": A_map_expr,
        "sliding_exponent": n_glen,
        "sliding_coefficient": K_map_expr,
    }
    L = (
        model.minimization.viscous_power(**flds, **rh)
        + model.minimization.friction_power(**flds, **rh)
        + model.minimization.momentum_balance(**flds)
    )
    if use_calving_terminus:
        L += model.minimization.calving_terminus(**flds, outflow_ids=calving_ids)

    PETSc.Sys.Print("  Continuation...")
    stop_manager()
    prob = NonlinearVariationalProblem(
        derivative(L, z_adj), z_adj, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)
    for exp in np.linspace(1.0, glen_flow_law, 5):
        n_glen.assign(exp)
        slvr.solve()

    u_ref = z_adj.subfunctions[0]
    PETSc.Sys.Print(
        f"  u_max = {global_max(Function(Q).interpolate(sqrt(inner(u_ref, u_ref)))):.0f} m/yr"
    )

    rho_ice = 917.0
    flux_ref_total = (
        abs(float(assemble(h_ctrl * inner(u_ref, grad(grounded_indicator)) * dx)))
        * rho_ice
        / 1e12
    )
    PETSc.Sys.Print(f"  Reference GL flux: {flux_ref_total:.2f} Gt/yr")

    # ══════════════════════════════════════════════════════════════════════
    # Adjoint: dJ/dH
    # ══════════════════════════════════════════════════════════════════════
    PETSc.Sys.Print("\n=== Adjoint solve ===")
    reset_manager()
    start_manager()
    clear_caches()

    # Surface on tape
    s_tape = max_value(b + h_ctrl, (Constant(1.0) - rho_ratio) * h_ctrl)

    u_s, M_s, tau_s = split(z_adj)
    flds2 = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": h_ctrl,
        "surface": s_tape,
    }
    rh2 = {
        "flow_law_exponent": n_glen,
        "flow_law_coefficient": A_map_expr,
        "sliding_exponent": n_glen,
        "sliding_coefficient": K_map_expr,
    }
    L2 = (
        model.minimization.viscous_power(**flds2, **rh2)
        + model.minimization.friction_power(**flds2, **rh2)
        + model.minimization.momentum_balance(**flds2)
    )
    if use_calving_terminus:
        L2 += model.minimization.calving_terminus(**flds2, outflow_ids=calving_ids)
    F2 = derivative(L2, z_adj)

    PETSc.Sys.Print("  Annotated solve...")
    for exp in np.linspace(1.0, glen_flow_law, 5):
        n_glen.assign(exp)
        EquationSolver(
            F2 == 0,
            z_adj,
            solver_parameters=sparams,
            form_compiler_parameters=fc_params,
        ).solve()

    u_sol, _, _ = split(z_adj)
    J = Functional(name="J_gl_flux")
    J.assign(h_ctrl * inner(u_sol, grad(grounded_indicator)) * dx)

    stop_manager()
    PETSc.Sys.Print("  Computing dJ/dH...")
    dJ_dH = compute_gradient(J, h_ctrl)
    dJ_dH_func = Function(Q, name="dJ_dH")
    dJ_dH_func.dat.data[:] = dJ_dH.dat.data_ro
    _dj_lo, _dj_hi = global_range(dJ_dH_func)
    PETSc.Sys.Print(f"  dJ/dH range: [{_dj_lo:.4e}, {_dj_hi:.4e}]")
    dJ_arr = dJ_dH_func.dat.data_ro

    # ══════════════════════════════════════════════════════════════════════
    # GL thinning sensitivity
    #
    # Perturbation: ΔH(x) = -ε · φ(d_GL(x)) on grounded ice near GL
    # where φ(d) = exp(-d/w) concentrates thinning near the GL.
    #
    # The adjoint prediction: δJ = Σ_i (dJ/dH)_i · ΔH_i
    # Normalized by ε gives: dJ/dε = Σ_i (dJ/dH)_i · (-φ_i)
    # ══════════════════════════════════════════════════════════════════════
    PETSc.Sys.Print(f"\n{'='*60}")
    PETSc.Sys.Print("GL thinning sensitivity (adjoint)")
    PETSc.Sys.Print(f"{'='*60}")

    rows = []
    for w_km in THINNING_WIDTHS_KM:
        w_m = w_km * 1e3
        PETSc.Sys.Print(f"\n--- Thinning width w = {w_km} km ---")

        # Thinning profile: φ(d) = exp(-d/w) on grounded ice
        phi_profile = np.exp(-d_gl_arr / w_m) * grounded_arr.astype(float)

        for shelf_name in sorted(
            shelf_mask_arrays,
            key=lambda n: -global_sum(shelf_mask_fns.get(n, Function(Q))),
        ):
            in_shelf = shelf_mask_arrays[shelf_name]
            # Thinning on grounded ice near this shelf's GL
            phi_shelf = phi_profile * in_shelf.astype(float)

            n_affected = int((phi_shelf > 0.01).sum())
            if n_affected == 0:
                continue

            # Adjoint prediction: dJ/dε_k = Σ_i dJ/dH_i · (-φ_i)
            dJ_deps = -np.dot(dJ_arr, phi_shelf)
            # Convert to Gt/yr per m of thinning
            dJ_deps_gt = dJ_deps * rho_ice / 1e12

            PETSc.Sys.Print(
                f"  {shelf_name:20s}: {n_affected:5d} nodes, "
                f"dJ/dε = {dJ_deps_gt:+.4f} Gt/yr per m thinning"
            )

            rows.append(
                {
                    "shelf": shelf_name,
                    "width_km": w_km,
                    "n_affected": n_affected,
                    "dJ_deps_gt_yr_per_m": dJ_deps_gt,
                }
            )

    # ══════════════════════════════════════════════════════════════════════
    # Forward FD validation for selected shelves at w=10km
    # ══════════════════════════════════════════════════════════════════════
    w_val = 10e3  # 10km width for validation
    eps_val = 10.0  # 10m thinning (large enough for signal)
    PETSc.Sys.Print(
        f"\n=== Forward FD validation (w={w_val/1e3:.0f}km, ε={eps_val}m) ==="
    )

    for val_shelf in VALIDATE_SHELVES:
        if val_shelf not in shelf_mask_arrays:
            continue
        in_shelf = shelf_mask_arrays[val_shelf]
        phi_profile = np.exp(-d_gl_arr / w_val) * grounded_arr.astype(float)
        phi_shelf = phi_profile * in_shelf.astype(float)

        if (phi_shelf > 0.01).sum() == 0:
            continue

        PETSc.Sys.Print(f"  {val_shelf}: forward solve...")

        # Perturbed thickness: H - ε·φ
        h_pert = Function(Q)
        h_pert.dat.data[:] = H_arr - eps_val * phi_shelf
        h_pert.interpolate(max_value(h_pert, Constant(h_clamp)))

        s_val = Function(Q).interpolate(
            max_value(b + h_pert, (Constant(1.0) - rho_ratio) * h_pert)
        )

        z_val = Function(Z)
        z_val.assign(z_adj)
        stop_manager()
        u_v, M_v, tau_v = split(z_val)
        flds_v = {
            "velocity": u_v,
            "membrane_stress": M_v,
            "basal_stress": tau_v,
            "thickness": h_pert,
            "surface": s_val,
        }
        rh_v = {
            "flow_law_exponent": n_glen,
            "flow_law_coefficient": A_map_expr,
            "sliding_exponent": n_glen,
            "sliding_coefficient": K_map_expr,
        }
        L_v = (
            model.minimization.viscous_power(**flds_v, **rh_v)
            + model.minimization.friction_power(**flds_v, **rh_v)
            + model.minimization.momentum_balance(**flds_v)
        )
        if use_calving_terminus:
            L_v += model.minimization.calving_terminus(
                **flds_v, outflow_ids=calving_ids
            )

        prob_v = NonlinearVariationalProblem(
            derivative(L_v, z_val), z_val, form_compiler_parameters=fc_params
        )
        slvr_v = NonlinearVariationalSolver(prob_v, solver_parameters=sparams)
        slvr_v.solve()

        u_val = z_val.subfunctions[0]
        flux_pert = (
            abs(float(assemble(h_pert * inner(u_val, grad(grounded_indicator)) * dx)))
            * rho_ice
            / 1e12
        )

        delta_fd = flux_pert - flux_ref_total
        # Per-meter: ΔF / ε
        delta_fd_per_m = delta_fd / eps_val

        # Adjoint prediction at this width
        adj_row = [
            r for r in rows if r["shelf"] == val_shelf and r["width_km"] == w_val / 1e3
        ]
        adj_val = adj_row[0]["dJ_deps_gt_yr_per_m"] if adj_row else np.nan

        PETSc.Sys.Print(f"    FD:      ΔF/ε = {delta_fd_per_m:+.6f} Gt/yr per m")
        PETSc.Sys.Print(f"    Adjoint: dJ/dε = {adj_val:+.6f} Gt/yr per m")
        if not np.isnan(adj_val) and abs(adj_val) > 1e-10:
            PETSc.Sys.Print(f"    Ratio FD/adj: {delta_fd_per_m / adj_val:.3f}")

    # ── Output ──
    df = pd.DataFrame(rows)
    csv_fn = os.path.join(RESULTS_DIR, f"gl_sensitivity_{lc}.csv")
    df.to_csv(csv_fn, index=False, float_format="%.6f")
    PETSc.Sys.Print(f"\nSaved: {csv_fn}")

    fields_fn = os.path.join(RESULTS_DIR, f"gl_sensitivity_fields_{lc}.h5")
    with fd.CheckpointFile(fields_fn, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(dJ_dH_func, name="dJ_dH")
        chk.save_function(d_gl, name="gl_distance")
        chk.save_function(H, name="thickness")
    PETSc.Sys.Print(f"Saved fields: {fields_fn}")


if __name__ == "__main__":
    main()
