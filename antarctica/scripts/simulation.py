#!/usr/bin/env python3
r"""Shared simulation engine for ISMIP7 Antarctic experiments."""

import numpy as np
import os, sys, glob
from time import perf_counter

# Wall clock from as close to process start as this module can observe, so
# ISMIP7_WALL_STOP_MIN counts the setup (mesh + MAP load) it has to pay for
# too, not just the time loop.
_T_PROCESS_START = perf_counter()

import firedrake as fd
from firedrake import (
    Constant,
    Function,
    max_value,
    sqrt,
    derivative,
    dx,
    dS,
    ds,
    split,
    assemble,
    Mesh,
    FunctionSpace,
    VectorFunctionSpace,
    TensorFunctionSpace,
    FiniteElement,
    LinearVariationalProblem,
    LinearVariationalSolver,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    exp,
)
from firedrake.petsc import PETSc

import rasterio, icepack
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
)

# SI densities for diagnostics (icepack2 constants are in MPa-m-yr units)
_RHO_I_SI = 917.0
_RHO_W_SI = 1024.0

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# BedMachine, MEaSUREs velocity, and RACMO may live outside the checkout on a
# workstation with a local data volume. Keep the repository layout as the
# default, but make the observational root explicit rather than requiring
# large files to be copied or symlinked into the tree.
MESH_DIR = os.path.join(_ROOT, "mesh")
RESULTS_DIR = os.path.join(_ROOT, "results")

# Repo root on the path for the shared dual-friction operator.
sys.path.insert(0, os.path.dirname(_ROOT))
from mesh_naming import mesh_filename

from icepack2_tools.transfer import interpolate_with_fill
from icepack2_tools.mpi_stats import (
    global_count,
    global_extreme_location,
    global_mean,
    global_range,
    global_size,
)
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.geometry import sample_to_geometry
from icepack2_tools.naming import map_basename
from icepack2_tools.front import (
    clamp_thickness, clear_reference_where_ice_free, retreat_slivers,
    unforced_cells, applied_forcing,
    facet_neighbours, front_connected, ocean_drag_cells,
    collapse_banner, collapse_cell_counts, collapse_csv_fields,
    COLLAPSE_CSV_COLUMNS, COLLAPSE_MARKER, FRONT_OWNER_MARKER,
)
from icepack2_tools.runconfig import (
    obs_data_root,
    BUDD_SHELF_GATE as _BUDD_SHELF_GATE,
    residual_stabilizers,
    friction as _friction, geometry_space as _geometry_space,
    mesh_override as _mesh_override,
    lc as _lc, lc_coarse as _lc_coarse, n_flow as _n_flow,
    TARGET_MESH_GEOMETRY_METHOD,
    calving_law as _calving_law, calving_sigma_max as _calving_sigma_max,
    fracture as _fracture_mode, ismip7_output as _ismip7_output,
    FRACTURE_MASK_MODES,
    # auto_resume is re-exported, not used here: every forward driver imports
    # it from this module alongside latest_checkpoint, so they resolve the
    # knob through one import rather than each reaching into runconfig.
    fixed_front as _fixed_front, auto_resume, apparent_mb_mode,  # noqa: F401
    front_hmin as _front_hmin,
    GEOMETRY_YEAR,
)
DATA_DIR = obs_data_root()
from icepack2_tools.solverconfig import (
    final_solve_bounds,
    continuation_steps,
    diagnostic_solver_label,
    diagnostic_solver_mode,
    diagnostic_solver_parameters,
    linearization_state,
    mass_residual_tol_gt,
    rescue_enabled,
    rescue_max_it,
    snes_atol_scale,
    snes_monitor_enabled,
    snes_restart_failure_atol_scale,
    solver_view_enabled,
    subcycles,
    transport_solver_parameters,
)

lc = _lc()
lc_coarse = _lc_coarse()
buffer_m = float(os.environ.get("ISMIP7_BUFFER_M", "20000"))

# Flow-law exponent for the composite viscous rheology: owned by
# icepack2_tools.runconfig, which the inversion that produced the MAP reads
# too. THIS BRANCH (antarctica-n3) runs STANDARD GLEN n=3: A0 =
# rate_factor(260 K) is already the n=3 fluidity, so the composite main term
# needs no prefactor rescale (A4_FACTOR_DEFAULT = 1). The n=4
# Goldsby-Kohlstedt composite (a4_factor ~ 10, so A_4 tau_c^4 ~ A_3 tau_c^3 at
# tau_c) lives on the `antarctica` branch. Override either per-run with
# ISMIP7_N_FLOW / ISMIP7_A4_FACTOR; the two must match between an inversion
# and the forward runs that load its MAP. Setting ISMIP7_N_FLOW=4 alone
# therefore carries the n=4 prefactor with it (a4_factor_default below), so
# the legacy untagged n=4 MAP is used with the factor it was inverted at.
A4_FACTOR_DEFAULT = "1.0"
A4_FACTOR_N4 = "10.0"


def a4_factor_default():
    r"""Prefactor default derived from the flow exponent: the n=4 composite
    needs A_4 = 10 A_3 so A_4 tau_c^4 ~ A_3 tau_c^3, while n=3 is plain Glen
    and needs none. ISMIP7_A4_FACTOR still overrides."""
    n = _n_flow()
    return A4_FACTOR_N4 if abs(n - 4.0) < 1e-9 else A4_FACTOR_DEFAULT


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def load_onto(chk, mesh, name, space):
    r"""Load a checkpoint field onto `space`, converting representation if the
    file stores it elsewhere. Returns ``(Function, converted)``.

    Same space is the normal path and copies dofs directly (the mesh came from
    this same checkpoint, so the ordering matches by construction). A different
    space means the file predates the DG0-geometry change; project and let the
    caller warn. Projection makes the run *executable*, not *consistent*: a
    legacy MAP's theta/phi were inferred against CG1 geometry with the biased
    calving front, so the controls still carry that bias and the MAP must be
    re-inverted before the results mean anything.
    """
    f = chk.load_function(mesh, name=name)
    if f.function_space().ufl_element() == space.ufl_element():
        out = Function(space, name=name)
        out.dat.data[:] = f.dat.data_ro
        return out, False
    return Function(space, name=name).project(f), True


def latest_checkpoint(experiment_name, lc_val=None):
    r"""Newest self-contained state checkpoint for an experiment, or None.

    Scans RESULTS_DIR for `{experiment_name}_{lc}_final.h5` and the periodic
    `{experiment_name}_{lc}_t<year>.h5` files and returns the path with the
    largest saved `t_yr` (final wins ties). Enables unattended auto-resume:
    a rebooted run continues from where it left off with no manual bookkeeping.
    """
    lc_val = lc if lc_val is None else lc_val
    cands = []
    final_fn = os.path.join(RESULTS_DIR, f"{experiment_name}_{lc_val}_final.h5")
    if os.path.exists(final_fn):
        cands.append(final_fn)
    cands += glob.glob(
        os.path.join(RESULTS_DIR, f"{experiment_name}_{lc_val}_t*.h5")
    )
    best, best_t = None, None
    for fn in cands:
        try:
            with fd.CheckpointFile(fn, "r") as chk:
                t = (float(chk.get_attr("/", "t_yr"))
                     if chk.has_attr("/", "t_yr") else None)
        except Exception:
            t = None
        if t is None:
            continue
        if (best_t is None or t > best_t
                or (t == best_t and fn.endswith("_final.h5"))):
            best, best_t = fn, t
    return best


def historical_endpoint(esm_tag, tag_sfx, t_branch, lc_val=None):
    r"""The historical endpoint a control or projection branches from, or None
    when there is none.

    A historical chain rewrites its ``_final.h5`` at the end of every job, so
    a chain that stalled, or was stopped, leaves one that holds the year it
    reached rather than the handoff. A control or projection that branched
    from it would start its 2015 experiment on an earlier geometry and say
    nothing, so an endpoint short of ``t_branch`` is refused outright.
    """
    lc_val = lc if lc_val is None else lc_val
    path = os.path.join(RESULTS_DIR, f"hist_{esm_tag}{tag_sfx}_{lc_val}_final.h5")
    if not os.path.exists(path):
        return None
    with fd.CheckpointFile(path, "r") as chk:
        t = (float(chk.get_attr("/", "t_yr"))
             if chk.has_attr("/", "t_yr") else None)
    if t is None or t < t_branch - 1e-6:
        reached = "no t_yr" if t is None else f"t_yr={t:g}"
        raise RuntimeError(
            f"the historical endpoint {path} holds {reached}, short of the "
            f"{t_branch:g} handoff: its chain stopped early. Finish the "
            f"historical, or name a state with ISMIP7_RESTART.")
    return path


def setup_model(restart_from=None, *, allow_timing_cache_a_ref=False,
                backdate_years=0.0):
    r"""Load mesh, data, inversion fields, and build diagnostic solver.

    ``allow_timing_cache_a_ref`` is the narrow exception used after a timing
    manifest has been validated: it permits ``APPARENT_MB=div`` to be built
    from that pristine initial state. Evolved restarts must carry their frozen
    correction and cannot use this escape hatch.

    ``backdate_years`` dates a cold start before the 2015 geometry: that many
    years of the Smith mean dH/dt are undone on grounded ice after the friction
    anchors are built, so the 2015 friction and fluidity carry over unchanged
    and the run starts from the earlier ice (issue #117). A restart ignores it.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Friction law: "budd" (power-law dual, default) or "regularized_coulomb"
    # (RC residual driven by the inversion_icepack2_rc_<lc>.h5 MAP:
    # grounded-only theta, exact-zero shelf drag via the Coulomb cap c0*N).
    # Must mirror inversion_icepack2.py so theta/phi keep their meaning.
    # Resolved first because the compute mesh is loaded FROM this friction's
    # MAP checkpoint (below), not from a fresh Mesh(.msh).
    friction = _friction()
    # Reject a mistyped front configuration before the MAP load and the
    # initial Newton solve, which cost minutes to tens of minutes at 2500 m
    # on a detached launch.
    _calving_law()
    _calving_sigma_max()
    # Exact-zero-shelf residual laws (icepack2 dual, dual_friction.py):
    #   regularized_coulomb -> Coulomb cap tau_c=c0*N
    #   budd                -> tau_b ~ N_hat=N_eff/N_ref, PISM-delta grounded
    #                          floor, exact-zero shelf (matches gia ASE ase_model)
    # Both share the C_w0/He/composite structure and need a _rc/_budd MAP.
    # "budd_legacy" keeps the old phi_eff action Budd (residual shelf drag).
    use_residual = friction in ("regularized_coulomb", "budd")
    use_rc = use_residual  # geometry/alpha/h_clamp handling is shared
    drag_mask = None       # set in the residual branch below
    drag_rule = None       # ditto: rewrites drag_mask for a new thickness

    is_restart = restart_from is not None
    # Prefer the MAP inverted under this geometry space. Falling back to the
    # untagged legacy (CG1) MAP keeps the transition runnable, but the load
    # then projects the geometry and warns loudly - the controls still carry
    # the CG1 front bias, so such a run is a smoke test, not a result.
    # ISMIP7_INVERSION names the MAP explicitly (the IU timing matrix and
    # A/B forwards of differently-regularised MAPs on one mesh). It must
    # still be a MAP of THIS friction/n/geometry; the name is not checked.
    inv_override = os.environ.get("ISMIP7_INVERSION")
    inv_fn = inv_override or os.path.join(MESH_DIR, map_basename(friction, lc))
    if inv_override and not is_restart and not os.path.exists(inv_fn):
        raise FileNotFoundError(f"ISMIP7_INVERSION={inv_override} does not exist")
    if not inv_override and not is_restart and not os.path.exists(inv_fn):
        legacy_fn = os.path.join(
            MESH_DIR, map_basename(friction, lc, geometry=False)
        )
        if os.path.exists(legacy_fn):
            PETSc.Sys.Print(
                f"  No {os.path.basename(inv_fn)}; falling back to "
                f"{os.path.basename(legacy_fn)} (inverted under a different "
                f"geometry space)"
            )
            inv_fn = legacy_fn
    # Take the mesh + reference fields from the checkpoint we start from: the
    # MAP on a cold start, or the self-contained restart checkpoint on a warm
    # start. A fresh Mesh(.msh) repartitions with a different dof ordering,
    # so a raw copy of the saved theta/phi/geometry scrambles at any rank
    # count != the one the checkpoint was written on (the -n4 tripwire crash).
    # The checkpoint mesh carries the full boundary-marker set, so the
    # calving BC is preserved.
    source_chk = restart_from if is_restart else inv_fn
    PETSc.Sys.Print(f"Loading mesh + reference state: {source_chk}")
    with fd.CheckpointFile(source_chk, "r") as _chk:
        source_mesh = _chk.load_mesh()
        # A cold start normally binds itself to a MAP of the right
        # configuration through map_basename, which encodes friction, n and
        # geometry space. ISMIP7_INVERSION bypasses the name, so check the
        # MAP's own record instead: theta means a different thing under each
        # friction law, and driving RC controls through the Budd block runs
        # cleanly and returns wrong velocities.
        if inv_override and not is_restart:
            _want = {"friction": friction, "n_flow": float(_n_flow()),
                     "geometry_space": _geometry_space()}
            _missing = [k for k in _want if not _chk.has_attr("/", k)]
            if _missing:
                PETSc.Sys.Print(
                    f"  WARNING: ISMIP7_INVERSION={inv_override} records no "
                    f"{'/'.join(_missing)}; it predates the attribute and "
                    f"cannot be checked against friction={friction}, "
                    f"n={_n_flow():g}, geometry={_geometry_space()}. Confirm "
                    f"it was inverted under those."
                )
            for _k, _v in _want.items():
                if _k in _missing:
                    continue
                _got = _chk.get_attr("/", _k)
                _got = float(_got) if _k == "n_flow" else str(_got)
                if _got != _v:
                    raise RuntimeError(
                        f"ISMIP7_INVERSION={inv_override} was inverted with "
                        f"{_k}={_got!r} but this run resolves {_k}={_v!r}. "
                        f"The controls only mean anything under the "
                        f"configuration they were inverted for; point at a "
                        f"matching MAP or change the run's configuration."
                    )
        checkpoint_metadata = {}
        for _key in (
            "timing_cache_schema_version",
            "timing_cache_role",
            "source_inversion",
            "source_inversion_sha256",
            "source_mesh_sha256",
            "diagnostic_solver_mode",
            "solver_configuration",
            "solver_configuration_fingerprint",
            "geometry_space",
            "n_flow",
            "a4_factor",
            "geometry_source",
            "geometry_source_method",
            # Shelf gate the Budd state was solved under (runconfig
            # .BUDD_SHELF_GATE); the timing lane checks it against the manifest.
            "friction_gate",
            # Residual of the saved mixed state under its writer's F; the
            # restart fast path trusts the state only within a factor of it.
            "full_state_residual",
        ):
            if _chk.has_attr("/", _key):
                checkpoint_metadata[_key] = _chk.get_attr("/", _key)
        # The mesh this checkpoint was built on, recorded by the inversion and
        # carried through every restart. A CheckpointFile mesh is named
        # "firedrake_default", so this attribute is the only way the run can
        # name its own .msh and pick the matching per-mesh sidecar.
        mesh_basename = (
            str(_chk.get_attr("/", "mesh_basename"))
            if _chk.has_attr("/", "mesh_basename") else ""
        )

        # The collaborators additionally stamp the mesh PARAMETERS. Keep both:
        # the basename is direct and also covers meshes outside the standard
        # naming pattern (e.g. the 500 m aniso mesh), while lc_coarse/buffer_m
        # let bndids_filename() reconstruct the name parametrically. They
        # cross-check each other, and either alone is enough to resolve the
        # sidecar from the CHECKPOINT rather than from the live environment --
        # which is the point: ISMIP7_BUFFER_M/ISMIP7_LC_COARSE can drift, and
        # a mismatched sidecar puts the calving BC on the wrong facets, where
        # ds(absent_id) integrates to zero: silently wrong physics, no crash.
        chk_lc_coarse = (int(_chk.get_attr("/", "lc_coarse"))
                         if _chk.has_attr("/", "lc_coarse") else None)
        chk_buffer_m = (float(_chk.get_attr("/", "buffer_m"))
                        if _chk.has_attr("/", "buffer_m") else None)
        # Raster sampling the MAP was inverted with. MAPs older than the
        # attribute were all vertex-sampled. The MAP wins over the
        # environment: the controls absorbed that bed, so a different one
        # here would be silently inconsistent.
        chk_raster_sample = (str(_chk.get_attr("/", "raster_sample"))
                             if _chk.has_attr("/", "raster_sample") else "vertex")
        _env_rs = os.environ.get("ISMIP7_RASTER_SAMPLE")
        if _env_rs and _env_rs.lower() != chk_raster_sample:
            PETSc.Sys.Print(
                f"  WARNING: ISMIP7_RASTER_SAMPLE={_env_rs} but the MAP was "
                f"inverted with {chk_raster_sample}; using the MAP's."
            )
        # The FINE resolution is stamped too, and every component of the
        # reconstructed name must come from the checkpoint: falling back to the
        # live ISMIP7_LC here would reintroduce exactly the environment drift
        # the parametric scheme exists to remove. Older MAPs that predate the
        # lc attribute have no recorded resolution, so they keep the live value.
        chk_lc = (int(_chk.get_attr("/", "lc"))
                  if _chk.has_attr("/", "lc") else None)
        if not mesh_basename and chk_lc_coarse is not None \
                and chk_buffer_m is not None:
            mesh_basename = os.path.basename(
                mesh_filename(chk_lc_coarse,
                              chk_lc if chk_lc is not None else lc,
                              chk_buffer_m))

    source_mesh_basename = mesh_basename
    # None for unset, empty and the sentinel `checkpoint`: solve on the mesh
    # the checkpoint carries.
    mesh_fn = _mesh_override()
    if mesh_fn:
        # The timing matrix deliberately solves on a mesh different from the
        # MAP mesh. Keep the checkpoint mesh as the interpolation source and
        # use the requested mesh only for the target finite-element spaces.
        mesh = Mesh(mesh_fn)
        target_lc_coarse = lc_coarse
        target_buffer_m = buffer_m
        mesh_basename = os.path.basename(mesh_fn)
        chk_lc = lc
        chk_lc_coarse = target_lc_coarse
        chk_buffer_m = target_buffer_m
        PETSc.Sys.Print(
            f"  Compute mesh override: {mesh_fn} "
            f"({global_size(mesh.coordinates)} vertices, "
            f"{mesh.comm.allreduce(mesh.cell_set.size)} cells)"
        )
    else:
        mesh = source_mesh
        target_lc_coarse = chk_lc_coarse
        target_buffer_m = chk_buffer_m
    # num_vertices()/num_cells() count this rank's plex, halo included; the
    # coordinate dofs and the owned cell set are reduced to global totals.
    PETSc.Sys.Print(
        f"  {global_size(mesh.coordinates)} vertices, "
        f"{mesh.comm.allreduce(mesh.cell_set.size)} cells"
    )
    # The recorded basename is the provenance and WINS. ISMIP7_MESH is only a
    # fallback for legacy checkpoints that carry no attribute: it names the
    # mesh the caller intends to build, which is not necessarily the one this
    # checkpoint was written on, so it must not override the record - and only
    # the recorded value is carried into the checkpoints this run writes.
    if not mesh_basename:
        _env_mesh = _mesh_override() or ""
        mesh_basename = os.path.basename(_env_mesh) if _env_mesh else ""
        if mesh_basename:
            PETSc.Sys.Print(
                f"  No mesh_basename in {os.path.basename(source_chk)}; using "
                f"ISMIP7_MESH ({mesh_basename}) to resolve the boundary sidecar"
            )

    # Boundary-id sidecar: resolved (per-mesh preferred) and hard-checked
    # against this mesh by the shared helper, so the inversion, the forward and
    # every auxiliary solver cannot disagree about the ice front.
    use_calving_terminus = os.environ.get("ISMIP7_NO_CALVING_TERMINUS") is None
    bnd_ids, calving_ids, bndids_fn = load_boundary_ids(
        mesh, MESH_DIR, mesh_hint=mesh_basename,
        print_coverage=use_calving_terminus,
    )

    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    Q_dg = FunctionSpace(mesh, "DG", 0)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    # GEOMETRY SPACE (h, s, b and everything derived from them), separate from
    # the CONTROL space Q. Controls stay CG1: the Whittle-Matern prior
    # regularizes |grad theta|, which needs a differentiable field.
    #
    # DG0 (default) makes the geometry the SAME field the mass transport
    # evolves, which is the point. Under the old CG1 geometry the transport
    # carried a DG0 thickness and the momentum solve a CG1 one, bridged by a
    # lumped-mass lift. That lift conserves volume exactly but is one-sided at
    # the domain boundary, so it pulled calving-front nodes up toward interior
    # values: measured at 32 km it DOUBLED the front thickness (105 -> 200 m,
    # against a BedMachine ice-front mean of 145 m), and since the terminus
    # traction goes as h^2 that tripled the terminus force and multiplied the
    # outflux by ~4.7. Force and mass also disagreed by ~22% because they
    # integrated different fields. With DG0 geometry there is one thickness:
    # the terminus back-pressure and the transport flux integrate the same
    # cell value, and no lift exists to bias the front.
    #
    # The cost, measured on a smooth manufactured problem: the driving stress
    # drops from 2nd-order (CG1 cell term) to 1st-order (DG0 facet-jump term),
    # ~1.4% vs ~0.08% error at 32x32. That is well inside the bed/thickness
    # data error at Antarctic resolutions, and it buys an unbiased front.
    #
    # ISMIP7_GEOMETRY_SPACE=cg1 restores the old behaviour for A/B comparison.
    geometry_space = _geometry_space()
    geom_dg = geometry_space == "dg0"
    Q_g = FunctionSpace(mesh, "DG", 0) if geom_dg else Q
    PETSc.Sys.Print(
        f"  Geometry space: {geometry_space.upper()}"
        + (" (h, s, b cell-wise; one thickness for force and mass)" if geom_dg
           else " (LEGACY: CG1 geometry, lumped lift, front thickness biased high)")
    )

    PETSc.Sys.Print("Loading data...")
    # Two clamps:
    #   - h_clamp_init: floor on the *initial* thickness so the first
    #     diagnostic solve is well-posed everywhere (default 10 m; 0 in RC
    #     mode, where the MAP was inverted against the true h=0 geometry and
    #     h_visc_floor handles the ice-free buffer).
    #   - h_clamp: floor in the advection step. Default 0 so melt can drive
    #     cells to zero at the calving front; the composite rheology keeps the
    #     diagnostic SNES nonsingular at h=0.
    h_clamp_init = float(
        os.environ.get("ISMIP7_H_CLAMP_INIT", "0.0" if use_rc else "10.0")
    )
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "0.0"))
    rho_ratio = Constant(917.0 / 1024.0)

    # Velocity obs (raster; interpolation is mesh-agnostic, so it is valid on
    # the checkpoint mesh). Needed for the sliding scale u_c and, on a cold
    # start, the Weertman anchor + initial misfit; incidental on restart.
    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V,
        fillvalue=0.0,
    )

    geometry_source = checkpoint_metadata.get("geometry_source")
    geometry_source_method = checkpoint_metadata.get(
        "geometry_source_method"
    )
    if not is_restart:
        # Cold start: geometry from BedMachine (RC/Budd overwrites it with the
        # inversion-time geometry from the MAP only when no target-mesh
        # override is active; target timing meshes retain this cell average).
        bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
        geometry_source = os.path.realpath(bm_fn)
        geometry_source_method = (
            TARGET_MESH_GEOMETRY_METHOD if geom_dg
            else "target-native-bedmachine-nodal-v1"
        )
        # Cell average onto the geometry space, NOT a centroid point sample --
        # see geometry.sample_to_geometry for the measurements behind that.
        PETSc.Sys.Print(f"  Raster sampling onto geometry cells: {chk_raster_sample}")
        b = sample_to_geometry(
            rasterio.open(f"netcdf:{bm_fn}:bed"), Q_g, Q,
            method=chk_raster_sample)
        b.rename("bed")
        H = sample_to_geometry(
            rasterio.open(f"netcdf:{bm_fn}:thickness"),
            Q_g, Q, floor=h_clamp_init, method=chk_raster_sample)
        H.rename("thickness")
        s = Function(Q_g, name="surface").interpolate(
            max_value(b + H, (Constant(1.0) - rho_ratio) * H)
        )

    # Reference log-adjustments (theta=log_friction, phi=log_fluidity) plus,
    # for RC/Budd or any restart, the geometry and frozen anchors - all read
    # onto the mesh we already took from THIS checkpoint, so the dof order
    # matches by construction (no fragile .msh-vs-checkpoint node comparison).
    C_w0 = None
    N_ref = None
    H_init = None
    phi_eff = None
    u_guess = None
    M_guess = None
    tau_guess = None
    a_ref_mb = None
    phys_div = None
    h_dg_state = None
    t_restart = None
    A_prior_f = None
    ismip7_resume = None
    # The fluidity baseline, here because the loader below fills the prior
    # with it where the compute mesh reaches past the MAP's.
    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    a4_factor = float(os.environ.get("ISMIP7_A4_FACTOR", a4_factor_default()))
    A_prior_baseline = float(A0) * a4_factor
    # One entry per field loaded across meshes: how many target dofs lay
    # outside the source mesh and what they were filled with. Printed below
    # and carried in the context for the caches and the map checks.
    transfer_fill = {}

    def load_checkpoint_field(chk, name, space, optional=False,
                              fill=0.0, fill_label="0"):
        """Load a checkpoint field, interpolating it onto the compute mesh.

        With ISMIP7_MESH naming another mesh, a target dof outside the source
        mesh takes ``fill``: a float, or a Function on ``space`` (the raster
        sample for velocity_obs). The buffered production mesh reaches 20 km
        past a buffer-0 MAP, so the whole ring is filled; zero there is the
        prior for theta and phi and a singular block for the fluidity prior
        (icepack2_tools.transfer). Same-mesh loads miss nothing.
        """
        try:
            source_field = chk.load_function(source_mesh, name=name)
        except (KeyError, RuntimeError, ValueError):
            if optional:
                return None
            raise
        target_field = Function(space, name=name)
        n_missing, n_total, n_clamped = interpolate_with_fill(
            target_field, source_field, fill, mesh.comm
        )
        transfer_fill[name] = {
            "missing": n_missing, "total": n_total, "fill": fill_label,
            "clamped": n_clamped,
        }
        return target_field

    with fd.CheckpointFile(source_chk, "r") as chk:
        theta_f = load_checkpoint_field(chk, "log_friction", Q)
        theta_f.rename("theta")
        phi_f = load_checkpoint_field(chk, "log_fluidity", Q)
        phi_f.rename("phi")
        # Fluidity prior mean (physical thermomechanical field): the fluidity
        # control is phi = log(A / A_prior), so the forward must reconstruct
        # A = A_prior * exp(phi) with the SAME A_prior the inversion used. New
        # MAPs and restart checkpoints carry it; older ones (constant-baseline
        # MAPs) don't, and A4_base falls back to A0*a4_factor below.
        A_prior_f = load_checkpoint_field(
            chk, "fluidity_prior", Q, optional=True,
            fill=A_prior_baseline,
            fill_label=f"the constant baseline A0*a4_factor = {A_prior_baseline:.3g}",
        )
        if is_restart:
            # Self-contained restart: evolved geometry, frozen anchors, time.
            b = load_checkpoint_field(chk, "bed", Q_g)
            H = load_checkpoint_field(chk, "thickness", Q_g)
            s = load_checkpoint_field(chk, "surface", Q_g)
            H_init = load_checkpoint_field(chk, "H_init", Q_g)
            phi_eff = load_checkpoint_field(chk, "phi_eff", Q_g)
            u_guess = load_checkpoint_field(chk, "velocity", V)
            # Preserve the exact observation field used to construct the
            # inversion-time friction anchor and report its initial misfit.
            # A fresh raster interpolation is only a compatibility fallback
            # for ordinary checkpoints written before this field was saved;
            # current timing caches require velocity_obs and fail loudly if it
            # is absent.
            cached_u_obs = load_checkpoint_field(
                chk, "velocity_obs", V, optional=True,
                fill=u_obs, fill_label="the raster-sampled velocity_obs",
            )
            if cached_u_obs is not None:
                u_obs.assign(cached_u_obs)
            elif checkpoint_metadata.get("timing_cache_role") == \
                    "timing-initial-state":
                raise RuntimeError(
                    "Timing cache is missing required velocity_obs"
                )
            # Stress components (newer checkpoints): restoring them makes the
            # resume Newton start from the full converged state instead of
            # (u, 0, 0), which needed a fresh continuation ramp.
            try:
                M_guess = load_checkpoint_field(chk, "membrane_stress", Sigma)
                tau_guess = load_checkpoint_field(chk, "basal_stress", T)
            except (KeyError, RuntimeError, ValueError):
                M_guess = tau_guess = None
            # Frozen apparent-MB reference (present iff the run used
            # ISMIP7_APPARENT_MB): restarts must reuse the ORIGINAL t=0
            # correction, never recompute it from the evolved state.
            a_ref_mb = load_checkpoint_field(
                chk, "a_ref_mb", Q_dg, optional=True
            )
            # An adapted checkpoint carries the transferred PHYSICAL divergence
            # instead of a_ref (icepack2_tools.adapt_mesh); a_ref is rebuilt
            # below with this mesh's own operator.
            phys_div = load_checkpoint_field(
                chk, "phys_div", Q_dg, optional=True
            )
            # Separate DG0 prognostic thickness (CG1-geometry runs only, where
            # the stored CG h was its lumped lift). Under DG0 geometry the
            # thickness IS the transport state, so there is nothing to restore.
            if not geom_dg:
                try:
                    h_dg_state = load_checkpoint_field(
                        chk, "thickness_dg", Q_dg
                    )
                except (KeyError, RuntimeError, ValueError):
                    h_dg_state = None
            if use_residual:
                C_w0 = load_checkpoint_field(chk, "C_w0", Q_g)
            if friction == "budd":
                N_ref = load_checkpoint_field(chk, "N_ref", Q_g)
            if chk.has_attr("/", "t_yr"):
                t_restart = float(chk.get_attr("/", "t_yr"))
            # The ISMIP7 year in progress, so a link that stopped mid-year
            # continues the same year's flux means instead of losing them.
            from icepack2_tools.ismip7_output import AnnualOutput as _AnnualOutput
            ismip7_resume = _AnnualOutput.read_state(chk, mesh)
            # Guard the resume environment against the checkpoint's recorded
            # state: a silently mismatched friction law or dropped apparent-MB
            # correction runs cleanly but produces wrong physics.
            if chk.has_attr("/", "friction"):
                chk_friction = str(chk.get_attr("/", "friction"))
                if chk_friction != friction:
                    raise RuntimeError(
                        f"Restart checkpoint {source_chk} was written with "
                        f"friction='{chk_friction}' but the environment "
                        f"resolves friction='{friction}'; set "
                        f"ISMIP7_FRICTION={chk_friction} to resume."
                    )
            amb_env = apparent_mb_mode()
            if (a_ref_mb is not None or phys_div is not None) and amb_env is None:
                # phys_div counts as the same evidence: adapt_mesh writes it
                # only in place of an a_ref_mb it found on the source, so an
                # adapted checkpoint that carries it came from an apparent-MB
                # run even though the a_ref itself was replaced.
                carried = "a_ref_mb" if a_ref_mb is not None else "phys_div"
                raise RuntimeError(
                    f"Restart checkpoint {source_chk} carries a frozen "
                    f"{carried} (the run used ISMIP7_APPARENT_MB) but "
                    f"ISMIP7_APPARENT_MB is off here; set it to resume with "
                    f"the same mass-balance correction."
                )
            _adapted_t0 = bool(chk.has_attr("/", "adapted_initial")
                               and int(chk.get_attr("/", "adapted_initial")))
            if a_ref_mb is None and amb_env is not None and phys_div is not None:
                PETSc.Sys.Print("  Apparent MB: adapted checkpoint, a_ref will be "
                                "rebuilt from the transferred physical divergence")
            elif a_ref_mb is None and amb_env is not None and _adapted_t0:
                # adapt_mesh.py --rebuild-aref: a t=0 state moved onto an
                # adapted mesh. Building a fresh a_ref here is legitimate
                # (nothing has evolved) and is the only way the correction can
                # cancel the NEW mesh's discrete divergence exactly.
                PETSc.Sys.Print("  Apparent MB: adapted t=0 state, a_ref will be "
                                "rebuilt on this mesh")
            elif a_ref_mb is None and amb_env is not None:
                is_timing_initial_state = (
                    checkpoint_metadata.get("timing_cache_role")
                    == "timing-initial-state"
                )
                if not (
                    allow_timing_cache_a_ref
                    and is_timing_initial_state
                    and amb_env == "div"
                ):
                    raise RuntimeError(
                        f"ISMIP7_APPARENT_MB is set but restart checkpoint "
                        f"{source_chk} has no a_ref_mb; a fresh a_ref cannot "
                        f"be built from an evolved state. Unset "
                        f"ISMIP7_APPARENT_MB or restart from a checkpoint "
                        f"that carries a_ref_mb."
                    )
                PETSc.Sys.Print(
                    "  Validated timing initial state: fresh div(h*u) "
                    "apparent-MB construction permitted"
                )
            # projection.sbatch parses t_yr out of this line to learn
            # the year the job STARTED at, which is how it tells a link that
            # advanced from one that spent its whole wall budget on setup.
            PETSc.Sys.Print(
                f"  Restart: evolved geometry + frozen anchors loaded "
                f"(t_yr={t_restart}, friction={friction})"
            )
        elif use_rc:
            # Cold RC/Budd normally uses the inversion-time MAP geometry. A
            # timing mesh override is different: interpolating a discontinuous
            # source DG0 field directly onto target DG0 samples one source cell
            # at each target centroid. The resulting aliasing is read as
            # driving stress because DG0 surface slope lives in facet jumps.
            # Retain the target-native BedMachine cell averages constructed
            # above; C_w0, N_ref, phi_eff, and H_init are then built from this
            # exact target geometry below. Continuous controls and u_obs still
            # come from the imported inversion.
            target_mesh_differs = (
                mesh_fn
                and geom_dg
                and (
                    not source_mesh_basename
                    or os.path.basename(mesh_fn) != source_mesh_basename
                )
            )
            if target_mesh_differs:
                PETSc.Sys.Print(
                    "  Target-mesh geometry: cell-averaged BedMachine "
                    f"({os.path.basename(geometry_source)})"
                )
            else:
                H = load_checkpoint_field(chk, "thickness", Q_g)
                b = load_checkpoint_field(chk, "bed", Q_g)
                s = load_checkpoint_field(chk, "surface", Q_g)
                geometry_source = os.path.realpath(source_chk)
                geometry_source_method = "checkpoint-native-v1"
            _uo = load_checkpoint_field(
                chk, "velocity_obs", V,
                fill=u_obs, fill_label="the raster-sampled velocity_obs",
            )
            u_obs.dat.data[:] = _uo.dat.data_ro
            if h_clamp_init > 0.0:
                H.interpolate(max_value(H, Constant(h_clamp_init)))
                s.interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))

    _filled = {k: v for k, v in transfer_fill.items()
               if v["missing"] or v["clamped"]}
    for _name, _info in _filled.items():
        PETSc.Sys.Print(
            f"  Transfer fill: {_name}: {_info['missing']} of {_info['total']} "
            f"dofs lie outside the source mesh; filled with {_info['fill']}; "
            f"{_info['clamped']} located dofs clamped to the source range "
            "(beyond it after strict location)"
        )
    if mesh_fn and not _filled:
        PETSc.Sys.Print(
            "  Transfer fill: none (every target dof lies inside the source mesh)"
        )

    # Clip the log-adjustments to a sane band. theta/phi are O(1) in a
    # converged MAP, so anything far beyond that is optimization noise from an
    # unconverged checkpoint (e.g. the rc_500 MAP snapshot carries ~700 nodes
    # with |.| up to 23 -> exp() ~1e10 local singularities the diagnostic SNES
    # cannot solve through). Default 6 for n=4 MAPs; the physical n=3 controls
    # legitimately reach ~8, so default 10 otherwise. Set 0 to disable.
    _n_flow_env = _n_flow()
    _map_clip_default = "6.0" if _n_flow_env == 4.0 else "10.0"
    map_clip = float(os.environ.get("ISMIP7_MAP_CLIP", _map_clip_default))
    if map_clip > 0.0:
        n_clip = 0
        for fld in (theta_f, phi_f):
            d = fld.dat.data
            n_clip += int((np.abs(d) > map_clip).sum())
            np.clip(d, -map_clip, map_clip, out=d)
        if n_clip:
            PETSc.Sys.Print(
                f"  MAP clip: bounded {n_clip} theta/phi node(s) to "
                f"|.|<={map_clip:.0f} (unconverged-checkpoint outliers)"
            )

    # Composite flow exponent (must match the inversion that produced the
    # MAP file we load above). This branch: n=3 standard Glen (A4_FACTOR=1).
    n_flow_val = _n_flow()
    m_slide_val = float(os.environ.get("ISMIP7_M_SLIDE", "3.0"))
    n_flow = Constant(n_flow_val)
    m_slide = Constant(m_slide_val)
    tau_c = Constant(0.1)
    # global_mean: .dat.data_ro.mean() is the rank-local owned slice (see
    # icepack2_tools/mpi_stats). Legacy action path only; build_rc_residual
    # anchors on C_w0.
    u_c = Constant(global_mean(Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    )))

    # Phi_eff (effective-pressure fraction). Uses a small floor on H so it
    # is well-defined where the original BedMachine thickness is 0. Loaded
    # (frozen) from the checkpoint on a restart; computed here on a cold start.
    _H_FLOOR_PHI = Constant(1.0)
    if phi_eff is None:
        phi_eff = Function(Q_g, name="phi_eff").interpolate(
            max_value(
                Constant(1.0)
                - rho_W * g * max_value(Constant(0.0), -b)
                  / (rho_I * g * max_value(H, _H_FLOOR_PHI)),
                Constant(0.01),
            )
        )
    # A4_base is the fluidity prior mean: the loaded thermomechanical A_prior
    # (phi = log(A/A_prior)), or the legacy constant A0*a4_factor for MAPs that
    # predate the physical prior. Must match the inversion that made the MAP.
    if A_prior_f is not None:
        A4_base = A_prior_f
        A_prior_lo, A_prior_hi = global_range(A_prior_f)
        PETSc.Sys.Print(
            f"  Fluidity prior A_prior loaded "
            f"[{A_prior_lo:.2f}, {A_prior_hi:.2f}]"
        )
        if not A_prior_lo > 0.0:
            raise RuntimeError(
                f"fluidity_prior minimum {A_prior_lo:g} is not positive after "
                "loading: A_eff = A_prior*exp(phi) vanishes there, which zeroes "
                "the dislocation term, the lin_reg regularizer and the alpha_gl "
                "collar in dual_friction.build_rc_residual at once (a singular "
                "membrane block). A transferred MAP fills the dofs outside its "
                "mesh with the constant baseline; a MAP carrying zeros is not "
                "usable as shipped."
            )
    else:
        A_prior_f = Function(Q, name="fluidity_prior").interpolate(
            Constant(A_prior_baseline)
        )
        A4_base = A_prior_f
        PETSc.Sys.Print(
            f"  Fluidity prior: checkpoint has no fluidity_prior; using LEGACY "
            f"constant baseline A0*a4_factor = {A_prior_baseline:.2f}"
        )
    A_map = A4_base * exp(phi_f)
    K_base = u_c / (phi_eff * tau_c) ** m_slide
    K_map = K_base * exp(-m_slide * theta_f)

    # Composite rheology: dislocation creep (n=n_flow, this branch n=3)
    # + α · linear regularizer (n=1) with a CONSTANT reference thickness
    # H_ref. This pins M where h → 0 so the SNES Jacobian stays
    # nonsingular at the calving front and h is allowed to reach zero.
    # Mirrors icepack2 dome_test.py.
    # RC MAP was inverted under alpha=1e-2 (SNES robustness across many
    # forward solves); keep the forward diagnostic consistent with it.
    alpha_reg = Constant(float(
        os.environ.get("ISMIP7_COMPOSITE_ALPHA", "1e-2" if use_rc else "1e-4")
    ))
    H_ref = Constant(float(os.environ.get("ISMIP7_H_REF", "100.0")))
    A_linear = A_map * tau_c ** (n_flow_val - 1)   # linearized at tau_c
    K_linear = u_c / (phi_eff * tau_c) * exp(-theta_f)  # linearized at tau_c

    # C_w0/N_ref were initialized (and, on a restart, loaded frozen) in the
    # reference-load block above; do NOT re-init to None here or a restart
    # would clobber the loaded anchors and pass None into build_rc_residual.
    if use_residual:
        from icepack2_tools.dual_friction import (
            build_rc_residual, weertman_anchor, effective_pressure,
        )
        c0_rc = float(os.environ.get("ISMIP7_RC_C0", "0.5"))
        rc_hvisc_floor = float(os.environ.get("ISMIP7_RC_HVISC_FLOOR", "10.0"))
        rc_cw0_floor = float(os.environ.get("ISMIP7_RC_CW0_FLOOR", "0.0"))
        rc_eps_tauc = float(os.environ.get("ISMIP7_RC_EPS_TAUC", "0.0"))
        # Budd N_hat knobs (used only for fric_law="budd"):
        #   PISM-delta grounded floor (Bueler & van Pelt 2015, ~0.02 of local
        #   overburden) removes the frictionless-GL degeneracy; N_hat cap
        #   (Joughin reduceNearGLBeta) bounds it; alpha_gl is the STRONG
        #   GL-gated viscosity coercivity that buys back the damping the
        #   exact-zero shelf removes (RC lacked it -> blew up). See Jul 2026
        #   RC-forward blow-up diagnosis + gia ase_model.momentum_F.
        budd_nhat_floor = float(os.environ.get("ISMIP7_BUDD_DELTA", "0.02"))
        budd_nhat_cap = float(os.environ.get("ISMIP7_BUDD_NHAT_CAP", "3.0"))
        # Floor-cell coercivity drag (gia ase_model ocean_drag): frictionless
        # ice-free buffer cells otherwise settle each diagnostic solve at a
        # velocity-runaway equilibrium (the Jul 2026 blow-up: ~2x/step outflux
        # growth, dt-independent, alpha_gl-immune). Linear drag ramping to
        # zero at h_ocean; ON by default. The knobs are owned by
        # runconfig.residual_stabilizers and the inversion passes the SAME
        # ones, so an inverted mixed state is a solution of this F too --
        # until 2026-09-14 it was not (||F|| 1e1 there vs 1e10 here), and
        # the restart fast path below trusted it anyway. The gia soft speed
        # limiter (u_lim/k_lim) is OFF by default: its max() kink at u_lim
        # breaks the nleqerr continuation (isolated Jul 18 2026); gia only
        # tolerates it under newtontr with dt-retry.
        _stabilizers = residual_stabilizers()
        ocean_drag = _stabilizers["ocean_drag"]
        h_ocean = _stabilizers["h_ocean"]
        # DG0 gate on the ocean drag, a live Function in the residual so the
        # gate changes without re-assembly. The drag never touches ice,
        # floating or grounded: it acts only in open water outside the t=0
        # extent that shares no facet with a cell holding ice now
        # (front.ocean_drag_cells, which has the measurement of the front it
        # used to pin). The t=0 extent is H_init on a restart and the
        # starting thickness on a cold start, so a restart applies the same
        # rule; run_simulation re-applies it after every transport advance,
        # and a level-set front, when one runs, writes its own gate instead.
        drag_mask = Function(FunctionSpace(mesh, "DG", 0), name="drag_mask")
        _Q0 = drag_mask.function_space()
        _drag_neighbours_of = facet_neighbours(_Q0)
        _drag_hmin = _front_hmin()
        _drag_extent0 = Function(_Q0).project(
            H if H_init is None else H_init).dat.data_ro >= _drag_hmin

        def _drag_rule(h_cells):
            r"""Write the ocean-drag gate for the per-cell thickness ``h_cells``."""
            drag_mask.dat.data[:] = ocean_drag_cells(
                h_cells >= _drag_hmin, _drag_neighbours_of, _drag_extent0)

        drag_rule = _drag_rule
        drag_rule(Function(_Q0).project(H).dat.data_ro)
        _n_drag = mesh.comm.allreduce(int(drag_mask.dat.data_ro.sum()))
        _n_cells = mesh.comm.allreduce(int(drag_mask.dat.data_ro.size))
        PETSc.Sys.Print(
            f"  Ocean drag gate: {_n_drag} of {_n_cells} cells, open water "
            f"outside the t=0 extent a cell away from any ice (h < {_drag_hmin:g} m)"
        )
        # Speed limiter: structurally present (threshold u_lim > 0) but INERT
        # by default - k_lim is a live Constant at 0 (term vanishes
        # identically; the cold continuation is unaffected, unlike a built-in
        # limiter, which breaks it). The run loop's rescue ladder raises
        # k_lim to ISMIP7_K_LIM for trust-region rescue solves at wall
        # geometries (runaway front nodes), then zeroes it again.
        u_lim = _stabilizers["u_lim"]
        k_lim = Constant(0.0)
        k_lim_rescue = float(os.environ.get("ISMIP7_K_LIM", "1e-3"))
        # GL-gated coercivity only for Budd (RC keeps its established form).
        alpha_gl = (float(os.environ.get("ISMIP7_ALPHA_GL", "0.5"))
                    if friction == "budd" else 0.0)
        # Anchor + reference effective pressure from the inversion-time
        # geometry (BEFORE any restart overwrites s): theta is a log-adjustment
        # on this C_w0, and N_ref pins N_hat=1 at the reference so the inverted
        # friction is reproduced at t=0. On a restart these frozen anchors are
        # loaded from the checkpoint above, never recomputed from evolved h.
        if not is_restart:
            # weertman_anchor needs |grad s|; under DG0 geometry it takes that
            # from a CG1 reconstruction internally (see geometry.surface_slope)
            # since a cell-wise surface has no cell gradient.
            C_w0 = weertman_anchor(H, s, u_obs, m_slide_val, Q_g)
            if friction == "budd":
                N_ref = Function(Q_g, name="N_ref").interpolate(
                    max_value(effective_pressure(H, s), Constant(0.0))
                )
        # ISMIP7_BUDD_NREF=none reproduces the inversion's own call, which
        # passes N_ref=None so N_hat = N/N is 1 wherever the gate is open. The
        # forward otherwise divides by the N_ref above, computed on a cold start
        # from its own geometry via effective_pressure(H, s) and loaded from the
        # checkpoint only on a restart. The two should agree cell-wise under DG0
        # geometry; this knob is how that is checked.
        if friction == "budd" and os.environ.get("ISMIP7_BUDD_NREF", "").lower() == "none":
            N_ref = None
        if friction == "budd":
            PETSc.Sys.Print(
                f"  Friction: Budd N_hat (N_ref={'none' if N_ref is None else 'reference'}; exact-zero shelf; delta="
                f"{budd_nhat_floor:.3f}, N_hat_cap={budd_nhat_cap:.1f}, "
                f"alpha_gl={alpha_gl:.2f}, h_visc_floor={rc_hvisc_floor:.0f}m, "
                f"ocean_drag={ocean_drag:.0e}@h<{h_ocean:.0f}m, "
                f"u_lim={u_lim:.0e})"
            )
        else:
            PETSc.Sys.Print(
                f"  Friction: regularized Coulomb (c0={c0_rc}, "
                f"h_visc_floor={rc_hvisc_floor:.0f}m, cw0_floor={rc_cw0_floor:.1e}, "
                f"eps_tauc={rc_eps_tauc:.1e} MPa, alpha={float(alpha_reg):.1e})"
            )

    # Issue #117: a cold start dated before the 2015 geometry starts from that
    # geometry with the Smith et al. (2020) mean thinning undone on grounded
    # ice. It comes after the anchors above, so C_w0 and N_ref stay those of
    # the 2015 geometry the MAP was inverted on, and before the prognostic
    # thickness is copied from H below, so the initial solve, the transport,
    # H_init and the apparent-MB reference all start from the earlier ice.
    if backdate_years > 0.0 and not is_restart:
        if not geom_dg:
            raise RuntimeError(
                "backdating the geometry (issue #117) needs the DG0 geometry "
                "(ISMIP7_GEOMETRY_SPACE=dg0)")
        from icepack2_tools.obs_dhdt import load_dhdt_obs, backdate_thickness
        _dhdt, _observed = load_dhdt_obs(Q_g)
        _h_obs = H.dat.data_ro.copy()
        _h_new, _changed = backdate_thickness(
            _h_obs, b.dat.data_ro, _dhdt.dat.data_ro, _observed.dat.data_ro,
            backdate_years, float(rho_ratio))
        H.dat.data[:] = _h_new
        s.interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))
        _area = assemble(fd.TestFunction(Q_g) * dx).dat.data_ro
        _dm_gt = mesh.comm.allreduce(
            float(((_h_new - _h_obs) * _area).sum())) * 917.0 / 1e12
        PETSc.Sys.Print(
            f"  Geometry backdated {backdate_years:g} yr (issue #117): "
            f"{mesh.comm.allreduce(int(_changed.sum()))} grounded cells with "
            f"Smith dH/dt coverage, ice mass {_dm_gt:+.0f} Gt; floating ice and "
            f"unobserved cells keep their {GEOMETRY_YEAR:g} thickness; friction "
            f"anchors from the {GEOMETRY_YEAR:g} geometry")

    # ISMIP7_SNES_TYPE=newtontr switches the diagnostic Newton to trust
    # region (gia COUPLED_SOLVER's choice: more robust than line search at
    # stiff melt-driven GL-retreat geometries, where nleqerr hit walls
    # ~9 yr into the 32 km ssp585 run). Line search stays the default.
    # gia: hard-era Budd steps converge LINEARLY (~2%/iter under active trust
    # region) and were being executed by the cap while still descending -
    # patience beats retries. 200 suffices for newtonls eras; raise via env for
    # newtontr pushes through hard geometry.
    sparams = diagnostic_solver_parameters()
    linear_solver = diagnostic_solver_mode()
    PETSc.Sys.Print(
        f"  Linear solver: {linear_solver} "
        f"({diagnostic_solver_label(linear_solver)})"
    )
    if linear_solver == "scpc_gamg":
        # A tuning rung differs from the next only in these; name them in
        # the job's own log, which outlives an overwritten record.
        PETSc.Sys.Print("  Condensed solve: " + " ".join(
            f"{key[len('condensed_field_'):]}={value}"
            for key, value in sparams.items()
            if key.startswith("condensed_field_") and value is not None
        ))
    # Optional SNES/KSP convergence monitoring (ISMIP7_SNES_MONITOR=1).
    # ISMIP7_SNES_LOG routes the output to a file (per run, so concurrent
    # debug runs don't interleave); otherwise it goes to stdout.
    _solver_log = None
    _solver_monitor = snes_monitor_enabled()
    _solver_view = solver_view_enabled()
    if _solver_monitor or _solver_view:
        _solver_log = os.environ.get("ISMIP7_SNES_LOG")
        if _solver_log:
            os.makedirs(
                os.path.dirname(os.path.abspath(_solver_log)), exist_ok=True
            )
        # The fourth ASCII-viewer field is the PETSc file mode.  Append mode
        # keeps our per-solve headers when a monitor opens the same file.
        _viewer = f"ascii:{_solver_log}::append" if _solver_log else None
        if _solver_monitor:
            sparams.update({
                "snes_monitor": _viewer,
                "snes_converged_reason": _viewer,
                "snes_linesearch_monitor": _viewer,
                "ksp_monitor_short": _viewer,
                "ksp_monitor_true_residual": _viewer,
                "ksp_converged_reason": _viewer,
            })
        if _solver_monitor and linear_solver.startswith("scpc_"):
            # The top-level FGMRES is still useful, but SCPC's condensed KSP
            # identifies whether the velocity solve or the exact local
            # elimination is responsible for a failure.
            sparams.update({
                "condensed_field_ksp_monitor_short": _viewer,
                "condensed_field_ksp_monitor_true_residual": _viewer,
                "condensed_field_ksp_converged_reason": _viewer,
            })
        if _solver_view:
            sparams.update({
                "snes_view": _viewer,
                "ksp_view": _viewer,
            })
            if linear_solver.startswith("scpc_"):
                sparams["condensed_field_ksp_view"] = _viewer
        PETSc.Sys.Print(
            f"  Solver diagnostics log: {_solver_log}" if _solver_log
            else "  Solver diagnostics: stdout"
        )
    fc_params = {"quadrature_degree": 4}

    z = Function(Z)
    z.sub(0).interpolate(Constant(0.1) * u_obs)
    if u_guess is not None:
        # Warm start: seed the diagnostic solve with the checkpoint velocity
        # (H/s/phi_eff/anchors were already restored in the reference block).
        z.sub(0).dat.data[:] = u_guess.dat.data_ro
        if M_guess is not None:
            z.sub(1).dat.data[:] = M_guess.dat.data_ro
            z.sub(2).dat.data[:] = tau_guess.dat.data_ro

    h = H.copy(deepcopy=True)
    h.rename("thickness")

    u_s, M_s, tau_s = split(z)
    fields = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": h,
        "surface": s,
    }
    rheo_glen = {
        "flow_law_exponent": n_flow,
        "flow_law_coefficient": A_map,
        "sliding_exponent": m_slide,
        "sliding_coefficient": K_map,
    }
    rheo_linear = {
        "flow_law_exponent": Constant(1.0),
        "flow_law_coefficient": A_linear,
        "sliding_exponent": Constant(1.0),
        "sliding_coefficient": K_linear,
    }
    # Composite-rheology fields: use the constant reference thickness for the
    # linear regularization terms so they stay positive-definite at h=0.
    fields_reg = dict(fields)
    fields_reg["thickness"] = H_ref

    if use_residual:
        # Residual closure on the LIVE prognostic fields (h, s): the grounded
        # gate, effective pressure, and driving stress all track the evolving
        # geometry, so the grounding line migrates freely with exact-zero
        # shelf drag. N_ref (Budd) is the fixed reference effective pressure.
        def _build_F(theta_c=None, phi_c=None, h_c=None, s_c=None, z_c=None):
            r"""The dual residual for a given set of controls and geometry.

            ``setup_model`` builds one residual on its own fields; a
            time-dependent assimilation needs a residual per window step,
            because each step has its own thickness and surface and an
            adjoint has to see them as distinct variables rather than as one
            mutated Function. Every argument defaults to this context's own
            field, so ``_build_F()`` reproduces the solver's residual
            exactly. Used by ``scripts/inversion_td``."""
            return build_rc_residual(
                z_c if z_c is not None else z,
                theta_c if theta_c is not None else theta_f,
                phi_c if phi_c is not None else phi_f,
                H=h_c if h_c is not None else h,
                s=s_c if s_c is not None else s,
                b=b, C_w0=C_w0,
                A4_base=A4_base, n_flow=n_flow, n_flow_val=n_flow_val,
                m_slide=m_slide_val, tau_c=tau_c, alpha=alpha_reg, H_ref=H_ref,
                fric_law=friction, N_ref=N_ref,
                nhat_floor=budd_nhat_floor, nhat_cap=budd_nhat_cap,
                alpha_gl=alpha_gl,
                c0=c0_rc, eps_tauc=rc_eps_tauc,
                c_w0_floor=rc_cw0_floor, h_visc_floor=rc_hvisc_floor,
                ocean_drag=ocean_drag, h_ocean=h_ocean, u_lim=u_lim,
                k_lim=k_lim, drag_mask=drag_mask,
                calving_ids=calving_ids if use_calving_terminus else None,
            )

        # The closure above is the single definition of this residual: the
        # forward solves exactly what a time-dependent assimilation rebuilds.
        F = _build_F()
    else:
        L = (
            model.minimization.viscous_power(**fields, **rheo_glen)
            + alpha_reg * model.minimization.viscous_power(**fields_reg, **rheo_linear)
            + model.minimization.friction_power(**fields, **rheo_glen)
            + alpha_reg * model.minimization.friction_power(**fields, **rheo_linear)
            + model.minimization.momentum_balance(**fields)
        )
        if use_calving_terminus:
            L += model.minimization.calving_terminus(**fields, outflow_ids=calving_ids)
        F = derivative(L, z)

    if linear_solver.startswith("scpc_"):
        # Firedrake's three-field SCPC expects both off-diagonal entries of the
        # eliminated (M, tau) block to be present in split_form.  These fields
        # are physically uncoupled, so UFL otherwise omits both structural-zero
        # blocks and SCPC raises KeyError before assembly.  A runtime Constant
        # preserves the block metadata while contributing exactly zero to the
        # residual and Jacobian.  Do not replace it with the literal 0: UFL
        # simplifies that away and recreates the missing-block failure.
        scpc_structural_zero = Constant(0.0)
        F += derivative(
            scpc_structural_zero * M_s[0, 0] * tau_s[0] * dx, z
        )

    # A matrix-free Jacobian follows the state Function, which the line
    # search's residual evaluations overwrite with its trial point; hold it at
    # the Newton iterate SCPC assembled its condensed system at.
    jacobian, pre_jacobian = None, None
    linearization = linearization_state(linear_solver)
    if linearization == "frozen":
        from icepack2_tools.preconditioners import frozen_linearization
        jacobian, pre_jacobian = frozen_linearization(F, z)
    PETSc.Sys.Print(f"  Jacobian linearization state: {linearization}")

    prob = NonlinearVariationalProblem(
        F, z, J=jacobian, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(
        prob,
        solver_parameters=sparams,
        options_prefix="ismip7_diagnostic_",
        pre_jacobian_callback=pre_jacobian,
    )

    _solve_count = 0
    solver_stats = []
    _header_viewer = None
    if _solver_log and mesh.comm.rank == 0:
        _header_viewer = PETSc.Viewer().createASCII(
            _solver_log,
            mode=PETSc.Viewer.FileMode.APPEND,
            comm=PETSc.COMM_SELF,
        )

    def _write_solve_header(label, **metadata):
        nonlocal _solve_count
        _solve_count += 1
        details = " ".join(
            f"{key}={value}" for key, value in metadata.items()
        )
        header = (
            f"\n=== DIAGNOSTIC SOLVE {_solve_count:04d} | {label} | "
            f"n_flow={float(n_flow):.6g} m_slide={float(m_slide):.6g} "
            f"linear={linear_solver} snes={slvr.snes.getType()} "
            f"ksp={sparams['ksp_type']} "
            f"ksp_rtol={sparams.get('ksp_rtol', 'direct')} "
            f"ksp_max_it={sparams.get('ksp_max_it', 'direct')}"
            f"{(' ' + details) if details else ''} ===\n"
        )
        if _solver_log:
            if mesh.comm.rank == 0:
                _header_viewer.printfASCII(header)
                _header_viewer.flush()
            mesh.comm.barrier()
        else:
            PETSc.Sys.Print(header.rstrip())

    def _condensed_work():
        """``(solves, iterations)`` SCPC has made on the condensed system so
        far, or None under a solver without one. Unlike SNES's
        ``linear_iterations`` it includes the line search's solves."""
        if not linear_solver.startswith("scpc_"):
            return None
        pc = slvr.snes.ksp.pc
        context = pc.getPythonContext() if pc.getType() == "python" else None
        # Zero until the first solve has built the PC.
        return (
            getattr(context, "condensed_solves", 0),
            getattr(context, "condensed_iterations", 0),
        )

    def solve_diagnostic(label, **metadata):
        """Execute one solve and always emit one compact convergence record."""
        _write_solve_header(label, **metadata)
        work_before = _condensed_work()
        t0_solve = perf_counter()
        try:
            return slvr.solve()
        finally:
            elapsed = perf_counter() - t0_solve
            reason = slvr.snes.getConvergedReason()
            reason_name = getattr(reason, "name", str(int(reason)))
            stat = {
                "solve": _solve_count,
                "label": label,
                "snes_reason": reason_name,
                "snes_iterations": slvr.snes.getIterationNumber(),
                "linear_iterations": slvr.snes.getLinearSolveIterations(),
                "function_norm": slvr.snes.getFunctionNorm(),
                "seconds": elapsed,
            }
            condensed = ""
            if work_before is not None:
                solves, iterations = _condensed_work()
                stat["condensed_solves"] = solves - work_before[0]
                stat["condensed_iterations"] = iterations - work_before[1]
                condensed = (
                    f"condensed_solves={stat['condensed_solves']} "
                    f"condensed_its={stat['condensed_iterations']} "
                )
            solver_stats.append(stat)
            PETSc.Sys.Print(
                "=== DIAGNOSTIC RESULT "
                f"{_solve_count:04d} | {label} | reason={reason_name} "
                f"snes_its={stat['snes_iterations']} "
                f"linear_its={stat['linear_iterations']} "
                f"{condensed}"
                f"fnorm={stat['function_norm']:.6e} "
                f"seconds={elapsed:.3f} ==="
            )

    # Adaptive n/m continuation for the cold-start diagnostic solve. On a
    # fine mesh with a rough (mid-optimization) MAP the n=1→n_flow_val jump
    # can outrun Newton (DIVERGED_MAX_IT); restore the initial guess and
    # re-ramp with more, smaller steps rather than crashing. Escalates
    # ISMIP7_CONTINUATION_STEPS (default 8) → 2× → 4×.
    base_steps = continuation_steps()
    z_init = z.copy(deepcopy=True)

    # Restart fast path: a checkpoint carrying the full mixed state (u, M,
    # tau) is trusted WITHOUT a setup solve only when its residual under THIS
    # F is at the level its writer recorded (full_state_residual, stamped by
    # the inversion and by save_model_state):
    #     ||F(z_loaded)|| <= snes_atol_scale x recorded.
    # Solving an unchanged converged state again is redundant and pathological
    # at the residual floor (the 2015.2 qualification cache ran 200 Newton
    # iterations at ||F||=4e-5), which is why the fast path exists. Until
    # 2026-09-14 it accepted ANY finite residual and scaled the run atol to
    # it: the forward's own reloaded prepare caches (||F||~1e7-1e8) and
    # inversion MAPs solved without the ocean_drag term (1.5e10) were
    # "accepted" with run atol 1e1-1e4, and every strict scout ran away
    # within four steps. A state above the limit, or without a record, is
    # re-solved from the loaded guess with the step-size exit live and the
    # iteration count bounded (final_solve_bounds: a floor-level start exits
    # at iteration 0-1, a far one is an ordinary Newton solve), and the run
    # atol then follows the achieved norm exactly as after a cold
    # continuation. Velocity-only checkpoints take the same re-solve path
    # because (u, 0, 0) is not a cached mixed solution.
    #
    # The accepted-state run tolerance stays TIGHT (1e-6 x ||F(z_loaded)||);
    # never use the acceptance residual itself as a loose persistent atol:
    # that let later steps "converge" at iteration zero and silently freeze
    # the velocity (the bug that invalidated the first 1873->2014 resume).
    restart_solved = False
    if is_restart and u_guess is not None:
        n_flow.assign(n_flow_val)
        m_slide.assign(m_slide_val)
        with assemble(F, form_compiler_parameters=fc_params).dat.vec_ro as _rv:
            fnorm0 = _rv.norm()
        if not np.isfinite(fnorm0):
            raise RuntimeError(
                "Restart checkpoint has a non-finite full-state residual"
            )
        restart_atol_scale = snes_restart_failure_atol_scale()
        restart_atol = restart_atol_scale * fnorm0
        recorded = checkpoint_metadata.get("full_state_residual")
        try:
            recorded = float(recorded) if recorded is not None else None
        except (TypeError, ValueError):
            recorded = None
        if recorded is not None and not (np.isfinite(recorded) and recorded > 0.0):
            recorded = None
        accept_scale = snes_atol_scale()
        accept_limit = accept_scale * recorded if recorded is not None else None
        have_mixed = M_guess is not None and tau_guess is not None

        if have_mixed and accept_limit is not None and fnorm0 <= accept_limit:
            restart_solved = True
            if fnorm0 > 0.0:
                slvr.snes.setTolerances(atol=restart_atol)
            PETSc.Sys.Print(
                "Restart full mixed state accepted without a setup solve "
                f"(||F||={fnorm0:.2e} <= {accept_limit:.2e} = {accept_scale:g} x "
                f"recorded {recorded:.2e}; run atol={restart_atol:.2e})"
            )
        else:
            if not have_mixed:
                why = "velocity-only checkpoint"
            elif accept_limit is None:
                why = "checkpoint records no full_state_residual"
            else:
                why = (
                    f"||F||={fnorm0:.2e} exceeds {accept_limit:.2e} = "
                    f"{accept_scale:g} x recorded {recorded:.2e}"
                )
            PETSc.Sys.Print(f"Restart state must be re-solved: {why}")
            _rtol0, _atol0, _stol0, _max_it0 = slvr.snes.getTolerances()
            _bounds = final_solve_bounds()
            slvr.snes.setTolerances(
                atol=accept_limit if accept_limit is not None else _atol0,
                stol=_bounds["snes_stol"],
                max_it=_bounds["snes_max_it"],
            )
            try:
                solve_diagnostic(
                    "restart-loaded-state",
                    loaded_fnorm=f"{fnorm0:.3e}",
                    recorded=("none" if recorded is None else f"{recorded:.3e}"),
                )
                restart_solved = True
                fnorm_conv = slvr.snes.getFunctionNorm()
                run_atol = (
                    snes_atol_scale() * fnorm_conv if fnorm_conv > 0.0 else _atol0
                )
                slvr.snes.setTolerances(atol=run_atol, stol=_stol0, max_it=_max_it0)
                PETSc.Sys.Print(
                    f"Restart state re-solved (||F|| {fnorm0:.2e} -> "
                    f"{fnorm_conv:.2e}; run atol={run_atol:.2e})"
                )
            except fd.ConvergenceError:
                z.assign(z_init)
                restart_solved = True
                slvr.snes.setTolerances(
                    atol=restart_atol if fnorm0 > 0.0 else _atol0,
                    stol=_stol0,
                    max_it=_max_it0,
                )
                PETSc.Sys.Print(
                    "Restart state re-solve did not converge; keeping the "
                    "loaded state and handing the step to the rescue ladder "
                    f"(atol={restart_atol:.2e})"
                )

    if not restart_solved:
        PETSc.Sys.Print(
            f"Initial diagnostic solve (continuation n_flow 1→{n_flow_val:.1f}, "
            f"m_slide 1→{m_slide_val:.1f})..."
        )
        _run_continuation = True
    else:
        _run_continuation = False
    for attempt, steps in enumerate(
        (base_steps, 2 * base_steps, 4 * base_steps) if _run_continuation else ()
    ):
        try:
            for step, t in enumerate(np.linspace(0.0, 1.0, steps), 1):
                n_flow.assign(1.0 + t * (n_flow_val - 1.0))
                m_slide.assign(1.0 + t * (m_slide_val - 1.0))
                solve_diagnostic(
                    "initial-continuation",
                    attempt=attempt + 1,
                    step=f"{step}/{steps}",
                    t=f"{t:.6g}",
                )
            PETSc.Sys.Print(f"  Done ({steps} continuation steps)")
            # Self-scaled absolute tolerance: a solve that STARTS at the
            # converged state (restart step 1: geometry unchanged since this
            # continuation) has ||F|| at the rounding floor already, and the
            # nleqerr linesearch then fails on a residual it cannot reduce.
            # Accepting anything within 100x of the achieved converged norm
            # makes such solves report converged at iteration 0. Geometry
            # changes during stepping push ||F_0|| far above this, so the
            # usual rtol path is untouched.
            fnorm_conv = slvr.snes.getFunctionNorm()
            if fnorm_conv > 0.0:
                atol_scale = snes_atol_scale()
                slvr.snes.setTolerances(atol=atol_scale * fnorm_conv)
                PETSc.Sys.Print(
                    f"  snes_atol set to {atol_scale * fnorm_conv:.2e} "
                    f"({atol_scale:g}x converged residual norm)"
                )
            break
        except fd.ConvergenceError:
            if attempt == 2:
                PETSc.Sys.Print(
                    f"  Continuation diverged at {steps} steps - giving up."
                )
                raise
            PETSc.Sys.Print(
                f"  Continuation diverged at {steps} steps; "
                f"restarting with {2 * steps}..."
            )
            z.assign(z_init)
            n_flow.assign(1.0)
            m_slide.assign(1.0)

    u0 = z.subfunctions[0]
    area = assemble(Constant(1.0) * dx(mesh))
    misfit0 = float(assemble(
        0.5 / area * ((u0[0] - u_obs[0]) ** 2 + (u0[1] - u_obs[1]) ** 2) * dx
    ))
    PETSc.Sys.Print(f"  Initial velocity misfit vs obs: {misfit0:.6e}")

    # Forcing fields live on the GEOMETRY space: the melt parameterization is
    # evaluated from the local draft (s - h) and must be cell-wise wherever the
    # geometry is, and the transport cell-averages (accum - ocean_melt) anyway.
    accum = Function(Q_g, name="accumulation").assign(0.0)
    ocean_melt = Function(Q_g, name="ocean_melt").assign(0.0)

    # Coordinates of the geometry dofs, for forcing callbacks that assign into
    # .dat.data directly (they cannot assume mesh vertices any more).
    _xy = Function(VectorFunctionSpace(mesh, Q_g.ufl_element())).interpolate(
        fd.SpatialCoordinate(mesh)
    ).dat.data_ro
    geom_xy = (_xy[:, 0].copy(), _xy[:, 1].copy())

    # t=0 fixed-front anchor: on a cold start it is the initial (BedMachine/
    # inversion) thickness; on a restart it was loaded from the checkpoint, so
    # it stays the ORIGINAL observed extent rather than the evolved geometry.
    if H_init is None:
        H_init = Function(Q_g, name="H_init")
        H_init.assign(H)

    return {
        "mesh": mesh,
        "Q": Q,
        "Q_g": Q_g,
        "geom_dg": geom_dg,
        "geom_xy": geom_xy,
        "mesh_basename": mesh_basename,
        "geometry_source": geometry_source,
        "geometry_source_method": geometry_source_method,
        "transfer_fill": transfer_fill,
        "initial_misfit": misfit0,
        "V": V,
        "Z": Z,
        "z": z,
        "h": h,
        "s": s,
        "b": b,
        "slvr": slvr,
        # The momentum residual and its form-compiler parameters, so a state
        # checkpoint can record ||F|| under the F its readers will assemble.
        "F": F,
        "fc_params": fc_params,
        "solve_diagnostic": solve_diagnostic,
        "solver_stats": solver_stats,
        "n_flow": n_flow,
        "n_flow_val": n_flow_val,
        "m_slide": m_slide,
        "m_slide_val": m_slide_val,
        "accum": accum,
        "ocean_melt": ocean_melt,
        "phi_eff": phi_eff,
        "rho_ratio": rho_ratio,
        "h_clamp": h_clamp,
        # The floor a cold start put under the initial thickness (0 on a
        # restart, which reads its geometry). A floor turns ice-free cells
        # into floating ice the melt calibration never fitted, so the melt
        # contract refuses it (forcing.check_melt_contract).
        "thickness_floor": 0.0 if is_restart else h_clamp_init,
        "calving_ids": calving_ids,
        "u_obs": u_obs,
        "friction": friction,
        # Level-set calving front support: the drag gate the residual was
        # built with. The front itself is not carried here: run_simulation
        # rebuilds it from the current thickness every step, and the `fixed`
        # law anchors on H_init below, so a restart needs no saved front.
        "drag_mask": drag_mask,
        # Rewrites drag_mask from a per-cell thickness under the ocean-drag
        # rule (front.ocean_drag_cells); None without a residual.
        "drag_rule": drag_rule,
        # Residual builder for a time-dependent assimilation (None for the
        # legacy action formulation, which has no residual to rebuild).
        "build_F": _build_F if use_residual else None,
        # Solver options the time-dependent assimilation
        # (antarctica/scripts/inversion_td, not yet committed) reuses so its
        # forward solves match this one, alongside build_F and ISMIP7_INVERSION.
        "sparams": sparams,
        "A_map": A_map,
        # Mesh provenance from the source checkpoint (re-stamped into every
        # state checkpoint so warm restarts stay self-describing).
        "lc": chk_lc,
        "lc_coarse": chk_lc_coarse,
        "buffer_m": chk_buffer_m,
        "raster_sample": chk_raster_sample,
        # Rescue speed limiter (residual laws): live Constant, 0 = inert.
        "k_lim": k_lim if use_residual else None,
        "k_lim_rescue": k_lim_rescue if use_residual else 0.0,
        # Reference/frozen fields persisted into every checkpoint so a restart
        # is self-contained and rank-count-robust (no recompute from evolved h).
        "theta": theta_f,
        "phi": phi_f,
        "C_w0": C_w0,
        "N_ref": N_ref,
        "A_prior": A_prior_f,
        "H_init": H_init,
        # Frozen t=0 apparent-MB correction (restart only; else None).
        "a_ref_mb": a_ref_mb,
        "phys_div": phys_div,
        # DG0 prognostic thickness state (restart only; else None).
        "h_dg_state": h_dg_state,
        # Resume time (None on a cold start); run_simulation continues the
        # timeline from here instead of the caller's t_start.
        "t_restart": t_restart,
        # Timing-cache identity, if this is a prepared timing restart.
        # Ordinary production checkpoints legitimately omit these fields.
        "checkpoint_metadata": checkpoint_metadata,
        # The partly accumulated ISMIP7 year carried by the restart (or None).
        "ismip7_resume": ismip7_resume,
    }


class LiveCalvingState:
    r"""The fields a calving law reads, taken live from the forward's state.

    hoffmaao/calving's ``laws.Law.rate(model, t)`` reads seven fields from
    its model: the dual solution ``u``, ``M``, ``tau``, the DG0 thickness
    ``h``, the height above flotation ``haf`` and the grounded indicator
    ``chi_gr`` on the cells, and the outward front normal ``nfront``. Here
    they are the forward's own: ``(u, M, tau)`` are the subfunctions of the
    mixed solution, ``haf`` and ``chi_gr`` are UFL on the cells so they
    follow the geometry without an update, and ``nfront`` is the unit
    gradient of the level set the forward advances, the same object
    ``calving/antarctic.py`` builds when a law is tuned against the Greene
    fronts, so the tuned threshold means the same thing here.

    Densities follow that tuning harness (CalvingMIP's 917 / 1028) rather
    than the forward's 1024, as ``antarctic.AntarcticState`` does: a
    threshold fitted there is applied under the same flotation test.
    """
    RHO_I = 917.0
    RHO_W = 1028.0

    def __init__(self, z, h_dg, b, level_set):
        from firedrake import conditional, gt
        self.u, self.M, self.tau = z.subfunctions
        self.h = h_dg
        self.b = b
        self.Q0 = h_dg.function_space()
        self.haf = self.h - Constant(self.RHO_W / self.RHO_I) * max_value(
            -self.b, Constant(0.0))
        self.chi_gr = conditional(gt(self.haf, Constant(0.0)),
                                  Constant(1.0), Constant(0.0))
        self.levelset = level_set
        self.nfront = level_set.ghat
        self.front_len = level_set.front_len


def save_model_state(ctx, final_path, t_now, extra_attrs=None):
    r"""Atomically save one self-contained mixed state.

    Ordinary simulation checkpoints and timing caches share this writer so a
    cache cannot silently omit one of the frozen fields required on restart.
    """
    mesh = ctx["mesh"]
    z = ctx["z"]
    h = ctx["h"]
    # Residual of the state being written, under the same F a restart will
    # assemble: the restart fast path accepts the mixed state without a solve
    # only within snes_atol_scale x this value. Measured at full n/m, which
    # is where every writer calls this (after the cold continuation, or after
    # a step's converged diagnostic solve).
    full_state_residual = None
    if ctx.get("F") is not None:
        with assemble(
            ctx["F"], form_compiler_parameters=ctx.get("fc_params")
        ).dat.vec_ro as _rv:
            full_state_residual = float(_rv.norm())
    tmp = final_path + ".tmp"
    with fd.CheckpointFile(tmp, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(ctx["theta"], name="log_friction")
        chk.save_function(ctx["phi"], name="log_fluidity")
        chk.save_function(ctx["u_obs"], name="velocity_obs")
        chk.save_function(ctx["b"], name="bed")
        chk.save_function(h, name="thickness")
        chk.save_function(ctx["s"], name="surface")
        chk.save_function(z.subfunctions[0], name="velocity")
        chk.save_function(z.subfunctions[1], name="membrane_stress")
        chk.save_function(z.subfunctions[2], name="basal_stress")
        chk.save_function(ctx.get("H_init", h), name="H_init")
        chk.save_function(ctx["phi_eff"], name="phi_eff")
        if ctx.get("C_w0") is not None:
            chk.save_function(ctx["C_w0"], name="C_w0")
        if ctx.get("N_ref") is not None:
            chk.save_function(ctx["N_ref"], name="N_ref")
        if ctx.get("A_prior") is not None:
            chk.save_function(ctx["A_prior"], name="fluidity_prior")
        if ctx.get("a_ref_mb") is not None:
            chk.save_function(ctx["a_ref_mb"], name="a_ref_mb")
        if ctx.get("level_set") is not None:
            chk.save_function(ctx["level_set"].phi, name="levelset")
        if not ctx.get("geom_dg", False):
            h_dg = ctx.get("h_dg_state")
            if h_dg is None:
                raise RuntimeError(
                    "CG1 state checkpoint requested before h_dg was prepared"
                )
            chk.save_function(h_dg, name="thickness_dg")
        # The ISMIP7 year in progress (ismip7_output.AnnualOutput), so a
        # chained resume continues the same year instead of losing it.
        annual = ctx.get("annual")
        if annual is not None:
            for _name, _f in annual.state_fields().items():
                chk.save_function(_f, name=_name)
            for _key, _val in annual.state_attrs().items():
                chk.set_attr("/", _key, _val)

        chk.set_attr("/", "t_yr", float(t_now))
        if ctx.get("calving_law") is not None:
            chk.set_attr("/", "calving_law", str(ctx["calving_law"].describe()))
        chk.set_attr("/", "friction", str(ctx.get("friction", "budd")))
        if str(ctx.get("friction", "budd")) == "budd":
            # Provenance of the shelf gate this state was solved under
            # (runconfig.BUDD_SHELF_GATE); timing-cache manifests require it.
            chk.set_attr("/", "friction_gate", _BUDD_SHELF_GATE)
        chk.set_attr(
            "/", "geometry_space", "dg0" if ctx.get("geom_dg") else "cg1"
        )
        if ctx.get("mesh_basename"):
            chk.set_attr("/", "mesh_basename", str(ctx["mesh_basename"]))
        for name in ("geometry_source", "geometry_source_method"):
            if ctx.get(name):
                chk.set_attr("/", name, str(ctx[name]))
        for name in ("lc", "lc_coarse"):
            if ctx.get(name) is not None:
                chk.set_attr("/", name, int(ctx[name]))
        if ctx.get("buffer_m") is not None:
            chk.set_attr("/", "buffer_m", float(ctx["buffer_m"]))
        if ctx.get("raster_sample"):
            chk.set_attr("/", "raster_sample", str(ctx["raster_sample"]))
        if full_state_residual is not None:
            chk.set_attr("/", "full_state_residual", full_state_residual)
        for name, value in (extra_attrs or {}).items():
            if value is not None:
                chk.set_attr("/", name, value)

    mesh.comm.barrier()
    if mesh.comm.rank == 0:
        os.replace(tmp, final_path)
    mesh.comm.barrier()


def _global_argmax_with_payload(values, xy, payload, comm):
    r"""Global argmax of ``values`` over the owned cells of every rank, with
    the cell centre and the ``payload`` columns (name -> per-cell array) at
    that cell. Non-finite values are ignored, one candidate per rank is
    communicated and ties resolve by rank, so every rank sees the same
    answer. Returns ``(value, (x, y), {name: value})``; the location and
    fields are empty when no rank holds a finite value."""
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if finite.any():
        index = int(np.argmax(np.where(finite, values, -np.inf)))
        candidate = (
            float(values[index]),
            int(comm.rank),
            tuple(float(v) for v in np.asarray(xy)[index]),
            {name: float(np.asarray(col)[index]) for name, col in payload.items()},
        )
    else:
        candidate = (-np.inf, int(comm.rank), (), {})
    candidates = comm.allgather(candidate)
    value, _rank, location, fields = max(
        candidates, key=lambda item: (item[0], -item[1])
    )
    return value, location, fields


def run_simulation(
    ctx,
    experiment_name,
    t_start,
    t_end,
    dt=1.0,
    output_interval=10,
    checkpoint_interval=100,
    forcing_callback=None,
):
    r"""Run the split diagnostic-prognostic time-stepping loop."""
    mesh = ctx["mesh"]
    Q = ctx["Q"]
    Q_g = ctx.get("Q_g", Q)
    z = ctx["z"]
    h = ctx["h"]
    s = ctx["s"]
    b = ctx["b"]
    slvr = ctx["slvr"]
    solve_diagnostic = ctx["solve_diagnostic"]
    n_flow = ctx["n_flow"]
    n_flow_val = ctx["n_flow_val"]
    m_slide = ctx["m_slide"]
    m_slide_val = ctx["m_slide_val"]
    accum = ctx["accum"]
    ocean_melt = ctx["ocean_melt"]
    phi_eff = ctx["phi_eff"]
    rho_ratio = ctx["rho_ratio"]
    h_clamp = ctx["h_clamp"]
    # Warm restart: continue the timeline from the checkpoint's saved year so
    # time-varying forcing (SSP projections) is applied at the correct year.
    t_restart = ctx.get("t_restart")
    if t_restart is not None:
        PETSc.Sys.Print(
            f"  Resuming timeline at t={t_restart:.2f} "
            f"(caller t_start={t_start} overridden)"
        )
        t_start = t_restart

    # round, don't truncate: int((2300-2015)/0.1) = 2849 loses the last step
    nsteps = int(round((t_end - t_start) / dt))
    dt_c = Constant(dt)
    PETSc.Sys.Print(
        f"\nTime-stepping: {t_start}->{t_end}, dt={dt}yr, {nsteps} steps"
    )

    # Checkpoint cadence in YEARS (default 5) so a reboot loses bounded wall
    # time regardless of dt; keep only the last few (plus _final.h5) to bound
    # disk. The step-count `checkpoint_interval` arg is the fallback.
    ckpt_every_yr = float(os.environ.get("ISMIP7_CHECKPOINT_EVERY_YR", "5.0"))
    ckpt_steps = (max(1, int(round(ckpt_every_yr / dt)))
                  if ckpt_every_yr > 0 else checkpoint_interval)
    keep_ckpts = int(os.environ.get("ISMIP7_KEEP_CHECKPOINTS", "3"))

    Q_dg = FunctionSpace(mesh, "DG", 0)
    geom_dg = ctx.get("geom_dg", False)
    # Under DG0 geometry the transport state IS the geometry - the same
    # Function object, not a copy. That is the whole point: the terminus
    # back-pressure and the boundary flux then integrate one thickness, so
    # they cannot disagree. Under CG1 geometry h_dg is a separate DG0 carrier
    # bridged by the lumped lift below.
    h_dg = h if geom_dg else Function(Q_dg, name="h_dg")
    ctx["h_dg_state"] = h_dg
    h_dg_old = Function(Q_dg)
    phi_dg = fd.TestFunction(Q_dg)
    h_dg_trial = fd.TrialFunction(Q_dg)

    n_facet = fd.FacetNormal(mesh)

    s_float = Function(Q_g).interpolate(
        b + (rho_W / rho_I) * max_value(-b, Constant(0.0))
    )

    rho_gt = _RHO_I_SI / 1e12  # m^3 ice -> Gt

    # The t=0 ice extent, needed by EVERY front mechanism (the legacy
    # ISMIP7_FIXED_FRONT flag and every ISMIP7_CALVING law): cells outside it
    # held no ice at t=0, so no apparent-MB reference belongs there, and under
    # a pinned front whatever flows into them is removed each step and tallied
    # as calving flux. Without a front mechanism a buffered mesh has NO calving
    # sink (~1300 Gt/yr in reality) and the sheet must gain mass. Only
    # meaningful when the initial state is the true BedMachine geometry
    # (RC mode / h_clamp_init=0) - with a clamped initial state every cell has
    # ice and the mask is empty.
    fixed_front = _fixed_front()
    front_hmin = _front_hmin()
    calving = _calving_law()
    # An external calving law (ctx["calving_law"], any object with
    # rate(model, t) -> UFL and describe(); hoffmaao/calving's laws.Law is
    # the reference) drives the shared level set through its "prescribed"
    # law, evaluated on the live dual state each transport advance.
    calving_law_obj = ctx.get("calving_law")
    if calving_law_obj is not None:
        if calving != "none":
            raise ValueError(
                f"ISMIP7_CALVING={calving} and an external calving law were "
                f"both requested; leave ISMIP7_CALVING=none for the law object")
        calving = "prescribed"
    beyond_front = None
    n_beyond = 0
    cell_area = assemble(fd.TestFunction(Q_dg) * dx).dat.data_ro.copy()
    # Runaway tripwire (ISMIP7_TRIPWIRE_U_MAX / _H_MAX / _DH_RATE / _HMIN; off
    # unless set). A lane that is running away used to be reported only when
    # the transport budget or Newton finally failed, three steps and half an
    # hour after the fact; failing at the first step that exceeds a bound,
    # naming the cell, lets the ladder experiments answer in minutes.
    # U_MAX and H_MAX are absolute bounds (no Antarctic cell moves faster
    # than 2e4 m/yr or is thicker than 5 km). DH_RATE bounds the relative
    # thickening rate (dh/h)/dt [1/yr] of cells that entered the step at
    # least HMIN thick. It is a rate, not a per-step fraction, so the same
    # physics scores the same on every rung of a dt ladder: the 2026-09-15
    # ladder was first stopped by a 0.3 m buffer cell filling at 70 m/yr
    # (dh/h = 1.6 under a max(h, 10 m) floor) and then, with a 0.5 per-step
    # bound, by a 168 m floating Amundsen cell fed at ~700 m/yr that scored
    # 0.51 at dt = 0.125 and 0.27 at dt = 0.0625 -- the run at 0.0625 went
    # on to complete 1.25 yr with flat speed and thickness. The passing
    # run's largest relative rate was 6.8/yr (a 107 m buffer cell); the
    # dt = 0.25 pile-up thickened 500 -> 3000 m within 0.25 yr (>= 20/yr at
    # onset, ~100/yr later) with speed_max already at 3e4 m/yr.
    _trip = os.environ.get("ISMIP7_TRIPWIRE_U_MAX", "").strip()
    tripwire_u_max = float(_trip) if _trip else None
    _trip = os.environ.get("ISMIP7_TRIPWIRE_H_MAX", "").strip()
    tripwire_h_max = float(_trip) if _trip else None
    _trip = os.environ.get("ISMIP7_TRIPWIRE_DH_RATE", "").strip()
    tripwire_dh_rate = float(_trip) if _trip else None
    tripwire_hmin = float(os.environ.get("ISMIP7_TRIPWIRE_HMIN", "100.0"))
    tripwire_on = any(
        bound is not None
        for bound in (tripwire_u_max, tripwire_h_max, tripwire_dh_rate)
    )
    tripwire_cells = None
    if tripwire_on:
        # Per-cell context for the step report: the t=0 extent (buffer flag,
        # same definition as the fixed front) and the bed (flotation flag).
        tripwire_cells = {
            "extent": Function(Q_dg).project(
                ctx.get("H_init", h)
            ).dat.data_ro.copy(),
            "bed": Function(Q_dg).project(b).dat.data_ro.copy(),
        }
    if fixed_front or calving != "none":
        # Mask from the t=0 observed extent (ctx["H_init"]), not the
        # current h: a restarted run must not re-mask cells that
        # legitimately retreated mid-run inside the observed extent.
        # Use a scratch Function, NOT h_dg: under DG0 geometry h_dg IS the
        # live thickness and borrowing it here would overwrite the geometry
        # with H_init before the run even starts.
        _extent = Function(Q_dg).project(ctx.get("H_init", h))
        beyond_front = _extent.dat.data_ro < front_hmin
        n_beyond = mesh.comm.allreduce(int(beyond_front.sum()))
        PETSc.Sys.Print(
            f"  Initial ice extent: {n_beyond} initially ice-free cells "
            f"masked (h < {front_hmin} m)"
        )
    ctx["fixed_front"] = fixed_front
    ctx["front_hmin"] = front_hmin
    ctx["fixed_front_cells"] = n_beyond

    # Which mechanism owns the REMOVAL of ice past the t=0 extent. A configured
    # level-set law owns the front outright: the ISMIP7 control holds calving
    # at end-of-2014 conditions and gets that from `fixed`, while a projection
    # law (vonmises) must never be silently pinned, and the legacy mask would
    # pin it - the front could retreat but never advance past the 2015 outline,
    # with the inflow mis-tallied as calving. So the legacy sink applies only
    # when no law is configured at all; run_core_matrix.sh exports
    # ISMIP7_FIXED_FRONT=1 unconditionally, which is why the flag may not
    # override an explicit ISMIP7_CALVING choice.
    legacy_front_sink = fixed_front and calving == "none"
    # A free law moves the front, so the frozen a_ref must follow the live
    # extent; `fixed` and the legacy flag pin it on purpose and keep the
    # t=0-only mask.
    free_front = calving not in ("none", "fixed")
    if calving_law_obj is not None:
        front_owner = (
            f"level-set prescribed law (external: {calving_law_obj.describe()})"
            + ("; ISMIP7_FIXED_FRONT is set but ignored for removal"
               if fixed_front else ""))
    elif calving != "none":
        front_owner = f"level-set {calving} law (ISMIP7_CALVING={calving})" + (
            "; ISMIP7_FIXED_FRONT is set but ignored for removal"
            if fixed_front else ""
        )
    elif fixed_front:
        front_owner = "legacy fixed-front mask (ISMIP7_FIXED_FRONT)"
    else:
        front_owner = "none (no calving sink)"
    PETSc.Sys.Print(f"  {FRONT_OWNER_MARKER} {front_owner}")

    # ISMIP7_LEGACY_TRANSPORT=1 restores the pre-Jul-2026 scheme: the
    # -h*div(u*phi) volume term (non-conservative for DG0: it adds
    # spurious h*div(u) mass at thickness jumps) and the L2 DG0->CG1
    # projection (overshoots negative at fronts; the h floor then
    # injects mass). The default is the exactly-conservative FV form
    # plus a lumped-mass projection (convex combination of adjacent
    # cell values: bounded and integral-preserving).
    legacy_transport = os.environ.get("ISMIP7_LEGACY_TRANSPORT") is not None
    if legacy_transport and geom_dg:
        raise ValueError(
            "ISMIP7_LEGACY_TRANSPORT is incompatible with DG0 geometry: the "
            "legacy scheme re-projects CG1 h <-> DG0 h_dg every step, but "
            "under DG0 geometry they are the same Function and the projection "
            "would be self-referential. Use ISMIP7_GEOMETRY_SPACE=cg1 to run "
            "the legacy transport."
        )
    if legacy_transport:
        PETSc.Sys.Print("  LEGACY transport: non-conservative volume term + L2 projection")
    m_lump = assemble(fd.TestFunction(Q) * dx)
    proj_rhs = fd.Cofunction(Q.dual())
    src_dg = Function(Q_dg, name="mass_source")
    src_cof = fd.Cofunction(Q_dg.dual())
    # Where the surface and ocean forcing act: 1 on cells that can hold ice,
    # 0 on open ocean and on cells a front rule holds ice-free
    # (front.unforced_cells). Refreshed at the start of every advance, and at
    # t=0 for the balancing reference, so the reference, the transport source,
    # the budget and the ISMIP7 fields all count the same forcing.
    forced = Function(Q_dg, name="forced")
    bed_cell = Function(Q_dg).project(b).dat.data_ro.copy()

    # h_dg is the PERSISTENT prognostic state (DG0). The old scheme
    # re-projected CG1 h -> DG0 every step; that roundtrip (L2 project +
    # lumped lift) is a per-step smoother, and its meter-scale perturbation
    # at the steep PIG grounding-zone thickness gradient re-triggered the
    # hypersensitive velocity response even from an exactly balanced (a_ref)
    # state. Transport now evolves h_dg directly; the CG1 h is derived
    # (lumped lift, for the diagnostic geometry / VAF), one-way. On a cold
    # start the initial CG h is NOT the lift of its own DG projection, so we
    # lift once here and re-solve the diagnostic on the lifted geometry -
    # otherwise step 1 applies that perturbation mid-run. Restarts skip all
    # of this: checkpoints carry h_dg (thickness_dg), and the stored CG h is
    # its lift by construction.
    #
    # NONE of that applies under DG0 geometry: h_dg IS h, there is no second
    # representation to reconcile, and the one-time lift is skipped entirely.
    # Skipping it is not an optimization - applying it would be the bug. At
    # 32 km it moved the calving front from 105 m to 200 m thick and, through
    # the h^2 terminus traction, multiplied the outflux by ~4.7.
    if geom_dg:
        PETSc.Sys.Print(
            "  DG0 geometry: transport state is the geometry (no lift)"
        )
    elif not legacy_transport:
        if ctx.get("h_dg_state") is not None:
            h_dg.dat.data[:] = ctx["h_dg_state"].dat.data_ro
            PETSc.Sys.Print("  DG thickness state restored from checkpoint")
        else:
            h_dg.project(h)
            _h_lift = proj_rhs  # reuse the cofunction as scratch
            assemble(fd.TestFunction(Q) * h_dg * dx, tensor=_h_lift)
            _new = _h_lift.dat.data_ro / m_lump.dat.data_ro
            _dh = float(np.abs(_new - h.dat.data_ro).max()) if _new.size else 0.0
            from mpi4py import MPI as _MPI4
            _dh = mesh.comm.allreduce(_dh, op=_MPI4.MAX)
            h.dat.data[:] = np.maximum(_new, 0.0)
            s.interpolate(max_value(b + h, (Constant(1.0) - rho_ratio) * h))
            phi_eff.interpolate(
                max_value(
                    Constant(1.0)
                    - rho_W * g * max_value(Constant(0.0), -b)
                    / (rho_I * g * max_value(h, Constant(1.0))),
                    Constant(0.01),
                )
            )
            PETSc.Sys.Print(
                f"  DG-consistent thickness lift (one-time, max |dh|={_dh:.2f} m); "
                f"re-solving diagnostic..."
            )
            try:
                solve_diagnostic("geometry-lift")
            except fd.ConvergenceError:
                PETSc.Sys.Print("    warm-start solve failed; re-ramping n...")
                for step, _t in enumerate(np.linspace(0.0, 1.0, 10), 1):
                    n_flow.assign(1.0 + _t * (n_flow_val - 1.0))
                    m_slide.assign(1.0 + _t * (m_slide_val - 1.0))
                    solve_diagnostic(
                        "geometry-lift-continuation",
                        step=f"{step}/10",
                        t=f"{_t:.6g}",
                    )

    # Apparent-mass-balance reference (ISMIP7_APPARENT_MB=1): a frozen DG0
    # correction equal to the DISCRETE FV flux divergence of the initial
    # (h0, u0), added to the mass source. With it, the t=0 thickness tendency
    # is EXACTLY the forcing (SMB - melt): the init-state flux-divergence
    # spikes (inversion u not flux-consistent with BedMachine h; ~1000 m/yr
    # locally at the PIG grounding zone at 32 km) otherwise dig a surface
    # depression in one step whose driving-stress response runs away - the
    # Jul 2026 forward blow-up, reproduced with NO forcing at all. Same cure
    # as gia forward_monolithic's smoothed a_ref, but built with THIS
    # scheme's own upwind operator so the cancellation is exact. Frozen in
    # time (a fixed MB correction, initMIP-style), saved in checkpoints, and
    # tallied as its own budget column. ISMIP7_AMB_CAP=<m/yr> optionally
    # clips it to [-cap, 5*cap] (gia's asymmetric clip); default uncapped.
    # Modes, resolved by runconfig.apparent_mb_mode: "div" cancels only the
    # flux divergence (t=0 tendency = SMB-melt, gia-style); "balance" (from
    # "1" or "balance") also subtracts the initial forcing, so the t=0
    # tendency is EXACTLY ZERO - a balanced control in the ISMIP6 ctrl_proj
    # sense, against which projections difference cleanly. Both freeze the
    # correction at t=0. None is off.
    amb_mode = apparent_mb_mode()
    apparent_mb = amb_mode is not None
    a_ref = None
    if apparent_mb:
        a_ref = Function(Q_dg, name="a_ref_mb")
        if ctx.get("a_ref_mb") is not None:
            a_ref.dat.data[:] = ctx["a_ref_mb"].dat.data_ro
            PETSc.Sys.Print("  Apparent MB: frozen a_ref loaded from checkpoint")
        else:
            # h_dg already holds the (lifted) state; u0 is the diagnostic
            # velocity solved ON that state - the pair the loop will see.
            u0 = z.subfunctions[0]
            if legacy_transport:
                h_dg.project(h)
            un0 = fd.dot(u0, n_facet)
            un0p = (un0 + abs(un0)) / 2
            flux0 = assemble(
                (un0p("+") * h_dg("+") - un0p("-") * h_dg("-"))
                * fd.jump(phi_dg) * dS
                + un0p * h_dg * phi_dg * ds
            )
            a_ref.dat.data[:] = flux0.dat.data_ro / cell_area
            if ctx.get("phys_div") is not None:
                # Remesh: cancel this mesh's discrete divergence net of the
                # physical divergence the run carried over (adapt_mesh.py).
                a_ref.dat.data[:] -= ctx["phys_div"].dat.data_ro
                PETSc.Sys.Print("  Apparent MB: a_ref rebuilt on the adapted mesh "
                                "from the transferred physical divergence")
            elif amb_mode != "div":
                # balanced control: evaluate the t=0 forcing and fold it in
                if forcing_callback is not None:
                    forcing_callback(ctx, t_start + dt)
                # the forcing the loop will apply, so the t=0 tendency is
                # zero where it acts and no reference is left where it
                # does not (open ocean would otherwise be handed a source
                # equal to the melt it never receives)
                forced.dat.data[:] = np.where(
                    unforced_cells(h_dg.dat.data_ro, bed_cell, beyond_front),
                    0.0, 1.0)
                b_smb = assemble(forced * (accum - ocean_melt) * phi_dg * dx)
                a_ref.dat.data[:] -= b_smb.dat.data_ro / cell_area
            if beyond_front is not None:
                # No ice existed outside the t=0 extent, so no balancing
                # reference belongs there, for ANY front law. Under a pinned
                # front the front tally stays the sink for flux into those
                # cells and must not be absorbed into a_ref; under a free law
                # (vonmises) a frozen sink here would re-empty every cell the
                # front advances into, pinning it at t=0 with no error.
                a_ref.dat.data[beyond_front] = 0.0
            # The same rule against the LIVE extent, both directions:
            #   advance - a frozen SINK outside the extent re-empties the
            #     cells a free front advances into (the t=0 mask above);
            #   retreat - a frozen SOURCE inside the t=0 extent regrows the
            #     cells a free front has calved, since a_ref at a t=0 front
            #     cell is the terminus outflow and is large and positive.
            # A free law therefore re-masks a_ref each step (below) against
            # the cells the level set reports ice-free. This changes free-law
            # projection numbers under ISMIP7_APPARENT_MB; they were wrong
            # before. The pinned laws keep the t=0-only mask: there the front
            # is meant to stay put.
            amb_cap = float(os.environ.get("ISMIP7_AMB_CAP", "0"))
            if amb_cap > 0.0:
                np.clip(a_ref.dat.data, -amb_cap, 5.0 * amb_cap,
                        out=a_ref.dat.data)
            from mpi4py import MPI as _MPI
            _ad = a_ref.dat.data_ro
            _lo = mesh.comm.allreduce(float(_ad.min() if _ad.size else 0), op=_MPI.MIN)
            _hi = mesh.comm.allreduce(float(_ad.max() if _ad.size else 0), op=_MPI.MAX)
            _net = float(assemble(a_ref * dx)) * rho_gt
            PETSc.Sys.Print(
                f"  Apparent MB: a_ref in [{_lo:+.1f}, {_hi:+.1f}] m/yr, "
                f"net {_net:+.1f} Gt/yr"
                + (f" (cap [-{amb_cap:.0f}, +{5*amb_cap:.0f}])" if amb_cap > 0 else "")
            )
        ctx["a_ref_mb"] = a_ref

    # The transport operator has fixed topology; only its Function/Constant
    # coefficients change.  Reuse one solver so every substep does not rebuild
    # the variational problem and PETSc objects.  The stable prefix also makes
    # transport-only PETSc inspection possible without touching the diagnostic
    # solve.
    u_transport = z.subfunctions[0]
    un_transport = fd.dot(u_transport, n_facet)
    un_transport_plus = (un_transport + abs(un_transport)) / 2
    F_transport = (
        (h_dg_trial - h_dg_old) / dt_c * phi_dg * dx
        + (
            un_transport_plus("+") * h_dg_trial("+")
            - un_transport_plus("-") * h_dg_trial("-")
        )
        * fd.jump(phi_dg)
        * dS
        + un_transport_plus * h_dg_trial * phi_dg * ds
        - src_dg * phi_dg * dx
    )
    if legacy_transport:
        F_transport += -h_dg_trial * fd.div(
            u_transport * phi_dg
        ) * dx
    transport_problem = LinearVariationalProblem(
        fd.lhs(F_transport), fd.rhs(F_transport), h_dg
    )
    transport_solver = LinearVariationalSolver(
        transport_problem,
        solver_parameters=transport_solver_parameters(),
        options_prefix="ismip7_transport_",
    )

    # Transport and field telemetry is retained in ctx so timing jobs can
    # persist it beside the diagnostic-solver summary.  The extrema use owned
    # dofs plus collective reductions; a rank-0 .dat statistic would depend on
    # the mesh partition and can miss the cell where a runaway starts.
    transport_stats = []
    field_stats = []
    ctx["transport_stats"] = transport_stats
    ctx["field_stats"] = field_stats
    transport_solve_count = 0
    mass_tol_gt = mass_residual_tol_gt()
    h_diag_xy = Function(
        VectorFunctionSpace(mesh, Q_dg.ufl_element())
    ).interpolate(fd.SpatialCoordinate(mesh))
    speed_diag = Function(Q, name="speed_diagnostic")
    speed_diag_xy = Function(
        VectorFunctionSpace(mesh, Q.ufl_element())
    ).interpolate(fd.SpatialCoordinate(mesh))

    def _field_diagnostics(label):
        r"""Log global h and |u| extrema, including their locations."""
        speed_diag.interpolate(sqrt(fd.dot(u_transport, u_transport)))
        h_min, h_min_xy = global_extreme_location(
            h_dg, h_diag_xy, mode="min"
        )
        h_max, h_max_xy = global_extreme_location(
            h_dg, h_diag_xy, mode="max"
        )
        u_min, u_min_xy = global_extreme_location(
            speed_diag, speed_diag_xy, mode="min"
        )
        u_max, u_max_xy = global_extreme_location(
            speed_diag, speed_diag_xy, mode="max"
        )
        stat = {
            "label": label,
            "thickness_min": h_min,
            "thickness_min_xy": h_min_xy,
            "thickness_max": h_max,
            "thickness_max_xy": h_max_xy,
            "speed_min": u_min,
            "speed_min_xy": u_min_xy,
            "speed_max": u_max,
            "speed_max_xy": u_max_xy,
        }
        field_stats.append(stat)

        def _xy(location):
            return "(" + ", ".join(f"{value:.0f}" for value in location) + ")"

        PETSc.Sys.Print(
            f"=== FIELD RANGE | {label} | "
            f"h=[{h_min:.6e} at {_xy(h_min_xy)}, "
            f"{h_max:.6e} at {_xy(h_max_xy)}] m "
            f"speed=[{u_min:.6e} at {_xy(u_min_xy)}, "
            f"{u_max:.6e} at {_xy(u_max_xy)}] m/yr ==="
        )
        return stat

    def _record_transport(label, elapsed, mass_residual_gt=None):
        nonlocal transport_solve_count
        transport_solve_count += 1
        ksp = transport_solver.snes.getKSP()
        reason = ksp.getConvergedReason()
        reason_name = getattr(reason, "name", str(int(reason)))
        stat = {
            "solve": transport_solve_count,
            "label": label,
            "ksp_reason": reason_name,
            "ksp_iterations": ksp.getIterationNumber(),
            "residual_norm": ksp.getResidualNorm(),
            "mass_residual_gt": mass_residual_gt,
            "seconds": elapsed,
        }
        transport_stats.append(stat)
        mass_text = (
            "unavailable" if mass_residual_gt is None
            else f"{mass_residual_gt:+.6e}"
        )
        PETSc.Sys.Print(
            "=== TRANSPORT RESULT "
            f"{transport_solve_count:04d} | {label} | reason={reason_name} "
            f"ksp_its={stat['ksp_iterations']} "
            f"rnorm={stat['residual_norm']:.6e} "
            f"mass_resid_gt={mass_text} seconds={elapsed:.3f} ==="
        )
        return stat, int(reason)

    mass_prev = float(assemble(h * dx)) * rho_gt

    def _prune_checkpoints():
        r"""Keep only the `keep_ckpts` most recently WRITTEN periodic state
        checkpoints. Recency, not the largest year: a re-run rewinds (every
        projection branches from the historical endpoint), so ranking by year
        would delete each checkpoint this run writes in favour of higher-year
        ones left by a previous run - exactly the states a resume needs."""
        if mesh.comm.rank != 0 or keep_ckpts <= 0:
            return
        import re
        pat = re.compile(re.escape(f"{experiment_name}_{lc}_t")
                         + r"(-?\d+\.\d+)\.h5$")
        found = []
        for fn in glob.glob(
            os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_t*.h5")
        ):
            if not pat.search(os.path.basename(fn)):
                continue
            try:
                found.append((os.path.getmtime(fn), fn))
            except OSError:
                pass
        for _, fn in sorted(found)[:-keep_ckpts]:
            try:
                os.remove(fn)
            except OSError:
                pass

    # ISMIP7 output (ISMIP7_OUTPUT=1): yearly state snapshots and flux means
    # on the model mesh, one checkpoint per year beside this stem, converted
    # and regridded afterwards by antarctica/scripts/write_ismip7_output.py.
    # Off by default.
    #
    # BEFORE the timeseries rewrite below: constructing this is what
    # refuses a run that would overwrite a banked submission series, and
    # that refusal has to happen while the run's own record is still
    # intact, not after the rewrite has dropped the rows past t_start.
    annual = None
    if _ismip7_output():
        from icepack2_tools.ismip7_output import AnnualOutput
        annual = AnnualOutput(
            mesh, Q_dg,
            os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_ismip7_annual.h5"),
            os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_ismip7_scalars.csv"),
            first_year=t_start, rho_ratio=float(rho_ratio), log=PETSc.Sys.Print,
            resume=ctx.get("ismip7_resume"))
        if annual.h_year_start is None:
            annual.start_year(h_dg)
        _stem, _ext = os.path.splitext(annual.out_path)
        PETSc.Sys.Print(
            f"  ISMIP7 output: one checkpoint per year -> {_stem}_<year>{_ext}")

    results = []
    ctx["results"] = results
    ctx["failure"] = None

    # Crash-safe timeseries: append each row and flush, so a reboot keeps the
    # budget-audit history (it used to be dumped only at completion). On a
    # warm restart, drop any rows at/after the resume year, then append.
    csv_fn = os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_timeseries.csv")
    csv_header = ("year,vaf_mm_sle,mass_gt,smb_gtyr,melt_gtyr,"
                  "outflux_gtyr,calv_gt,clamp_gt,resid_gt,amb_gtyr,"
                  + ",".join(COLLAPSE_CSV_COLUMNS) + "\n")
    # The header the file really carries: a series begun before the collapse
    # columns existed keeps its own on resume (collapse_csv_fields).
    csv_head = csv_header
    csv_f = None
    if mesh.comm.rank == 0:
        if t_restart is not None and os.path.exists(csv_fn):
            with open(csv_fn) as _cf:
                _lines = _cf.readlines()
            kept = ([_lines[0]] if _lines and _lines[0].startswith("year")
                    else [csv_header])
            csv_head = kept[0]
            for _ln in _lines[1:]:
                try:
                    if float(_ln.split(",", 1)[0]) <= t_start + 0.5 * dt:
                        kept.append(_ln)
                except (ValueError, IndexError):
                    pass
            with open(csv_fn, "w") as _cf:
                _cf.writelines(kept)
            csv_f = open(csv_fn, "a")
        else:
            csv_f = open(csv_fn, "w")
            csv_f.write(csv_header)
            csv_f.flush()
    # Rank 0 alone read the file, and the note below is printed by all ranks.
    csv_head = mesh.comm.bcast(csv_head, root=0)

    def _grounded_cells():
        return Function(Q_dg).interpolate(s - s_float).dat.data_ro > 0.0

    def _any_rank(local):
        from mpi4py import MPI as _MPI
        return mesh.comm.allreduce(bool(local), op=_MPI.LOR)

    # ISMIP7 ice-shelf collapse forcing (ISMIP7_FRACTURE=mask or mask_front):
    # the forcing callback fills ctx["collapse"] with the year's mask on the
    # geometry cells, and every transport advance removes FLOATING cells it
    # flags, booked as calving. Grounded ice is never touched (protocol
    # path C). `mask` takes every flagged floating cell; `mask_front` only
    # those open water has reached, so no hole opens behind a standing front
    # (discussion #30, icepack2_tools.front.front_connected). Off by default.
    collapse = None
    collapse_neighbours = None
    collapse_edge = None
    fracture_mode = _fracture_mode()
    if fracture_mode in FRACTURE_MASK_MODES:
        if not geom_dg:
            raise RuntimeError(f"ISMIP7_FRACTURE={fracture_mode} needs ISMIP7_GEOMETRY_SPACE=dg0 (cell-wise removal)")
        collapse = np.zeros(len(cell_area), dtype=bool)
        ctx["collapse"] = collapse
        if fracture_mode == "mask_front":
            collapse_neighbours = facet_neighbours(Q_dg)
            # Where the mesh ends at the calving front (the boundary-id
            # sidecar says which exterior facets those are) a cell faces open
            # water even though no ice-free buffer cell is there to say so.
            collapse_edge = assemble(
                phi_dg * ds(tuple(ctx["calving_ids"]))).dat.data_ro > 0.0
    # Under every mode, `none` included, so the mode of a run can be read off
    # its log (issue #10).
    PETSc.Sys.Print(f"  {collapse_banner(fracture_mode)}")
    if not collapse_csv_fields(csv_head, (0, 0, 0)):
        PETSc.Sys.Print(
            "  Timeseries header predates the collapse columns: this resume "
            "keeps it, and the cell counts go to the log alone")

    def _global_cells(mask):
        return global_count(mask, mesh.comm)

    def _write_csv_row(row, collapse_cells):
        if csv_f is None:
            return
        csv_f.write(
            f"{row[0]:.1f},{row[1]:.6f},{row[2]:.2f},"
            + ",".join(f"{v:.4f}" for v in row[3:])
            + collapse_csv_fields(csv_head, collapse_cells) + "\n"
        )
        csv_f.flush()

    # Rescue ladder state: the last CONVERGED mixed state, restored between
    # attempts so each rescue starts from a valid warm point instead of a
    # diverged Newton iterate.
    z_entry = z.copy(deepcopy=True)
    snes_type0 = slvr.snes.getType()
    allow_rescue = rescue_enabled()

    def _ramp(label):
        for step, t in enumerate(np.linspace(0.0, 1.0, 10), 1):
            n_flow.assign(1.0 + t * (n_flow_val - 1.0))
            m_slide.assign(1.0 + t * (m_slide_val - 1.0))
            solve_diagnostic(
                label,
                step=f"{step}/10",
                t=f"{t:.6g}",
            )

    def _solve_with_rescue(k):
        r"""Diagnostic solve with an escalation ladder for hard eras
        (evidence at the 32 km ssp585 2024.5 wall: nleqerr fails, a fresh
        continuation fails, but trust region converges the same step in one
        try). Sequence: direct -> re-continuation -> newtontr direct ->
        newtontr continuation. Restores the entry state between attempts;
        always restores the configured SNES type and full n/m on exit.
        Returns True on success."""
        try:
            solve_diagnostic(f"step-{k}-direct")
            return True
        except fd.ConvergenceError as exc:
            if not allow_rescue:
                _field_diagnostics(f"step-{k}-diagnostic-failed")
                PETSc.Sys.Print(
                    f"  Step {k}: direct diagnostic solve failed; "
                    "rescue disabled"
                )
                ctx["failure"] = {
                    "category": "diagnostic_convergence",
                    "phase": f"step-{k}-direct",
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                }
                raise
        k_lim_c = ctx.get("k_lim")
        k_rescue = ctx.get("k_lim_rescue", 0.0)
        # gia: hard-era steps under trust region converge LINEARLY and are
        # executed by the iteration cap while still descending - rescue
        # rungs get extra patience, restored afterwards.
        rescue_maxit = rescue_max_it()
        _rt, _at, _dt_, _mi = slvr.snes.getTolerances()
        attempts = [
            (
                "re-doing continuation",
                snes_type0,
                lambda: _ramp(f"step-{k}-rescue-continuation"),
                False,
            ),
            (
                "trust-region retry",
                "newtontr",
                lambda: solve_diagnostic(f"step-{k}-trust-region"),
                False,
            ),
            (
                "trust-region continuation",
                "newtontr",
                lambda: _ramp(f"step-{k}-trust-region-continuation"),
                False,
            ),
        ]
        if k_lim_c is not None and k_rescue > 0.0:
            # Deepest rungs: pin the runaway front nodes with the soft speed
            # limiter (only |u| > u_lim feels it) while trust region solves.
            attempts += [
                (
                    "trust-region + speed limiter",
                    "newtontr",
                    lambda: solve_diagnostic(
                        f"step-{k}-trust-region-speed-limiter"
                    ),
                    True,
                ),
                (
                    "trust-region + limiter continuation",
                    "newtontr",
                    lambda: _ramp(
                        f"step-{k}-trust-region-limiter-continuation"
                    ),
                    True,
                ),
            ]
        try:
            for label, stype, action, use_lim in attempts:
                PETSc.Sys.Print(f"  Step {k}: {label}...")
                z.assign(z_entry)
                n_flow.assign(n_flow_val)
                m_slide.assign(m_slide_val)
                slvr.snes.setType(stype)
                slvr.snes.setTolerances(max_it=rescue_maxit)
                if k_lim_c is not None:
                    k_lim_c.assign(k_rescue if use_lim else 0.0)
                try:
                    action()
                    if stype != snes_type0 or use_lim:
                        PETSc.Sys.Print(f"  Step {k}: recovered via {label}")
                    return True
                except fd.ConvergenceError:
                    continue
            return False
        finally:
            slvr.snes.setType(snes_type0)
            slvr.snes.setTolerances(max_it=_mi)
            n_flow.assign(n_flow_val)
            m_slide.assign(m_slide_val)
            if k_lim_c is not None:
                k_lim_c.assign(0.0)

    h_dg_entry = Function(Q_dg)

    def _lift_h():
        r"""Refresh the geometry derived from the thickness state (s, phi_eff).

        Under DG0 geometry h already IS h_dg, so there is nothing to lift and
        only the derived fields are recomputed. Under CG1 geometry this is the
        lumped-mass lift that bridges the two representations.
        """
        if geom_dg:
            # h is h_dg: same Function, already current AND already floored by
            # the transport advance. Re-applying the clamp here would refill
            # the cells the fixed-front mask just emptied, turning the calving
            # sink into an h_clamp source whenever ISMIP7_H_CLAMP > 0.
            pass
        else:
            if legacy_transport:
                h.project(h_dg)
            else:
                # Lumped-mass projection: h_i = ∫phi_i h_dg / ∫phi_i.
                assemble(fd.TestFunction(Q) * h_dg * dx, tensor=proj_rhs)
                h.dat.data[:] = proj_rhs.dat.data_ro / m_lump.dat.data_ro
            h.interpolate(max_value(h, Constant(h_clamp)))
        s.interpolate(max_value(b + h, (Constant(1.0) - rho_ratio) * h))
        phi_eff.interpolate(
            max_value(
                Constant(1.0)
                - rho_W * g * max_value(Constant(0.0), -b)
                / (rho_I * g * max_value(h, Constant(1.0))),
                Constant(0.01),
            )
        )

    # Level-set calving front (ISMIP7_CALVING != none): each step the level
    # set is the eikonal distance to the current ice extent, the calving
    # rate retreats it by normal flow, and the cells it leaves behind plus
    # the sub-cell mass the front cells shed go into the calving tally below
    # (icepack2_tools.levelset). The drag gate it writes is the Function the
    # momentum residual holds. Every law but `fixed` rebuilds the front from
    # the current thickness each step, so a restart needs no saved front: the
    # retreat is already carried in h by the sub-cell shed.
    level_set = None
    phi_entry = None
    last_c_mean = 0.0
    drag_rule = ctx.get("drag_rule")
    if calving != "none":
        from icepack2_tools.levelset import LevelSet, initial_distance
        sig_g, sig_f = _calving_sigma_max()
        # `fixed` holds the front at phi0, which LevelSet captures at
        # construction. On a warm restart h_dg is the RESTARTED extent, so
        # anchor phi0 on the t=0 thickness (ctx["H_init"], reloaded from every
        # checkpoint) instead: a resumed run must not re-freeze the front where
        # it had already retreated to, permanently barring cells a continuous
        # run of the same length would keep. The distance field comes from a
        # distance-only construction on that thickness, because the eikonal
        # solve reads the thickness of the object it belongs to. Use a scratch
        # Function, NOT h_dg: under DG0 geometry h_dg IS the live geometry.
        phi_init = None
        if calving == "fixed":
            _h0 = Function(Q_dg).project(ctx.get("H_init", h))
            phi_init = initial_distance(mesh, _h0, h_min=front_hmin)
        level_set = LevelSet(
            mesh, h_dg, law=calving, h_min=front_hmin,
            sigma_max_grounded=sig_g, sigma_max_floating=sig_f,
            drag_mask=ctx.get("drag_mask"), phi_init=phi_init,
        )
        phi_entry = Function(level_set.Q0)
    live_calving_state = None
    if calving_law_obj is not None:
        live_calving_state = LiveCalvingState(z, h_dg, b, level_set)
        PETSc.Sys.Print(
            f"  Calving law (external, on the live dual state): "
            f"{calving_law_obj.describe()}")
    # save_model_state writes the front and the ISMIP7 year in progress.
    ctx["level_set"] = level_set
    ctx["annual"] = annual
    a_ref_entry = Function(a_ref.function_space()) if a_ref is not None else None
    A_map = ctx.get("A_map")

    def _advance(dt_local, label):
        r"""One transport advance of dt_local with the CURRENT velocity
        (transport-first ordering: the velocity was solved at the current
        geometry). Mutates h_dg and the derived CG fields; returns the
        advance's mass tallies [Gt]."""
        nonlocal last_c_mean
        u_vel = z.subfunctions[0]
        if legacy_transport:
            h_dg.project(h)
        h_dg_old.assign(h_dg)

        # Front first, with the same velocity the transport is about to
        # use, so the cells emptied below are the ones the front left.
        # Whenever a level-set law is configured it is the sole authority on
        # removal, so the legacy t=0 mask contributes nothing here.
        beyond = beyond_front if legacy_front_sink else None
        calv_frac = None
        ls_ice_free = None
        if level_set is not None:
            ext_rate = (calving_law_obj.rate(live_calving_state, t_yr)
                        if calving_law_obj is not None else None)
            last_c_mean = level_set.advance(
                dt_local, u_vel, h_dg, b, A_map, n_flow_val, rate=ext_rate)
            lsb, calv_frac = level_set.calving_masks()
            beyond = lsb
            ls_ice_free = level_set.beyond_front()
            if a_ref is not None and free_front:
                clear_reference_where_ice_free(a_ref.dat.data, ls_ice_free)

        # No forcing where there can be no ice: open ocean at the start of
        # this advance and the cells the front rules hold ice-free. SMB there
        # would make ice the front removes and books as calving; melt there
        # has nothing to melt (front.unforced_cells has the measured cost).
        forced.dat.data[:] = np.where(
            unforced_cells(h_dg_old.dat.data_ro, bed_cell, beyond, ls_ice_free),
            0.0, 1.0)
        smb_f, melt_f, ref_f = applied_forcing(forced, accum, ocean_melt, a_ref)
        smb_gt = float(assemble(smb_f * dx)) * rho_gt * dt_local
        melt_gt = float(assemble(melt_f * dx)) * rho_gt * dt_local
        src = smb_f - melt_f
        amb_gt = 0.0
        if ref_f is not None:
            src = src + ref_f
            # The reference as APPLIED here: the live-extent mask above may
            # have zeroed cells since the step's entry measurement, and the
            # forcing mask withholds it from open ocean, so the budget and
            # the CSV must use this, not the entry value.
            amb_gt = float(assemble(ref_f * dx)) * rho_gt * dt_local
        # Cell-averaged DG0 source (exact for the DG0 test space) with a
        # positivity limit (gia a_step clamp): the net sink may not draw a
        # cell below h_clamp within one advance. With the limited source
        # the implicit upwind update is an M-matrix system with nonnegative
        # RHS, so h stays >= h_clamp and the post-solve floor is a no-op.
        # The withheld sink is tallied into the clamp budget column.
        # The bound is capped at zero so it can only ever CAP A SINK: for a
        # cell already below h_clamp (an emptied fixed-front cell) the raw
        # bound is positive and would make the limiter a mandatory SOURCE,
        # injecting h_clamp of ice per step for the mask to re-calve.
        assemble(src * phi_dg * dx, tensor=src_cof)
        src_dg.dat.data[:] = src_cof.dat.data_ro / cell_area
        _src_want = src_dg.dat.data_ro.copy()
        np.maximum(
            src_dg.dat.data,
            np.minimum(-(h_dg.dat.data_ro - h_clamp) / dt_local, 0.0),
            out=src_dg.dat.data,
        )
        limit_gt = mesh.comm.allreduce(float(
            ((src_dg.dat.data_ro - _src_want) * cell_area).sum()
        )) * rho_gt * dt_local
        m0 = float(assemble(h_dg_old * dx)) * rho_gt
        source_gt = float(assemble(src_dg * dx)) * rho_gt * dt_local
        dt_c.assign(dt_local)
        transport_t0 = perf_counter()
        try:
            transport_solver.solve()
        except Exception as exc:
            elapsed = perf_counter() - transport_t0
            _record_transport(label, elapsed)
            _field_diagnostics(f"{label}-transport-failed")
            ctx["failure"] = {
                "category": "transport_convergence",
                "phase": label,
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
            raise

        out_gt = float(assemble(
            un_transport_plus * h_dg * ds
        )) * rho_gt * dt_local
        m1 = float(assemble(h_dg * dx)) * rho_gt
        transport_resid_gt = None
        if not legacy_transport:
            transport_resid_gt = (m1 - m0) - (source_gt - out_gt)
        elapsed = perf_counter() - transport_t0
        _transport_stat, transport_reason = _record_transport(
            label, elapsed, transport_resid_gt
        )
        if transport_reason <= 0:
            _field_diagnostics(f"{label}-transport-diverged")
            ctx["failure"] = {
                "category": "transport_convergence",
                "phase": label,
                "exception_type": "RuntimeError",
                "message": (
                    f"Transport KSP diverged in {label}: "
                    f"reason={_transport_stat['ksp_reason']}"
                ),
            }
            raise RuntimeError(
                ctx["failure"]["message"]
            )
        if transport_resid_gt is not None and (
            not np.isfinite(transport_resid_gt)
            or abs(transport_resid_gt) > mass_tol_gt
        ):
            _field_diagnostics(f"{label}-transport-budget-failed")
            message = (
                f"Transport mass residual {transport_resid_gt:+.6e} Gt "
                f"exceeds {mass_tol_gt:.6e} Gt in {label}"
            )
            PETSc.Sys.Print(f"ERROR: {message}")
            ctx["failure"] = {
                "category": "transport_mass_budget",
                "phase": label,
                "exception_type": "RuntimeError",
                "message": message,
            }
            raise RuntimeError(message)
        # One mesh-wide interpolation for the whole advance: `s` is only
        # refreshed by _lift_h() at the end, so the booking, the floor
        # exemption and the collapse removal all mean the same grounding
        # state and must see it.
        grounded = _grounded_cells()
        if annual is not None:
            # the sources as the limiter left them, and only where there can
            # be ice: the withheld sink comes off the booked SMB, melt and
            # reference (ismip7_output)
            annual.book_advance(dt_local, smb_f, melt_f, ref_f,
                                h_dg, u_vel, grounded,
                                withheld=src_dg.dat.data_ro - _src_want)

        # Floor to h_clamp, EXCEPT in the cells the front rules report as
        # holding no ice: see clamp_thickness for why every such rule has to
        # name its cells here.
        collapsed = (collapse & ~grounded) if collapse is not None else None
        flagged = collapsed
        if collapsed is not None and collapse_neighbours is not None:
            # Open water as the advance found it: floating cells below the
            # ice-mask thickness (the ones emptied on an earlier advance
            # among them; after the transport they carry a trace of inflow,
            # so this advance's h cannot say) and the cells the front rules
            # name. Collective, so every rank sweeps or none does.
            open_water = ~grounded & (h_dg_old.dat.data_ro < front_hmin)
            for ice_free in (beyond, ls_ice_free):
                if ice_free is not None:
                    open_water |= ice_free & ~grounded
            collapsed = front_connected(
                flagged, open_water | (collapse_edge & ~grounded),
                collapse_neighbours, _any_rank)
        # Global counts, on every rank under either mask mode and on none
        # under `none`. `mask` holds nothing by construction, so its held
        # count is zero and says so.
        collapse_cells = collapse_cell_counts(flagged, collapsed, _global_cells)
        clamp_thickness(h_dg.dat.data, h_clamp, ls_ice_free, beyond, collapsed)
        m2 = float(assemble(h_dg * dx)) * rho_gt
        clamp_gt = m2 - m1                                       # Gt added by DG floor

        calv_gt = 0.0
        if collapse is not None and _any_rank(collapse.any()):
            data = h_dg.dat.data
            hit = collapsed & (data > 0.0)
            calv_gt += mesh.comm.allreduce(
                float((data[hit] * cell_area[hit]).sum())) * rho_gt
            if annual is not None:
                annual.book_removal(hit, data[hit])
            data[hit] = 0.0
        if beyond is not None:
            data = h_dg.dat.data
            calv_gt += mesh.comm.allreduce(
                float((data[beyond] * cell_area[beyond]).sum())
            ) * rho_gt
            if annual is not None:
                annual.book_removal(beyond, data[beyond])
            data[beyond] = 0.0
        if calv_frac is not None:
            # Sub-cell calving: the front cells shed the fraction of their
            # thickness that a front retreating at rate c removes in dt.
            data = h_dg.dat.data
            shed = data * calv_frac
            calv_gt += mesh.comm.allreduce(
                float((shed * cell_area).sum())) * rho_gt
            if annual is not None:
                annual.book_removal(slice(None), shed)
            data -= shed
        if level_set is not None:
            # Retreat slivers only: what the sub-cell shed and the melt leave
            # behind in a cell that HELD ice when the advance began is the
            # front's own loss, so it is removed and tallied as calving. A
            # cell that was ice-free keeps whatever the transport just put
            # there, however little: that sub-threshold inflow is how the
            # front ADVANCES, and zeroing it would pin the front wherever the
            # one-step influx is under front_hmin. So outside the level set's
            # extent the thickness is NOT exactly zero - it may hold inflow
            # accumulating toward the threshold. The level set's drag gate has
            # switched the ocean drag OFF in those cells, and what damps them
            # splits by whether they lie inside the t=0 extent:
            #   inside  - h_visc_floor, the basal friction law, and the
            #     alpha_gl collar. N > 0 on a thin cell (effective_pressure
            #     floors the overburden at 1 m while p_W uses the true
            #     thickness), and C_w0 > 0, so tau_b is nonzero.
            #   outside (the advance strip of a free law) - h_visc_floor and
            #     the alpha_gl collar ONLY. N > 0 for the same reason, but
            #     C_w0 is weertman_anchor at the t=0 geometry, where H = 0,
            #     and c_w0_floor defaults to 0, so tau_W = 0 and tau_b = 0
            #     under both laws whatever N is.
            data = h_dg.dat.data
            sliver = retreat_slivers(data, h_dg_old.dat.data_ro, front_hmin)
            calv_gt += mesh.comm.allreduce(
                float((data[sliver] * cell_area[sliver]).sum())) * rho_gt
            if annual is not None:
                annual.book_removal(sliver, data[sliver])
            data[sliver] = 0.0

        _lift_h()
        # The ocean-drag gate follows the ice: without a level set (which
        # writes its own) the next diagnostic solve sees drag only in water
        # outside the t=0 extent that no ice cell now touches.
        if level_set is None and drag_rule is not None:
            drag_rule(h_dg.dat.data_ro)
        m3 = m2 - calv_gt                                        # ∫h preserved by projection
        mass_now = float(assemble(h * dx)) * _RHO_I_SI / 1e12
        clamp_cg_gt = mass_now - m3          # Gt added by the CG floor after projection
        return {
            "out_gt": out_gt,
            "calv_gt": calv_gt,
            # clamp_gt keeps the full sum (the step budget identity below
            # depends on it); limit_gt is ALSO reported alone because it is
            # exactly the error the positivity limiter introduces into an
            # apparent-MB cancellation in thin converging cells.
            "clamp_gt": clamp_gt + clamp_cg_gt + limit_gt,
            "limit_gt": limit_gt,
            "amb_gt": amb_gt,
            "smb_gt": smb_gt,
            "melt_gt": melt_gt,
            "collapse_cells": collapse_cells,
        }

    # Time loop, transport-first: each step advances the geometry with the
    # velocity SOLVED AT it (the previous solve), then solves at the new
    # geometry - the same (G, u) sequence as the old solve-then-advance
    # ordering, but a hard step can now rewind ITS OWN advance and retry it
    # as dt/4 then dt/16 subcycles with rescue solves between (the 1873
    # front-cell-emptying eras: the ladder alone fails at dt=0.1 but a
    # dt=0.025 approach crosses them - smaller geometry increments let
    # Newton track the branch through the event). Checkpoints improve too:
    # saved (h, u) are now mutually consistent.
    SUBCYCLES = subcycles()

    # Timing drivers use this marker rather than timing setup_model or the
    # transport-operator construction above.  It is also available after an
    # exception, so a partial failure record retains the useful elapsed time.
    ctx["transient_t0"] = perf_counter()

    # Wall-clock budget, in minutes from process start. A batch job that runs
    # into its Slurm limit is killed mid-step and loses everything since the
    # last periodic checkpoint, and its chained successor never starts. With a
    # budget the run stops itself between steps, falls through to the final
    # save below, and exits cleanly with t_yr short of t_end, which is exactly
    # what the chain resubmits from. The runner derives the value from the
    # partition's own limit; 0 (the default) is no budget.
    wall_stop_min = float(os.environ.get("ISMIP7_WALL_STOP_MIN", "0"))
    if wall_stop_min > 0:
        PETSc.Sys.Print(
            f"  Wall-clock budget: stopping after {wall_stop_min:g} min "
            f"(ISMIP7_WALL_STOP_MIN)"
        )
    # Longest step so far, the estimate of what the NEXT one could cost. The
    # budget is checked before entering a step rather than after finishing
    # one: a step that goes through the rescue ladder and then 16 subcycles
    # runs far longer than a normal one, and a run that starts such a step
    # just under the budget overshoots it by that whole step and gets killed
    # mid-write, which is the outcome the budget exists to avoid.
    step_max_min = 0.0
    stalled = False
    collapse_cells = (0, 0, 0)
    collapse_held_peak = (0, t_start)      # most cells held at once, and when

    for k in range(1, nsteps + 1):
        if wall_stop_min > 0:
            mins = (perf_counter() - _T_PROCESS_START) / 60.0
            # Every rank must leave the loop together: the final save is
            # collective, so a rank-local decision would hang the job.
            if mesh.comm.allreduce(
                    int(mins + step_max_min >= wall_stop_min)) > 0:
                PETSc.Sys.Print(
                    f"  Wall-clock budget: {mins:.1f} min used of "
                    f"{wall_stop_min:g}, and the longest step so far took "
                    f"{step_max_min:.1f} min. Stopping at t_yr={t_start + (k - 1) * dt:.1f} "
                    f"with the budget intact, for the next job in the chain"
                )
                break
        t_step_start = perf_counter()
        t_yr = t_start + k * dt

        if forcing_callback is not None:
            forcing_callback(ctx, t_yr)


        z_entry.assign(z)
        h_dg_entry.assign(h_dg)
        if level_set is not None:
            phi_entry.assign(level_set.phi)
        # a_ref is mutated by the live-extent mask inside _advance, so it is
        # step state and must rewind with the rest: an abandoned attempt that
        # calved a cell the accepted trajectory keeps must not leave that
        # cell without its balancing reference.
        if a_ref_entry is not None:
            a_ref_entry.assign(a_ref)
        tallies = None
        for m in SUBCYCLES:
            if annual is not None:
                annual.begin_step()          # a rewound attempt must not double-count
            if m > 1:
                PETSc.Sys.Print(
                    f"  Step {k}: subcycling x{m} (dt={dt / m:.4g})..."
                )
                h_dg.assign(h_dg_entry)
                _lift_h()
                z.assign(z_entry)
                if a_ref_entry is not None:
                    a_ref.assign(a_ref_entry)
                if level_set is not None:
                    level_set.phi.assign(phi_entry)
                    level_set.update_cell_fields()
            acc = {"out_gt": 0.0, "calv_gt": 0.0, "clamp_gt": 0.0,
                   "limit_gt": 0.0, "amb_gt": 0.0,
                   "smb_gt": 0.0, "melt_gt": 0.0}
            ok = True
            for _j in range(m):
                sub = _advance(
                    dt / m, f"step-{k}-substep-{_j + 1}/{m}"
                )
                for key in acc:
                    acc[key] += sub[key]
                # A state where the tallies are sums, so the last advance of
                # the accepted attempt stands for the step.
                collapse_cells = sub["collapse_cells"]
                if not _solve_with_rescue(k):
                    ok = False
                    break
            if ok:
                if m > 1:
                    PETSc.Sys.Print(f"  Step {k}: completed via x{m} subcycle")
                tallies = acc
                if annual is not None:
                    annual.commit_step()
                break
        if tallies is None:
            ctx["failure"] = {
                "category": "diagnostic_convergence",
                "phase": f"step-{k}-rescue-exhausted",
                "exception_type": "ConvergenceError",
                "message": "diagnostic rescue ladder and subcycles exhausted",
            }
            PETSc.Sys.Print(
                f"  Step {k}: rescue ladder + subcycles exhausted, "
                f"saving and stopping"
            )
            stalled = True
            h_dg.assign(h_dg_entry)
            _lift_h()
            z.assign(z_entry)   # checkpoint the last converged pair
            if a_ref_entry is not None:
                a_ref.assign(a_ref_entry)
            if level_set is not None:
                level_set.phi.assign(phi_entry)
                level_set.update_cell_fields()
            break

        t_elapsed = perf_counter() - t_step_start

        haf = Function(Q_g).interpolate(max_value(s - s_float, Constant(0.0)))
        vaf = float(assemble(haf * dx)) * _RHO_I_SI / 1e12 / 362.5
        total_mass = float(assemble(h * dx)) * _RHO_I_SI / 1e12

        # SMB and melt as the advances applied them, where there was ice to
        # force, so the budget, the timeseries and the ISMIP7 fields agree.
        smb_rate = tallies["smb_gt"] / dt                        # Gt/yr
        melt_rate = tallies["melt_gt"] / dt                      # Gt/yr
        out_rate = tallies["out_gt"] / dt                        # Gt/yr
        calv_gt = tallies["calv_gt"]
        clamp_all = tallies["clamp_gt"]
        amb_rate = tallies["amb_gt"] / dt                        # Gt/yr
        dm = total_mass - mass_prev
        resid_gt = dm - (
            (smb_rate - melt_rate + amb_rate - out_rate) * dt
            + clamp_all - calv_gt
        )
        mass_prev = total_mass

        results.append((t_yr, vaf, total_mass, smb_rate, melt_rate,
                        out_rate, calv_gt, clamp_all, resid_gt,
                        amb_rate))
        _write_csv_row(results[-1], collapse_cells)
        (ctx["collapse_flagged_cells"], ctx["collapse_removed_cells"],
         ctx["collapse_held_cells"]) = collapse_cells
        if collapse_cells[2] > collapse_held_peak[0]:
            collapse_held_peak = (collapse_cells[2], t_yr)
        ctx.setdefault("step_budget", []).append({
            "step": k,
            "t_yr": float(t_yr),
            "out_gt": float(tallies["out_gt"]),
            "calv_gt": float(calv_gt),
            "clamp_gt": float(clamp_all),
            "limit_gt": float(tallies.get("limit_gt", 0.0)),
            "amb_gt_per_yr": float(amb_rate),
            "dm_gt": float(dm),
            "resid_gt": float(resid_gt),
        })

        step_stat = _field_diagnostics(f"step-{k}-t={t_yr:.6g}")

        if tripwire_on:
            _h_now = np.asarray(h_dg.dat.data_ro)
            _h_was = np.asarray(h_dg_entry.dat.data_ro)
            _dh = _h_now - _h_was
            _eligible = _h_was >= tripwire_hmin
            growth = np.divide(
                _dh, _h_was,
                out=np.full(_dh.shape, -np.inf), where=_eligible,
            )
            _payload = {
                "h_entry": _h_was,
                "h_now": _h_now,
                "dh": _dh,
                "extent": tripwire_cells["extent"],
                "bed": tripwire_cells["bed"],
            }
            _xy = np.asarray(h_diag_xy.dat.data_ro)
            growth_max, growth_xy, growth_cell = _global_argmax_with_payload(
                growth, _xy, _payload, mesh.comm
            )
            growth_rate_max = growth_max / dt          # (dh/h)/dt [1/yr]
            dh_abs_max, dh_xy, dh_cell = _global_argmax_with_payload(
                np.abs(_dh), _xy, _payload, mesh.comm
            )
            speed_max = float(step_stat["speed_max"])
            speed_xy = step_stat["speed_max_xy"]
            h_max = float(step_stat["thickness_max"])
            h_max_xy = step_stat["thickness_max_xy"]

            def _cell_txt(cell):
                if not cell:
                    return "-"
                flags = []
                if cell["extent"] < front_hmin:
                    flags.append("buffer")
                afloat = cell["h_now"] < (
                    max(0.0, -cell["bed"]) * _RHO_W_SI / _RHO_I_SI
                )
                flags.append("floating" if afloat else "grounded")
                return (
                    f"h {cell['h_entry']:.1f}->{cell['h_now']:.1f} m, "
                    + "/".join(flags)
                )

            def _xy_txt(location):
                if not location:
                    return "(-)"
                return "(" + ", ".join(f"{v:.0f}" for v in location) + ")"

            _growth_txt = (
                f"{growth_rate_max:+.2f}/yr (dh/h={growth_max:+.3f})"
                if np.isfinite(growth_max) else "n/a"
            )
            PETSc.Sys.Print(
                f"  tripwire step-{k}: max (dh/h)/dt={_growth_txt} "
                f"(cells h>={tripwire_hmin:g} m) at {_xy_txt(growth_xy)} "
                f"[{_cell_txt(growth_cell)}]; max |dh|="
                f"{dh_cell.get('dh', dh_abs_max):+.1f} m at {_xy_txt(dh_xy)} "
                f"[{_cell_txt(dh_cell)}]; speed_max={speed_max:.3e} m/yr; "
                f"h_max={h_max:.1f} m"
            )
            ctx.setdefault("tripwire", {
                "u_max": tripwire_u_max,
                "h_max": tripwire_h_max,
                "dh_rate": tripwire_dh_rate,
                "hmin": tripwire_hmin,
                "steps": [],
            })["steps"].append({
                "step": k,
                "speed_max": speed_max,
                "speed_max_xy": list(speed_xy),
                "thickness_max": h_max,
                "thickness_max_xy": list(h_max_xy),
                "growth_max": (
                    float(growth_max) if np.isfinite(growth_max) else None
                ),
                "growth_rate_max": (
                    float(growth_rate_max)
                    if np.isfinite(growth_rate_max) else None
                ),
                "growth_max_xy": list(growth_xy),
                "growth_max_cell": growth_cell,
                "dh_abs_max": float(dh_abs_max),
                "dh_abs_max_xy": list(dh_xy),
                "dh_abs_max_cell": dh_cell,
            })
            tripped = []
            if tripwire_u_max is not None and speed_max > tripwire_u_max:
                tripped.append(
                    f"speed_max={speed_max:.3e} m/yr at "
                    f"{_xy_txt(speed_xy)} > {tripwire_u_max:g}"
                )
            if tripwire_h_max is not None and h_max > tripwire_h_max:
                tripped.append(
                    f"thickness_max={h_max:.1f} m at {_xy_txt(h_max_xy)} "
                    f"> {tripwire_h_max:g}"
                )
            if (
                tripwire_dh_rate is not None
                and np.isfinite(growth_rate_max)
                and growth_rate_max > tripwire_dh_rate
            ):
                tripped.append(
                    f"(dh/h)/dt={growth_rate_max:.2f}/yr "
                    f"(dh/h={growth_max:.3f} in one step of {dt:g} yr) at "
                    f"{_xy_txt(growth_xy)} [{_cell_txt(growth_cell)}] "
                    f"> {tripwire_dh_rate:g}/yr"
                )
            if tripped:
                message = "; ".join(tripped)
                ctx["failure"] = {
                    "category": "runaway_tripwire",
                    "phase": f"step-{k}",
                    "exception_type": "RuntimeError",
                    "message": message,
                }
                PETSc.Sys.Print(f"RUNAWAY TRIPWIRE step-{k}: {message}")
                raise RuntimeError(f"runaway tripwire at step {k}: {message}")
        if annual is not None and abs(t_yr - round(t_yr)) < 1e-6:
            annual.year_end(h_dg, s, b, z.subfunctions[0], z.subfunctions[2],
                            _grounded_cells(), h_dg.dat.data_ro > 1.0)

        if k % output_interval == 0 or k == 1:
            _amb_txt = f"amb={amb_rate:+.0f} " if a_ref is not None else ""
            PETSc.Sys.Print(
                f"  t={t_yr:.1f}  VAF={vaf:.4f} mm SLE  "
                f"mass={total_mass:.1f} Gt  [{t_elapsed:.1f}s]\n"
                f"      budget [Gt/yr]: SMB={smb_rate:+.0f} melt={-melt_rate:+.0f} "
                f"{_amb_txt}"
                f"outflux={-out_rate:+.0f} calv={-calv_gt/dt:+.0f} "
                + (f"c_front={last_c_mean:.0f}m/yr " if level_set is not None else "")
                +
                f"clamp={clamp_all/dt:+.1f} "
                f"dM/dt={dm/dt:+.0f} resid={resid_gt/dt:+.2f}"
                + (f"\n      collapse [cells]: flagged={collapse_cells[0]} "
                   f"removed={collapse_cells[1]} held={collapse_cells[2]}"
                   if collapse is not None else "")
            )

        if not np.isfinite(resid_gt) or abs(resid_gt) > mass_tol_gt:
            message = (
                f"Step {k} mass residual {resid_gt:+.6e} Gt exceeds "
                f"{mass_tol_gt:.6e} Gt"
            )
            PETSc.Sys.Print(f"ERROR: {message}")
            ctx["failure"] = {
                "category": "step_mass_budget",
                "phase": f"step-{k}",
                "exception_type": "RuntimeError",
                "message": message,
            }
            raise RuntimeError(message)

        if k % ckpt_steps == 0:
            chk_fn = os.path.join(
                RESULTS_DIR, f"{experiment_name}_{lc}_t{t_yr:.1f}.h5"
            )
            save_model_state(ctx, chk_fn, t_yr)
            _prune_checkpoints()
            PETSc.Sys.Print(f"    [checkpoint: {os.path.basename(chk_fn)}]")

        step_max_min = max(step_max_min,
                           (perf_counter() - t_step_start) / 60.0)

    ctx["transient_seconds"] = perf_counter() - ctx["transient_t0"]
    PETSc.Sys.Print(f"\n{experiment_name} simulation complete.")
    if collapse is not None:
        # One line for the run record (core_report.py lifts the marker). The
        # peak covers the steps this process ran, so a chained run prints one
        # such line per link.
        _flagged, _removed, _held = collapse_cells
        PETSc.Sys.Print(
            f"{COLLAPSE_MARKER} ISMIP7_FRACTURE={fracture_mode} ended "
            f"t={results[-1][0] if results else t_start:.1f} with "
            f"flagged={_flagged} removed={_removed} held={_held} cells; "
            f"most held at once {collapse_held_peak[0]} cells at "
            f"t={collapse_held_peak[1]:.1f} (steps from t={t_start:.1f})"
        )

    # Final state is a self-contained checkpoint too (a valid restart source).
    # Save the ACTUAL last year so a resume after an early stop continues from
    # where it really stopped, not the nominal t_end.
    final_fn = os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_final.h5")
    last_t = results[-1][0] if results else t_start
    save_model_state(
        ctx, final_fn, last_t,
        # The driver's own verdict on why it stopped here, so a batch
        # chain does not have to infer it from log text or from the year
        # alone: 1 means the solver gave up, and resuming would re-attempt
        # the same years and give up again.
        extra_attrs={"stalled": int(bool(stalled))},
    )
    # Printed on every exit, early stop included. projection.sbatch reads
    # this exact "Saved: <...>_final.h5" line out of its own Slurm log to find
    # THIS job's checkpoint (the results directory is flat and shared, so the
    # newest file does not identify the run); keep the prefix and the path on
    # one line.
    PETSc.Sys.Print(f"Saved: {final_fn}")

    if csv_f is not None:
        csv_f.close()
    PETSc.Sys.Print(f"Saved: {csv_fn}")
    if annual is not None:
        annual.close()
        _yrs = annual.years_on_disk(annual.out_path)
        PETSc.Sys.Print(
            f"Saved: {len(_yrs)} ISMIP7 years"
            + (f" ({_yrs[0]}-{_yrs[-1]})" if _yrs else "")
            + f" beside {annual.out_path}")

    return results
