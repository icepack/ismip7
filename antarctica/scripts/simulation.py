#!/usr/bin/env python3
r"""Shared simulation engine for ISMIP7 Antarctic experiments."""

import numpy as np
import os, sys, glob
from time import perf_counter

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
DATA_DIR = os.environ.get("ISMIP7_OBS_DATA_ROOT", os.path.join(_ROOT, "data"))
MESH_DIR = os.path.join(_ROOT, "mesh")
RESULTS_DIR = os.path.join(_ROOT, "results")

# Repo root on the path for the shared dual-friction operator.
sys.path.insert(0, os.path.dirname(_ROOT))
from mesh_naming import mesh_filename

from icepack2_tools.mpi_stats import (
    global_extreme_location,
    global_mean,
    global_range,
)
from icepack2_tools.boundary import load_boundary_ids
from icepack2_tools.geometry import sample_to_geometry
from icepack2_tools.naming import map_basename
from icepack2_tools.runconfig import (
    friction as _friction, geometry_space as _geometry_space,
    lc as _lc, lc_coarse as _lc_coarse, n_flow as _n_flow,
)
from icepack2_tools.solverconfig import (
    continuation_steps,
    diagnostic_solver_label,
    diagnostic_solver_mode,
    diagnostic_solver_parameters,
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


def setup_model(restart_from=None):
    r"""Load mesh, data, inversion fields, and build diagnostic solver."""
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Friction law: "budd" (power-law dual, default) or "regularized_coulomb"
    # (RC residual driven by the inversion_icepack2_rc_<lc>.h5 MAP:
    # grounded-only theta, exact-zero shelf drag via the Coulomb cap c0*N).
    # Must mirror inversion_icepack2.py so theta/phi keep their meaning.
    # Resolved first because the compute mesh is loaded FROM this friction's
    # MAP checkpoint (below), not from a fresh Mesh(.msh).
    friction = _friction()
    # Exact-zero-shelf residual laws (icepack2 dual, dual_friction.py):
    #   regularized_coulomb -> Coulomb cap tau_c=c0*N
    #   budd                -> tau_b ~ N_hat=N_eff/N_ref, PISM-delta grounded
    #                          floor, exact-zero shelf (matches gia ASE ase_model)
    # Both share the C_w0/He/composite structure and need a _rc/_budd MAP.
    # "budd_legacy" keeps the old phi_eff action Budd (residual shelf drag).
    use_residual = friction in ("regularized_coulomb", "budd")
    use_rc = use_residual  # geometry/alpha/h_clamp handling is shared

    is_restart = restart_from is not None
    # Prefer the MAP inverted under this geometry space. Falling back to the
    # untagged legacy (CG1) MAP keeps the transition runnable, but the load
    # then projects the geometry and warns loudly - the controls still carry
    # the CG1 front bias, so such a run is a smoke test, not a result.
    inv_fn = os.environ.get(
        "ISMIP7_INVERSION",
        os.path.join(MESH_DIR, map_basename(friction, lc)),
    )
    if not is_restart and "ISMIP7_INVERSION" not in os.environ \
            and not os.path.exists(inv_fn):
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

        # Dan/David additionally stamp the mesh PARAMETERS. Keep both: the
        # basename is direct and also covers meshes outside the standard
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

    mesh_fn = os.environ.get("ISMIP7_MESH")
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
            f"({mesh.num_vertices()} vertices, {mesh.num_cells()} cells)"
        )
    else:
        mesh = source_mesh
        target_lc_coarse = chk_lc_coarse
        target_buffer_m = chk_buffer_m
    PETSc.Sys.Print(f"  {mesh.num_vertices()} vertices, {mesh.num_cells()} cells")
    # The recorded basename is the provenance and WINS. ISMIP7_MESH is only a
    # fallback for legacy checkpoints that carry no attribute: it names the
    # mesh the caller intends to build, which is not necessarily the one this
    # checkpoint was written on, so it must not override the record - and only
    # the recorded value is carried into the checkpoints this run writes.
    if not mesh_basename:
        _env_mesh = os.environ.get("ISMIP7_MESH", "")
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

    if not is_restart:
        # Cold start: geometry from BedMachine (RC/Budd overwrites it with the
        # inversion-time geometry from the MAP in the reference-load block).
        bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
        # Cell average onto the geometry space, NOT a centroid point sample --
        # see geometry.sample_to_geometry for the measurements behind that.
        b = sample_to_geometry(
            lambda sp: icepack.interpolate(
                rasterio.open(f"netcdf:{bm_fn}:bed"), sp), Q_g, Q)
        b.rename("bed")
        H = sample_to_geometry(
            lambda sp: icepack.interpolate(
                rasterio.open(f"netcdf:{bm_fn}:thickness"), sp),
            Q_g, Q, floor=h_clamp_init)
        H.rename("thickness")
        s = Function(Q_g, name="surface").interpolate(
            max_value(b + H, (Constant(1.0) - rho_ratio) * H)
        )

    # Reference log-adjustments (theta=log_friction, phi=log_fluidity) plus,
    # for RC/Budd or any restart, the geometry and frozen anchors — all read
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
    h_dg_state = None
    t_restart = None
    A_prior_f = None
    def load_checkpoint_field(chk, name, space, optional=False):
        """Load a checkpoint field, interpolating it for a timing mesh."""
        try:
            source_field = chk.load_function(source_mesh, name=name)
        except (KeyError, RuntimeError, ValueError):
            if optional:
                return None
            raise
        target_field = Function(space, name=name)
        target_field.interpolate(
            source_field,
            allow_missing_dofs=True,
            default_missing_val=0.0,
        )
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
            chk, "fluidity_prior", Q, optional=True
        )
        if is_restart:
            # Self-contained restart: evolved geometry, frozen anchors, time.
            b = load_checkpoint_field(chk, "bed", Q_g)
            H = load_checkpoint_field(chk, "thickness", Q_g)
            s = load_checkpoint_field(chk, "surface", Q_g)
            H_init = load_checkpoint_field(chk, "H_init", Q_g)
            phi_eff = load_checkpoint_field(chk, "phi_eff", Q_g)
            u_guess = load_checkpoint_field(chk, "velocity", V)
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
            amb_env = os.environ.get("ISMIP7_APPARENT_MB")
            if a_ref_mb is not None and amb_env is None:
                raise RuntimeError(
                    f"Restart checkpoint {source_chk} carries a frozen "
                    f"a_ref_mb (the run used ISMIP7_APPARENT_MB) but "
                    f"ISMIP7_APPARENT_MB is unset; set it to resume with "
                    f"the same mass-balance correction."
                )
            if a_ref_mb is None and amb_env is not None:
                raise RuntimeError(
                    f"ISMIP7_APPARENT_MB is set but restart checkpoint "
                    f"{source_chk} has no a_ref_mb; a fresh a_ref cannot be "
                    f"built from an evolved state. Unset ISMIP7_APPARENT_MB "
                    f"or restart from a checkpoint that carries a_ref_mb."
                )
            PETSc.Sys.Print(
                f"  Restart: evolved geometry + frozen anchors loaded "
                f"(t_yr={t_restart}, friction={friction})"
            )
        elif use_rc:
            # Cold RC/Budd: geometry + velocity_obs from the MAP so the
            # Weertman anchor C_w0 (hence the meaning of theta) is reproduced.
            H = load_checkpoint_field(chk, "thickness", Q_g)
            b = load_checkpoint_field(chk, "bed", Q_g)
            s = load_checkpoint_field(chk, "surface", Q_g)
            _uo = load_checkpoint_field(chk, "velocity_obs", V)
            u_obs.dat.data[:] = _uo.dat.data_ro
            if h_clamp_init > 0.0:
                H.interpolate(max_value(H, Constant(h_clamp_init)))
                s.interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))

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

    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    # Composite flow exponent (must match the inversion that produced the
    # MAP file we load above). This branch: n=3 standard Glen (A4_FACTOR=1).
    n_flow_val = _n_flow()
    m_slide_val = float(os.environ.get("ISMIP7_M_SLIDE", "3.0"))
    a4_factor = float(os.environ.get("ISMIP7_A4_FACTOR", a4_factor_default()))
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
    else:
        A_prior_f = Function(Q, name="fluidity_prior").interpolate(A0 * Constant(a4_factor))
        A4_base = A_prior_f
        PETSc.Sys.Print(
            f"  Fluidity prior: checkpoint has no fluidity_prior; using LEGACY "
            f"constant baseline A0*a4_factor = {float(A0) * a4_factor:.2f}"
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
        # zero at h_ocean; ON by default in the forward. The inversion
        # operator (which never passes it) is unchanged, so existing MAPs
        # stay consistent. The gia soft speed limiter (u_lim/k_lim) is OFF by
        # default: its max() kink at u_lim breaks the nleqerr continuation
        # (isolated Jul 18 2026); gia only tolerates it under newtontr with
        # dt-retry. Enable via ISMIP7_U_LIM if a mid-run runaway ever needs a
        # backstop.
        ocean_drag = float(os.environ.get("ISMIP7_OCEAN_DRAG", "1e-2"))
        h_ocean = float(os.environ.get("ISMIP7_H_OCEAN", "10.0"))
        # Speed limiter: structurally present (threshold u_lim > 0) but INERT
        # by default - k_lim is a live Constant at 0 (term vanishes
        # identically; the cold continuation is unaffected, unlike a built-in
        # limiter, which breaks it). The run loop's rescue ladder raises
        # k_lim to ISMIP7_K_LIM for trust-region rescue solves at wall
        # geometries (runaway front nodes), then zeroes it again.
        u_lim = float(os.environ.get("ISMIP7_U_LIM", "2e4"))
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
        if friction == "budd":
            PETSc.Sys.Print(
                f"  Friction: Budd N_hat (exact-zero shelf; delta="
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
        F = build_rc_residual(
            z, theta_f, phi_f, H=h, s=s, b=b, C_w0=C_w0,
            A4_base=A4_base, n_flow=n_flow, n_flow_val=n_flow_val,
            m_slide=m_slide_val, tau_c=tau_c, alpha=alpha_reg, H_ref=H_ref,
            fric_law=friction, N_ref=N_ref,
            nhat_floor=budd_nhat_floor, nhat_cap=budd_nhat_cap,
            alpha_gl=alpha_gl,
            c0=c0_rc, eps_tauc=rc_eps_tauc,
            c_w0_floor=rc_cw0_floor, h_visc_floor=rc_hvisc_floor,
            ocean_drag=ocean_drag, h_ocean=h_ocean, u_lim=u_lim, k_lim=k_lim,
            calving_ids=calving_ids if use_calving_terminus else None,
        )
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

    prob = NonlinearVariationalProblem(
        F, z, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(
        prob,
        solver_parameters=sparams,
        options_prefix="ismip7_diagnostic_",
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

    def solve_diagnostic(label, **metadata):
        """Execute one solve and always emit one compact convergence record."""
        _write_solve_header(label, **metadata)
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
            solver_stats.append(stat)
            PETSc.Sys.Print(
                "=== DIAGNOSTIC RESULT "
                f"{_solve_count:04d} | {label} | reason={reason_name} "
                f"snes_its={stat['snes_iterations']} "
                f"linear_its={stat['linear_iterations']} "
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

    # Restart fast path: a current checkpoint holds the last CONVERGED full
    # mixed state (u, M, tau) at its saved post-transport geometry. Trust that
    # cache after checking that its residual at FULL n/m is finite; solving the
    # unchanged state again is both redundant and pathological at the residual
    # floor. The 2015.2 qualification cache drove such a solve for 200 Newton
    # iterations even though it reached ||F||=4e-5, wasting eight minutes and
    # recording DIVERGED_MAX_IT before five otherwise-clean steps.
    #
    # Install a TIGHT run tolerance derived from the loaded-state residual
    # (1e-6 x ||F(z_loaded)||); geometry changes then use the normal relative
    # convergence path. Never use the cache-acceptance residual itself as a
    # loose persistent atol: that can let later steps "converge" at iteration
    # zero and silently freeze the velocity (the bug that invalidated the first
    # 1873->2014 resume). Older velocity-only checkpoints still take the direct
    # solve path because (u, 0, 0) is not a cached mixed solution.
    restart_solved = False
    if is_restart and u_guess is not None:
        n_flow.assign(n_flow_val)
        m_slide.assign(m_slide_val)
        with assemble(F).dat.vec_ro as _rv:
            fnorm0 = _rv.norm()
        if not np.isfinite(fnorm0):
            raise RuntimeError(
                "Restart checkpoint has a non-finite full-state residual"
            )
        restart_atol_scale = snes_restart_failure_atol_scale()
        restart_atol = restart_atol_scale * fnorm0
        if fnorm0 > 0.0:
            slvr.snes.setTolerances(atol=restart_atol)

        if M_guess is not None and tau_guess is not None:
            restart_solved = True
            PETSc.Sys.Print(
                "Restart full mixed state accepted without a setup solve "
                f"(||F||={fnorm0:.2e}, run atol={restart_atol:.2e})"
            )
        else:
            try:
                solve_diagnostic("restart-loaded-velocity")
                restart_solved = True
                fnorm_conv = slvr.snes.getFunctionNorm()
                if fnorm_conv > 0.0:
                    slvr.snes.setTolerances(
                        atol=snes_atol_scale() * fnorm_conv
                    )
                PETSc.Sys.Print(
                    "Restart velocity-only state re-solved "
                    f"(||F|| {fnorm0:.2e} -> {fnorm_conv:.2e})"
                )
            except fd.ConvergenceError:
                z.assign(z_init)
                restart_solved = True
                PETSc.Sys.Print(
                    "Restart velocity-only solve did not converge; keeping "
                    "the loaded state and handing the step to the rescue "
                    f"ladder (atol={restart_atol:.2e})"
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
                    f"  Continuation diverged at {steps} steps — giving up."
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
        "V": V,
        "Z": Z,
        "z": z,
        "h": h,
        "s": s,
        "b": b,
        "slvr": slvr,
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
        "calving_ids": calving_ids,
        "u_obs": u_obs,
        "friction": friction,
        # Mesh provenance from the source checkpoint (re-stamped into every
        # state checkpoint so warm restarts stay self-describing).
        "lc": chk_lc,
        "lc_coarse": chk_lc_coarse,
        "buffer_m": chk_buffer_m,
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
        # DG0 prognostic thickness state (restart only; else None).
        "h_dg_state": h_dg_state,
        # Resume time (None on a cold start); run_simulation continues the
        # timeline from here instead of the caller's t_start.
        "t_restart": t_restart,
        # Timing-cache identity, if this is a prepared timing restart.
        # Ordinary production checkpoints legitimately omit these fields.
        "checkpoint_metadata": checkpoint_metadata,
    }


def save_model_state(ctx, final_path, t_now, extra_attrs=None):
    r"""Atomically save one self-contained mixed state.

    Ordinary simulation checkpoints and timing caches share this writer so a
    cache cannot silently omit one of the frozen fields required on restart.
    """
    mesh = ctx["mesh"]
    z = ctx["z"]
    h = ctx["h"]
    tmp = final_path + ".tmp"
    with fd.CheckpointFile(tmp, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(ctx["theta"], name="log_friction")
        chk.save_function(ctx["phi"], name="log_fluidity")
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
        if not ctx.get("geom_dg", False):
            h_dg = ctx.get("h_dg_state")
            if h_dg is None:
                raise RuntimeError(
                    "CG1 state checkpoint requested before h_dg was prepared"
                )
            chk.save_function(h_dg, name="thickness_dg")

        chk.set_attr("/", "t_yr", float(t_now))
        chk.set_attr("/", "friction", str(ctx.get("friction", "budd")))
        chk.set_attr(
            "/", "geometry_space", "dg0" if ctx.get("geom_dg") else "cg1"
        )
        if ctx.get("mesh_basename"):
            chk.set_attr("/", "mesh_basename", str(ctx["mesh_basename"]))
        for name in ("lc", "lc_coarse"):
            if ctx.get(name) is not None:
                chk.set_attr("/", name, int(ctx[name]))
        if ctx.get("buffer_m") is not None:
            chk.set_attr("/", "buffer_m", float(ctx["buffer_m"]))
        for name, value in (extra_attrs or {}).items():
            if value is not None:
                chk.set_attr("/", name, value)

    mesh.comm.barrier()
    if mesh.comm.rank == 0:
        os.replace(tmp, final_path)
    mesh.comm.barrier()


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

    # Fixed calving front (ISMIP7_FIXED_FRONT=1): cells that are ice-free
    # in the initial state may not accumulate ice; whatever flows into
    # them is removed each step and tallied as calving flux. Without this
    # a buffered mesh has NO calving sink (~1300 Gt/yr in reality) and the
    # sheet must gain mass. Only meaningful when the initial state is the
    # true BedMachine geometry (RC mode / h_clamp_init=0) — with a clamped
    # initial state every cell has ice and the mask is empty.
    fixed_front = os.environ.get("ISMIP7_FIXED_FRONT") is not None
    front_hmin = float(os.environ.get("ISMIP7_FRONT_HMIN", "1.0"))
    beyond_front = None
    cell_area = assemble(fd.TestFunction(Q_dg) * dx).dat.data_ro.copy()
    if fixed_front:
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
            f"  Fixed calving front: {n_beyond} initially ice-free cells "
            f"masked (h < {front_hmin} m)"
        )

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
    # Modes: "div" cancels only the flux divergence (t=0 tendency = SMB-melt,
    # gia-style); any other value ("1"/"balance") also subtracts the initial
    # forcing, so the t=0 tendency is EXACTLY ZERO - a balanced control in
    # the ISMIP6 ctrl_proj sense, against which projections difference
    # cleanly. Both freeze the correction at t=0.
    amb_mode = os.environ.get("ISMIP7_APPARENT_MB")
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
            if amb_mode != "div":
                # balanced control: evaluate the t=0 forcing and fold it in
                if forcing_callback is not None:
                    forcing_callback(ctx, t_start + dt)
                b_smb = assemble((accum - ocean_melt) * phi_dg * dx)
                a_ref.dat.data[:] -= b_smb.dat.data_ro / cell_area
            if beyond_front is not None:
                # the fixed-front tally stays the sink for flux into the
                # initially ice-free cells; do not absorb it into a_ref
                a_ref.dat.data[beyond_front] = 0.0
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

    results = []
    ctx["results"] = results
    ctx["failure"] = None

    # Crash-safe timeseries: append each row and flush, so a reboot keeps the
    # budget-audit history (it used to be dumped only at completion). On a
    # warm restart, drop any rows at/after the resume year, then append.
    csv_fn = os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_timeseries.csv")
    csv_header = ("year,vaf_mm_sle,mass_gt,smb_gtyr,melt_gtyr,"
                  "outflux_gtyr,calv_gt,clamp_gt,resid_gt,amb_gtyr\n")
    csv_f = None
    if mesh.comm.rank == 0:
        if t_restart is not None and os.path.exists(csv_fn):
            with open(csv_fn) as _cf:
                _lines = _cf.readlines()
            kept = ([_lines[0]] if _lines and _lines[0].startswith("year")
                    else [csv_header])
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

    def _write_csv_row(row):
        if csv_f is None:
            return
        csv_f.write(
            f"{row[0]:.1f},{row[1]:.6f},{row[2]:.2f},"
            + ",".join(f"{v:.4f}" for v in row[3:]) + "\n"
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

    def _advance(dt_local, label):
        r"""One transport advance of dt_local with the CURRENT velocity
        (transport-first ordering: the velocity was solved at the current
        geometry). Mutates h_dg and the derived CG fields; returns the
        advance's mass tallies [Gt]."""
        if legacy_transport:
            h_dg.project(h)
        h_dg_old.assign(h_dg)

        src = accum - ocean_melt
        if a_ref is not None:
            src = src + a_ref
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

        # Floor to h_clamp, EXCEPT beyond the fixed front: those cells are
        # outside the ice domain, so flooring them would hand the mask below
        # h_clamp of fresh ice to re-calve every step and report as terminus
        # discharge. There the floor is zero and the front stays a pure sink.
        data = h_dg.dat.data
        floor = np.full_like(data, h_clamp)
        if fixed_front:
            floor[beyond_front] = 0.0
        np.maximum(data, floor, out=data)
        m2 = float(assemble(h_dg * dx)) * rho_gt
        clamp_gt = m2 - m1                                       # Gt added by DG floor

        calv_gt = 0.0
        if fixed_front:
            data = h_dg.dat.data
            calv_gt = mesh.comm.allreduce(
                float((data[beyond_front] * cell_area[beyond_front]).sum())
            ) * rho_gt
            data[beyond_front] = 0.0

        _lift_h()
        m3 = m2 - calv_gt                                        # ∫h preserved by projection
        mass_now = float(assemble(h * dx)) * _RHO_I_SI / 1e12
        clamp_cg_gt = mass_now - m3          # Gt added by the CG floor after projection
        return {
            "out_gt": out_gt,
            "calv_gt": calv_gt,
            "clamp_gt": clamp_gt + clamp_cg_gt + limit_gt,
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
    for k in range(1, nsteps + 1):
        t_step_start = perf_counter()
        t_yr = t_start + k * dt

        if forcing_callback is not None:
            forcing_callback(ctx, t_yr)

        # Forcing-field integrals are constant within the step.
        smb_rate = float(assemble(accum * dx)) * rho_gt          # Gt/yr
        melt_rate = float(assemble(ocean_melt * dx)) * rho_gt    # Gt/yr
        amb_rate = (float(assemble(a_ref * dx)) * rho_gt
                    if a_ref is not None else 0.0)               # Gt/yr

        z_entry.assign(z)
        h_dg_entry.assign(h_dg)
        tallies = None
        for m in SUBCYCLES:
            if m > 1:
                PETSc.Sys.Print(
                    f"  Step {k}: subcycling x{m} (dt={dt / m:.4g})..."
                )
                h_dg.assign(h_dg_entry)
                _lift_h()
                z.assign(z_entry)
            acc = {"out_gt": 0.0, "calv_gt": 0.0, "clamp_gt": 0.0}
            ok = True
            for _j in range(m):
                sub = _advance(
                    dt / m, f"step-{k}-substep-{_j + 1}/{m}"
                )
                for key in acc:
                    acc[key] += sub[key]
                if not _solve_with_rescue(k):
                    ok = False
                    break
            if ok:
                if m > 1:
                    PETSc.Sys.Print(f"  Step {k}: completed via x{m} subcycle")
                tallies = acc
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
            h_dg.assign(h_dg_entry)
            _lift_h()
            z.assign(z_entry)   # checkpoint the last converged pair
            break

        t_elapsed = perf_counter() - t_step_start

        haf = Function(Q_g).interpolate(max_value(s - s_float, Constant(0.0)))
        vaf = float(assemble(haf * dx)) * _RHO_I_SI / 1e12 / 362.5
        total_mass = float(assemble(h * dx)) * _RHO_I_SI / 1e12

        out_rate = tallies["out_gt"] / dt                        # Gt/yr
        calv_gt = tallies["calv_gt"]
        clamp_all = tallies["clamp_gt"]
        dm = total_mass - mass_prev
        resid_gt = dm - (
            (smb_rate - melt_rate + amb_rate - out_rate) * dt
            + clamp_all - calv_gt
        )
        mass_prev = total_mass

        results.append((t_yr, vaf, total_mass, smb_rate, melt_rate,
                        out_rate, calv_gt, clamp_all, resid_gt,
                        amb_rate))
        _write_csv_row(results[-1])

        _field_diagnostics(f"step-{k}-t={t_yr:.6g}")

        if k % output_interval == 0 or k == 1:
            _amb_txt = f"amb={amb_rate:+.0f} " if a_ref is not None else ""
            PETSc.Sys.Print(
                f"  t={t_yr:.1f}  VAF={vaf:.4f} mm SLE  "
                f"mass={total_mass:.1f} Gt  [{t_elapsed:.1f}s]\n"
                f"      budget [Gt/yr]: SMB={smb_rate:+.0f} melt={-melt_rate:+.0f} "
                f"{_amb_txt}"
                f"outflux={-out_rate:+.0f} calv={-calv_gt/dt:+.0f} "
                f"clamp={clamp_all/dt:+.1f} "
                f"dM/dt={dm/dt:+.0f} resid={resid_gt/dt:+.2f}"
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

    ctx["transient_seconds"] = perf_counter() - ctx["transient_t0"]
    PETSc.Sys.Print(f"\n{experiment_name} simulation complete.")

    # Final state is a self-contained checkpoint too (a valid restart source).
    # Save the ACTUAL last year so a resume after an early stop continues from
    # where it really stopped, not the nominal t_end.
    final_fn = os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_final.h5")
    last_t = results[-1][0] if results else t_start
    save_model_state(ctx, final_fn, last_t)
    PETSc.Sys.Print(f"Saved: {final_fn}")

    if csv_f is not None:
        csv_f.close()
    PETSc.Sys.Print(f"Saved: {csv_fn}")

    return results
