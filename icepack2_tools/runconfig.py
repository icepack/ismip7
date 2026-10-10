r"""Run-shaping knobs: the discretization and physics choices that decide
which MAP a run writes, and which MAP a forward is allowed to load.

The inversion, the forward, the mesh pipeline, the probes, the preflight gate
and the per-core report all need these. Each used to re-declare its own
literal default, which is how ISMIP7_LC came to mean 2500 in ``simulation.py``
and ``preflight.py``, 8000 in ``inversion_icepack2.py`` and 32000 in
``thermo_prior.py``: with the variable unset, the gate blessed
``inversion_icepack2_budd_n3_dg0_2500.h5`` while the inversion would have
built and written the 8000 MAP. That is exactly the mismatch ``naming.py``
exists to prevent, so the values it builds names from are owned here and every
reader calls the accessors below rather than repeating a literal. A script
that genuinely needs a different value passes it at the call site or exports
the variable, so the divergence is visible instead of hiding in a default.

A knob left at its default is also absent from ``os.environ``, so
``core_report.py`` resolves the run-env block through this module: the report
is the only committed record of a run, so it has to state the value the run
used, not only the ones that happened to be exported.

Pure Python: importable without Firedrake so the preflight and the report stay
fast.
"""

import hashlib
import json
import os

# 1000 m / 10 km is the production pair (``antarctica_10000_1000_buffered20000``),
# the submission mesh since 25 September 2026 (issue 20): the finest mesh the
# Quartz timing matrix carries through a 285-year run in two days at dt 0.05,
# under scpc_gamg on 64 ranks (antarctica/TIMING_MATRIX_QUARTZ_SCPC_GAMG.md).
# It is the pair the batch runners export (batch_runners/site_env.sh) and the
# README documents; until 2026-09-19 that was 2500 m / 64 km. The old 8000 and
# 32000 module-level defaults were dev-probe leftovers; a coarse probe exports
# ISMIP7_LC / ISMIP7_LC_COARSE instead of disagreeing with the gate about what
# "unset" means.
LC_DEFAULT = "1000"
LC_COARSE_DEFAULT = "10000"
# Metres the ice outline is pushed outward before meshing, so the calving
# front sits inside the domain beside ice-free cells (the `_buffered<N>` of
# mesh_naming.mesh_basename). 20 km is the production mesh's (issue 20); the
# 32 km probes export 0. The outline extraction used to default to 0 while the
# mesh names defaulted to 20000, so a bare call could build an unbuffered mesh
# under a buffered name.
BUFFER_M_DEFAULT = "20000"
# The ice edge whose marine front a buffered mesh's nodes and edges lie on
# (the `_front<edge>` of naming.mesh_basename): "bm", BedMachine's own, or
# "none" for a mesh that follows no front. Issue #167: on a mesh that does
# not follow it, the front crosses cells, and vertex sampling gives them a
# fraction of its thickness.
MESH_FRONTS = ("none", "bm")
MESH_FRONT_DEFAULT = "none"
# The production forward step [yr], chosen with the mesh (issue 20). The timing
# matrix's rule gives 0.05 at 1000 m. At 0.05 a 1 km control from a transferred
# 2 km Budd MAP diverged in 2016.1 at Rice, and on Quartz it grew a two-step
# grounded/floating oscillation at the Lambert confluence (MAP_CHECK.md). At
# 0.025 the same control ran five years, and Rice's 1 km historicals ran 78
# (CESM2-WACCM) and 66 (MRI-ESM2-0) model years with no rescue step.
# batch_runners/projection.sbatch exports the same value, and a test holds the
# two equal.
DT_DEFAULT = "0.025"
GEOMETRY_SPACE_DEFAULT = "dg0"
FRICTION_DEFAULT = "budd"
# The closed set friction() accepts; an unknown spelling is an error at
# startup, never a silent fall-through to another law's block.
FRICTION_LAWS = ("budd", "regularized_coulomb", "budd_legacy")
# THIS BRANCH (antarctica-n3) runs standard Glen n=3. An inversion and every
# forward that loads its MAP must agree on this.
N_FLOW_DEFAULT = "3.0"
# How the Budd law (dual_friction.budd_nhat) zeroes shelf friction. Before
# 2026-09-13 it tested the sign of N = max(p_I - p_W, 0), a roundoff residue
# on floating ice, and the delta floor lifted every roundoff-positive shelf
# cell to the friction cap (hoffmaao/antarctica e602705). "haf" gates on
# height above flotation. Stamped into every Budd state checkpoint and required
# by the timing-cache manifest, so a state solved under the old gate can never
# seed a fixed-law lane.
BUDD_SHELF_GATE = "haf"

# Exact-mesh timing caches initialize DG0 geometry from BedMachine on the
# TARGET mesh. The raster is first sampled into CG1 and then L2-projected to
# DG0, so each cell stores an average rather than one centroid pixel. Keep the
# method name stable: it is stamped into cache provenance and changing the
# construction must invalidate old caches.
TARGET_MESH_GEOMETRY_METHOD = "target-native-bedmachine-cell-average-v1"
# The same with BedMachine's front cells rebuilt (geometry.front_cells, issue
# #167); target_mesh_geometry_method picks one by raster sampling.
TARGET_MESH_GEOMETRY_METHOD_FRONT = "target-native-bedmachine-cell-average-front-v1"

GEOMETRY_SPACES = ("dg0", "cg1")

# How a raster (BedMachine) is put onto a DG0 geometry cell.
#   vertex    - icepack's bilinear sample at the three CG1 vertices, then the
#               L2 projection of that linear interpolant (= the mean of the
#               3 vertex values). The pre-Sep-2026 behaviour. A 20 km interior
#               cell sees 3 of its ~1600 BedMachine pixels.
#   cell_mean - the mean of the raster over the cell itself, sampled on an
#               equal-area sub-triangle lattice at pixel density
#               (geometry.raster_cell_mean).
#   vertex_front - vertex, except that BedMachine's own mask decides which
#               cells hold ice at the marine front, and the front cells take
#               BedMachine's mean thickness and bed over their ice
#               (geometry.front_cells, issue #167). Vertex sampling alone
#               gives a cell the front crosses a fraction of the front's
#               thickness: 40 m against 163 m on the 2 km buffered mesh, and
#               a vertex on the front of a _frontbm mesh samples a blend.
# MAPs record the method used; the forward reads it back from the MAP.
# vertex_front is the default since 8 October 2026 (IU, issue #167): IU's
# final MAPs are refitted under it, and it has its own melt calibration.
FRONT_RASTER_SAMPLES = ("vertex_front",)
RASTER_SAMPLES = ("vertex", "cell_mean") + FRONT_RASTER_SAMPLES
RASTER_SAMPLE_DEFAULT = "vertex_front"


# Floor-cell coercivity drags of the dual-friction residual (dual_friction.py
# ``ocean_drag`` / ``u_lim``). The forward has applied them since Jul 2026; the
# inversion never passed them, so a MAP that solved the inversion's F was not
# a solution of the forward's F in the ice-free buffer cells: the same mixed
# state measured ||F||=1.3e1 in the inversion and 1.5e10 in the forward
# (2026-09-14), and every strict scout that trusted it ran away within four
# steps. Both now read the knobs here. ``u_lim`` is only the threshold of the
# soft speed limiter; its gain ``k_lim`` stays a live Constant at 0 in the
# forward (raised for rescue solves only) and is passed as 0 by the inversion.
OCEAN_DRAG_DEFAULT = "1e-2"   # MPa yr/m, linear drag ramping to zero at H_OCEAN
H_OCEAN_DEFAULT = "10.0"      # m
U_LIM_DEFAULT = "2e4"         # m/yr, ~5x the fastest observed Antarctic flow


def residual_stabilizers():
    r"""``ocean_drag``, ``h_ocean`` and ``u_lim`` keyword arguments for
    ``build_rc_residual``; identical in the inversion and the forward."""
    return {
        "ocean_drag": float(os.environ.get("ISMIP7_OCEAN_DRAG", OCEAN_DRAG_DEFAULT)),
        "h_ocean": float(os.environ.get("ISMIP7_H_OCEAN", H_OCEAN_DEFAULT)),
        "u_lim": float(os.environ.get("ISMIP7_U_LIM", U_LIM_DEFAULT)),
    }


# ``ISMIP7_MESH=checkpoint`` names the mesh embedded in the MAP or restart
# file. site_env.sh always exports a derived .msh path, so a job submitted
# through it (submit.sh projection) has no other way to run MAP-native.
MESH_FROM_CHECKPOINT = "checkpoint"


def mesh_override():
    r"""``ISMIP7_MESH`` as a compute-mesh override, or None.

    Unset, empty and the sentinel ``checkpoint`` all mean: solve on the mesh
    the MAP or restart checkpoint carries. Anything else is the path of the
    mesh to solve on, with the checkpoint kept as the interpolation source
    (the timing matrix and the 2 km to 1 km transfer).
    """
    value = os.environ.get("ISMIP7_MESH", "").strip()
    if value in ("", MESH_FROM_CHECKPOINT):
        return None
    return value


def inversion_mesh_source(derived):
    r"""Where the inversion reads its mesh: ``(path, from_checkpoint)``.

    Unset or empty ``ISMIP7_MESH`` is ``derived``, the .msh the mesh knobs
    name. The sentinel ``checkpoint`` is the mesh inside ``ISMIP7_WARM_START``,
    so an inversion can continue a MAP on its own mesh when the .msh did not
    travel with it (Rice's 2 km MAPs are released without theirs, and another
    site's build of the same name is another triangulation). It needs a warm
    start to read the mesh from. Anything else is a .msh path.
    """
    value = os.environ.get("ISMIP7_MESH", "").strip()
    if not value:
        return derived, False
    if value != MESH_FROM_CHECKPOINT:
        return value, False
    warm = os.environ.get("ISMIP7_WARM_START", "").strip()
    if not warm:
        raise ValueError(
            f"ISMIP7_MESH={MESH_FROM_CHECKPOINT} reads the mesh from the warm "
            f"start, and ISMIP7_WARM_START is not set.")
    return warm, True


def eval_continuation():
    r"""``ISMIP7_EVAL_CONTINUATION`` with the direct forward disabled:
    whether every inversion evaluation ramps n_flow and m_slide from 1 in
    five annotated solves. ``0`` solves once at the full exponents from the
    previous evaluation's state, which the startup ramp has already brought
    there. The objective depends on the full-n solution alone, so it is the
    same either way; the cost differs. On unless ``0``."""
    return _int_flag("ISMIP7_EVAL_CONTINUATION", True)


TRANSFER_FILL_MODES = ("extend", "constant")


def transfer_fill():
    r"""``ISMIP7_TRANSFER_FILL``: what the controls and the fluidity prior
    take on the dofs of a compute mesh beyond the mesh they were read from (a
    MAP loaded by a forward, an inversion's warm start). ``extend`` (the
    default): theta, phi and alpha continue harmonically from the source
    outline, and the fluidity prior's logarithm does too
    (``icepack2_tools.transfer.harmonic_extension``). ``constant``: theta =
    phi = 0 and the constant baseline prior ``A0 * a4_factor``."""
    mode = os.environ.get("ISMIP7_TRANSFER_FILL", "extend").strip().lower()
    if mode not in TRANSFER_FILL_MODES:
        raise ValueError(
            f"ISMIP7_TRANSFER_FILL must be one of {', '.join(TRANSFER_FILL_MODES)}, "
            f"not {mode!r}")
    return mode


def lc():
    r"""Target edge length [m] in the refined region of the mesh."""
    return int(os.environ.get("ISMIP7_LC", LC_DEFAULT))


def lc_coarse():
    r"""Target edge length [m] in the coarse region of the mesh."""
    return int(os.environ.get("ISMIP7_LC_COARSE", LC_COARSE_DEFAULT))


def dt():
    r"""Forward time step [yr]."""
    return float(os.environ.get("ISMIP7_DT", DT_DEFAULT))


def buffer_m():
    r"""Outline buffer [m] a mesh is built with, and named by."""
    return float(os.environ.get("ISMIP7_BUFFER_M", BUFFER_M_DEFAULT))


def mesh_front():
    r"""``ISMIP7_MESH_FRONT``: the ice edge whose marine front the mesh
    follows (``bm``), or None for ``none``."""
    value = os.environ.get("ISMIP7_MESH_FRONT", MESH_FRONT_DEFAULT).strip().lower()
    if value not in MESH_FRONTS:
        raise ValueError(
            f"ISMIP7_MESH_FRONT must be one of {MESH_FRONTS}, not {value!r}")
    return None if value == "none" else value


def geometry_space():
    r"""Discretization of h/s/b, ``'dg0'`` or ``'cg1'``. Validated here so
    every reader rejects the same set."""
    value = os.environ.get(
        "ISMIP7_GEOMETRY_SPACE", GEOMETRY_SPACE_DEFAULT).lower()
    if value not in GEOMETRY_SPACES:
        raise ValueError(
            f"ISMIP7_GEOMETRY_SPACE must be 'dg0' or 'cg1', got {value!r}"
        )
    return value


def raster_sample():
    r"""How BedMachine is sampled onto a DG0 cell, one of RASTER_SAMPLES."""
    value = os.environ.get(
        "ISMIP7_RASTER_SAMPLE", RASTER_SAMPLE_DEFAULT).lower()
    if value not in RASTER_SAMPLES:
        raise ValueError(
            f"ISMIP7_RASTER_SAMPLE must be one of {RASTER_SAMPLES}, got {value!r}"
        )
    return value


def warm_start_geometry(*, same_mesh, same_lake, warm_sampling, run_sampling):
    r"""Whether an inversion takes its warm start's thickness, bed and
    surface. By default it does on the warm start's mesh when the warm start
    was built with this run's lake_ice_base and raster sampling;
    ``ISMIP7_WARM_START_GEOMETRY=0/1`` overrides that. A taken geometry keeps
    the sampling it was built with and the MAP records the run's
    (``raster_sample``), so taking one sampled another way is refused."""
    default = same_mesh and same_lake and warm_sampling == run_sampling
    take = os.environ.get(
        "ISMIP7_WARM_START_GEOMETRY", "1" if default else "0").strip() != "0"
    if take and warm_sampling != run_sampling:
        raise ValueError(
            f"ISMIP7_WARM_START_GEOMETRY takes the geometry of a warm start "
            f"that records raster_sample={warm_sampling}, and this run samples "
            f"with {run_sampling}. A taken geometry keeps the sampling it was "
            f"built with and the MAP records the run's: set "
            f"ISMIP7_RASTER_SAMPLE={warm_sampling}.")
    return take


def warm_start_fluidity():
    r"""``ISMIP7_WARM_START_FLUIDITY=<MAP>``: an inversion warm-started from
    ``ISMIP7_WARM_START`` takes its log fluidity, and the fluidity prior it
    is a deviation from, from this other MAP; theta and the state stay the
    warm start's. Fluidity does not depend on the friction law, so RC's refit
    under the front-cell rule can start from the fluidity Budd's refit fitted
    on the same geometry (issue #167). None when unset; a chain link resuming
    its own checkpoint drops it (inversion.sbatch)."""
    path = os.environ.get("ISMIP7_WARM_START_FLUIDITY", "").strip()
    if not path:
        return None
    if not os.path.isfile(path):
        raise FileNotFoundError(f"ISMIP7_WARM_START_FLUIDITY={path}: no such file")
    return path


def warm_start_front_extend():
    r"""``ISMIP7_WARM_START_FRONT_EXTEND=1``: a refit under a front sampling
    from a warm start sampled without it continues theta and phi harmonically
    over the cells the front rule rebuilt or emptied, and the inversion sets
    them afresh there (issue #167). Every other node keeps the warm start's
    value, so the band becomes a harmonic blend of the ice upstream and the
    nodes seaward of it; on a synthetic shelf it came out as the linear
    interpolation between the two sides. The warm start fitted them
    against its own front: RC's final MAP left a band stiffer and with more
    friction than the rest of the ice (log fluidity -1.33 against -0.17, log
    friction +1.06 against -0.02), which the rule makes four times thicker.
    A phi taken from a MAP fitted under a front sampling
    (`warm_start_fluidity`) is kept. For RC's refit the continuation left the
    start worse: from the ef2 state the residual began at ||F|| 4.1e14
    against 1.7e14 without it, and the first forward solve failed (job
    11884485, with the band blended between both sides as above); RC's
    refit starts from Budd's fluidity and state instead. The band is found
    by comparing this run's thickness with a vertex sample, so only a warm
    start sampled with vertex is continued; `front_band_extends` refuses any
    other sampling."""
    return os.environ.get("ISMIP7_WARM_START_FRONT_EXTEND", "0").strip() not in ("", "0")


def front_band_extends(*, same_mesh, geometry_taken, run_sampling, warm_sampling):
    r"""Whether an inversion continues the front band
    (`warm_start_front_extend`): the knob is on, this run's sampling is a
    front sampling, and the warm start is on this mesh with its geometry not
    taken. A warm start already sampled under a front sampling has no band.
    The band is the cells where this run's thickness differs from a vertex
    sample, so a warm start sampled any other way (cell_mean, say) would flag
    interior cells too, and it is refused."""
    if not (warm_start_front_extend() and raster_front(run_sampling)
            and same_mesh and not geometry_taken):
        return False
    if raster_front(warm_sampling):
        return False
    base = raster_base_method(run_sampling)
    if str(warm_sampling).lower() != base:
        raise ValueError(
            f"ISMIP7_WARM_START_FRONT_EXTEND=1 finds the front band by comparing "
            f"this run's thickness with a {base} sample, and the warm start "
            f"records raster_sample={warm_sampling}: that comparison would flag "
            f"interior cells. Unset the knob or warm start from a {base} MAP.")
    return True


def front_band_controls(fluidity_sampling=None):
    r"""The controls ``ISMIP7_WARM_START_FRONT_EXTEND`` continues over the
    front band: theta, and phi unless it came from a MAP fitted under a front
    sampling (``fluidity_sampling``, the raster sampling of the
    ``ISMIP7_WARM_START_FLUIDITY`` MAP; None without one), whose phi already
    belongs to the rule's front. Each is a harmonic blend over the band of
    its values on the ice upstream and on the nodes seaward of the band,
    which keep the warm start's values."""
    if fluidity_sampling is not None and raster_front(fluidity_sampling):
        return ("theta",)
    return ("theta", "phi")


def ramp_slide_fixed():
    r"""``ISMIP7_RAMP_SLIDE_FIXED=1``: the inversion's startup ramp climbs the
    flow exponent from 1 with the sliding exponent held at its target. RC's
    refit under the front-cell rule on the 2 km mesh could not ramp from
    n = m = 1: under full MUMPS (job 11869835) Newton ran 200 iterations at
    n = m = 1 in each of three rungs and diverged (||F|| 3.0e11 to 5.3e11),
    and under scpc_gamg and scpc_mumps (jobs 11883518 and 11883519) the
    first linear solve at n = m = 1 failed (1,000 Krylov iterations from
    ||F|| 8.9e9). IU tried this ramp there and stopped it as too costly:
    after 2 h 29 min (job 11884073) it stood at n = 1.13, its first step
    having taken three rungs (37, 34 and 30 min). The refit starts from the
    fluidity and state of Budd's refit instead (``ISMIP7_WARM_START_FLUIDITY``,
    ``ISMIP7_WARM_START_STATE=fluidity``; issue #167)."""
    return os.environ.get("ISMIP7_RAMP_SLIDE_FIXED", "0").strip() not in ("", "0")


def warm_start_state(*, geometry_taken, same_mesh):
    r"""Whether an inversion loads its warm start's mixed state, and whether
    it is only a first guess: ``(load, guess)``.

    The state comes with a taken geometry (`warm_start_geometry`), as the
    solution of the same equations. ``ISMIP7_WARM_START_STATE=1`` also loads
    it on the warm start's own mesh when the geometry is not taken, as the
    first guess on this run's geometry (a warm start sampled another way,
    issue #167). The first evaluation's forward then solves it at the full
    exponents. For RC's refit under vertex_front it did not work: from the
    ef2 state the run started at ||F|| 1.7e14 against the recorded 9.3e-3
    and the first forward solve failed (jobs 11883148 and 11883149). On
    another mesh the knob is refused. ``ISMIP7_WARM_START_STATE=fluidity``
    loads no state from the warm start: the first guess is the state of the
    MAP ``ISMIP7_WARM_START_FLUIDITY`` names (`warm_start_state_fluidity`),
    which is how RC's refit starts. A chain link resuming its own checkpoint
    drops the knob (inversion.sbatch)."""
    mode = os.environ.get("ISMIP7_WARM_START_STATE", "0").strip().lower()
    if mode == "fluidity":
        if warm_start_fluidity() is None:
            raise ValueError(
                "ISMIP7_WARM_START_STATE=fluidity takes the state of the MAP "
                "ISMIP7_WARM_START_FLUIDITY names, and that is unset")
        return False, False
    if geometry_taken:
        return True, False
    want = mode not in ("", "0")
    if want and not same_mesh:
        raise ValueError(
            "ISMIP7_WARM_START_STATE=1 loads the warm start's mixed state on "
            "its own mesh as a first guess, and this run's mesh differs: unset "
            "it.")
    return want, want


def warm_start_state_fluidity(*, same_mesh):
    r"""``ISMIP7_WARM_START_STATE=fluidity``: the first guess is the mixed
    state of the MAP ``ISMIP7_WARM_START_FLUIDITY`` names, the solution under
    the fluidity this run starts from, on that MAP's own mesh only
    (``same_mesh``). RC's refit under the front-cell rule takes the state of
    Budd's refit with its fluidity: RC's own ef2 state started at ||F||
    1.7e14 with RC's fluidity (job 11883148) and 2.5e21 with Budd's (job
    11884749), its stresses balanced against the other fluidity
    (issue #167)."""
    if os.environ.get("ISMIP7_WARM_START_STATE", "0").strip().lower() != "fluidity":
        return False
    if not same_mesh:
        raise ValueError(
            "ISMIP7_WARM_START_STATE=fluidity loads the state of the "
            "ISMIP7_WARM_START_FLUIDITY MAP on its own mesh, and this run's "
            "mesh differs")
    return True


def raster_front(method):
    r"""Whether the raster sampling ``method`` rebuilds the front cells."""
    return str(method).lower() in FRONT_RASTER_SAMPLES


def raster_base_method(method):
    r"""How a single raster is put on a cell under ``method``: a front
    sampling samples every raster by its vertices and changes only the front
    cells' thickness and bed afterwards."""
    return "vertex" if raster_front(method) else str(method).lower()


def target_mesh_geometry_method(method):
    r"""The provenance name of geometry rebuilt from BedMachine on a target
    mesh under the raster sampling ``method``: vertex keeps the name its
    caches carry, a front sampling gets its own, so a cache built one way is
    never read as the other."""
    if raster_front(method):
        return TARGET_MESH_GEOMETRY_METHOD_FRONT
    return TARGET_MESH_GEOMETRY_METHOD


def forward_raster_sample(recorded, transfer, source="the MAP"):
    r"""The raster sampling a forward builds its geometry with.

    On the MAP's own mesh the sampling is the one the MAP records (``vertex``
    for a MAP older than the record), which an explicitly set
    ``ISMIP7_RASTER_SAMPLE`` may repeat and may not change: the geometry is
    the MAP's, or for a relaxed MAP BedMachine's rebuilt with that sampling,
    and the controls and the melt calibration follow it. A
    forward on another mesh (``transfer``) rebuilds the geometry from
    BedMachine there, with ``ISMIP7_RASTER_SAMPLE`` or its default: a MAP's
    controls carry over, and the front its new mesh holds is that mesh's own
    (issue #167)."""
    if transfer:
        return raster_sample()
    recorded = _recorded(recorded)
    method = "vertex" if recorded is None else str(recorded).lower()
    env = os.environ.get("ISMIP7_RASTER_SAMPLE")
    if env and env.lower() != method:
        raise RuntimeError(
            f"ISMIP7_RASTER_SAMPLE={env} but {source} was sampled with "
            f"{method}: a forward on its MAP's mesh follows the MAP's "
            f"sampling, which its controls and melt calibration were fitted under")
    return method


def friction():
    r"""Friction law: ``budd``, ``regularized_coulomb`` or ``budd_legacy``.
    Validated here so every reader rejects the same set: theta means a
    different thing under each law, and a mistyped value used to run the
    legacy action branch without a word."""
    value = os.environ.get("ISMIP7_FRICTION", FRICTION_DEFAULT).strip().lower()
    if value not in FRICTION_LAWS:
        raise ValueError(
            f"ISMIP7_FRICTION must be one of {FRICTION_LAWS}, got {value!r}"
        )
    return value


def n_flow():
    r"""Glen flow-law exponent."""
    return float(os.environ.get("ISMIP7_N_FLOW", N_FLOW_DEFAULT))


# Calving front. ``none`` (the default) runs no level set: on a buffered mesh
# the front advances freely and never calves, and ISMIP7_FIXED_FRONT decides
# whether the ice that flows past the t=0 extent is removed. Any other value
# names a law in icepack_tools.calving, the one registry every project that
# runs a front selects from (ISMIP7_CALVING_MODULE registers a law from a
# file first), made with the parameters in ISMIP7_CALVING_PARAMS.
CALVING_DEFAULT = "none"
#: Knobs that named one law's parameters before the laws had one home. Their
#: values now go in ISMIP7_CALVING_PARAMS; set, they are refused rather than
#: silently ignored.
RETIRED_CALVING_KNOBS = {
    "ISMIP7_CALVING_SIGMA_MAX_GROUNDED": "sigma_max_gr",
    "ISMIP7_CALVING_SIGMA_MAX_FLOATING": "sigma_max_fl",
}


FRACTURE_MODES = ("none", "mask", "mask_front")
FRACTURE_MASK_MODES = ("mask", "mask_front")     # the modes that read the collapse mask
FRACTURE_DEFAULT = "none"


def fracture():
    r"""``ISMIP7_FRACTURE``: how the ISMIP7 ice-shelf collapse forcing is
    applied. ``none`` (default) loads nothing. Both mask modes act on FLOATING
    ice only and book what they remove as calving (protocol path C,
    discussions #30 and #33; no mask exists for historical or OCX, so those
    runs see nothing either way), and they are the two end-members the
    modelling groups arrived at in discussion #30:

    ``mask`` empties every floating cell the year's mask flags, wherever it
    is. The masks flag the Ross and Filchner-Ronne shelves near their
    grounding lines first, so this opens holes far behind the front which the
    momentum balance treats as open ocean.

    ``mask_front`` empties a flagged floating cell only once open water has
    reached it through other flagged cells, so a shelf collapses from its
    front and nothing happens until the flagged region touches it (see
    ``icepack2_tools.front.front_connected``).

    A stress-gated variant (Lai et al. 2020) is not implemented."""
    value = os.environ.get("ISMIP7_FRACTURE", FRACTURE_DEFAULT).lower()
    if value not in FRACTURE_MODES:
        raise ValueError(f"ISMIP7_FRACTURE must be one of {FRACTURE_MODES}, got {value!r}")
    return value


OCX_FORCING_MODES = ("protocol", "stopgap")
OCX_FORCING_DEFAULT = "protocol"


def ocx_forcing():
    r"""``ISMIP7_OCX_FORCING``: what core 11 runs on. ``protocol`` (default)
    is the ISMIP7 OCX product, RACMO2.3p2-ERA downscaled SMB and the
    expert-judgment ocean, and the run refuses to start without it.
    ``stopgap`` is what the core ran on before the product was readable here:
    RACMO2.4p1 actual-year SMB and the constant OI ocean climatology. It used
    to be the silent fallback, which is how a core ran on it for weeks with
    the real product on disk; it is now something a run has to ask for."""
    value = os.environ.get("ISMIP7_OCX_FORCING", OCX_FORCING_DEFAULT).lower()
    if value not in OCX_FORCING_MODES:
        raise ValueError(f"ISMIP7_OCX_FORCING must be one of {OCX_FORCING_MODES}, got {value!r}")
    return value


def ocx_ocean():
    r"""``ISMIP7_OCX_OCEAN``: which of the four expert-judgment OCX ocean
    scenarios to read. ``main`` (default) is the core experiment's; ``cold``,
    ``warm`` and ``vary`` are its sensitivity members."""
    from .forcing import OCX_OCEAN_VARIANTS
    value = os.environ.get("ISMIP7_OCX_OCEAN", "main").lower()
    if value not in OCX_OCEAN_VARIANTS:
        raise ValueError(f"ISMIP7_OCX_OCEAN must be one of {OCX_OCEAN_VARIANTS}, got {value!r}")
    return value


def ismip7_output():
    r"""``ISMIP7_OUTPUT``: record the ISMIP7 yearly fields and scalars.

    ``1`` enables it; ``0`` and the empty string disable it. The value set is
    closed, like ``fracture`` and ``apparent_mb_mode``: there is one spelling
    each way, and anything else raises rather than silently deciding whether
    a submission gets written.
    """
    value = (os.environ.get("ISMIP7_OUTPUT") or "").strip()
    if value in ("", "0"):
        return False
    if value == "1":
        return True
    raise ValueError(
        f"ISMIP7_OUTPUT must be 1 to enable or 0/empty to disable, "
        f"got {value!r}"
    )


# The SMB-elevation feedback from the ISMIP7 SMB gradient ``dacabfdz``
# (forcing.SMBElevationFeedback), decided on icepack/ismip7 issue 116. On by
# default here, so every driver and the report resolve the same value whether
# or not a runner exports it.
SMB_ELEVATION_FEEDBACK_DEFAULT = "1"


def smb_elevation_feedback():
    r"""``ISMIP7_SMB_ELEVATION_FEEDBACK``: add ``dacabfdz`` times the surface
    change since the chain's initial state to the SMB.

    Unset or ``1`` enables it; ``0`` and the empty string disable it. The
    value set is closed, like ``ismip7_output``: anything else raises rather
    than silently deciding whether a run carries the feedback.
    """
    value = os.environ.get("ISMIP7_SMB_ELEVATION_FEEDBACK",
                           SMB_ELEVATION_FEEDBACK_DEFAULT).strip()
    if value in ("", "0"):
        return False
    if value == "1":
        return True
    raise ValueError(
        f"ISMIP7_SMB_ELEVATION_FEEDBACK must be 1 to enable or 0/empty to "
        f"disable, got {value!r}"
    )


def calving_law():
    r"""``ISMIP7_CALVING``, lower-cased: ``none`` or the name of a law in
    :mod:`icepack_tools.calving`.

    Pure, like the rest of this module: it refuses parameters given for no
    law and, when a law is configured, the retired per-law knobs, but
    whether the name is a registered law, and whether its parameters are its
    own, only :func:`calving_law_object` can say.
    """
    value = os.environ.get("ISMIP7_CALVING", CALVING_DEFAULT).strip().lower()
    if value == "none":
        if calving_params():
            raise ValueError(
                "ISMIP7_CALVING_PARAMS is set but ISMIP7_CALVING is none: the "
                "parameters would configure no law")
        return value
    for knob, key in RETIRED_CALVING_KNOBS.items():
        if knob in os.environ:
            raise ValueError(
                f"{knob} is retired: a law's parameters live with the law "
                f"(icepack_tools.calving), so set "
                f"ISMIP7_CALVING_PARAMS={key}=<value> instead")
    return value


DRAG_GATES = ("facet", "vertex")
DRAG_GATE_DEFAULT = "vertex"
# What an inversion records as its drag gate when no cell was dragged (a mesh
# ending at the ice front, or ISMIP7_OCEAN_DRAG=0): its controls absorbed no
# gate, so a forward from it runs ISMIP7_DRAG_GATE.
DRAG_GATE_NONE = "none"


def drag_gate():
    r"""``ISMIP7_DRAG_GATE``: which water cells beside the ice the ocean drag
    skips (``front.ocean_drag_cells``). ``vertex`` (the default) skips every
    cell that touches the ice, at an edge or at a single vertex; ``facet``
    skips only the cells sharing an edge with ice, so the drag of a cell
    touching it at one vertex acts on the ice's own front node."""
    gate = os.environ.get("ISMIP7_DRAG_GATE", DRAG_GATE_DEFAULT).strip().lower()
    if gate not in DRAG_GATES:
        raise ValueError(
            f"ISMIP7_DRAG_GATE must be one of {', '.join(DRAG_GATES)}, not {gate!r}")
    return gate


# The membrane-only thickness floor [m] (dual_friction.build_rc_residual,
# h_visc_floor). An inversion runs HVISC_FLOOR_DEFAULT, 2.5 m since 6 October
# 2026 (GEOMETRY_DISCRETIZATION.md, issue #153); every MAP inverted before the
# floor was recorded ran HVISC_FLOOR_UNRECORDED.
HVISC_FLOOR_DEFAULT = "2.5"
HVISC_FLOOR_UNRECORDED = "10.0"


def hvisc_floor():
    r"""``ISMIP7_RC_HVISC_FLOOR`` [m], the inversion's membrane floor: a cell
    thinner than this carries this much ice in the membrane term alone. The
    first row of water cells beside the ice carries it too, coupling the ice
    front to the ocean drag one cell further out; 2.5 m keeps about a fifth
    of the 10 m coupling for 14 % more forward time at 1 km (issue #153)."""
    return float(os.environ.get("ISMIP7_RC_HVISC_FLOOR", HVISC_FLOOR_DEFAULT))


def _recorded(value):
    if isinstance(value, bytes):
        value = value.decode()
    return value


def forward_drag_gate(recorded, source="the MAP", restart=False):
    r"""The drag gate a forward runs: the gate its MAP records
    (``drag_gate``), which ``ISMIP7_DRAG_GATE`` may repeat and may not
    change. A MAP whose inversion dragged no cell (``none``) or that predates
    the record runs ``ISMIP7_DRAG_GATE``. A forward state from before the
    record (``restart``) ran ``facet`` unless that knob said otherwise, so
    its restart keeps ``facet`` unless the knob is set."""
    recorded = _recorded(recorded)
    if recorded is None and restart and "ISMIP7_DRAG_GATE" not in os.environ:
        return "facet"
    if recorded is None or str(recorded).strip().lower() == DRAG_GATE_NONE:
        return drag_gate()
    gate = str(recorded).strip().lower()
    if gate not in DRAG_GATES:
        raise RuntimeError(f"{source} records drag_gate={gate!r}; this code "
                           f"builds {', '.join(DRAG_GATES)}")
    if "ISMIP7_DRAG_GATE" in os.environ and drag_gate() != gate:
        raise RuntimeError(
            f"ISMIP7_DRAG_GATE={os.environ['ISMIP7_DRAG_GATE']} but {source} was "
            f"inverted under the {gate} drag gate: a forward follows its MAP")
    return gate


def forward_hvisc_floor(recorded, source="the MAP"):
    r"""The membrane floor [m] a forward runs: the floor its MAP records
    (``h_visc_floor``), which ``ISMIP7_RC_HVISC_FLOOR`` may repeat and may
    not change. A MAP that predates the record runs
    ``ISMIP7_RC_HVISC_FLOOR``, by default the 10 m every such MAP was
    inverted with."""
    env = os.environ.get("ISMIP7_RC_HVISC_FLOOR")
    recorded = _recorded(recorded)
    if recorded is None:
        return float(env if env is not None else HVISC_FLOOR_UNRECORDED)
    floor = float(recorded)
    if env is not None and float(env) != floor:
        raise RuntimeError(
            f"ISMIP7_RC_HVISC_FLOOR={env} but {source} was inverted with a "
            f"{floor:g} m membrane floor: a forward follows its MAP")
    return floor


# The forms of dual_friction.front_cliff_correction, recorded in the MAP as
# exact_front (0 off). A change to the push at the same geometry takes a new
# version.
#   1: the free-cliff push on every ice edge (60c0262, 5 Oct 2026).
#   2: the push of the face above the ice-free neighbour's bed (issue #166).
EXACT_FRONT_VERSIONS = (0, 1, 2)


def exact_front_version(value):
    r"""``ISMIP7_EXACT_FRONT`` or a MAP's ``exact_front`` record as one of
    :data:`EXACT_FRONT_VERSIONS`. A bool is the version 1 switch it was
    before version 2 existed; anything else is refused."""
    raw = _recorded(value)
    raw = raw.strip() if isinstance(raw, str) else raw
    try:
        version = int(raw)
    except (TypeError, ValueError):
        version = None
    if version is not None and not isinstance(raw, str) and version != raw:
        version = None                                   # 1.5 is no version
    if version not in EXACT_FRONT_VERSIONS:
        raise ValueError(
            f"exact_front must be one of {EXACT_FRONT_VERSIONS}, not {value!r}")
    return version


def forward_exact_front(recorded, source="the MAP"):
    r"""The cliff push a forward runs: the version its MAP records
    (``exact_front``), which ``ISMIP7_EXACT_FRONT`` may repeat and may not
    change. A MAP that predates the record was inverted without one, 0."""
    recorded = _recorded(recorded)
    version = exact_front_version(0 if recorded is None else recorded)
    env = os.environ.get("ISMIP7_EXACT_FRONT")
    if env is not None and exact_front_version(env) != version:
        raise RuntimeError(
            f"ISMIP7_EXACT_FRONT={env} but {source} was inverted with "
            f"exact_front={version}: a forward follows its MAP")
    return version


FRONT_HMIN_DEFAULT = "1.0"    # m


def front_hmin():
    r"""``ISMIP7_FRONT_HMIN`` [m]: the thickness below which a cell holds no
    ice. It draws the fixed front's t=0 extent and decides where the ocean
    drag may act (``front.ocean_drag_cells``), in the forward and the
    inversion alike."""
    return float(os.environ.get("ISMIP7_FRONT_HMIN", FRONT_HMIN_DEFAULT))


ANCHOR_LENGTH_DEFAULT = "0"    # m: local balance


def anchor_length():
    r"""``ISMIP7_ANCHOR_LENGTH`` [m]: the reach of the driving stress in the
    inversion's friction anchor (``dual_friction.weertman_anchor``). ``0``
    keeps the local balance, which vanishes with the surface slope at ice
    divides; a positive length averages the grounded driving stress over about
    that distance. The inversion records it in the MAP, and a forward takes it
    from there, never from this variable."""
    length = float(os.environ.get("ISMIP7_ANCHOR_LENGTH", ANCHOR_LENGTH_DEFAULT))
    if length < 0.0:
        raise ValueError(f"ISMIP7_ANCHOR_LENGTH must be >= 0 m, got {length}")
    return length


def lake_ice_base():
    r"""``ISMIP7_LAKE_ICE_BASE``: under BedMachine's subglacial-lake mask
    (``mask == 4``, Lake Vostok) raise the bed to the ice base, so the model
    surface ``b + H`` is BedMachine's surface rather than a bowl the depth of
    the lake's water column. On unless set to ``"0"``. The inversion records
    it in the MAP, and a forward takes it from there, so a MAP inverted without
    it keeps the geometry it was inverted on."""
    return os.environ.get("ISMIP7_LAKE_ICE_BASE", "1").strip() != "0"


def calving_params():
    r"""``ISMIP7_CALVING_PARAMS``: the law's parameters as ``key=value``
    items, from a comma-separated list (``sigma_max_fl=0.2,sigma_max_gr=1``).
    Parsing and validation are the law's (:func:`calving_law_object`)."""
    raw = os.environ.get("ISMIP7_CALVING_PARAMS", "")
    return [item.strip() for item in raw.split(",") if item.strip()]


#: Law files already registered in this process, so a second call does not
#: register the same law twice.
_CALVING_MODULES_LOADED = set()


def calving_law_object():
    r"""The configured calving law, made with its parameters, or ``None`` for
    ``none``.

    Raises for an unknown law, a misspelt or impossible parameter, a retired
    knob or a missing ``ISMIP7_CALVING_MODULE``, so a mistyped front fails at
    startup rather than after the MAP load. Imports
    :mod:`icepack_tools.calving`, and with it Firedrake, only when a law is
    configured: under the default ``none`` this stays pure.
    """
    name = calving_law()
    if name == "none":
        return None
    from icepack_tools import calving
    module = os.environ.get("ISMIP7_CALVING_MODULE")
    if module:
        path = os.path.realpath(module)
        if not os.path.isfile(path):
            raise ValueError(f"ISMIP7_CALVING_MODULE={module!r} is not a file")
        if path not in _CALVING_MODULES_LOADED:
            calving.load_module(path)
            _CALVING_MODULES_LOADED.add(path)
    params = calving.parse_params(calving_params())
    try:
        return calving.make(name, **params)
    except ValueError as err:
        raise ValueError(f"ISMIP7_CALVING={name}: {err}") from None


def fixed_front():
    r"""``ISMIP7_FIXED_FRONT``: the legacy pinned front.

    On when the variable is set to anything but the exact string ``"0"``:
    ``run_core_matrix.sh`` exports it unconditionally, so ``=0`` has to be the
    way to turn it off from there.
    """
    return os.environ.get("ISMIP7_FIXED_FRONT") not in (None, "0")


FRONT_ADVANCE_MODES = ("free", "none")


def front_advance():
    r"""``ISMIP7_FRONT_ADVANCE``: what a level-set law's front may do besides
    retreat.

    ``free`` (default): the transport advances the front wherever ice reaches
    an empty cell that the law does not remove (the extent anchor follows the
    thickness). ``none``: nothing advances past the t=0 extent; ice reaching
    a cell beyond it is removed each step and booked as calving, while the
    law still retreats the front inside it. That is the retreat-only front of
    most ISMIP6 models (Seroussi et al. 2020) and the ISMIP7 submission's
    (Andrew, 26 Sep 2026). Needs a law: without one ``ISMIP7_FIXED_FRONT``
    already holds the front.
    """
    value = os.environ.get("ISMIP7_FRONT_ADVANCE", "free").strip().lower()
    if value not in FRONT_ADVANCE_MODES:
        raise ValueError(
            f"ISMIP7_FRONT_ADVANCE must be one of {FRONT_ADVANCE_MODES}, got {value!r}")
    return value


# The year the initial geometry is dated (BedMachine v4.1's nominal year, the
# year the MAPs are inverted at), and the window of the Smith et al. (2020)
# mean dH/dt that dates a cold start before it (issue #117).
GEOMETRY_YEAR = 2015.0
DHDT_WINDOW = (2003.0, 2019.0)


def geometry_backdate_years(t_start):
    r"""Years of the Smith mean dH/dt a cold start at ``t_start`` subtracts from
    the 2015 geometry (issue #117, 25 September 2026 meeting: "2003, subtracting
    the change map").

    ``ISMIP7_GEOMETRY_BACKDATE`` names the years outright (``0`` turns it off).
    Unset, a start from 2003 up to 2015 subtracts ``2015 - t_start`` years, a
    start at or after 2015 none, and a start before 2003 is refused: the
    2003 to 2019 mean would be extrapolated past its own window.
    """
    named = os.environ.get("ISMIP7_GEOMETRY_BACKDATE", "").strip()
    if named:
        return float(named)
    years = GEOMETRY_YEAR - float(t_start)
    if years <= 0.0:
        return 0.0
    if float(t_start) < DHDT_WINDOW[0]:
        raise ValueError(
            f"a cold start at {t_start:g} would subtract {years:g} years of the "
            f"{DHDT_WINDOW[0]:g} to {DHDT_WINDOW[1]:g} Smith mean dH/dt from the "
            f"{GEOMETRY_YEAR:g} geometry, past the window it was measured over; set "
            f"ISMIP7_GEOMETRY_BACKDATE=0 to start from the {GEOMETRY_YEAR:g} "
            f"geometry as it is, or to the years to subtract")
    return years


def apparent_mb_mode():
    r"""``ISMIP7_APPARENT_MB``: the apparent-mass-balance init, or None for off.

    ``"div"`` cancels only the flux divergence, so the t=0 tendency is
    SMB minus melt (gia-style). ``"1"`` or ``"balance"`` also subtracts the
    initial forcing, so the t=0 tendency is exactly zero: a balanced control
    in the ISMIP6 ctrl_proj sense.

    ``0``, ``off``, ``none`` and the empty string mean OFF. The batch runners
    export this unconditionally and ``sbatch --export=ALL,VAR=...`` cannot
    unset a variable, so there has to be an off value; without one a run asked
    to drop the correction would silently get the full balanced one.

    The value set is closed, like ``calving_law`` and ``raster_sample``:
    anything else raises. ``no`` and ``false`` are not off spellings, and
    ``divergence`` is not ``div``, so accepting them would hand back the
    balanced control, which differs from both by the whole t=0 forcing.
    """
    value = (os.environ.get("ISMIP7_APPARENT_MB") or "").strip().lower()
    if value in ("", "0", "off", "none"):
        return None
    if value in ("1", "balance"):
        return "balance"
    if value == "div":
        return "div"
    raise ValueError(
        f"ISMIP7_APPARENT_MB must be 1 or balance (balanced control), div "
        f"(divergence only), or 0/off/none/empty to disable; got {value!r}"
    )


def auto_resume():
    r"""``ISMIP7_AUTO_RESUME``: continue unattended from this experiment's own
    newest checkpoint when no explicit restart is given.

    An integer flag, so ``=0`` turns it OFF. The batch runners export it
    unconditionally and ``sbatch --export=ALL,VAR=...`` gives no way to unset a
    variable, so ``0`` has to be the off switch; testing the string for mere
    presence would silently resume a run the user asked to start clean.
    """
    value = (os.environ.get("ISMIP7_AUTO_RESUME") or "").strip()
    if not value:
        return False
    try:
        return int(value) != 0
    except ValueError:
        raise ValueError(
            f"ISMIP7_AUTO_RESUME must be an integer flag (0 to disable), "
            f"got {value!r}"
        ) from None


def _int_flag(name, default):
    r"""An integer flag read like ``auto_resume``: unset or empty is
    ``default``, ``0`` is off, any other integer on, anything else an error."""
    value = (os.environ.get(name) or "").strip()
    if not value:
        return default
    try:
        return int(value) != 0
    except ValueError:
        raise ValueError(f"{name} must be an integer flag, got {value!r}") from None


def mesh_build_check():
    r"""``ISMIP7_MESH_BUILD_CHECK``: refuse a MAP or restart whose mesh has the
    name of the ``ISMIP7_MESH`` file and a different triangulation, as two
    sites' builds of the production mesh do. On unless ``0``."""
    return _int_flag("ISMIP7_MESH_BUILD_CHECK", True)


# ── Roots a second checkout does not carry ──────────────────────────────
# The code moves with the invocation; the large gitignored artifacts do not.
# site_env.sh names them for the shell half of a run; these are the Python
# half, so a driver started by hand resolves the same paths.
_ANTARCTICA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "antarctica"
)


def obs_data_root():
    r"""``ISMIP7_OBS_DATA_ROOT``: BedMachine, MEaSUREs, RACMO and the dH/dt
    cache rasters. Falls back to this checkout's ``antarctica/data``, which is
    right wherever the data sits beside the code and empty where it does not."""
    return os.environ.get("ISMIP7_OBS_DATA_ROOT",
                          os.path.join(_ANTARCTICA, "data"))


# ── The melt calibration every run reads ─────────────────────────────────
# Issue 26, decided on 25 September 2026: K50 of IU's rule-based toolbox
# selection (scripts/select_melt_parameters.py, run record
# calibration-melt-toolbox-1km-rule), one K with a thermal-forcing offset per
# IMBIE basin, fitted through the forward's own DG0 melt path on the
# 1000 m / 10 km production mesh. The file is tracked, so a fresh clone melts
# with it and nothing is copied or rerun. Its sidecar (<name>.source.json)
# records the file's sha256 and the conventions the fit holds under, which
# the forward checks (forcing.check_melt_contract). The offsets were fitted
# again at that K under the front-cell rule (issue #167, run record
# calibration-melt-1km-vertex-front), whose melt-receiving area is 2.4 %
# smaller on the 1 km mesh; that file is the default, the calibration of the
# default raster sampling.
MELT_CALIBRATION_VERTEX = os.path.join(
    _ANTARCTICA, "calibration", "deltaT_per_basin_1000_K6.500e-05.npz")
MELT_CALIBRATION_DEFAULT = os.path.join(
    _ANTARCTICA, "calibration", "deltaT_per_basin_1000_K6.500e-05_vertex_front.npz")
# The tracked calibration of each raster sampling. A calibration sums melt
# over the cells its sampling builds, so a run takes the one fitted under its
# own (forcing.check_melt_contract refuses another); a sampling with no entry
# needs ISMIP7_DELTAT_PER_BASIN_NPZ.
MELT_CALIBRATIONS = {"vertex": MELT_CALIBRATION_VERTEX,
                     "vertex_front": MELT_CALIBRATION_DEFAULT}

# Knobs that no longer shape a run. They are refused rather than ignored, so
# a job script written before the change fails at startup.
REMOVED_MELT_KNOBS = {
    "ISMIP7_K_MELT": "the scalar K fallback went with the tracked melt "
                     "calibration; a run takes K from its calibration file",
}


def file_sha256(path):
    r"""sha256 of a file, as hex."""
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def melt_calibration_sidecar(npz_path):
    r"""The ``.source.json`` beside a melt calibration npz."""
    return os.path.splitext(npz_path)[0] + ".source.json"


def melt_calibration_contract(npz_path):
    r"""The sidecar of ``npz_path`` as a dict, or None where there is none.

    It records what the npz does not: the file's sha256, the mesh it was
    fitted on, the raster sampling of its geometry and the observation table.
    calibrate_deltaT.py and select_melt_parameters.py write one beside every
    offsets file they fit (`write_melt_calibration_sidecar`), and the tracked
    default carries one. A file named with ISMIP7_DELTAT_PER_BASIN_NPZ that
    has none, written before the fits wrote sidecars, is read as given."""
    sidecar = melt_calibration_sidecar(npz_path)
    if not os.path.exists(sidecar):
        return None
    with open(sidecar) as f:
        return json.load(f)


# What every melt calibration sidecar records: the keys the forward reads
# (sha256, mesh, vertices, raster_sample) and the ones a reader needs to tell
# one fit from another. The fits write more; the tracked default also carries
# the decision fields a person adds when a file is promoted (decided,
# decision, selected_as, selection, run_record, mesh_build).
MELT_CALIBRATION_REQUIRED = (
    "file", "sha256", "K", "mesh", "vertices", "raster_sample", "melt_slope",
    "geometry_space", "obs_table", "job", "code",
)


def refuse_tracked_calibration_out(directory):
    r"""Refuse to write a fit into the tracked calibration's directory.

    A fit writes its offsets and their sidecar together, so a fit written
    there would replace the tracked file and its hand-kept record with a
    consistent pair that has lost the decision fields, and every run would
    read it as the default without a word. Promoting a fit takes a reviewed
    change: write it elsewhere, copy both files in, add the decision fields
    and a run record."""
    tracked = os.path.dirname(MELT_CALIBRATION_DEFAULT)
    if os.path.realpath(directory) == os.path.realpath(tracked):
        raise ValueError(
            f"{directory} holds the tracked melt calibration. Write the fit to "
            f"another directory; promoting it means copying the npz and its "
            f".source.json here, adding the decision fields and writing a run "
            f"record.")


def _json_value(value):
    r"""A numpy scalar or array as the plain value JSON writes (``np.int64``
    and ``np.bool_`` are not JSON serialisable)."""
    if hasattr(value, "tolist"):
        return value.tolist()
    raise TypeError(f"{type(value).__name__} is not JSON serialisable")


def _strings(value):
    r"""Every string inside a record, however deeply nested."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for v in value.values():
            yield from _strings(v)
    elif isinstance(value, (list, tuple)):
        for v in value:
            yield from _strings(v)


def _names_a_path(text):
    r"""Whether a word of ``text`` is an absolute or home-directory path. A
    lone slash is prose ("1000 m / 10 km"), which the tracked sidecar has."""
    for word in text.split():
        word = word.strip("()[]{}<>,;:'\"")
        if len(word) > 1 and (os.path.isabs(word) or word.startswith("~")):
            return True
    return False


def check_melt_calibration_record(record):
    r"""Refuse a sidecar record, and return it as the JSON a sidecar holds.

    Refused: a `MELT_CALIBRATION_REQUIRED` key missing (``file`` and
    ``sha256`` excepted, which the writer takes from the npz) or empty
    (``job`` excepted, since a fit run outside Slurm has none); an absolute
    path anywhere, since a promoted sidecar is tracked and a path would carry
    a home directory into git; a NaN or an infinity, which JSON cannot hold.
    Pure, so every rank of a parallel fit checks the record it is about to
    write and a bad one stops them all together."""
    problems = []
    absent = [k for k in MELT_CALIBRATION_REQUIRED
              if k not in ("file", "sha256") and k not in record]
    if absent:
        problems.append(f"lacks {', '.join(absent)}")
    empty = [k for k in MELT_CALIBRATION_REQUIRED
             if k != "job" and k in record and record[k] is None]
    if empty:
        problems.append(f"leaves {', '.join(empty)} empty")
    paths = [s for s in _strings(record) if _names_a_path(s)]
    if paths:
        problems.append(f"names a path ({paths[0]}); name files by basename "
                        f"and sha256")
    if problems:
        raise ValueError("a melt calibration sidecar record "
                         + "; ".join(problems))
    return json.dumps(record, indent=2, allow_nan=False, default=_json_value)


def write_melt_calibration_sidecar(npz_path, record):
    r"""Write the sidecar of the offsets file ``npz_path`` from ``record``,
    and return its path.

    ``file`` and ``sha256`` are taken from the npz as written, the bytes the
    forward hashes. The record is checked first
    (`check_melt_calibration_record`), and the sidecar is replaced in one
    step, so a job that dies mid-write leaves the old sidecar, whose hash then
    refuses a changed npz."""
    refuse_tracked_calibration_out(os.path.dirname(os.path.abspath(npz_path)))
    contract = {"file": os.path.basename(npz_path),
                "sha256": file_sha256(npz_path)}
    contract.update((k, v) for k, v in record.items() if k not in contract)
    text = check_melt_calibration_record(contract)
    sidecar = melt_calibration_sidecar(npz_path)
    with open(sidecar + ".tmp", "w") as f:
        f.write(text + "\n")
    os.replace(sidecar + ".tmp", sidecar)
    return sidecar


def deltat_per_basin_npz(raster_sample_of_run=None):
    r"""The per-basin thermal-forcing offsets a run melts with, or None on the
    legacy per-basin K path.

    The ISMIP7 ocean-forcing recommendation calibrates ONE dimensionless K
    from the 4-term toolbox and then a thermal-forcing offset deltaT_b per
    IMBIE basin at that K (``optimise_deltaT``); a per-basin K is not part of
    it. The ocean callbacks add the offset to TF before the melt law and melt
    with the file's K everywhere.

    ``ISMIP7_DELTAT_PER_BASIN_NPZ`` names the file; unset, it is the tracked
    calibration of the run's raster sampling (``MELT_CALIBRATIONS``;
    ``raster_sample_of_run``, else ``ISMIP7_RASTER_SAMPLE``), which for the
    default sampling is ``MELT_CALIBRATION_DEFAULT``. ``ISMIP7_K_PER_BASIN_NPZ`` selects the legacy
    per-basin K path instead (``k_per_basin_npz``), and then this returns
    None. Naming both is refused.

    Checked where it is read, so a driver that calls this before its model
    setup fails before the MAP is loaded: the file must exist, it must match
    the sha256 its sidecar records, and ``ISMIP7_K_SCALE`` must be 1, because
    the offsets were fitted at the file's K and a scaled K invalidates them."""
    for knob, why in REMOVED_MELT_KNOBS.items():
        if os.environ.get(knob):
            raise ValueError(f"{knob} is no longer read: {why}. Unset it.")
    path = os.environ.get("ISMIP7_DELTAT_PER_BASIN_NPZ") or None
    legacy = os.environ.get("ISMIP7_K_PER_BASIN_NPZ") or None
    if path is not None and legacy is not None:
        raise ValueError(
            f"ISMIP7_DELTAT_PER_BASIN_NPZ={path} and "
            f"ISMIP7_K_PER_BASIN_NPZ={legacy} are both set. A run melts with "
            f"one calibration: unset one of them.")
    if legacy is not None:
        return None
    if path is None:
        sampling = (raster_sample() if raster_sample_of_run is None
                    else str(raster_sample_of_run).lower())
        if sampling not in MELT_CALIBRATIONS:
            raise FileNotFoundError(
                f"No tracked melt calibration was fitted on geometry sampled "
                f"with {sampling}. Fit one (calibrate_deltaT.py under "
                f"ISMIP7_RASTER_SAMPLE={sampling}) and name it with "
                f"ISMIP7_DELTAT_PER_BASIN_NPZ.")
        path = MELT_CALIBRATIONS[sampling]
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"The tracked melt calibration {path} is missing from this "
                f"checkout. Restore it with git checkout, or name another "
                f"file with ISMIP7_DELTAT_PER_BASIN_NPZ.")
    elif not os.path.exists(path):
        raise FileNotFoundError(
            f"ISMIP7_DELTAT_PER_BASIN_NPZ={path} does not exist. Write it with "
            f"antarctica/scripts/calibrate_deltaT.py, or unset the knob for "
            f"the tracked calibration.")
    contract = melt_calibration_contract(path)
    if contract is not None and contract.get("sha256") != file_sha256(path):
        raise ValueError(
            f"{path} does not match the sha256 its sidecar "
            f"{melt_calibration_sidecar(path)} records. The file changed "
            f"without its record: restore it with git checkout, or rewrite "
            f"the sidecar with the calibration that produced it.")
    k_scale = float(os.environ.get("ISMIP7_K_SCALE", "1.0"))
    if k_scale != 1.0:
        raise ValueError(
            f"The melt calibration {path} and ISMIP7_K_SCALE={k_scale:g} "
            f"are both in effect. The offsets were fitted at the file's K, so "
            f"a scaled K invalidates them: unset ISMIP7_K_SCALE, or refit with "
            f"calibrate_deltaT.py --K at the K you want and name the file "
            f"with ISMIP7_DELTAT_PER_BASIN_NPZ.")
    return path


def k_per_basin_npz():
    r"""``ISMIP7_K_PER_BASIN_NPZ``: a per-basin K file from calibrate_melt.py,
    the calibration the runs used before issue 26 was decided, or None.

    It is read only when named: no directory is searched, so a run on a
    machine that happens to hold an old file still melts with the tracked
    calibration. A named file that does not exist is refused."""
    path = os.environ.get("ISMIP7_K_PER_BASIN_NPZ") or None
    if path is not None and not os.path.exists(path):
        raise FileNotFoundError(
            f"ISMIP7_K_PER_BASIN_NPZ={path} does not exist. Unset it for the "
            f"tracked melt calibration.")
    return path
