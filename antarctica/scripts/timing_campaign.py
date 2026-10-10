#!/usr/bin/env python3
"""Shared policy and validation for the cached transient timing campaign.

This module deliberately has no Firedrake dependency.  The Makefile helpers,
Slurm wrappers, report builder, and unit tests all import the same lane and
record policy so the submission matrix cannot drift away from the report.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import socket
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.fspath(Path(__file__).resolve().parents[2]))

from icepack2_tools.runconfig import (
    BUDD_SHELF_GATE, TARGET_MESH_GEOMETRY_METHOD, U_LIM_DEFAULT,
)


RECORD_SCHEMA_VERSION = 3
CACHE_SCHEMA_VERSION = 3
# Physics generation of the caches and lanes, spelled into CACHE_TAG and every
# campaign tag. v3: the 09-14 caches (explicit field inventory, shared residual
# stabilizers). v4 (2026-09-16): the Budd shelf gate tests height above
# flotation instead of the sign of roundoff N (runconfig.BUDD_SHELF_GATE,
# ported from hoffmaao/antarctica e602705); nothing solved under the old gate
# is reused. Bump it whenever the physics behind a prepared state changes.
CAMPAIGN_VERSION = 4
CACHE_ROLE = "timing-initial-state"
CACHE_REQUIRED_FIELDS = (
    "log_friction",
    "log_fluidity",
    "velocity_obs",
    "thickness",
    "bed",
    "surface",
    "fluidity_prior",
    "velocity",
    "membrane_stress",
    "basal_stress",
    "H_init",
    "phi_eff",
    "C_w0",
    "N_ref",
)


def cache_required_fields(friction="budd"):
    """The fields a published initial-state cache has to carry. N_ref is
    Budd's frozen reference effective pressure (simulation.setup_model builds
    it under that law alone), so a regularized Coulomb state never has it
    (found on Quartz 10569348, whose publish refused a complete RC cache)."""
    if friction == "budd":
        return CACHE_REQUIRED_FIELDS
    return tuple(name for name in CACHE_REQUIRED_FIELDS if name != "N_ref")


LCS = (500, 1000, 2000, 2500, 5000)
RATIOS = (10, 20)
DISPLAY_CORES = (16, 32, 64)
CORES_BY_LC = {
    500: (32, 64),
    1000: (16, 32, 64),
    2000: (16, 32),
    2500: (16, 32),
    5000: (16,),
}

# Two solvers, deliberately separate. Every initial-state cache is prepared
# (and every per-mesh invert publishes) under CACHE_SOLVER_MODE; it is spelled
# into CACHE_TAG and so into the cache filenames. The solver a lane TIMES is a
# campaign parameter like the interval: it leads the campaign tag, so records
# under different solvers never mix, and a lane under any of them starts from
# the same prepared state -- the linear solver is the only thing that differs
# between two campaigns' lanes.
CACHE_SOLVER_MODE = "scpc_mumps"
LANE_SOLVER_MODES = ("scpc_mumps", "scpc_gamg")
LANE_SOLVER_DEFAULT = CACHE_SOLVER_MODE
SOURCE_TAG = "dg0_logvelnet"
# Campaign source MAP. Promoted on 2026-09-17 from the imported old-gate MAP
# (inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5, 64000/2500 mesh) to
# its 250-iteration re-inversion on the 2500/25000 campaign mesh under the
# HAF-gated Budd law (quartz job 10489636: L-BFGS converged at 219 iterations,
# velocity chi2 2.29e4 -> 695, published ||F|| 2.46e3, friction_gate=haf).
# The old-gate controls run away at step 3 on every fixed-law lane (13.3 % of
# floating cells had sat at the friction cap), so lanes must descend from a
# fixed-law MAP; one source for all meshes keeps the 500 m meshes, which are
# never re-inverted, consistent with the rest. By construction this equals
# mesh_inversion_basename(2500, 25000): that mesh's invert-published cache IS
# the source-published cache, and the manager never re-inverts or re-prepares
# the source mesh. Must agree with the Makefile's TIMING_INVERSION.
SOURCE_INVERSION_BASENAME = (
    "inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5"
)
# Matrix interval. Step count and the 2.5 km timestep are PARAMETERS (the dt
# ladder of 2026-09-15): both are written into the campaign tag, so records
# from different rungs can never mix, and validate_timing_record derives the
# expected interval from a record's own tag rather than from these constants.
MATRIX_T_START = 2015.0
# 2026-09-16 dt ladder on 2500/25000 (strict contract, pristine prepare state):
# 5 x 0.25 ran away at step 2, 10 x 0.125 and 20 x 0.0625 completed 1.25 yr
# with flat speed and thickness; 10 x 0.125 is the cheapest rung that passes.
# Must match the Makefile's MATRIX_STEPS / MATRIX_DT_2500 defaults.
MATRIX_STEPS_DEFAULT = 10
MATRIX_DT_2500_DEFAULT = 0.125
# The archived legacy campaign (tag scpc_mumps_5step_dt0p25at2500, no
# parseable interval) was always 5 x 0.25; judge its records by that.
LEGACY_MATRIX_STEPS = 5
LEGACY_MATRIX_DT_2500 = 0.25
MATRIX_REFERENCE_LC = 2500.0
# The timestep scales with LC only up to an absolute cap. 2026-09-17: both
# 5000 m scouts (dt 0.25 by the LC/2500 rule) ran away exactly like the
# 5 x 0.25 probe on 2500/25000 -- a cell within ~15 m of flotation flips
# grounded <-> floating every step with growing dh -- while 2500/50000
# (dt 0.125) and 2000/40000 (dt 0.1) passed from the same source MAP. The
# lagged thickness/velocity coupling through the HAF friction gate sets a
# stable dt that does not grow with the mesh size. A constant, not a campaign
# parameter (it is not in the tag); campaigns before v4 were never capped and
# are still judged by the plain LC/2500 rule. A rung whose own 2.5 km step
# exceeds the cap runs what its tag says. Must match the Makefile's
# MATRIX_DT_MAX.
MATRIX_DT_MAX = 0.125
MATRIX_DT_MAX_SINCE_VERSION = 4


def matrix_steps(override=None):
    if override is not None:
        return int(override)
    return int(os.environ.get("ISMIP7_MATRIX_STEPS", str(MATRIX_STEPS_DEFAULT)))


def matrix_dt_2500(override=None):
    if override is not None:
        return float(override)
    return float(
        os.environ.get("ISMIP7_MATRIX_DT_2500", str(MATRIX_DT_2500_DEFAULT))
    )


def dt_tag(dt_2500):
    """``0.25 -> '0p25'``: the Makefile's ``$(subst .,p,$(MATRIX_DT_2500))``."""
    return f"{float(dt_2500):g}".replace(".", "p")


def lane_solver(override=None):
    """The diagnostic solver a campaign's lanes time (not the caches')."""
    name = (
        override
        if override is not None
        else os.environ.get("ISMIP7_TIMING_SOLVER", LANE_SOLVER_DEFAULT)
    )
    if name not in LANE_SOLVER_MODES:
        raise ValueError(
            f"timing solver must be one of {LANE_SOLVER_MODES}, not {name!r}"
        )
    return name


# Lane physics contracts. ``strict`` is the original v3 contract (no apparent
# mass balance, no calving sink at the 2015 front). The others add the
# production closure (``ISMIP7_APPARENT_MB=div``, uncapped) and/or the fixed
# calving front (``ISMIP7_FIXED_FRONT``); a non-strict contract suffixes the
# campaign tag so its records never masquerade as strict ones. Rescue stays
# off and subcycles stay [1] under every contract: the benchmark times the
# bare solver.
CONTRACTS = {
    "strict": {
        "apparent_mb_mode": None,
        "apparent_mb_cap_m_per_yr": 0.0,
        "fixed_front": False,
    },
    "divfront": {
        "apparent_mb_mode": "div",
        "apparent_mb_cap_m_per_yr": 0.0,
        "fixed_front": True,
    },
    "front": {
        "apparent_mb_mode": None,
        "apparent_mb_cap_m_per_yr": 0.0,
        "fixed_front": True,
    },
    "div": {
        "apparent_mb_mode": "div",
        "apparent_mb_cap_m_per_yr": 0.0,
        "fixed_front": False,
    },
}
CONTRACT_DEFAULT = "strict"


def contract_name(override=None):
    name = (
        override
        if override is not None
        else os.environ.get("ISMIP7_TIMING_CONTRACT", CONTRACT_DEFAULT)
    )
    if name not in CONTRACTS:
        raise ValueError(
            f"timing contract must be one of {tuple(CONTRACTS)}, not {name!r}"
        )
    return name


def contract_exports(name=None):
    """Environment a lane needs to run under ``name``'s contract."""
    contract = CONTRACTS[contract_name(name)]
    exports = {"ISMIP7_AMB_CAP": f"{contract['apparent_mb_cap_m_per_yr']:g}"}
    if contract["apparent_mb_mode"]:
        exports["ISMIP7_APPARENT_MB"] = contract["apparent_mb_mode"]
    if contract["fixed_front"]:
        exports["ISMIP7_FIXED_FRONT"] = "1"
    return exports


# Runaway tripwire (simulation.run_simulation): fail a lane at the first step
# whose speed or per-step thickness growth exceeds these, naming the cell.
TRIPWIRE_DEFAULTS = {
    "ISMIP7_TRIPWIRE_U_MAX": U_LIM_DEFAULT,
    # Absolute thickness cap: no Antarctic cell is 5 km thick (the 2500 m
    # mesh's maximum is 4.25 km); a runaway pile-up reaches 1e4-3e4 m.
    "ISMIP7_TRIPWIRE_H_MAX": "5000",
    # Relative thickening rate (dh/h)/dt [1/yr], tested only on cells that
    # entered the step at least HMIN thick. A rate, not a per-step fraction,
    # so every rung of a dt ladder is judged by the same physics: the lane
    # that completed 1.25 yr at dt = 0.0625 peaked at 6.8/yr (a 107 m
    # buffer cell being fed), the dt = 0.125 rung scored 4.1/yr in the same
    # Amundsen cell that a 0.5 per-step bound had tripped on, and the
    # dt = 0.25 pile-up thickened 500 -> 3000 m within 0.25 yr (>= 20/yr at
    # onset) with speed_max already past 2e4. Thin cells (a 0.3 m Beardmore
    # buffer cell filling at 70 m/yr) are reported, never tripped.
    "ISMIP7_TRIPWIRE_DH_RATE": "20",
    "ISMIP7_TRIPWIRE_HMIN": "100.0",
}


def campaign_tag(steps=None, dt_2500=None, contract=None, solver=None):
    name = contract_name(contract)
    tag = (
        f"{lane_solver(solver)}_{matrix_steps(steps)}step_"
        f"dt{dt_tag(matrix_dt_2500(dt_2500))}at2500"
        f"_dg0_logvelnet_cached_strict_v{CAMPAIGN_VERSION}"
    )
    if name != "strict":
        tag += f"_{name}"
    return tag


# Lanes from the transferred prepare state are THE campaign (TIMING_INITIAL_STATE
# =prepare, the default); lanes from a per-mesh invert-published cache carry
# their own suffix so the two initial states never mix in one report.
LANE_INITIAL_STATES = ("prepare", "invert")
LANE_INITIAL_STATE_DEFAULT = "prepare"


def lane_tag(initial_state=LANE_INITIAL_STATE_DEFAULT, campaign=None):
    if initial_state not in LANE_INITIAL_STATES:
        raise ValueError(
            f"initial state must be one of {LANE_INITIAL_STATES}, not {initial_state!r}"
        )
    campaign = campaign_tag() if campaign is None else campaign
    return f"{campaign}_reinverted" if initial_state == "invert" else campaign


def probe_tag(steps=None, dt_2500=None, contract=None, solver=None):
    """Probe lanes (``cache_probe`` kind) run one contract on one mesh without
    touching campaign records; they are the experiment ladder."""
    return f"{campaign_tag(steps, dt_2500, contract, solver)}_probe"


_NONSTRICT_CONTRACTS = "|".join(
    sorted((name for name in CONTRACTS if name != "strict"), key=len, reverse=True)
)
_CAMPAIGN_TAG_RE = re.compile(
    rf"^(?P<solver>{'|'.join(LANE_SOLVER_MODES)})"
    r"_(?P<steps>\d+)step_dt(?P<dt>\d+(?:p\d+)?)at2500"
    r"_dg0_logvelnet_cached_strict_v(?P<version>\d+)"
    rf"(?:_(?P<contract>{_NONSTRICT_CONTRACTS}))?"
    r"(?:_(?P<lane>reinverted|probe))?$"
)


def parse_campaign_tag(tag):
    """``{solver, steps, dt_2500, contract, lane, version}`` encoded in a lane tag.

    Older versions still parse (``make matrix TIMING_TAG=<v3 tag>`` renders
    the archive); only the current CAMPAIGN_VERSION is ever submitted."""
    match = _CAMPAIGN_TAG_RE.match(str(tag))
    if match is None:
        raise ValueError(f"tag {tag!r} does not name a cached-strict lane")
    return {
        "solver": match["solver"],
        "steps": int(match["steps"]),
        "dt_2500": float(match["dt"].replace("p", ".")),
        "contract": match["contract"] or "strict",
        "lane": match["lane"],
        "version": int(match["version"]),
    }


# Import-time defaults (env may still override via the helpers above).
CAMPAIGN_TAG = campaign_tag()
REINVERTED_TAG = lane_tag("invert")
PROBE_TAG = probe_tag()
CACHE_TAG = f"{CACHE_SOLVER_MODE}_dg0_logvelnet_v{CAMPAIGN_VERSION}"
# Per-mesh short invert length. Override with ISMIP7_TIMING_INVERSION_MAXITER
# or `make timing-inversion TIMING_INVERSION_MAXITER=5` for a debug pass.
INVERSION_MAXITER_DEFAULT = 250


def inversion_maxiter(override=None):
    if override is not None:
        return int(override)
    return int(
        os.environ.get(
            "ISMIP7_TIMING_INVERSION_MAXITER",
            str(INVERSION_MAXITER_DEFAULT),
        )
    )


def inversion_tag(maxiter=None):
    return f"{CACHE_TAG}_{inversion_maxiter(maxiter)}iter"


# Import-time defaults (env may still override via the helpers above).
INVERSION_MAXITER = inversion_maxiter()
INVERSION_TAG = inversion_tag(INVERSION_MAXITER)

# Import-time snapshots for the legacy (cold-start, 30-lane) matrix code.
MATRIX_STEPS = matrix_steps()
MATRIX_DT_2500 = matrix_dt_2500()
BUFFER_M = 20000
MASS_RESIDUAL_TOL_GT = 5.0e-5

MEMORY_BY_LC = {
    500: "240G",
    1000: "96G",
    2000: "80G",
    2500: "64G",
    5000: "32G",
}


def inversion_cores(lc):
    """Ranks for the per-mesh short logvelnet inversion."""
    return 32 if int(lc) < 2500 else 16


# The request is sized for the per-mesh invert's full_mumps solve and its
# tlm_adjoint tape, which leaves room under scpc_gamg (the default since
# 3 October, about half the memory), so its footprint differs from the
# scpc_mumps transient's (MEMORY_BY_LC). Measured on Quartz: 2500/25000 on 16 ranks peaks
# near 2 GB/rank; both 500 m meshes exceeded 13 GB/rank on 32 ranks and were
# OOM-killed at 240G before their first evaluation, so 500 m is not inverted.
INVERSION_MEMORY_BY_LC = {
    1000: "192G",
    2000: "160G",
    2500: "64G",
    5000: "32G",
}
INVERSION_SKIP_LCS = (500,)


def inversion_required(lc):
    """Whether the campaign re-inverts on this mesh.

    Meshes in INVERSION_SKIP_LCS start their lanes from the prepared cache,
    i.e. the transferred 2.5 km MAP; TIMING_MATRIX.md records which initial
    state each mesh used so the two are never mixed silently.
    """
    return int(lc) not in INVERSION_SKIP_LCS


def inversion_memory(lc):
    if not inversion_required(lc):
        raise ValueError(f"no per-mesh inversion is run at LC={int(lc)}")
    return INVERSION_MEMORY_BY_LC[int(lc)]


def mesh_inversion_basename(lc, lc_coarse, maxiter=None):
    n = inversion_maxiter(maxiter)
    return (
        f"inversion_icepack2_budd_n3_dg0_logvelnet_"
        f"{int(lc)}_{int(lc_coarse)}_{n}iter.h5"
    )


_MESH_INVERSION_RE = re.compile(
    r"^inversion_icepack2_budd_n3_dg0_logvelnet_"
    r"(?P<lc>\d+)_(?P<lc_coarse>\d+)_(?P<maxiter>\d+)iter\.h5$"
)


def mesh_inversion_source_mesh(basename):
    """``(lc, lc_coarse)`` when ``basename`` is a per-mesh invert MAP, else None.

    The campaign source is such a MAP since 2026-09-17; the manager uses this
    to recognise the source mesh at any --maxiter (its cache path does not
    depend on maxiter, so any invert there would republish the source cache).
    """
    match = _MESH_INVERSION_RE.match(os.path.basename(os.fspath(basename)))
    if not match:
        return None
    return int(match["lc"]), int(match["lc_coarse"])


def mesh_inversion_map_path(root, lc, lc_coarse, maxiter=None):
    return (
        Path(root)
        / "results"
        / "timing"
        / "inversion"
        / mesh_inversion_basename(lc, lc_coarse, maxiter=maxiter)
    )


def mesh_inversion_timing_json_path(root, lc, lc_coarse, ncores=None, maxiter=None):
    if ncores is None:
        ncores = inversion_cores(lc)
    return (
        Path(root)
        / "results"
        / "timing"
        / (
            f"inversion_timing_{inversion_tag(maxiter)}"
            f"_{int(lc)}_{int(lc_coarse)}_{int(ncores)}.json"
        )
    )


def mesh_inversion_status_path(root, lc, lc_coarse, maxiter=None):
    return (
        Path(root)
        / "results"
        / "timing"
        / (
            f"status_inversion_{inversion_tag(maxiter)}"
            f"_{int(lc)}_{int(lc_coarse)}.txt"
        )
    )


def mesh_rows(lcs=LCS, ratios=RATIOS):
    """Return unique ``(lc, lc_coarse)`` rows in display order."""
    return tuple(sorted(
        {(int(lc), int(lc) * int(ratio))
         for lc in lcs for ratio in ratios},
        key=lambda item: (item[0], item[1]),
    ))


def planned_lanes():
    """Return the selected 20 ``(lc, lc_coarse, ncores)`` lanes."""
    return tuple(
        (lc, lc_coarse, ncores)
        for lc, lc_coarse in mesh_rows()
        for ncores in CORES_BY_LC[lc]
    )


def scout_lanes():
    """Return one lowest-retained-core scout for every target mesh."""
    return tuple(
        (lc, lc_coarse, CORES_BY_LC[lc][0])
        for lc, lc_coarse in mesh_rows()
    )


def scaling_lanes():
    scouts = set(scout_lanes())
    return tuple(lane for lane in planned_lanes() if lane not in scouts)


def expected_dt(lc, dt_2500=None, version=None):
    """A lane's timestep: the 2.5 km step scaled by LC/2500, capped at
    MATRIX_DT_MAX for campaigns that have the cap (``version`` defaults to
    the current one)."""
    base = matrix_dt_2500(dt_2500)
    dt = base * float(lc) / MATRIX_REFERENCE_LC
    version = CAMPAIGN_VERSION if version is None else int(version)
    if version < MATRIX_DT_MAX_SINCE_VERSION:
        return dt
    return min(dt, max(MATRIX_DT_MAX, base))


def expected_t_end(lc, steps=None, dt_2500=None, version=None):
    return MATRIX_T_START + matrix_steps(steps) * expected_dt(lc, dt_2500, version)


def mesh_basename(lc, lc_coarse, buffer_m=BUFFER_M):
    return f"antarctica_{int(lc_coarse)}_{int(lc)}_buffered{int(buffer_m)}.msh"


def cache_stem(lc, lc_coarse, buffer_m=BUFFER_M):
    return (
        f"initial_state_{CACHE_TAG}_{int(lc)}_{int(lc_coarse)}"
        f"_buffered{int(buffer_m)}"
    )


def cache_paths(cache_dir, lc, lc_coarse, buffer_m=BUFFER_M):
    stem = cache_stem(lc, lc_coarse, buffer_m)
    root = Path(cache_dir)
    return root / f"{stem}.h5", root / f"{stem}.json"


def pristine_cache_paths(cache_dir, lc, lc_coarse, buffer_m=BUFFER_M):
    """``(checkpoint, manifest)`` of the untouched copy of the prepared cache.

    The per-mesh invert warm-starts from it and then republishes the cache
    path above from its own MAP, so a re-run must never warm-start from that
    output; the copied manifest keeps naming the imported source MAP so the
    warm start can still be validated after the published one has moved on.
    """
    stem = f"{cache_stem(lc, lc_coarse, buffer_m)}.prepare"
    root = Path(cache_dir)
    return root / f"{stem}.h5", root / f"{stem}.json"


def pristine_sibling(cache_path):
    """The pristine copy's path for a published cache path (same stem)."""
    text = os.fspath(cache_path)
    if text.endswith(".prepare.h5"):
        return text
    return text[:-3] + ".prepare.h5" if text.endswith(".h5") else text


# --- MAP checks (make map-check) ------------------------------------------
# One released MAP under one law, on its own mesh and transferred onto the
# production mesh (antarctica/MAP_CHECK.md). Its caches, tags and records are
# named apart from the campaign's: the cache role differs, the cache stem
# never starts with ``initial_state_``, and the tag never parses as a
# campaign tag, so neither manager can claim the other's artifacts.
MAP_CHECK_CACHE_ROLE = "map-check-initial-state"
MAP_CHECK_KIND = "map_check"
MAP_CHECK_STEPS = 10
MAP_CHECK_LAW_TAGS = {"budd": "budd", "regularized_coulomb": "rc"}


def map_check_law_tag(friction):
    try:
        return MAP_CHECK_LAW_TAGS[friction]
    except KeyError:
        raise ValueError(
            f"a MAP check runs under one of {tuple(MAP_CHECK_LAW_TAGS)}, "
            f"not {friction!r}"
        ) from None


def map_check_stem(map_path):
    """The MAP's basename without its extension: the per-MAP directory name."""
    return Path(os.fspath(map_path)).stem


def map_check_tag(stem, friction, solver):
    return (
        f"mapcheck_{map_check_law_tag(friction)}_{lane_solver(solver)}"
        f"_{MAP_CHECK_STEPS}step_{stem}"
    )


def map_check_cache_stem(stem, friction, lc, lc_coarse, buffer_m):
    return (
        f"map_check_state_{map_check_law_tag(friction)}_{stem}"
        f"_{int(lc)}_{int(lc_coarse)}_buffered{int(buffer_m)}"
    )


def map_check_run_tag(stem, friction, mesh_label):
    """``ISMIP7_RUN_TAG`` of a control: short, since the driver's experiment
    name and its checkpoint names carry it."""
    match = re.search(r"snap\d+_\d+", stem)
    snap = match.group(0) if match else stem[:24]
    return f"mapcheck_{map_check_law_tag(friction)}_{snap}_{mesh_label}"


def timing_record_basename(tag, lc, lc_coarse, ncores):
    return f"timing_{tag}_{int(lc)}_{int(lc_coarse)}_{int(ncores)}.json"


def timing_status_basename(tag, lc, lc_coarse, ncores):
    return f"status_{tag}_{int(lc)}_{int(lc_coarse)}_{int(ncores)}.txt"


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


# MUMPS's analysis and ordering (ICNTL 28 and 29) follow the PETSc build and
# ISMIP7_MUMPS_ANALYSIS. They change how the factorization is computed and
# leave the converged state alone, so a MUMPS-factored block is fingerprinted
# with the values every cache was prepared under.
_MUMPS_BLOCK_PREFIXES = ("condensed_field_", "fieldsplit_1_")
_MUMPS_ANALYSIS_FINGERPRINT = {"mat_mumps_icntl_28": 2, "mat_mumps_icntl_29": 1}


def solver_configuration_fingerprint(configuration):
    """Fingerprint setup-relevant solver settings, excluding recovery policy.

    Cache construction deliberately enables adaptive continuation while strict
    transient lanes disable rescue/subcycling.  Those recovery settings must
    differ, so the cache identity covers the nonlinear/linear operators and
    tolerances that define the prepared state, not the transient recovery
    policy, nor the MUMPS analysis that factored them.
    """
    petsc_options = configuration.get("diagnostic_petsc_options")
    if petsc_options is not None:
        petsc_options = dict(petsc_options)
        for prefix in _MUMPS_BLOCK_PREFIXES:
            if petsc_options.get(f"{prefix}pc_factor_mat_solver_type") == "mumps":
                petsc_options.update({
                    f"{prefix}{key}": value
                    for key, value in _MUMPS_ANALYSIS_FINGERPRINT.items()
                })
    selected = {
        "diagnostic_mode": configuration.get("diagnostic_mode"),
        "diagnostic_petsc_options": petsc_options,
        "transport_petsc_options": configuration.get(
            "transport_petsc_options"
        ),
        "snes_atol_policy": configuration.get("snes_atol_policy"),
    }
    encoded = json.dumps(
        selected, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def atomic_write_json(path, payload):
    """Write JSON by rename so interrupted jobs never publish half a record."""
    path = os.fspath(path)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + ".tmp.",
                               dir=directory)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise


def host_provenance(environ=None):
    """Where a timing record was measured, for the record only.

    The matrix is run at several sites, and ``seconds_per_step`` has no
    warm-up excluded, so a record says which site and node measured it and
    whether the kernel cache it started from was empty (a cold first step is
    then the compiler, not the solver).  ``jit_cache_dir`` is None when the
    job left the cache at Firedrake's default location, whose state is not
    inspected.  Call this before anything compiles.

    None of this enters a cache manifest or a solver fingerprint: the same
    cache is valid wherever it is read.
    """
    environ = os.environ if environ is None else environ
    cache_dir = environ.get("PYOP2_CACHE_DIR") or None
    was_empty = None
    if cache_dir is not None:
        try:
            was_empty = not any(os.scandir(cache_dir))
        except OSError:
            was_empty = True
    return {
        "site": environ.get("ISMIP7_SITE") or None,
        "hostname": socket.gethostname(),
        "slurm_job_id": environ.get("SLURM_JOB_ID") or None,
        "slurm_partition": environ.get("SLURM_JOB_PARTITION") or None,
        "slurm_nodelist": environ.get("SLURM_JOB_NODELIST") or None,
        "container": environ.get("ISMIP7_CONTAINER") or None,
        "jit_cache_dir": cache_dir,
        "jit_cache_was_empty": was_empty,
    }


def atomic_write_status(path, state, **fields):
    tokens = [state]
    for key, value in fields.items():
        if value is None or value == "":
            continue
        text = str(value).replace(" ", "_").replace("\n", "_")
        tokens.append(f"{key}={text}")
    path = os.fspath(path)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + ".tmp.",
                               dir=directory)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(" ".join(tokens) + "\n")
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise


def reason_diverged(reason):
    text = str(reason)
    if text.startswith("DIVERGED"):
        return True
    try:
        return int(text) < 0
    except ValueError:
        return False


def diverged_reasons(record, summary_name):
    reasons = record.get(summary_name, {}).get("reason_counts", {})
    return {
        str(reason): count
        for reason, count in reasons.items()
        if reason_diverged(reason)
    }


def validate_cache_manifest(
    manifest,
    *,
    lc,
    lc_coarse,
    cache_path=None,
    source_sha256=None,
    mesh_sha256=None,
    solver_fingerprint=None,
):
    """Return ``(valid, detail)`` for a target-mesh initial-state cache."""
    expected = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_role": CACHE_ROLE,
        "lc": int(lc),
        "lc_coarse": int(lc_coarse),
        "buffer_m": BUFFER_M,
        "diagnostic_solver_mode": CACHE_SOLVER_MODE,
        "friction": "budd",
        "friction_gate": BUDD_SHELF_GATE,
        "geometry_space": "dg0",
        "n_flow": 3.0,
        "a4_factor": 1.0,
        "t_yr": MATRIX_T_START,
        "mesh_basename": mesh_basename(lc, lc_coarse),
        "geometry_source_method": TARGET_MESH_GEOMETRY_METHOD,
    }
    allowed_sources = {
        SOURCE_INVERSION_BASENAME,
        mesh_inversion_basename(lc, lc_coarse),
    }
    return _validate_manifest(
        manifest, expected, allowed_sources,
        cache_path=cache_path, source_sha256=source_sha256,
        mesh_sha256=mesh_sha256, solver_fingerprint=solver_fingerprint,
    )


def validate_map_check_manifest(
    manifest,
    *,
    lc,
    lc_coarse,
    buffer_m,
    friction,
    source_basename,
    mesh_name=None,
    cache_path=None,
    source_sha256=None,
    mesh_sha256=None,
    solver_fingerprint=None,
):
    """Return ``(valid, detail)`` for a map-check initial-state cache: one
    released MAP transferred onto one target mesh under one law. The role
    keeps it apart from the campaign's caches, and the source is that MAP."""
    expected = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_role": MAP_CHECK_CACHE_ROLE,
        "lc": int(lc),
        "lc_coarse": int(lc_coarse),
        "buffer_m": int(buffer_m),
        "diagnostic_solver_mode": CACHE_SOLVER_MODE,
        "friction": friction,
        "friction_gate": BUDD_SHELF_GATE if friction == "budd" else None,
        "geometry_space": "dg0",
        "n_flow": 3.0,
        "a4_factor": 1.0,
        "t_yr": MATRIX_T_START,
        "mesh_basename": mesh_name or mesh_basename(lc, lc_coarse, buffer_m),
        "geometry_source_method": TARGET_MESH_GEOMETRY_METHOD,
    }
    return _validate_manifest(
        manifest, expected, {os.path.basename(os.fspath(source_basename))},
        cache_path=cache_path, source_sha256=source_sha256,
        mesh_sha256=mesh_sha256, solver_fingerprint=solver_fingerprint,
    )


def _validate_manifest(
    manifest,
    expected,
    allowed_sources,
    *,
    cache_path=None,
    source_sha256=None,
    mesh_sha256=None,
    solver_fingerprint=None,
):
    for key, value in expected.items():
        actual = manifest.get(key)
        if isinstance(value, float):
            try:
                matches = math.isclose(float(actual), value, abs_tol=1e-12)
            except (TypeError, ValueError):
                matches = False
        else:
            matches = actual == value
        if not matches:
            return False, f"cache {key}={actual!r}; expected {value!r}"
    source_basename = manifest.get("source_inversion_basename")
    if source_basename not in allowed_sources:
        return False, (
            f"cache source_inversion_basename={source_basename!r}; "
            f"expected one of {sorted(allowed_sources)}"
        )
    for key in ("source_inversion_sha256", "source_mesh_sha256"):
        if not manifest.get(key):
            return False, f"cache {key} is missing"
    for key in ("geometry_source", "geometry_source_basename"):
        if not manifest.get(key):
            return False, f"cache {key} is missing"
    if os.path.basename(os.fspath(manifest["geometry_source"])) != manifest[
        "geometry_source_basename"
    ]:
        return False, "cache geometry-source basename is inconsistent"
    checkpoint_fields = manifest.get("checkpoint_fields")
    if not isinstance(checkpoint_fields, list):
        return False, "cache checkpoint_fields is missing"
    missing_fields = sorted(
        set(cache_required_fields(expected.get("friction", "budd")))
        - set(checkpoint_fields)
    )
    if missing_fields:
        return False, (
            "cache checkpoint is incomplete; missing fields: "
            + ", ".join(missing_fields)
        )
    if cache_path is not None:
        actual = os.path.realpath(os.fspath(manifest.get("cache_path", "")))
        expected_path = os.path.realpath(os.fspath(cache_path))
        # A lane may restart from the pristine prepare copy, whose manifest
        # (a hard link of the original) still names the published path.
        if actual != expected_path and pristine_sibling(actual) != expected_path:
            return False, f"cache path={actual!r}; expected {expected_path!r}"
        if not os.path.isfile(expected_path):
            return False, f"cache file is missing: {expected_path}"
    if source_sha256 is not None \
            and manifest.get("source_inversion_sha256") != source_sha256:
        return False, "cache source-inversion checksum is stale"
    if mesh_sha256 is not None \
            and manifest.get("source_mesh_sha256") != mesh_sha256:
        return False, "cache source-mesh checksum is stale"
    if solver_fingerprint is not None and manifest.get(
        "solver_configuration_fingerprint"
    ) != solver_fingerprint:
        return False, "cache solver configuration is stale"
    return True, "cache provenance matches"


def validate_timing_record(
    record,
    *,
    lc=None,
    lc_coarse=None,
    ncores=None,
    require_success=True,
    timing_kind="matrix",
    timing_tag=None,
):
    """Validate one lane record against the contract its tag names.

    Step count, the 2.5 km timestep and the physics contract are all read
    from ``timing_tag`` (default: this process's campaign tag), never from
    module constants, so a report built for one ladder rung cannot accept
    records from another. Returns ``(valid, detail)``.
    """
    if timing_tag is None:
        timing_tag = campaign_tag()
    if record.get("record_schema_version") != RECORD_SCHEMA_VERSION:
        return False, "record schema is not the cached-strict schema"
    if require_success and record.get("run_status") != "success":
        failure = record.get("failure") or {}
        detail = failure.get("category", record.get("run_status", "failed"))
        phase = failure.get("phase")
        return False, detail + (f" at {phase}" if phase else "")
    if record.get("timing_kind") != timing_kind:
        return False, f"timing_kind={record.get('timing_kind')!r}"
    if record.get("timing_tag") != timing_tag:
        return False, f"timing_tag={record.get('timing_tag')!r}"
    try:
        spec = parse_campaign_tag(timing_tag)
    except ValueError as exc:
        return False, str(exc)
    expected_kind = "cache_probe" if spec["lane"] == "probe" else "matrix"
    if timing_kind != expected_kind:
        return False, (
            f"tag {timing_tag!r} names a {expected_kind} lane, not {timing_kind}"
        )
    contract = CONTRACTS[spec["contract"]]
    if record.get("apparent_mb_mode") != contract["apparent_mb_mode"]:
        return False, (
            f"apparent_mb_mode={record.get('apparent_mb_mode')!r}; contract "
            f"{spec['contract']!r} expects {contract['apparent_mb_mode']!r}"
        )
    try:
        cap = float(record.get("apparent_mb_cap_m_per_yr", 0.0))
    except (TypeError, ValueError):
        cap = math.nan
    if not math.isclose(
        cap, contract["apparent_mb_cap_m_per_yr"], abs_tol=1e-12
    ):
        return False, (
            f"apparent MB cap={cap!r}; contract expects "
            f"{contract['apparent_mb_cap_m_per_yr']!r}"
        )
    fixed = bool(record.get("fixed_front", False))
    if fixed != contract["fixed_front"]:
        return False, (
            "lane fixed the calving front but the contract does not"
            if fixed
            else "contract fixes the calving front but the lane did not"
        )
    if fixed:
        try:
            masked = int(record.get("fixed_front_cells") or 0)
        except (TypeError, ValueError):
            masked = 0
        if masked <= 0:
            return False, "fixed front masked no cells (clamped initial state?)"
    if record.get("diagnostic_solver_mode") != spec["solver"]:
        return False, (
            "diagnostic_solver_mode="
            f"{record.get('diagnostic_solver_mode')!r}; tag names "
            f"{spec['solver']!r}"
        )

    checks = (("lc", lc), ("lc_coarse", lc_coarse), ("ncores", ncores))
    for key, expected in checks:
        if expected is not None and record.get(key) != int(expected):
            return False, f"record {key}={record.get(key)!r}; expected {expected}"

    # A current-campaign lane must descend from the campaign source MAP or
    # from its own mesh's invert (the _reinverted lanes). A record left by an
    # earlier source (the old-gate MAP before the 2026-09-17 promotion) is
    # never a pass; archived campaigns are judged by their own tag only.
    source = record.get("initial_state_source")
    if source and spec["version"] == CAMPAIGN_VERSION:
        allowed = {SOURCE_INVERSION_BASENAME}
        try:
            allowed.add(
                mesh_inversion_basename(record["lc"], record["lc_coarse"])
            )
        except (KeyError, TypeError, ValueError):
            pass
        if source not in allowed:
            return False, (
                f"record initial state descends from {source!r}; expected "
                f"one of {sorted(allowed)}"
            )

    steps = spec["steps"]
    try:
        record_lc = int(record["lc"])
        dt = expected_dt(record_lc, spec["dt_2500"], spec["version"])
        t_end = expected_t_end(record_lc, steps, spec["dt_2500"], spec["version"])
        interval_ok = (
            math.isclose(float(record["t_start"]), MATRIX_T_START,
                         abs_tol=1e-12)
            and math.isclose(float(record["dt"]), dt, abs_tol=1e-12)
            and int(record["nsteps"]) == steps
            and int(record["completed_steps"]) == steps
            and math.isclose(float(record["t_end"]), t_end, abs_tol=1e-12)
            and math.isclose(float(record["t_final"]), t_end, abs_tol=1e-12)
        )
    except (KeyError, TypeError, ValueError):
        interval_ok = False
    if not interval_ok:
        return False, (
            f"record did not complete the required {steps}-step interval "
            f"(dt {expected_dt(int(record.get('lc', 2500)), spec['dt_2500'], spec['version']):g} yr)"
        )

    diagnostic = diverged_reasons(record, "diagnostic_solve_summary")
    if diagnostic:
        return False, f"diagnostic divergence: {diagnostic}"
    labels = [stat.get("label", "")
              for stat in record.get("diagnostic_solves", [])]
    if len(labels) < steps:
        return False, f"record has fewer than {steps} diagnostic solves"
    if any("rescue" in label or "trust-region" in label for label in labels):
        return False, "record used the rescue ladder"
    if any(label != f"step-{index}-direct"
           for index, label in enumerate(labels, 1)):
        return False, f"unexpected diagnostic solve sequence: {labels}"

    transport = diverged_reasons(record, "transport_solve_summary")
    if transport:
        return False, f"transport divergence: {transport}"
    transport_summary = record.get("transport_solve_summary", {})
    if int(transport_summary.get("count", 0)) != steps:
        return False, f"record does not contain exactly {steps} transport solves"
    transport_labels = [
        stat.get("label", "") for stat in record.get("transport_solves", [])
    ]
    expected_transport_labels = [
        f"step-{index}-substep-1/1"
        for index in range(1, steps + 1)
    ]
    if transport_labels != expected_transport_labels:
        return False, f"unexpected transport solve sequence: {transport_labels}"
    residual = transport_summary.get("mass_residual_gt_max")
    try:
        residual_ok = (
            residual is not None
            and math.isfinite(float(residual))
            and float(residual) <= MASS_RESIDUAL_TOL_GT
        )
    except (TypeError, ValueError):
        residual_ok = False
    if not residual_ok:
        return False, f"transport mass residual {residual!r} exceeds tolerance"
    step_residual = record.get("step_mass_residual_gt_max")
    try:
        step_residual_ok = (
            step_residual is not None
            and math.isfinite(float(step_residual))
            and float(step_residual) <= MASS_RESIDUAL_TOL_GT
        )
    except (TypeError, ValueError):
        step_residual_ok = False
    if not step_residual_ok:
        return False, f"step mass residual {step_residual!r} exceeds tolerance"

    cache = record.get("cache_validation", {})
    if cache.get("status") != "valid":
        return False, f"cache validation={cache.get('status')!r}"
    if record.get("rescue_enabled") is not False:
        return False, "strict timing record did not disable rescue"
    if record.get("solver_configuration", {}).get("subcycles") != [1]:
        return False, "strict timing record did not restrict subcycles to [1]"
    if record.get("timing_scope") != "transient_loop_only":
        return False, f"timing_scope={record.get('timing_scope')!r}"
    what = "timing run" if timing_kind == "matrix" else "probe"
    return True, (
        f"completed strict cached {steps}-step {what} "
        f"({spec['contract']} contract)"
    )


def synthetic_record(lc, lc_coarse, ncores, timing_tag=None, timing_kind="matrix"):
    """A minimal PASSING record for ``timing_tag`` (self-tests, dry runs)."""
    timing_tag = campaign_tag() if timing_tag is None else timing_tag
    spec = parse_campaign_tag(timing_tag)
    contract = CONTRACTS[spec["contract"]]
    steps = spec["steps"]
    dt = expected_dt(lc, spec["dt_2500"], spec["version"])
    t_end = expected_t_end(lc, steps, spec["dt_2500"], spec["version"])
    return {
        "record_schema_version": RECORD_SCHEMA_VERSION,
        "run_status": "success",
        "failure": None,
        "lc": int(lc),
        "lc_coarse": int(lc_coarse),
        "ncores": int(ncores),
        "timing_kind": timing_kind,
        "timing_tag": timing_tag,
        "apparent_mb_mode": contract["apparent_mb_mode"],
        "apparent_mb_cap_m_per_yr": contract["apparent_mb_cap_m_per_yr"],
        "fixed_front": contract["fixed_front"],
        "fixed_front_cells": 48843 if contract["fixed_front"] else 0,
        "diagnostic_solver_mode": spec["solver"],
        "t_start": MATRIX_T_START,
        "dt": dt,
        "nsteps": steps,
        "completed_steps": steps,
        "t_end": t_end,
        "t_final": t_end,
        "diagnostic_solve_summary": {"reason_counts": {"2": steps}},
        "diagnostic_solves": [
            {"label": f"step-{index}-direct"} for index in range(1, steps + 1)
        ],
        "transport_solve_summary": {
            "count": steps,
            "reason_counts": {"2": steps},
            "mass_residual_gt_max": 1e-8,
        },
        "transport_solves": [
            {"label": f"step-{index}-substep-1/1"}
            for index in range(1, steps + 1)
        ],
        "step_mass_residual_gt_max": 1e-8,
        "cache_validation": {"status": "valid"},
        "rescue_enabled": False,
        "solver_configuration": {"subcycles": [1]},
        "timing_scope": "transient_loop_only",
    }


def selftest():
    """Pure-Python checks of the tag/contract/record machinery."""
    for solver in LANE_SOLVER_MODES:
        for steps, dt in ((5, 0.25), (10, 0.125), (20, 0.0625)):
            for contract in CONTRACTS:
                tag = campaign_tag(steps, dt, contract, solver)
                spec = parse_campaign_tag(tag)
                assert spec == {
                    "solver": solver, "steps": steps, "dt_2500": dt,
                    "contract": contract, "lane": None,
                    "version": CAMPAIGN_VERSION,
                }, (tag, spec)
                assert parse_campaign_tag(lane_tag("invert", tag))["lane"] == "reinverted"
                assert parse_campaign_tag(
                    probe_tag(steps, dt, contract, solver)
                ) == dict(spec, lane="probe")
    assert campaign_tag(5, 0.25, "strict") == (
        "scpc_mumps_5step_dt0p25at2500_dg0_logvelnet_cached_strict_v4"
    )
    assert campaign_tag(10, 0.125, "divfront").endswith("dt0p125at2500_dg0_logvelnet_cached_strict_v4_divfront")
    # The old-law archive still parses, and is never the current campaign.
    archived = parse_campaign_tag(
        "scpc_mumps_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v3"
    )
    assert archived["version"] == 3 and archived["version"] != CAMPAIGN_VERSION
    assert CACHE_TAG.endswith(f"_v{CAMPAIGN_VERSION}")
    # The lane solver leads the tag; the caches stay the prepare solver's, so
    # a GAMG campaign reads the cache files the MUMPS campaign prepared.
    assert campaign_tag(10, 0.125, "strict", "scpc_gamg") == (
        "scpc_gamg_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v4"
    )
    assert CACHE_TAG == f"scpc_mumps_dg0_logvelnet_v{CAMPAIGN_VERSION}"
    assert cache_stem(2500, 25000).startswith("initial_state_scpc_mumps_")
    for unqualified in ("full_mumps", "schur_gamg", "iterative"):
        try:
            campaign_tag(10, 0.125, "strict", unqualified)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{unqualified} must not name a campaign")
    # The promoted source is the 2500/25000 invert MAP; the manager must
    # recognise that mesh from the basename alone, at any maxiter.
    assert SOURCE_INVERSION_BASENAME == mesh_inversion_basename(2500, 25000, 250)
    assert mesh_inversion_source_mesh(SOURCE_INVERSION_BASENAME) == (2500, 25000)
    assert mesh_inversion_source_mesh(
        "/x/inversion_icepack2_budd_n3_dg0_logvelnet_2000_20000_5iter.h5"
    ) == (2000, 20000)
    assert mesh_inversion_source_mesh(
        "inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5"
    ) is None
    assert contract_exports("strict") == {"ISMIP7_AMB_CAP": "0"}
    assert contract_exports("divfront") == {
        "ISMIP7_AMB_CAP": "0", "ISMIP7_APPARENT_MB": "div", "ISMIP7_FIXED_FRONT": "1"
    }
    assert contract_exports("front") == {"ISMIP7_AMB_CAP": "0", "ISMIP7_FIXED_FRONT": "1"}

    # The LC-scaled step stops at MATRIX_DT_MAX from v4 on; a rung above the
    # cap keeps its own step, and archived campaigns keep the plain rule.
    assert [expected_dt(lc, 0.125) for lc in LCS] == [0.025, 0.05, 0.1, 0.125, 0.125]
    assert expected_t_end(5000, 10, 0.125) == MATRIX_T_START + 1.25
    assert expected_dt(5000, 0.0625) == 0.125 and expected_dt(2500, 0.25) == 0.25
    assert expected_dt(5000, 0.25) == 0.25
    assert expected_dt(5000, 0.125, version=3) == 0.25

    def check(record, expected_valid, fragment=None, **kwargs):
        valid, detail = validate_timing_record(record, **kwargs)
        assert valid is expected_valid, (expected_valid, detail, kwargs)
        if fragment is not None:
            assert fragment in detail, (fragment, detail)

    tag10 = campaign_tag(10, 0.125, "strict")
    good = synthetic_record(2500, 25000, 16, tag10)
    check(good, True, "10-step", lc=2500, lc_coarse=25000, ncores=16, timing_tag=tag10)
    # Initial-state provenance: the campaign source, or the mesh's own invert.
    check(dict(good, initial_state_source=SOURCE_INVERSION_BASENAME), True,
          timing_tag=tag10)
    good2000 = synthetic_record(2000, 20000, 16, tag10)
    check(dict(good2000, initial_state_source=mesh_inversion_basename(2000, 20000)),
          True, timing_tag=tag10, lc=2000, lc_coarse=20000)
    check(dict(good2000, initial_state_source=mesh_inversion_basename(2500, 50000)),
          False, "descends from", timing_tag=tag10)
    check(dict(good, initial_state_source=(
        "inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5")),
          False, "descends from", timing_tag=tag10)
    good5000 = synthetic_record(5000, 50000, 16, tag10)
    assert good5000["dt"] == 0.125
    check(good5000, True, "10-step", lc=5000, timing_tag=tag10)
    uncapped = dict(good5000, dt=0.25, t_end=MATRIX_T_START + 2.5,
                    t_final=MATRIX_T_START + 2.5)
    check(uncapped, False, "(dt 0.125 yr)", timing_tag=tag10)
    check(good, False, "timing_tag=", timing_tag=campaign_tag(5, 0.25))
    short = dict(good, completed_steps=9, t_final=good["t_end"] - 0.125)
    check(short, False, "10-step interval", timing_tag=tag10)
    check(dict(good, apparent_mb_mode="div"), False, "apparent_mb_mode", timing_tag=tag10)
    check(dict(good, fixed_front=True), False, "contract does not", timing_tag=tag10)
    tagdf = campaign_tag(10, 0.125, "divfront")
    gooddf = synthetic_record(2500, 25000, 16, tagdf)
    check(gooddf, True, "divfront", timing_tag=tagdf)
    check(dict(gooddf, fixed_front=False), False, "lane did not", timing_tag=tagdf)
    check(dict(gooddf, fixed_front_cells=0), False, "masked no cells", timing_tag=tagdf)
    check(dict(gooddf, apparent_mb_cap_m_per_yr=5.0), False, "cap=", timing_tag=tagdf)
    ptag = probe_tag(10, 0.125, "div")
    probe = synthetic_record(2500, 25000, 16, ptag, timing_kind="cache_probe")
    check(probe, True, "probe", timing_kind="cache_probe", timing_tag=ptag)
    check(dict(probe, timing_kind="matrix"), False, "names a cache_probe lane",
          timing_tag=ptag)
    # A record passes only under the solver its own tag names.
    gtag = campaign_tag(10, 0.125, "strict", "scpc_gamg")
    goodg = synthetic_record(2500, 25000, 16, gtag)
    assert goodg["diagnostic_solver_mode"] == "scpc_gamg"
    check(goodg, True, "10-step", lc=2500, lc_coarse=25000, ncores=16, timing_tag=gtag)
    check(dict(goodg, diagnostic_solver_mode="scpc_mumps"), False,
          "tag names 'scpc_gamg'", timing_tag=gtag)
    check(dict(good, diagnostic_solver_mode="scpc_gamg"), False,
          "tag names 'scpc_mumps'", timing_tag=tag10)
    check(goodg, False, "timing_tag=", timing_tag=tag10)
    old = dict(good, record_schema_version=RECORD_SCHEMA_VERSION - 1)
    check(old, False, "schema", timing_tag=tag10)
    assert pristine_sibling("/x/initial_state_a.h5") == "/x/initial_state_a.prepare.h5"
    assert pristine_sibling("/x/initial_state_a.prepare.h5") == "/x/initial_state_a.prepare.h5"

    # MAP checks: named apart from the campaign at every level.
    stem = map_check_stem(
        "/m/inversion_icepack2_rc_n3_dg0_logvelnet_2000_int5000_bilap_snap20260922_0948.h5"
    )
    mtag = map_check_tag(stem, "regularized_coulomb", "scpc_gamg")
    assert mtag.startswith("mapcheck_rc_scpc_gamg_10step_inversion_icepack2_rc")
    try:
        parse_campaign_tag(mtag)
    except ValueError:
        pass
    else:
        raise AssertionError("a map-check tag must never parse as a campaign tag")
    assert not map_check_cache_stem(stem, "budd", 1000, 10000, 20000).startswith(
        "initial_state_"
    )
    assert map_check_run_tag(stem, "budd", "1000") == "mapcheck_budd_snap20260922_0948_1000"
    assert mesh_inversion_source_mesh(stem + ".h5") is None
    for bad in ("budd_legacy", "rc"):
        try:
            map_check_tag(stem, bad, "scpc_gamg")
        except ValueError:
            pass
        else:
            raise AssertionError(f"{bad} must not name a map check")
    campaign_manifest = {
        "cache_schema_version": CACHE_SCHEMA_VERSION, "cache_role": CACHE_ROLE,
        "lc": 1000, "lc_coarse": 10000, "buffer_m": BUFFER_M,
        "diagnostic_solver_mode": CACHE_SOLVER_MODE, "friction": "budd",
        "friction_gate": BUDD_SHELF_GATE, "geometry_space": "dg0",
        "n_flow": 3.0, "a4_factor": 1.0, "t_yr": MATRIX_T_START,
        "mesh_basename": mesh_basename(1000, 10000),
        "geometry_source_method": TARGET_MESH_GEOMETRY_METHOD,
        "source_inversion_basename": SOURCE_INVERSION_BASENAME,
        "source_inversion_sha256": "a", "source_mesh_sha256": "b",
        "geometry_source": "/d/bedmachine.nc", "geometry_source_basename": "bedmachine.nc",
        "checkpoint_fields": list(CACHE_REQUIRED_FIELDS),
    }
    assert validate_cache_manifest(campaign_manifest, lc=1000, lc_coarse=10000)[0]
    valid, detail = validate_map_check_manifest(
        campaign_manifest, lc=1000, lc_coarse=10000, buffer_m=BUFFER_M,
        friction="budd", source_basename=SOURCE_INVERSION_BASENAME,
    )
    assert not valid and "cache_role" in detail, detail
    rc_manifest = dict(
        campaign_manifest, cache_role=MAP_CHECK_CACHE_ROLE,
        friction="regularized_coulomb", friction_gate=None,
        source_inversion_basename=stem + ".h5",
    )
    valid, detail = validate_map_check_manifest(
        rc_manifest, lc=1000, lc_coarse=10000, buffer_m=BUFFER_M,
        friction="regularized_coulomb", source_basename="/x/" + stem + ".h5",
    )
    assert valid, detail
    valid, detail = validate_cache_manifest(rc_manifest, lc=1000, lc_coarse=10000)
    assert not valid and "cache_role" in detail, detail
    valid, detail = validate_map_check_manifest(
        rc_manifest, lc=1000, lc_coarse=10000, buffer_m=BUFFER_M,
        friction="budd", source_basename=stem + ".h5",
    )
    assert not valid and "friction" in detail, detail
    # A regularized Coulomb state has no N_ref and publishes without it; a
    # Budd state without one is incomplete.
    assert "N_ref" not in cache_required_fields("regularized_coulomb")
    assert "N_ref" in cache_required_fields("budd")
    rc_no_nref = dict(
        rc_manifest,
        checkpoint_fields=list(cache_required_fields("regularized_coulomb")),
    )
    valid, detail = validate_map_check_manifest(
        rc_no_nref, lc=1000, lc_coarse=10000, buffer_m=BUFFER_M,
        friction="regularized_coulomb", source_basename=stem + ".h5",
    )
    assert valid, detail
    budd_no_nref = dict(
        campaign_manifest,
        checkpoint_fields=list(cache_required_fields("regularized_coulomb")),
    )
    valid, detail = validate_cache_manifest(budd_no_nref, lc=1000, lc_coarse=10000)
    assert not valid and "N_ref" in detail, detail
    print("timing_campaign selftest OK")


assert len(mesh_rows()) == 10
assert len(scout_lanes()) == 10
assert len(scaling_lanes()) == 10
assert len(planned_lanes()) == 20


def _print_lanes(lanes):
    for lane in lanes:
        print(" ".join(str(value) for value in lane))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("rows", "planned", "scout", "scale", "selftest", "tag"):
        subparsers.add_parser(command)
    validate = subparsers.add_parser("validate-cache")
    validate.add_argument("--manifest", required=True)
    validate.add_argument("--cache", required=True)
    validate.add_argument("--lc", required=True, type=int)
    validate.add_argument("--lc-coarse", required=True, type=int)
    validate.add_argument("--source-sha256")
    validate.add_argument("--mesh-sha256")
    validate.add_argument("--solver-fingerprint")
    args = parser.parse_args()

    if args.command == "rows":
        _print_lanes(mesh_rows())
    elif args.command == "planned":
        _print_lanes(planned_lanes())
    elif args.command == "scout":
        _print_lanes(scout_lanes())
    elif args.command == "scale":
        _print_lanes(scaling_lanes())
    elif args.command == "selftest":
        selftest()
    elif args.command == "tag":
        print(campaign_tag())
    elif args.command == "validate-cache":
        try:
            with open(args.manifest) as stream:
                manifest = json.load(stream)
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"INVALID: {exc}") from exc
        valid, detail = validate_cache_manifest(
            manifest,
            lc=args.lc,
            lc_coarse=args.lc_coarse,
            cache_path=args.cache,
            source_sha256=args.source_sha256,
            mesh_sha256=args.mesh_sha256,
            solver_fingerprint=args.solver_fingerprint,
        )
        if not valid:
            raise SystemExit(f"INVALID: {detail}")
        print(f"VALID: {detail}")


if __name__ == "__main__":
    main()
