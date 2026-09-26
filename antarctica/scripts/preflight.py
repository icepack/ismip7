#!/usr/bin/env python3
r"""Data preflight for the ISMIP7 core experiments.

Answers "which experiments can run on this machine right now?" in a few
seconds, checking every input each core needs: mesh, MAP inversion,
boundary ids, per-basin K, RACMO, the OI climatology (core 11's stopgap),
and the (ESM, scenario) atmosphere/ocean trees over the run period, the
control's `ctrl` ocean among them, and under the SMB-elevation feedback the
`dacabfdz` gradient each core reads. Honors the same
environment knobs as the runs (ISMIP7_LC, ISMIP7_FRICTION,
ISMIP7_OI_VERSION, ISMIP7_SMB_ELEVATION_FEEDBACK, ...).

Usage:
    python scripts/preflight.py
    ISMIP7_LC=500 ISMIP7_FRICTION=regularized_coulomb python scripts/preflight.py
"""
import os, sys, glob

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_ANT = os.path.dirname(_SCRIPTS)
_PROJECT = os.path.dirname(_ANT)
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from icepack2_tools.forcing import (
    ISMIP7Atmosphere, ISMIP7Ocean, OCX, OCX_ATMOSPHERE_SOURCE, SMB_GRADIENT,
    _oi_climatology_path, _find_ismip7_data,
)
from icepack2_tools.boundary import sidecar_path
from icepack2_tools.naming import map_basename
from icepack2_tools.climatology import clim_start, clim_end, clim_scenario
from icepack2_tools.runconfig import (
    obs_data_root,
    calving_law as _calving_law, calving_sigma_max as _calving_sigma_max,
    friction as _friction, geometry_space as _geometry_space, lc as _lc,
    lc_coarse as _lc_coarse, ocx_forcing as _ocx_forcing, ocx_ocean as _ocx_ocean,
    mesh_override, deltat_per_basin_npz as _deltat_per_basin_npz,
    k_per_basin_candidates, smb_elevation_feedback as _smb_elevation_feedback,
)
DATA_DIR = obs_data_root()
from mesh_naming import get_buffer_m, mesh_filename

MESH_DIR = os.path.join(_ANT, "mesh")
RESULTS_DIR = os.path.join(_ANT, "results")

lc = _lc()
lc_coarse = _lc_coarse()
friction = _friction()
oi_version = os.environ.get("ISMIP7_OI_VERSION", "30_sep")
CLIM_START = clim_start()
CLIM_END = clim_end()
CLIM_SCENARIO = clim_scenario()
root = _find_ismip7_data()

CORES = [
    (1, "hist CESM2-WACCM", "CESM2-WACCM", "historical", 1850, 2014),
    (2, "hist MRI-ESM2-0", "MRI-ESM2-0", "historical", 1850, 2014),
    (3, "ssp370 CESM2-WACCM", "CESM2-WACCM", "ssp370", 2015, 2100),
    (4, "ssp370 MRI-ESM2-0", "MRI-ESM2-0", "ssp370", 2015, 2100),
    (5, "ssp126 CESM2-WACCM", "CESM2-WACCM", "ssp126", 2015, 2300),
    (6, "ssp126 MRI-ESM2-0", "MRI-ESM2-0", "ssp126", 2015, 2300),
    (7, "ssp585 CESM2-WACCM", "CESM2-WACCM", "ssp585", 2015, 2300),
    (8, "ssp585 MRI-ESM2-0", "MRI-ESM2-0", "ssp585", 2015, 2300),
    (9, "CTRL2015 (CESM2-WACCM clim)", "CESM2-WACCM", None, 2015, 2300),
    (10, "CTRL2015 (MRI-ESM2-0 clim)", "MRI-ESM2-0", None, 2015, 2300),
    (11, "OCX obs-constrained", None, None, 1979, 2025),
]


def atm_years(esm, scenario, var="acabf-anomaly"):
    r"""The years the reader can open, by the reader's own match: a gate that
    counts files the loader would not find clears a run that then fails."""
    return ISMIP7Atmosphere(esm=esm, scenario=scenario).available_years(var)


def clim_pool_years(esm, var):
    r"""Years the runs actually pool to build the reference climatology:
    historical + CLIM_SCENARIO, restricted to the CLIM_START-CLIM_END window.

    The window filter matters - both the control's ``compute_climatology`` and
    the projections' ``smb_scheme`` keep only years inside it, so an ESM with
    plenty of files outside the window still yields an empty pool. Deduplicated
    because a year present in both scenarios is one year of coverage.
    """
    return sorted({y for y in (atm_years(esm, "historical", var)
                               + atm_years(esm, CLIM_SCENARIO, var))
                   if CLIM_START <= y <= CLIM_END})


def clim_pool_gaps(esm, var, years=None):
    r"""Years of the climatology window the pool is MISSING.

    A partial pool is a failure, not a warning: the pooled mean is the baseline
    the aSMB anomalies are re-referenced to, so a core built from 15 of the 30
    years is referenced to a different mean than a sibling core with full
    coverage, while both are differenced against the same CTRL. That is the
    baseline mismatch this gate exists to catch, so it is held to the same
    standard as the run-scenario coverage check below.
    """
    have = clim_pool_years(esm, var) if years is None else years
    return sorted(set(range(CLIM_START, CLIM_END + 1)) - set(have))


def pool_status(esm, var, what_empty, what_partial):
    r"""``(bucket, detail)`` for the climatology pool: bucket is ``"ok"``,
    ``"partial"`` or ``"empty"``.

    The three are graded because the RUNTIME grades them: an empty pool makes
    the run take a different code path (the control raises, a projection
    silently falls back to full acabf(t)), while a partial pool runs and warns.
    A gate that blocked what the runtime happily runs would train people to
    ignore it, so a partial pool reports PARTIAL, not BLOCKED.

    The caller passes both consequences and the bucket picks between them, so
    the emptiness test lives in one place and the year list is read from disk
    once per (esm, var) rather than once per caller condition.
    """
    have = clim_pool_years(esm, var)
    gaps = clim_pool_gaps(esm, var, have)
    if not gaps:
        return "ok", None
    span = f"covers {have[0]}-{have[-1]}" if have else "no years"
    return ("partial" if have else "empty",
            f"{esm} {var} climatology pool (historical+{CLIM_SCENARIO}) "
            f"{span}, {len(gaps)} of {CLIM_START}-{CLIM_END} missing "
            f"({gaps[0]}..{gaps[-1]}): "
            f"{what_partial if have else what_empty}")


def ocx_atm_years():
    return ISMIP7Atmosphere(esm=OCX_ATMOSPHERE_SOURCE, scenario=OCX).available_years("acabf")


def ocx_ocean_cover():
    r"""The years both OCX ocean variables cover, or None."""
    ocean = ISMIP7Ocean(scenario=OCX, variant=_ocx_ocean())
    covers = [ocean.coverage(v) for v in ("tf", "so")]
    if None in covers:
        return None
    return max(c[0] for c in covers), min(c[1] for c in covers)


def ocean_cover(esm, scenario):
    r"""``(first, last)`` ocean forcing year on disk, as the reader sees it."""
    return ISMIP7Ocean(esm=esm, scenario=scenario).coverage("tf")


def gradient_problem(esm, scenario, y0, y1):
    r"""What keeps the SMB gradient the feedback reads from serving ``y0`` to
    ``y1`` at an accepted version, or None. The reader's own rule
    (``coverage_problem``, ``version_problem``), the one
    ``forcing.SMBElevationFeedback.check`` refuses a run on."""
    atm = ISMIP7Atmosphere(esm=esm, scenario=scenario)
    problem = (atm.coverage_problem(y0, y1, SMB_GRADIENT)
               or atm.version_problem(SMB_GRADIENT))
    if problem is None:
        return None
    return (f"{problem} (the SMB-elevation feedback reads it; "
            f"ISMIP7_SMB_ELEVATION_FEEDBACK=0 runs without it)")


def shared_missing(warn=None):
    r"""Missing shared inputs. Non-fatal caveats are appended to ``warn``."""
    miss = []
    warn = warn if warn is not None else []
    # Front configuration: a mistyped law would otherwise surface only after
    # the forward's MAP load and initial solve.
    try:
        _calving_law()
        _calving_sigma_max()
    except ValueError as e:
        miss.append(str(e))
    # ISMIP7_MESH=checkpoint means the mesh inside the MAP; there is no file
    # to look for and the MAP check below covers it.
    mesh_fn = mesh_override() or mesh_filename(lc_coarse, lc, get_buffer_m())
    if not os.path.exists(mesh_fn):
        miss.append(f"mesh ({os.path.basename(mesh_fn)})")
    # The MAP the forward will actually load: ISMIP7_INVERSION if it names one
    # explicitly, else the one tagged with this geometry space, else the legacy
    # untagged (CG1) MAP it falls back to with a warning. A legacy MAP runs, but
    # its controls carry the CG1 front bias, so the run is a smoke test rather
    # than a result. An explicit override deliberately bypasses that lookup, so
    # setup_model raises on a missing file rather than falling back: report it
    # as a hard miss, exactly as the forward would.
    inv_override = os.environ.get("ISMIP7_INVERSION")
    inv = inv_override or os.path.join(MESH_DIR, map_basename(friction, lc))
    legacy = os.path.join(MESH_DIR, map_basename(friction, lc, geometry=False))
    if not os.path.exists(inv):
        if not inv_override and os.path.exists(legacy):
            warn.append(
                f"no {os.path.basename(inv)}; the forward would fall back to "
                f"{os.path.basename(legacy)} (inverted under a different "
                f"geometry space — smoke test only)"
            )
        else:
            miss.append(f"MAP ({os.path.basename(inv)})")
    # Same sidecar-name rule the solvers use: per-mesh preferred, shared file
    # as the fallback.
    bnd = sidecar_path(MESH_DIR, mesh_hint=mesh_fn)
    if not os.path.exists(bnd):
        miss.append(
            f"{os.path.basename(bnd)} (untracked — restore or regenerate "
            f"with make_boundary_ids.py)"
        )
    # The melt calibration the forward reads: per-basin deltaT at one K when
    # ISMIP7_DELTAT_PER_BASIN_NPZ names a file (the protocol's knob, and then
    # no driver reads a per-basin K), else a per-basin K file, found where the
    # drivers look for it.
    try:
        dT_npz = _deltat_per_basin_npz()
    except (FileNotFoundError, ValueError) as err:
        miss.append(str(err))
        dT_npz = err
    if dT_npz is None and not any(
            os.path.exists(c) for c in k_per_basin_candidates(RESULTS_DIR, lc)):
        miss.append("per-basin K npz (or ISMIP7_DELTAT_PER_BASIN_NPZ)")
    for d, pat, what in [
        (os.path.join(DATA_DIR, "bedmachine"), "*.nc", "BedMachine"),
        (os.path.join(DATA_DIR, "velocity"), "*.nc", "MEaSUREs velocity"),
    ]:
        if not glob.glob(os.path.join(d, pat)):
            miss.append(what)
    return miss


def racmo_ok():
    return os.path.exists(os.path.join(
        DATA_DIR, "racmo",
        "smbgl_monthlyS_ANT11_RACMO2.4p1_ERA5_197901_202312.nc",
    ))


def oi_ok():
    if root is None:
        return False
    return all(os.path.exists(_oi_climatology_path(root, v, oi_version))
               for v in ("tf", "so"))


def main():
    geom = _geometry_space()
    feedback = _smb_elevation_feedback()
    print(f"Preflight: lc={lc}, friction={friction}, geometry={geom}, "
          f"OI={oi_version}, climatology=historical+{CLIM_SCENARIO} "
          f"{CLIM_START}-{CLIM_END}")
    print(f"  SMB-elevation feedback: "
          + (f"on, every core needs its {SMB_GRADIENT}" if feedback
             else "off (ISMIP7_SMB_ELEVATION_FEEDBACK=0)"))
    warn = []
    base_missing = shared_missing(warn)
    if base_missing:
        print(f"  SHARED inputs missing: {', '.join(base_missing)}")
    else:
        print("  Shared inputs (mesh, MAP, bndids, K, BedMachine, velocity): OK")
    for w in warn:
        print(f"  WARNING: {w}")
    print(f"  RACMO baseline: {'OK' if racmo_ok() else 'MISSING (acabf fallback)'}")
    print("  Status: READY = every input present; PARTIAL = the run proceeds "
          "but warns and its provenance is degraded; BLOCKED = missing input")
    print()

    for core, title, esm, scenario, y0, y1 in CORES:
        miss = list(base_missing)
        degraded = []
        notes = []
        if core == 11 and _ocx_forcing() == "protocol":
            # projections/ocx.py refuses to start on anything less, so this
            # gate asks the same readers the same question.
            yrs = ocx_atm_years()
            gaps = sorted(set(range(y0, y1 + 1)) - set(yrs))
            if gaps:
                miss.append(f"OCX atmosphere acabf ({OCX_ATMOSPHERE_SOURCE}): "
                            + (f"{len(gaps)} of {y0}-{y1} missing" if yrs else "absent"))
            oc = ocx_ocean_cover()
            if oc is None or oc[0] > y0 or oc[1] < y1:
                miss.append(f"OCX ocean '{_ocx_ocean()}' tf/so"
                            + (f" covers {oc[0]}-{oc[1]}, need {y0}-{y1}" if oc else ": absent"))
            notes.append("K is fitted to the OI climatology, not to the OCX ocean: "
                         "read check_melt_bound.py --ocx first (discussion #48)")
        elif core == 11:
            notes.append("ISMIP7_OCX_FORCING=stopgap: RACMO2.4p1 + OI climatology, "
                         "not the ISMIP7 OCX product")
            if not racmo_ok():
                miss.append("RACMO (OCX SMB)")
            if not oi_ok():
                miss.append(f"OI climatology ({oi_version})")
        elif scenario is None:  # CTRL
            # the ESM's own ctrl ocean, which control/run.py refuses to start
            # without (icepack/ismip7#107)
            oc = ocean_cover(esm, "ctrl")
            if oc is None:
                miss.append(f"{esm}/ctrl ocean tf/so")
            elif oc[0] > y0 or oc[1] < y1 - 1:
                miss.append(f"ctrl ocean covers {oc[0]}-{oc[1]}, need {y0}-{y1}")
            if not racmo_ok():
                bucket, detail = pool_status(
                    esm, "acabf",
                    "the control would refuse to run",
                    "the constant SMB climatology would be a mean over part "
                    "of the window; the run warns and proceeds")
                if bucket == "empty":
                    miss.append(detail)
                elif bucket == "partial":
                    degraded.append(detail)
        else:
            # experiment.py re-references aSMB against the historical +
            # CLIM_SCENARIO acabf-anomaly pool and, on FileNotFoundError,
            # silently degrades to the full acabf(t) field - a different SMB
            # scheme from the one the CTRL this core is differenced against
            # uses. Only meaningful when RACMO is present: without it every
            # core degrades together, which the shared RACMO line reports.
            if racmo_ok():
                bucket, detail = pool_status(
                    esm, "acabf-anomaly",
                    "the run would silently fall back to full acabf(t), a "
                    "different SMB scheme than the CTRL",
                    "the aSMB re-reference baseline would differ from a "
                    "full-window sibling's; the run warns and proceeds")
                if bucket == "empty":
                    miss.append(detail)
                elif bucket == "partial":
                    degraded.append(detail)
            yrs = atm_years(esm, scenario) or atm_years(esm, scenario, "acabf")
            gaps = sorted(set(range(y0, y1 + 1)) - set(yrs))
            # The last year of a series may be absent: the reader persists
            # the last year on disk one year past the end (a CESM2-WACCM tree
            # fetched while the empty 2300 files were withdrawn stops at 2299,
            # discussion #8), so that is a note, not a missing input. Any
            # other gap is an error the reader raises on, so it blocks.
            bridged = gaps == [y1] and yrs and yrs[-1] == y1 - 1
            if not yrs:
                miss.append(f"{esm}/{scenario} atmosphere")
            elif bridged:
                notes.append(
                    f"atmosphere covers {yrs[0]}-{yrs[-1]}; {y1} is absent "
                    f"and the reader persists {y1 - 1} for it"
                )
            elif gaps:
                miss.append(
                    f"atmosphere covers {yrs[0]}-{yrs[-1]} with "
                    f"{len(gaps)} of {y0}-{y1} missing "
                    f"({gaps[0]}..{gaps[-1]})"
                )
            oc = ocean_cover(esm, scenario)
            if oc is None:
                miss.append(f"{esm}/{scenario} ocean tf/so")
            elif oc[0] > y0 or oc[1] < y1 - 1:
                miss.append(f"ocean covers {oc[0]}-{oc[1]}, need {y0}-{y1}")
        if feedback:
            # the gradient each driver reads: the OCX product's in either OCX
            # mode, the ESM's ctrl for the control, else the core's own
            tree = ((OCX_ATMOSPHERE_SOURCE, OCX) if core == 11
                    else (esm, "ctrl") if scenario is None else (esm, scenario))
            problem = gradient_problem(*tree, y0, y1)
            if problem:
                miss.append(problem)

        status = "BLOCKED" if miss else "PARTIAL" if degraded else "READY  "
        shown = miss + degraded + notes
        detail = "" if not shown else "  <- " + "; ".join(shown)
        print(f"  core {core:2d}  {status}  {title}{detail}")


if __name__ == "__main__":
    main()
