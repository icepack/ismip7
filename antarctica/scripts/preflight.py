#!/usr/bin/env python3
r"""Data preflight for the ISMIP7 core experiments.

Answers "which experiments can run on this machine right now?" in a few
seconds, checking every input each core needs: mesh, MAP inversion,
boundary ids, per-basin K, RACMO, OI climatology, and the (ESM,
scenario) atmosphere/ocean trees over the run period. Honors the same
environment knobs as the runs (ISMIP7_LC, ISMIP7_FRICTION,
ISMIP7_OI_VERSION, ...).

Usage:
    python scripts/preflight.py
    ISMIP7_LC=500 ISMIP7_FRICTION=regularized_coulomb python scripts/preflight.py
"""
import os, sys, glob, re

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_ANT = os.path.dirname(_SCRIPTS)
_PROJECT = os.path.dirname(_ANT)
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from icepack2_tools.forcing import (
    atmosphere_path, ocean_path, _oi_climatology_path, _find_ismip7_data,
)
from icepack2_tools.boundary import sidecar_path
from icepack2_tools.naming import map_basename
from icepack2_tools.climatology import clim_start, clim_end, clim_scenario
from icepack2_tools.runconfig import (
    friction as _friction, geometry_space as _geometry_space, lc as _lc,
    lc_coarse as _lc_coarse,
)
from mesh_naming import get_buffer_m, mesh_filename

MESH_DIR = os.path.join(_ANT, "mesh")
RESULTS_DIR = os.path.join(_ANT, "results")
DATA_DIR = os.path.join(_ANT, "data")

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
    (11, "OCX obs-constrained", None, None, 1990, 2025),
]


def atm_years(esm, scenario, var="acabf-anomaly"):
    d = atmosphere_path(scenario, esm, var)
    if d is None or not os.path.isdir(d):
        return []
    yrs = []
    for fn in os.listdir(d):
        m = re.search(r"_(\d{4})\.nc$", fn)
        if m:
            yrs.append(int(m.group(1)))
    return sorted(yrs)


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


def ocean_cover(esm, scenario):
    d = ocean_path(scenario, esm, "tf")
    if d is None or not os.path.isdir(d):
        return None
    spans = []
    for fn in os.listdir(d):
        m = re.search(r"_(\d{4})-(\d{4})\.nc$", fn)
        if m:
            spans.append((int(m.group(1)), int(m.group(2))))
    if not spans:
        return None
    return min(s[0] for s in spans), max(s[1] for s in spans)


def shared_missing(warn=None):
    r"""Missing shared inputs. Non-fatal caveats are appended to ``warn``."""
    miss = []
    warn = warn if warn is not None else []
    mesh_fn = os.environ.get(
        "ISMIP7_MESH", mesh_filename(lc_coarse, lc, get_buffer_m())
    )
    if not os.path.exists(mesh_fn):
        miss.append(f"mesh ({os.path.basename(mesh_fn)})")
    # The MAP the forward will actually load: the one tagged with this
    # geometry space, else the legacy untagged (CG1) MAP it falls back to with
    # a warning. A legacy MAP runs, but its controls carry the CG1 front bias,
    # so the run is a smoke test rather than a result.
    inv = os.path.join(MESH_DIR, map_basename(friction, lc))
    legacy = os.path.join(MESH_DIR, map_basename(friction, lc, geometry=False))
    if not os.path.exists(inv):
        if os.path.exists(legacy):
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
    if not (glob.glob(os.path.join(RESULTS_DIR, f"calibrated_K_per_basin_{lc}.npz"))
            or glob.glob(os.path.join(RESULTS_DIR, "calibrated_K_per_basin_2500.npz"))):
        miss.append("per-basin K npz")
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
    print(f"Preflight: lc={lc}, friction={friction}, geometry={geom}, "
          f"OI={oi_version}, climatology=historical+{CLIM_SCENARIO} "
          f"{CLIM_START}-{CLIM_END}")
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
        if core == 11:
            if not racmo_ok():
                miss.append("RACMO (OCX SMB)")
            if not oi_ok():
                miss.append(f"OI climatology ({oi_version})")
        elif scenario is None:  # CTRL
            if not oi_ok():
                miss.append(f"OI climatology ({oi_version})")
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
            if not yrs:
                miss.append(f"{esm}/{scenario} atmosphere")
            elif gaps:
                # get_field silently returns zeros for a missing year, so
                # interior gaps corrupt a run just like missing endpoints
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

        status = "BLOCKED" if miss else "PARTIAL" if degraded else "READY  "
        notes = miss + degraded
        detail = "" if not notes else "  <- " + "; ".join(notes)
        print(f"  core {core:2d}  {status}  {title}{detail}")


if __name__ == "__main__":
    main()
