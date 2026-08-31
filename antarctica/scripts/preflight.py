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

from icepack2_tools.forcing import (
    atmosphere_path, ocean_path, _oi_climatology_path, _find_ismip7_data,
)

MESH_DIR = os.path.join(_ANT, "mesh")
RESULTS_DIR = os.path.join(_ANT, "results")
DATA_DIR = os.path.join(_ANT, "data")

N_FLOW_DEFAULT = "3.0"


def map_n_tag():
    r"""`_n<N>` filename tag so MAPs at different flow exponents coexist;
    n=4 keeps the legacy untagged name. Mirrors simulation.map_n_tag."""
    n = float(os.environ.get("ISMIP7_N_FLOW", N_FLOW_DEFAULT))
    return "" if abs(n - 4.0) < 1e-9 else f"_n{int(round(n))}"


lc = int(os.environ.get("ISMIP7_LC", "2500"))
lc_coarse = int(os.environ.get("ISMIP7_LC_COARSE", "64000"))
friction = os.environ.get("ISMIP7_FRICTION", "budd")
map_tag = ({"regularized_coulomb": "_rc", "budd": "_budd"}.get(friction, "")
           + map_n_tag())
oi_version = os.environ.get("ISMIP7_OI_VERSION", "30_sep")
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


def shared_missing():
    miss = []
    mesh_fn = os.environ.get(
        "ISMIP7_MESH",
        os.path.join(MESH_DIR, f"antarctica_{lc_coarse}_{lc}.msh"),
    )
    if not os.path.exists(mesh_fn):
        miss.append(f"mesh ({os.path.basename(mesh_fn)})")
    inv = os.path.join(MESH_DIR, f"inversion_icepack2{map_tag}_{lc}.h5")
    if not os.path.exists(inv):
        miss.append(f"MAP ({os.path.basename(inv)})")
    bnd = os.environ.get(
        "ISMIP7_BNDIDS", os.path.join(MESH_DIR, "boundary_ids.json")
    )
    if not os.path.exists(bnd):
        miss.append("boundary_ids.json (untracked — restore or regenerate)")
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
    print(f"Preflight: lc={lc}, friction={friction}, OI={oi_version}")
    base_missing = shared_missing()
    if base_missing:
        print(f"  SHARED inputs missing: {', '.join(base_missing)}")
    else:
        print("  Shared inputs (mesh, MAP, bndids, K, BedMachine, velocity): OK")
    print(f"  RACMO baseline: {'OK' if racmo_ok() else 'MISSING (acabf fallback)'}")
    print()

    for core, title, esm, scenario, y0, y1 in CORES:
        miss = list(base_missing)
        if core == 11:
            if not racmo_ok():
                miss.append("RACMO (OCX SMB)")
            if not oi_ok():
                miss.append(f"OI climatology ({oi_version})")
        elif scenario is None:  # CTRL
            if not oi_ok():
                miss.append(f"OI climatology ({oi_version})")
            clim_years = (atm_years(esm, "historical", "acabf")
                          + atm_years(esm, os.environ.get(
                              "ISMIP7_CLIM_SCENARIO", "ssp585"), "acabf"))
            if not racmo_ok() and not clim_years:
                miss.append(f"{esm} acabf (no RACMO and no ESM climatology)")
        else:
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

        status = "READY  " if not miss else "BLOCKED"
        detail = "" if not miss else "  <- " + "; ".join(miss)
        print(f"  core {core:2d}  {status}  {title}{detail}")


if __name__ == "__main__":
    main()
