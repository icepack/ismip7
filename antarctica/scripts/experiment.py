#!/usr/bin/env python3
r"""Shared driver for the ISMIP7 core experiments.

Every ESM-forced core experiment (historical 1/2, projections 3-8) is the
same run with different (esm, scenario, period), so the drivers in
scripts/historical/ and scripts/projections/ are thin shims over
:func:`run_core_experiment`. The forcing scheme is the one the CTRL uses,
so projection minus control is a clean forced signal:

    SMB(t)  = RACMO_clim + [aSMB(t) - mean(aSMB | reference window)]
    ocean   = ISMIP7 bias-corrected tf/so at draft (Burgard quadratic
              mixed-slope, calibrated per-basin K)

The aSMB re-reference pool is historical + ISMIP7_CLIM_SCENARIO (default
ssp585) over ISMIP7_CLIM_START..END (default 2000-2029) — the SAME pool
for every experiment, so all cores share one baseline and the historical
-> projection handoff at 2014/2015 is seamless. Without RACMO the run
falls back to the full acabf(t) field; with no acabf data at all it
refuses to run (ISMIP7_ALLOW_ZERO_SMB=1 to override).

Fracture / shelf-collapse masks are loaded when present but NOT yet
applied by make_forcing_callback — wiring the collapse mask into the
thickness update is an open protocol item.
"""

import os, sys
import numpy as np

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from simulation import setup_model, run_simulation, RESULTS_DIR, PETSc, lc
from icepack2_tools.forcing import (
    ISMIP7Atmosphere, ISMIP7Ocean, ISMIP7Fracture,
    make_forcing_callback, load_racmo_smb_climatology,
)

CLIM_START = int(os.environ.get("ISMIP7_CLIM_START", "2000"))
CLIM_END = int(os.environ.get("ISMIP7_CLIM_END", "2029"))
CLIM_SCENARIO = os.environ.get("ISMIP7_CLIM_SCENARIO", "ssp585")


def find_k_npz():
    r"""Calibrated per-basin K npz: this mesh's calibration, else the 2500 m
    one (16 basin scalars remapped through the IMBIE2 8 km grid —
    mesh-independent), else None (scalar ISMIP7_K_MELT)."""
    override = os.environ.get("ISMIP7_K_PER_BASIN_NPZ")
    candidates = [override] if override else [
        os.path.join(RESULTS_DIR, f"calibrated_K_per_basin_{lc}.npz"),
        os.path.join(RESULTS_DIR, "calibrated_K_per_basin_2500.npz"),
    ]
    for c in candidates:
        if c and os.path.exists(c):
            return c
    return None


def smb_scheme(ctx, esm):
    r"""(smb_anomaly, smb_baseline, description) for make_forcing_callback.

    RACMO climatology + aSMB re-referenced over the shared pool; falls back
    to the full acabf(t) field when RACMO (or the pool) is unavailable.
    """
    mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
    mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]
    ref_atms = [
        ISMIP7Atmosphere(esm=esm, scenario="historical"),
        ISMIP7Atmosphere(esm=esm, scenario=CLIM_SCENARIO),
    ]
    try:
        pool = []
        for a in ref_atms:
            pool += [(a, y) for y in a.available_years("acabf-anomaly")
                     if CLIM_START <= y <= CLIM_END]
        if not pool:
            raise FileNotFoundError(
                f"no acabf-anomaly years in {CLIM_START}-{CLIM_END}"
            )
        racmo = load_racmo_smb_climatology(ctx["Q"], CLIM_START, CLIM_END)
        ref = np.zeros(len(mesh_x))
        for a, y in pool:
            ref += a.get_smb(y, mesh_x, mesh_y, anomaly=True)
        ref /= len(pool)
        years = sorted(y for _, y in pool)
        return True, racmo.dat.data_ro - ref, (
            f"RACMO2.4p1 baseline + aSMB re-referenced to "
            f"{years[0]}-{years[-1]} ({len(pool)} yr pooled)"
        )
    except FileNotFoundError as e:
        return False, None, f"full acabf(t) — no RACMO baseline ({e})"


def run_core_experiment(*, core, title, name, esm, scenario,
                        t_start_default, t_end_default,
                        restart_from_hist=True):
    r"""Run one ESM-forced ISMIP7 core experiment end to end."""
    t_start = float(os.environ.get("ISMIP7_T_START", str(t_start_default)))
    t_end = float(os.environ.get("ISMIP7_T_END", str(t_end_default)))
    dt = float(os.environ.get("ISMIP7_DT", "0.1"))
    output_interval = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

    esm_tag = esm.lower().replace("-", "_")
    restart = os.environ.get("ISMIP7_RESTART")
    if restart is None and restart_from_hist:
        cand = os.path.join(RESULTS_DIR, f"hist_{esm_tag}_{lc}_final.h5")
        restart = cand if os.path.exists(cand) else None
    if restart and not os.path.exists(restart):
        raise FileNotFoundError(f"ISMIP7_RESTART not found: {restart}")

    # ── Forcing availability BEFORE the expensive model setup ──
    atm = ISMIP7Atmosphere(esm=esm, scenario=scenario)
    yrs = atm.available_years("acabf-anomaly") or atm.available_years("acabf")
    if not yrs:
        if os.environ.get("ISMIP7_ALLOW_ZERO_SMB"):
            PETSc.Sys.Print("  WARNING: no atmosphere data, zero SMB")
            atm = None
        else:
            raise FileNotFoundError(
                f"No acabf data for {esm}/{scenario}. Download the "
                f"atmosphere tree first (ISMIP7_ALLOW_ZERO_SMB=1 to force "
                f"a zero-SMB run)."
            )
    else:
        PETSc.Sys.Print(
            f"  Atmosphere forcing: {len(yrs)} years ({yrs[0]}-{yrs[-1]})"
        )
    if not (ISMIP7Ocean(esm=esm, scenario=scenario)
            ._year_field("tf", t_start, nan_fill=0.0)):
        raise FileNotFoundError(
            f"No ocean tf data for {esm}/{scenario} covering {t_start:.0f}. "
            f"Download the ocean tree first."
        )

    if restart:
        PETSc.Sys.Print(f"  Restart: {restart}")
    else:
        PETSc.Sys.Print("  Cold start from BedMachine/inversion initial state")

    ctx = setup_model(restart_from=restart)

    smb_anomaly, smb_baseline = False, None
    if atm is not None:
        smb_anomaly, smb_baseline, desc = smb_scheme(ctx, esm)
        PETSc.Sys.Print(f"  SMB: {desc}")

    # ── Ocean + fracture ──
    ocean = ISMIP7Ocean(esm=esm, scenario=scenario)
    fracture = ISMIP7Fracture(esm=esm, scenario=scenario)
    try:
        fracture.load()
    except Exception:
        fracture = None

    K_npz = find_k_npz()
    K_melt = float(os.environ.get("ISMIP7_K_MELT", "1.15e-4"))
    if K_npz is not None:
        PETSc.Sys.Print(f"  Ocean melt: calibrated per-basin K from {K_npz}")
    else:
        PETSc.Sys.Print(
            f"  WARNING: no per-basin K calibration; scalar K={K_melt:.2e}"
        )

    callback = make_forcing_callback(
        atm=atm, ocean=ocean, fracture=fracture,
        K=K_melt, K_per_basin_npz=K_npz,
        smb_anomaly=smb_anomaly, smb_baseline=smb_baseline,
    )

    PETSc.Sys.Print(f"\nCore Experiment {core}: {title}")
    PETSc.Sys.Print(f"  Period: {t_start}-{t_end}, dt={dt}")

    results = run_simulation(
        ctx,
        experiment_name=name,
        t_start=t_start,
        t_end=t_end,
        dt=dt,
        output_interval=output_interval,
        forcing_callback=callback,
    )

    if ocean is not None:
        ocean.close()
    if fracture is not None:
        fracture.close()
    return results
