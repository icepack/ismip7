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
ssp126, per protocol) over ISMIP7_CLIM_START..END (default 2000-2029), the
SAME pool for every experiment, so all cores share one baseline and the
historical -> projection handoff at 2014/2015 is seamless. Without RACMO the run
falls back to the full acabf(t) field; with no acabf data at all it
refuses to run (ISMIP7_ALLOW_ZERO_SMB=1 to override).

Fracture / shelf-collapse masks are loaded when present, and
make_forcing_callback publishes the year's mask as ctx["collapse"].
Whether the run acts on it is ISMIP7_FRACTURE: under `mask` the transport
empties every FLOATING cell the mask flags and books it as calving
(protocol path C), under `mask_front` only those open water has reached
(discussion #30); grounded ice is never touched, and under the default
`none` the mask is loaded but unused. No mask exists for historical or OCX.
The stress-gated variant (Lai et al. 2020) is not implemented.
"""

import math, os, sys
import numpy as np

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from simulation import (setup_model, run_simulation, latest_checkpoint,
                        historical_endpoint,
                        auto_resume, PETSc)
from icepack2_tools.forcing import (
    ISMIP7Atmosphere, ISMIP7Ocean, ISMIP7Fracture,
    make_forcing_callback, load_racmo_smb_climatology, forcing_coords,
    describe_forcing_provenance, describe_melt_calibration, forcing_year,
)
from icepack2_tools.climatology import (
    clim_start, clim_end, clim_scenario, clim_pool_missing, describe_clim_pool,
)
from icepack2_tools.runconfig import (
    geometry_backdate_years,
    FRACTURE_MASK_MODES, fracture as fracture_mode, k_per_basin_npz,
    deltat_per_basin_npz,
)

# Owned by icepack2_tools.climatology: this pool must match the CONTROL's
# climatology, or the projections are re-referenced against a different
# baseline than the control they are differenced from.
CLIM_START = clim_start()
CLIM_END = clim_end()
CLIM_SCENARIO = clim_scenario()


def smb_scheme(ctx, esm):
    r"""(smb_anomaly, smb_baseline, description) for make_forcing_callback.

    RACMO climatology + aSMB re-referenced over the shared pool; falls back
    to the full acabf(t) field when RACMO (or the pool) is unavailable.
    """
    # Sample forcing at the geometry dofs, not the mesh vertices: under
    # DG0 geometry those are cell centroids (see forcing.forcing_coords).
    mesh_x, mesh_y = forcing_coords(ctx)
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
        racmo = load_racmo_smb_climatology(ctx["Q_g"], CLIM_START, CLIM_END)
        ref = np.zeros(len(mesh_x))
        for a, y in pool:
            ref += a.get_smb(y, mesh_x, mesh_y, anomaly=True)
        ref /= len(pool)
        years = sorted(y for _, y in pool)
        PETSc.Sys.Print(f"  {describe_clim_pool(years, 'acabf-anomaly')}")
        missing = clim_pool_missing(years)
        if missing:
            PETSc.Sys.Print(
                f"  WARNING: aSMB re-referenced over {len(set(years))} of the "
                f"{CLIM_END - CLIM_START + 1} window years ({len(missing)} "
                f"missing, {missing[0]}..{missing[-1]}). This core's baseline "
                f"differs from a full-window sibling's while both are "
                f"differenced against the same CTRL. Proceeding: the pool is "
                f"recorded above and in the per-core report."
            )
        return True, racmo.dat.data_ro - ref, (
            f"RACMO2.4p1 baseline + aSMB re-referenced to "
            f"{years[0]}-{years[-1]} ({len(pool)} yr pooled)"
        )
    except FileNotFoundError as e:
        return False, None, f"full acabf(t), no RACMO baseline ({e})"


def run_core_experiment(*, core, title, name, esm, scenario,
                        t_start_default, t_end_default,
                        restart_from_hist=True):
    r"""Run one ESM-forced ISMIP7 core experiment end to end."""
    t_start = float(os.environ.get("ISMIP7_T_START", str(t_start_default)))
    t_end = float(os.environ.get("ISMIP7_T_END", str(t_end_default)))
    dt = float(os.environ.get("ISMIP7_DT", "0.1"))
    output_interval = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

    esm_tag = esm.lower().replace("-", "_")
    # Optional run tag (ISMIP7_RUN_TAG / --tag) suffixed on the experiment name
    # so parallel method lines (e.g. n=3 vs n=4) write distinct output files
    # instead of clobbering each other. It also selects the tagged historical
    # to restart from, keeping the hist->projection chain within one line.
    tag = os.environ.get("ISMIP7_RUN_TAG", "")
    tag_sfx = f"_{tag}" if tag else ""
    experiment_name = f"{name}{tag_sfx}"
    restart = os.environ.get("ISMIP7_RESTART")
    # Unattended auto-resume (ISMIP7_AUTO_RESUME=1), the same lookup the
    # control driver does: with no explicit restart, continue from this
    # experiment's own newest checkpoint. A chained batch job depends on it,
    # and it takes precedence over the historical endpoint below, which is
    # only where the FIRST link of a projection starts.
    if restart is None and auto_resume():
        restart = latest_checkpoint(experiment_name)
        PETSc.Sys.Print(
            f"Auto-resume: {restart}" if restart
            else "Auto-resume: no prior checkpoint"
        )
    if restart is None and restart_from_hist:
        restart = historical_endpoint(esm_tag, tag_sfx, t_start)
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
    # The reader serves the years its chunk files hold and the single year
    # after them, and raises on anything else, so ask before the model setup
    # rather than find out at the first step, or at the last one.
    first, last = int(math.floor(t_start + 1e-9)), forcing_year(t_end)
    cover = ISMIP7Ocean(esm=esm, scenario=scenario).require_years(first, last)
    PETSc.Sys.Print(f"  Ocean forcing: tf, so cover {cover[0]}-{cover[1]}")

    if restart:
        PETSc.Sys.Print(f"  Restart: {restart}")
    else:
        PETSc.Sys.Print("  Cold start from BedMachine/inversion initial state")
    dT_npz = deltat_per_basin_npz()

    # A cold start before the 2015 geometry starts from it with the observed
    # thinning undone (issue #117); a restart carries its own geometry.
    ctx = setup_model(
        restart_from=restart,
        backdate_years=0.0 if restart else geometry_backdate_years(t_start))

    smb_anomaly, smb_baseline = False, None
    if atm is not None:
        smb_anomaly, smb_baseline, desc = smb_scheme(ctx, esm)
        PETSc.Sys.Print(f"  SMB: {desc}")

    # ── Ocean + fracture ──
    ocean = ISMIP7Ocean(esm=esm, scenario=scenario)
    fracture = ISMIP7Fracture(esm=esm, scenario=scenario)
    try:
        fracture.load()
    except Exception as e:
        # Under the default `none` the mask is never read, so an unreadable
        # fracture tree must not abort a run that does not want it; under a
        # mask mode the run asked for exactly this file, so it sees the failure.
        if fracture_mode() in FRACTURE_MASK_MODES:
            raise
        PETSc.Sys.Print(f"  Fracture tree not readable, ignored: {e}")
    if fracture_mode() in FRACTURE_MASK_MODES and not fracture.has_collapse_mask():
        raise FileNotFoundError(
            f"ISMIP7_FRACTURE={fracture_mode()} but no ice-shelf collapse mask was found "
            f"for {esm}/{scenario} under {fracture.fracture_dir()}. The masks "
            f"exist for the SSP scenarios only, so this is a configuration "
            f"error: download the fracture tree, or run with "
            f"ISMIP7_FRACTURE=none."
        )
    if fracture_mode() in FRACTURE_MASK_MODES:
        # the highest version on disk can be one replaced on Globus before the
        # mirror caught up; refuse it rather than run on it
        fracture.check_min_version()

    # The melt calibration: the tracked one or a named offsets file, else a
    # legacy per-basin K named with ISMIP7_K_PER_BASIN_NPZ.
    K_npz = None if dT_npz is not None else k_per_basin_npz()
    # What this run opens, for the committed report: a collapse mask only
    # counts when the run reads it, and the melt calibration by its sha256.
    for line in (describe_forcing_provenance(
            atm, ocean, fracture if fracture_mode() != "none" else None)
            + describe_melt_calibration(dT_npz, K_npz, ctx.get("mesh_basename"))):
        PETSc.Sys.Print(f"  {line}")

    if dT_npz is not None:
        PETSc.Sys.Print(f"  Ocean melt: per-basin deltaT at one K from {dT_npz}")
    else:
        PETSc.Sys.Print(f"  Ocean melt: legacy per-basin K from {K_npz}")

    callback = make_forcing_callback(
        atm=atm, ocean=ocean, fracture=fracture, K_per_basin_npz=K_npz,
        smb_anomaly=smb_anomaly, smb_baseline=smb_baseline,
    )

    PETSc.Sys.Print(f"\nCore Experiment {core}: {title}")
    PETSc.Sys.Print(f"  Period: {t_start}-{t_end}, dt={dt}")

    results = run_simulation(
        ctx,
        experiment_name=experiment_name,
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
