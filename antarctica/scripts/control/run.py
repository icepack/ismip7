#!/usr/bin/env python3
r"""ISMIP7 Core Experiments 9/10: CTRL2015 -- constant 2015 climate.

Atmosphere: RACMO2.4p1 2000-2029 SMB climatology held fixed (falls back
            to the pooled ISMIP7 acabf climatology if RACMO is absent).
Ocean:      the ESM's own `ctrl` tree (tf + so, v3), the organisers'
            2000-2029 mean of historical and ssp126, identical every year
            (forum threads 15 and 28, icepack/ismip7#107). It is read and
            melted as the projections read and melt theirs
            (`experiment.py`): the Burgard quadratic-mixed-slope melt,
            recomputed each step from the evolving geometry, with the melt
            calibration: one K and a thermal-forcing offset per IMBIE basin,
            fitted to the Zhou `30_sep` climatology (the tracked file,
            `runconfig.MELT_CALIBRATION_DEFAULT`, unless another is named).

Usage:
    mpiexec -n 12 python scripts/control/run.py
    ISMIP7_ESM=MRI-ESM2-0 mpiexec -n 12 python scripts/control/run.py
"""

import os, sys, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)
))))
sys.path.insert(0, _PROJECT)

from firedrake import assemble, dx, Constant
from simulation import (setup_model, run_simulation, latest_checkpoint,
                        auto_resume, RESULTS_DIR, PETSc, lc)
from icepack2_tools.forcing import (
    ISMIP7Atmosphere,
    ISMIP7Ocean,
    describe_forcing_provenance,
    describe_observational_forcing,
    load_racmo_smb_climatology,
    make_forcing_callback,
    reject_collapse_mask,
    compute_sin_alpha,
    quadratic_mixed_slope,
    forcing_coords,
    forcing_year,
    melt_receiving,
    describe_melt_calibration,
    _K_DEFAULT,
)
from icepack2_tools.runconfig import deltat_per_basin_npz, k_per_basin_npz
from icepack2_tools.climatology import (
    clim_start, clim_end, clim_scenario, clim_pool_missing, describe_clim_pool,
)

T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2301"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

ESM = os.environ.get("ISMIP7_ESM", "CESM2-WACCM")
# Reference-climate window and pool scenario for the constant SMB, owned by
# icepack2_tools.climatology so the control, the projections, the preflight
# and the report cannot drift apart. RACMO2.4p1 (1979-2023) is the primary
# baseline, so 2000-2029 yields its 2000-2023 mean. The acabf fallback pools
# historical + the protocol scenario (ssp126) and uses whatever subset of the
# window exists locally; the acabf-anomaly files are referenced to 1960-1989,
# so an anomaly-based baseline is NOT constructible. CLIM_SCENARIO is only
# reached when RACMO is unavailable, but the fallback should still be the
# protocol pool.
CLIM_START = clim_start()
CLIM_END = clim_end()
CLIM_SCENARIO = clim_scenario()

# The organisers' per-ESM control ocean: ISMIP7/AIS/<ESM>/ctrl/ocean/{tf,so}/.
CTRL_SCENARIO = "ctrl"


def area_weighted_mean(field, mesh):
    r"""Domain-area-weighted mean of a field (mesh refinement varies spatially)."""
    return assemble(field * dx) / assemble(Constant(1.0) * dx(domain=mesh))


def compute_climatology(atms, mesh_x, mesh_y):
    r"""Mean full-field (acabf) SMB over [CLIM_START, CLIM_END], pooling
    years across the given atmospheres (historical + projection scenario).

    Full field only: a mean of `acabf-anomaly` is an anomaly wrt the ESM's
    1960-1989 climatology, not a baseline, so no anomaly fallback exists.
    """
    smb_sum = np.zeros(len(mesh_x))
    n = 0
    pooled = []
    for atm in atms:
        years = [y for y in atm.available_years("acabf")
                 if CLIM_START <= y <= CLIM_END]
        for yr in years:
            smb_sum += atm.get_smb(yr, mesh_x, mesh_y, anomaly=False)
        n += len(years)
        pooled += years
        if years:
            PETSc.Sys.Print(
                f"  {atm.scenario}: {len(years)} acabf years "
                f"({years[0]}-{years[-1]}) in climatology window"
            )
    if n == 0:
        return None
    PETSc.Sys.Print(f"  {describe_clim_pool(pooled, 'acabf full field')}")
    missing = clim_pool_missing(pooled)
    if missing:
        PETSc.Sys.Print(
            f"  WARNING: the constant SMB climatology is a mean over "
            f"{len(set(pooled))} of the {CLIM_END - CLIM_START + 1} window "
            f"years ({len(missing)} missing, {missing[0]}..{missing[-1]}). "
            f"This control's baseline is NOT the full-window one, so a "
            f"projection re-referenced over the full window is differenced "
            f"against a different baseline. Proceeding: the pool is recorded "
            f"above and in the per-core report."
        )
    return smb_sum / n


def make_synthetic_ocean_callback(tf_max=1.5, depth_ref=1000.0, K=_K_DEFAULT):
    r"""STOPGAP ocean-melt callback that needs no ocean data.

    Prescribes a depth-dependent thermal forcing
    TF(draft) = clip(tf_max * |draft| / depth_ref, 0, tf_max) and runs it
    through the Burgard quadratic_mixed_slope melt with constant salinity and a
    scalar K. Physically structured (deeper ice melts more) but NOT calibrated --
    a development substitute for the `ctrl` ocean + per-basin K path while the
    ocean data is unavailable. Enable with ISMIP7_SYNTHETIC_MELT=1.
    """
    PETSc.Sys.Print(
        f"  SYNTHETIC ocean melt (uncalibrated stopgap): "
        f"TF=clip({tf_max:.2f}K * |draft|/{depth_ref:.0f}m), "
        f"K={K:.2e}, salinity=34.5 PSU"
    )

    def callback(ctx, t_yr):
        h = ctx["h"].dat.data_ro
        b = ctx["b"].dat.data_ro
        s = ctx["s"].dat.data_ro
        draft = np.minimum(s - h, 0.0)
        tf = np.clip(tf_max * np.abs(draft) / depth_ref, 0.0, tf_max)
        sal = np.full_like(tf, 34.5)
        sin_a = compute_sin_alpha(ctx)
        melt = quadratic_mixed_slope(tf, sal, sin_a, K=K)
        ctx["ocean_melt"].dat.data[:] = np.where(melt_receiving(s, b, h), melt, 0.0)

    return callback


def main():
    parser = argparse.ArgumentParser(description="ISMIP7 CTRL2015 control run")
    parser.add_argument("--snes-monitor", action="store_true",
                        help="enable PETSc SNES/KSP convergence monitoring")
    parser.add_argument("--snes-log", default=None,
                        help="route SNES monitor output to this file (default: stdout)")
    parser.add_argument("--restart", default=os.environ.get("ISMIP7_RESTART"),
                        help="restart from a checkpoint .h5 (default: branch from "
                             "the historical endpoint, else cold start)")
    parser.add_argument("--tag", default=os.environ.get("ISMIP7_RUN_TAG", ""),
                        help="suffix on experiment_name so output files are distinct")
    parser.add_argument("--checkpoint-interval", type=int,
                        default=int(os.environ.get("ISMIP7_CHECKPOINT_INTERVAL", "100")),
                        help="save a thickness+velocity checkpoint every N steps")
    args, _ = parser.parse_known_args()
    if args.snes_monitor or args.snes_log:
        os.environ["ISMIP7_SNES_MONITOR"] = "1"
    if args.snes_log:
        os.environ["ISMIP7_SNES_LOG"] = args.snes_log

    esm_tag = ESM.lower().replace("-", "_")
    # --tag lands in experiment_name BEFORE the auto-resume lookup so a tagged
    # run resumes its own checkpoints, not the untagged experiment's.
    experiment_name = f"ctrl2015_{esm_tag}" + (f"_{args.tag}" if args.tag else "")

    # Unattended auto-resume (ISMIP7_AUTO_RESUME=1): with no explicit restart,
    # continue from the newest self-contained checkpoint for this experiment.
    # A rebooted long run picks up where it left off; the mesh + frozen anchors
    # + timeline year all come from that checkpoint.
    restart_from = args.restart
    if restart_from is None and auto_resume():
        restart_from = latest_checkpoint(experiment_name)
        PETSc.Sys.Print(
            f"Auto-resume: {restart_from}" if restart_from
            else "Auto-resume: no prior checkpoint"
        )
    # ISMIP6/ISMIP7 ctrl_proj convention: the control and the scenario
    # projections MUST branch from the SAME initial state at the SAME time, so
    # their shared spin-up/relaxation drift (and the identical frozen a_ref)
    # cancels in the projection-minus-control difference and leaves only the
    # forced response. experiment.py restarts every projection from
    # hist_<esm>_<lc>_final.h5; the CTRL therefore does the same instead of
    # cold-starting from the 2015 inversion. Cold-starting from the pristine
    # inversion gives a DIFFERENT initial geometry than the projections, so the
    # projection's own historical-endpoint relaxation does NOT cancel and the
    # forced signal is mis-estimated. NB this is a correctness fix, not a
    # magnitude fix: at 32 km the hist-branched control drifts +37 mm SLE over
    # 2015-2100 vs the cold-start control's +8 mm, so the (correct) same-state
    # difference is LARGER, not smaller (ssp126 +161 vs +132 mm at n=4). a_ref is
    # a t=0 BALANCING correction (zeroes the initial tendency), not a net sink,
    # so it does not add a standalone SLE trend; an earlier note claiming a
    # ~160 mm spurious a_ref sink was wrong. The ISMIP6 overshoot is a real
    # forced-response bias, not a differencing artifact. Falls back to a cold
    # start (mis-matched control) only when the historical endpoint is absent.
    if restart_from is None:
        tag_sfx = f"_{args.tag}" if args.tag else ""
        hist = os.path.join(RESULTS_DIR, f"hist_{esm_tag}{tag_sfx}_{lc}_final.h5")
        if os.path.exists(hist):
            restart_from = hist
            PETSc.Sys.Print(
                f"  Branching CTRL from the historical endpoint (same initial "
                f"state as the projections; shared drift cancels in proj-CTRL): "
                f"{os.path.basename(hist)}"
            )
        else:
            PETSc.Sys.Print(
                f"  WARNING: no historical endpoint {os.path.basename(hist)}; "
                f"cold-starting from the inversion. The CTRL will start from a "
                f"DIFFERENT geometry than the hist-branched projections, so "
                f"projection-minus-CTRL will not cleanly isolate the forced "
                f"response. Run the historical first, or pass --restart."
            )
    if restart_from and not os.path.exists(restart_from):
        raise FileNotFoundError(f"CTRL restart not found: {restart_from}")
    # The melt calibration: per-basin deltaT at one K, the tracked file
    # unless another is named, else a legacy per-basin K named with
    # ISMIP7_K_PER_BASIN_NPZ. Both are checked here, before the setup.
    dT_npz = deltat_per_basin_npz()
    K_npz = None if dT_npz is not None else k_per_basin_npz()
    reject_collapse_mask("the control experiment")

    # The ocean is read year by year as in a projection, and a tree with no
    # files reads as zero thermal forcing, so the inputs are checked here,
    # before the expensive setup.
    ocean = None
    if not os.environ.get("ISMIP7_SYNTHETIC_MELT"):
        ocean = ISMIP7Ocean(esm=ESM, scenario=CTRL_SCENARIO)
        cover = ocean.require_years(int(T_START), forcing_year(T_END))
        PETSc.Sys.Print(f"  Ocean forcing: {ESM}/{CTRL_SCENARIO} tf, so cover "
                        f"{cover[0]}-{cover[1]}")

    ctx = setup_model(restart_from=restart_from)

    # Forcing fields live on the GEOMETRY space (DG0 by default), whose
    # dofs are cell centroids, not mesh vertices. Sampling climatologies
    # at vertices would write a wrong-length array into ctx["accum"].
    mesh_x, mesh_y = forcing_coords(ctx)

    # Atmosphere: RACMO climatology baseline; fall back to the pooled ISMIP7
    # acabf full-field climatology; refuse to run zero-SMB unless forced.
    try:
        ctx["accum"].assign(
            load_racmo_smb_climatology(ctx["Q_g"], CLIM_START, CLIM_END)
        )
        mean_smb = area_weighted_mean(ctx["accum"], ctx["mesh"])
        PETSc.Sys.Print(
            f"  Climatological SMB from RACMO2.4p1 ({CLIM_START}-{CLIM_END}): "
            f"area-weighted mean={mean_smb:.4f} m/yr"
        )
        for line in describe_observational_forcing(
                smb=f"RACMO2.4p1 SMB climatology {CLIM_START}-{CLIM_END}"):
            PETSc.Sys.Print(f"  {line}")
    except FileNotFoundError:
        PETSc.Sys.Print("  No RACMO data; falling back to ISMIP7 acabf climatology")
        atms = [
            ISMIP7Atmosphere(esm=ESM, scenario="historical"),
            ISMIP7Atmosphere(esm=ESM, scenario=CLIM_SCENARIO),
        ]
        clim_smb = compute_climatology(atms, mesh_x, mesh_y)
        # the pool is a mean of the FULL field, so that is all it opens
        for line in describe_forcing_provenance(
                *atms, variables={"atmosphere": ("acabf",)}):
            PETSc.Sys.Print(f"  {line}")
        if clim_smb is None:
            if os.environ.get("ISMIP7_ALLOW_ZERO_SMB"):
                PETSc.Sys.Print("  WARNING: no climatology data, using zero SMB")
                ctx["accum"].assign(0.0)
            else:
                raise FileNotFoundError(
                    f"No RACMO SMB and no acabf data for {ESM} in "
                    f"{CLIM_START}-{CLIM_END} (historical or {CLIM_SCENARIO}). "
                    f"A zero-SMB control is almost certainly not what you "
                    f"want; set ISMIP7_ALLOW_ZERO_SMB=1 to force it."
                )
        else:
            ctx["accum"].dat.data[:] = clim_smb
            mean_smb = area_weighted_mean(ctx["accum"], ctx["mesh"])
            PETSc.Sys.Print(
                f"  Climatological SMB: area-weighted mean={mean_smb:.4f} m/yr"
            )

    # Ocean melt: synthetic stopgap (no data) or the ESM's ctrl ocean.
    if ocean is None:
        tf_max = float(os.environ.get("ISMIP7_SYNTH_TF_MAX", "1.5"))
        depth_ref = float(os.environ.get("ISMIP7_SYNTH_DEPTH_REF", "1000.0"))
        callback = make_synthetic_ocean_callback(tf_max, depth_ref)
        melt_desc = "Synthetic ocean melt stopgap (ISMIP7_SYNTHETIC_MELT)"
    else:
        for line in (describe_forcing_provenance(ocean)
                     + describe_melt_calibration(ctx.get("mesh_basename"))):
            PETSc.Sys.Print(f"  {line}")
        # The projections' callback (experiment.py) with no atmosphere, so the
        # SMB assigned above stays: the offsets file's TF shift at its one K,
        # else the legacy per-basin K times ISMIP7_K_SCALE.
        callback = make_forcing_callback(ocean=ocean, K_per_basin_npz=K_npz)
        source = (f"per-basin deltaT from {dT_npz}" if dT_npz is not None
                  else f"legacy per-basin K from {K_npz}")
        melt_desc = f"Constant {ESM} {CTRL_SCENARIO} ocean + {source}"

    PETSc.Sys.Print(f"\nControl experiment: {ESM}")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}")
    PETSc.Sys.Print(f"  Constant {CLIM_START}-{CLIM_END} SMB climatology")
    PETSc.Sys.Print(f"  {melt_desc}")

    run_simulation(
        ctx,
        experiment_name=experiment_name,
        t_start=T_START,
        t_end=T_END,
        dt=DT,
        output_interval=OUTPUT_INTERVAL,
        checkpoint_interval=args.checkpoint_interval,
        forcing_callback=callback,
    )
    if ocean is not None:
        ocean.close()


if __name__ == "__main__":
    main()
