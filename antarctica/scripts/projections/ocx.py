#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 11: OCX observationally constrained (1979-2025).

Observation-forced run over the satellite era, for validating the initialized
model against the observed record, and independent of CMIP by design
(discussion #32). ``ISMIP7_OCX_FORCING`` says what it runs on:

``protocol`` (default), the ISMIP7 OCX product:
    SMB:   RACMO2.3p2-ERA, statistically downscaled (SDBN1, 8 km), the full
           ``acabf`` field year by year, 1979-2025.
    Ocean: the expert-judgment thermal forcing and salinity at draft,
           1950-2025, with the melt calibration (one K and a thermal-forcing
           offset per IMBIE basin). ``ISMIP7_OCX_OCEAN``
           picks the scenario: ``main`` (the core one), ``cold``, ``warm`` or
           ``vary``. The Antarctic OCX ocean cites no observational source
           (discussion #41).
    The run refuses to start if the product is not on disk.

``stopgap``, what this core ran on before the product was readable here:
    SMB:   RACMO2.4p1 actual-year fields (1979-2023; end years held at the
           last available RACMO year).
    Ocean: constant OI-climatology TF/so at draft with the melt
           calibration, the climatology it was fitted against.

OPEN, discussion #48 (17 September 2026): the OCX ``main`` thermal forcing
differs strongly from the Zhou climatology around Mertz, halving that
region's melt against a calibration made on the climatology, and may have
been built from an older extrapolated climatology. Every K here is fitted to
the climatology, so run ``check_melt_bound.py --ocx`` and read its per-basin
table before trusting a protocol-forced core 11.

No fracture forcing exists for OCX (discussion #33). The initial state is the
~2015 BedMachine/MAP geometry, so a 1979 start is anachronistic by
construction: treat the early years as relaxation and the 2000s-2025 as the
validation window.

Usage:
    mpiexec -n 24 python scripts/projections/ocx.py
"""
import os, sys

_SCRIPTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from simulation import (setup_model, run_simulation, latest_checkpoint,
                        auto_resume, PETSc)
import math

from icepack2_tools.runconfig import (
    ocx_forcing, ocx_ocean, deltat_per_basin_npz, k_per_basin_npz,
)
from icepack2_tools.forcing import (
    OCX, OCX_ATMOSPHERE_SOURCE,
    ISMIP7Atmosphere, ISMIP7Ocean, make_forcing_callback,
    make_climatology_ocean_callback, load_racmo_smb_climatology,
    load_K_per_basin, forcing_coords, reject_collapse_mask, forcing_year,
    describe_forcing_provenance, describe_observational_forcing,
    describe_melt_calibration,
)

T_START = float(os.environ.get("ISMIP7_T_START", "1979"))
# 1 January of the year AFTER the last one covered, the convention every core
# driver uses: years 1979 through 2025 run and 2025 is the last banked year.
T_END = float(os.environ.get("ISMIP7_T_END", "2026"))
DT = float(os.environ.get("ISMIP7_DT", "0.1"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))
RACMO_LAST = 2023  # smbgl_monthlyS_ANT11_RACMO2.4p1_ERA5_197901_202312


def protocol_forcing(t_start, t_end):
    r"""The OCX atmosphere and ocean readers, once both cover the run.

    Checked before the model is set up, and refused rather than degraded: the
    readers raise on a year they do not hold, and a core that quietly ran on
    something else is what this replaces.
    """
    atm = ISMIP7Atmosphere(esm=OCX_ATMOSPHERE_SOURCE, scenario=OCX)
    ocean = ISMIP7Ocean(scenario=OCX, variant=ocx_ocean())
    first, last = int(math.floor(t_start + 1e-9)), forcing_year(t_end)
    have = atm.available_years("acabf")
    # the reader holds the last year on disk for one year past it, no more
    short = sorted(set(range(first, last)) - set(have))
    if not have or short or last > have[-1] + 1:
        problem = (f"atmosphere acabf covers {have[0]}-{have[-1]} with {len(short)} of "
                   f"{first}-{last} missing" if have else "no atmosphere acabf")
    else:
        problem = None
        for var in ("tf", "so"):
            cover = ocean.coverage(var)
            if cover is None or cover[0] > first or cover[1] + 1 < last:
                problem = (f"ocean {var} covers {cover[0]}-{cover[1]}, need {first}-{last}"
                           if cover else f"no ocean {var} ({ocean.variant})")
                break
    if problem:
        raise FileNotFoundError(
            f"The ISMIP7 OCX product does not cover this run: {problem}. Fetch it with\n"
            f"  python antarctica/scripts/download_mirror.py "
            f"data/OCX/{OCX_ATMOSPHERE_SOURCE}/SDBN1-8000m/acabf/ data/OCX/ocean/{ocean.variant}/\n"
            f"or run the pre-product forcing on purpose with ISMIP7_OCX_FORCING=stopgap."
        )
    return atm, ocean


def main():
    forcing = ocx_forcing()
    experiment_name = "ocx"
    if forcing == "protocol" and ocx_ocean() != "main":
        experiment_name += f"_{ocx_ocean()}"          # a sensitivity member is not core 11
    if os.environ.get("ISMIP7_RUN_TAG"):
        experiment_name += f"_{os.environ['ISMIP7_RUN_TAG']}"
    reject_collapse_mask("the OCX experiment")
    readers = protocol_forcing(T_START, T_END) if forcing == "protocol" else None
    # Explicit restart, else unattended auto-resume from this experiment's own
    # newest checkpoint (ISMIP7_AUTO_RESUME=1), the same lookup the control
    # driver does. A chained batch job depends on it: without it every link
    # cold-starts and the chain never advances.
    restart = os.environ.get("ISMIP7_RESTART")
    if restart is None and auto_resume():
        restart = latest_checkpoint(experiment_name)
        PETSc.Sys.Print(
            f"Auto-resume: {restart}" if restart
            else "Auto-resume: no prior checkpoint"
        )
    dT_npz = deltat_per_basin_npz()
    ctx = setup_model(restart_from=restart)
    # Sample forcing at the geometry dofs, not the mesh vertices: under
    # DG0 geometry those are cell centroids (see forcing.forcing_coords).
    mesh_x, mesh_y = forcing_coords(ctx)

    # The melt calibration: per-basin deltaT at one K, the tracked file
    # unless another is named, else a legacy per-basin K named with
    # ISMIP7_K_PER_BASIN_NPZ, with its optional global scale.
    if dT_npz is not None:
        K_npz = None
        K_field = None
        melt_what = f"per-basin deltaT at one K ({dT_npz})"
    else:
        K_npz = k_per_basin_npz()
        K_field = load_K_per_basin(K_npz, mesh_x, mesh_y, fill=0.0)
        K_scale = float(os.environ.get("ISMIP7_K_SCALE", "1.0"))
        if K_scale != 1.0:
            K_field = K_field * K_scale
            PETSc.Sys.Print(f"  K scaled by ISMIP7_K_SCALE={K_scale:.3f}")
        melt_what = f"legacy per-basin K ({K_npz})"
    if readers is not None:
        atm, ocean = readers
        PETSc.Sys.Print(f"  Atmosphere: ISMIP7 OCX, {OCX_ATMOSPHERE_SOURCE} SDBN1 acabf")
        PETSc.Sys.Print(f"  Ocean melt: ISMIP7 OCX '{ocean.variant}' tf/so + {melt_what}")
        callback = make_forcing_callback(
            atm=atm, ocean=ocean, K_per_basin_npz=K_npz, smb_anomaly=False,
        )
        provenance = describe_forcing_provenance(
            atm, ocean, variables={"atmosphere": ("acabf",)})
    else:
        PETSc.Sys.Print(f"  Ocean melt: OI climatology + {melt_what}")
        PETSc.Sys.Print("  Atmosphere: RACMO2.4p1 actual-year SMB (ISMIP7_OCX_FORCING=stopgap)")
        provenance = describe_observational_forcing(
            smb=f"RACMO2.4p1 actual-year SMB, {RACMO_LAST} held after it",
            ocean=True)
        ocean = None
        oi_melt = make_climatology_ocean_callback(K_field)
        racmo_cache = {}

        def callback(ctx_, t_yr):
            yr = min(forcing_year(t_yr), RACMO_LAST)
            if yr not in racmo_cache:
                racmo_cache[yr] = load_racmo_smb_climatology(
                    ctx_["Q_g"], yr, yr
                ).dat.data_ro.copy()
                while len(racmo_cache) > 3:
                    racmo_cache.pop(next(iter(racmo_cache)))
            ctx_["accum"].dat.data[:] = racmo_cache[yr]
            oi_melt(ctx_, t_yr)

    for line in provenance + describe_melt_calibration(ctx.get("mesh_basename")):
        PETSc.Sys.Print(f"  {line}")

    PETSc.Sys.Print("\nCore Experiment 11: OCX (observationally constrained)")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}, dt={DT}")

    run_simulation(
        ctx,
        experiment_name=experiment_name,
        t_start=T_START,
        t_end=T_END,
        dt=DT,
        output_interval=OUTPUT_INTERVAL,
        forcing_callback=callback,
    )
    if ocean is not None:
        ocean.close()


if __name__ == "__main__":
    main()
