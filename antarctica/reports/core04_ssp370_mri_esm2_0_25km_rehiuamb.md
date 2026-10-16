# Core 4: ssp370_mri_esm2_0_rehiuamb (25 km)

- date: 2026-10-09
- git: 288382a
- log: `../logs/ismip7_fwd_11890013.out`
- timeseries: `results/ssp370_mri_esm2_0_rehiuamb_25000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly MRI-ESM2-0 ssp370 GEMB-SDBN1-8000m v2
- Forcing provenance: atmosphere acabf MRI-ESM2-0 ssp370 GEMB-SDBN1-8000m v2
- Forcing provenance: ocean tf MRI-ESM2-0 ssp370 ocean v3
- Forcing provenance: ocean so MRI-ESM2-0 ssp370 ocean v3
- Forcing provenance: ocean melt calibration deltaT_per_basin_25000_K6.500e-05.npz sha256 385c1f479ec5dfa97414597fe1dba5111f72410797f53d38ff08a306945c97a4 (named with ISMIP7_DELTAT_PER_BASIN_NPZ): K 6.500e-05 (K50) with a thermal-forcing offset per basin, fitted on antarctica_250000_25000_buffered20000; this run's mesh is antarctica_250000_25000_buffered20000
- Forcing provenance: atmosphere dacabfdz MRI-ESM2-0 ssp370 GEMB-SDBN1-8000m v2
- SMB-elevation feedback: dacabfdz from MRI-ESM2-0 ssp370 GEMB-SDBN1-8000m v2, surface change from the chain's initial state (H_init)
- Free-surface stabilization: theta 1 on the lagged thickness-velocity coupling (ISMIP7_FSSA_THETA)
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)
- ISMIP6 ensemble: not judged (exit 2: missing data or bad input) (pool: exp01,exp02,exp03,exp04,exp05)

25 km rehearsal evidence for IU's Budd line (the issue 138 recipe, from IU's non-relaxed vf MAP), not a submission result.

## Run environment

```
ISMIP7_APPARENT_MB=1
ISMIP7_AUTO_RESUME=1
ISMIP7_BNDIDS=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km_iu/mesh/boundary_ids_antarctica_250000_25000_buffered20000.json
ISMIP7_BUFFER_M=20000
ISMIP7_CHECKPOINT_EVERY_YR=5
ISMIP7_CLIM_END=2029    # default (not exported)
ISMIP7_CLIM_SCENARIO=ssp126    # default (not exported)
ISMIP7_CLIM_START=2000    # default (not exported)
ISMIP7_CONTINUATION_STEPS=8    # default (not exported)
ISMIP7_DATA_ROOT=/N/project/ice_rheology/ISMIP7/ISMIP7/AIS
ISMIP7_DELTAT_PER_BASIN_NPZ=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km_iu/calibration/deltaT_per_basin_25000_K6.500e-05.npz
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_mumps
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER_CANONICAL=scpc_mumps    # resolved canonical mode
ISMIP7_DT=0.025
ISMIP7_ESM=MRI-ESM2-0
ISMIP7_EXPERIMENT=ssp370_mri_esm2
ISMIP7_FIXED_FRONT=1
ISMIP7_FRACTURE=none
ISMIP7_FRICTION=budd
ISMIP7_FSSA_REFERENCE=auto    # default (not exported)
ISMIP7_FSSA_THETA=1    # default (not exported)
ISMIP7_GEOMETRY_SPACE=dg0
ISMIP7_INVERSION=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km_iu/maps/inversion_icepack2_budd_n3_dg0_logvel_25000_int250000_b20k_rho7500_floating_ef2_vf_wsiu.h5
ISMIP7_KEEP_CHECKPOINTS=3
ISMIP7_KSP_MAXIT=1000    # default (not exported)
ISMIP7_KSP_RTOL=1e-6    # default (not exported)
ISMIP7_LAKE_ICE_BASE=1
ISMIP7_LC=25000
ISMIP7_LC_COARSE=250000
ISMIP7_MAP_CLIP=0
ISMIP7_MAP_DEFAULT=/N/scratch/dlilien/ismip7_rehearsal25iu/run/antarctica/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_25000.h5
ISMIP7_MASS_RESIDUAL_TOL_GT=5e-5    # default (not exported)
ISMIP7_MELT_SLOPE=ant    # default (not exported)
ISMIP7_MESH=checkpoint
ISMIP7_N_FLOW=3.0
ISMIP7_OBS_DATA_ROOT=/N/project/ice_rheology/ISMIP7/antarctica/data
ISMIP7_OUTPUT=1
ISMIP7_RESCUE_ENABLED=1    # default (not exported)
ISMIP7_RESCUE_MAXIT=600    # default (not exported)
ISMIP7_RUN_TAG=rehiuamb
ISMIP7_SIN_ALPHA_ANT=0.005115    # default (not exported)
ISMIP7_SMB_ELEVATION_FEEDBACK=1
ISMIP7_SNES_ATOL=1e-50    # default (not exported)
ISMIP7_SNES_ATOL_SCALE=100    # default (not exported)
ISMIP7_SNES_DIVERGENCE_TOL=-3    # default (not exported)
ISMIP7_SNES_KSP_EW=0    # default (not exported)
ISMIP7_SNES_LINESEARCH=nleqerr    # default (not exported)
ISMIP7_SNES_LOG=stdout    # default (not exported)
ISMIP7_SNES_MAXIT=200    # default (not exported)
ISMIP7_SNES_MONITOR=0    # default (not exported)
ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE=1e-6    # default (not exported)
ISMIP7_SNES_RTOL=1e-8    # default (not exported)
ISMIP7_SNES_STOL=0    # default (not exported)
ISMIP7_SNES_TYPE=newtonls    # default (not exported)
ISMIP7_SOLVER_VIEW=0    # default (not exported)
ISMIP7_SUBCYCLES=1,4,16    # default (not exported)
ISMIP7_SUBSTEP_ADAPT=0    # default (not exported)
ISMIP7_SUBSTEP_HMIN=10    # default (not exported)
ISMIP7_SUBSTEP_INIT=1    # default (not exported)
ISMIP7_SUBSTEP_MAX=64    # default (not exported)
ISMIP7_SUBSTEP_QUIET=20    # default (not exported)
ISMIP7_SUBSTEP_TOL=1.0    # default (not exported)
ISMIP7_TRANSPORT_KSP_MAXIT=500    # default (not exported)
ISMIP7_TRANSPORT_KSP_RTOL=1e-10    # default (not exported)
OMP_NUM_THREADS=1
```

## Solver configuration

```json
{
  "continuation_steps": 8,
  "diagnostic_label": "exact-local-condensation-mumps-ptscotch",
  "diagnostic_mode": "scpc_mumps",
  "diagnostic_mode_requested": "scpc_mumps",
  "diagnostic_petsc_options": {
    "condensed_field_ksp_type": "preonly",
    "condensed_field_mat_mumps_icntl_28": 2,
    "condensed_field_mat_mumps_icntl_29": 1,
    "condensed_field_mat_type": "aij",
    "condensed_field_pc_factor_mat_solver_type": "mumps",
    "condensed_field_pc_type": "lu",
    "ksp_max_it": 1000,
    "ksp_rtol": 1e-06,
    "ksp_type": "fgmres",
    "mat_type": "matfree",
    "pc_python_type": "icepack2_tools.preconditioners.ISMIP7SCPC",
    "pc_sc_eliminate_fields": "1,2",
    "pc_type": "python",
    "pmat_type": "matfree",
    "snes_atol": 1e-50,
    "snes_divergence_tolerance": -3.0,
    "snes_linesearch_type": "nleqerr",
    "snes_max_it": 200,
    "snes_rtol": 1e-08,
    "snes_stol": 0.0,
    "snes_type": "newtonls"
  },
  "fssa_reference": "auto",
  "fssa_theta": 1.0,
  "linearization_state": "frozen",
  "mass_residual_tolerance_gt": 5e-05,
  "monitoring": {
    "destination": "stdout",
    "enabled": false,
    "view": false
  },
  "options_prefixes": {
    "diagnostic": "ismip7_diagnostic_",
    "transport": "ismip7_transport_"
  },
  "rescue_enabled": true,
  "rescue_max_it": 600,
  "snes_atol_policy": {
    "initial": 1e-50,
    "post_convergence_scale": 100.0,
    "restart_failure_residual_scale": 1e-06
  },
  "subcycles": [
    1,
    4,
    16
  ],
  "substep_adapt": null,
  "transport_petsc_options": {
    "ksp_error_if_not_converged": null,
    "ksp_max_it": 500,
    "ksp_rtol": 1e-10,
    "ksp_type": "gmres",
    "pc_type": "bjacobi",
    "sub_ksp_type": "preonly",
    "sub_pc_type": "ilu"
  }
}
```

## Budget at marker years

year | vaf_mm_sle | mass_gt | smb_gtyr | melt_gtyr | outflux_gtyr | calv_gt | clamp_gt | resid_gt | amb_gtyr | collapse_flagged_cells | collapse_removed_cells | collapse_held_cells
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
2015.025000 | 57448.257053 | 23909007.07 | 2436.4750 | 1118.5587 | 4.1996 | 31.3342 | -0.0000 | 0.0000 | -246.0581 | 0 | 0 | 0
2027.300000 | 57448.026576 | 23906703.72 | 2275.5831 | 1028.5794 | 4.2246 | 30.7636 | 0.0009 | 0.0000 | -246.0581 | 0 | 0 | 0
2039.600000 | 57451.223458 | 23905980.54 | 2658.6423 | 1273.3448 | 4.2547 | 30.6713 | 0.0000 | -0.0000 | -246.0581 | 0 | 0 | 0
2051.875000 | 57455.918930 | 23904588.36 | 2711.7563 | 1097.2361 | 4.2708 | 30.4749 | 0.0004 | 0.0000 | -246.0581 | 0 | 0 | 0
2064.150000 | 57462.776215 | 23902166.04 | 2540.7671 | 1241.2109 | 4.2649 | 30.3103 | 0.0005 | -0.0000 | -246.0581 | 0 | 0 | 0
2076.425000 | 57470.843895 | 23901243.96 | 2686.8962 | 1279.6135 | 4.2970 | 30.0930 | 0.0000 | 0.0000 | -246.0581 | 0 | 0 | 0
2088.725000 | 57481.122285 | 23898953.22 | 2805.3285 | 1528.1498 | 4.3668 | 29.7777 | 0.0196 | -0.0000 | -246.0581 | 0 | 0 | 0
2101.000000 | 57494.611934 | 23893205.96 | 3006.3584 | 1938.4248 | 4.4266 | 29.2351 | 0.0368 | 0.0000 | -246.0581 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ssp370_mri_esm2_0_rehiuamb_25000_timeseries.csv
  3440 steps, 2015.0->2101.0, dt=0.025 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2636.3   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1357.3   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1216.4   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)          -183.8   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           0.5   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1257.6   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```

## ISMIP6 ensemble overlay

```
no ISMIP6 members found under /N/u/dlilien/Quartz/data/ismip6 for ['exp01', 'exp02', 'exp03', 'exp04', 'exp05']
```
