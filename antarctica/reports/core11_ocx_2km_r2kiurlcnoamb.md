# Core 11: ocx_r2kiurlcnoamb (2 km)

- date: 2026-10-14
- git: c4c9eb8
- log: `/N/scratch/dlilien/ismip7_prod2kiu/logs/ismip7_fwd_11912190.out`
- timeseries: `/N/project/ice_rheology/ISMIP7/antarctica/results/prod_2km_iu/results/ocx_r2kiurlcnoamb_2000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: none reported in `/N/scratch/dlilien/ismip7_prod2kiu/logs/ismip7_fwd_11912190.out` (expected for a CTRL running on the RACMO climatology, which builds no ESM pool)
- Forcing provenance: atmosphere acabf RACMO2.3p2-ERA OCX SDBN1-8000m v1
- Forcing provenance: ocean tf expert-judgment OCX ocean/main v2
- Forcing provenance: ocean so expert-judgment OCX ocean/main v2
- Forcing provenance: ocean melt calibration deltaT_per_basin_2000_K6.500e-05.npz sha256 3d2e5348b207090e8f8ad39a44cb152e6c482f550f6041be1a311238bab1a7e4 (named with ISMIP7_DELTAT_PER_BASIN_NPZ): K 6.500e-05 (K50) with a thermal-forcing offset per basin, fitted on antarctica_5000_2000_buffered20000; this run's mesh is antarctica_5000_2000_buffered20000
- Forcing provenance: atmosphere dacabfdz RACMO2.3p2-ERA OCX SDBN1-8000m v2
- SMB-elevation feedback: dacabfdz from RACMO2.3p2-ERA OCX SDBN1-8000m v2, surface change from the chain's initial state (H_init)
- Free-surface stabilization: measured from the velocity of the last advance, with that advance's tendency as a load (ISMIP7_FSSA_REFERENCE=step)
- Free-surface stabilization: theta 1 on the lagged thickness-velocity coupling (ISMIP7_FSSA_THETA)
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)

IU's 2 km Budd production suite (tag r2kiurlcnoamb): the relaxed MAP's controls on BedMachine geometry, no apparent mass-balance reference. Production grade; the group chooses the submission line on 20 or 21 October.

## Run environment

```
ISMIP7_APPARENT_MB=0
ISMIP7_AUTO_RESUME=1
ISMIP7_BNDIDS=/N/scratch/dlilien/ismip7_rehearsal25iu/run2kp/antarctica/mesh/boundary_ids_antarctica_5000_2000_buffered20000.json
ISMIP7_BUFFER_M=20000
ISMIP7_CHECKPOINT_EVERY_YR=5
ISMIP7_CLIM_END=2029    # default (not exported)
ISMIP7_CLIM_SCENARIO=ssp126    # default (not exported)
ISMIP7_CLIM_START=2000    # default (not exported)
ISMIP7_CONTINUATION_STEPS=8    # default (not exported)
ISMIP7_DATA_ROOT=/N/project/ice_rheology/ISMIP7/ISMIP7/AIS
ISMIP7_DELTAT_PER_BASIN_NPZ=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_2km_iu/calibration/deltaT_per_basin_2000_K6.500e-05.npz
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_mumps
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER_CANONICAL=scpc_mumps    # resolved canonical mode
ISMIP7_DT=0.025
ISMIP7_EXPERIMENT=ocx
ISMIP7_FIXED_FRONT=1
ISMIP7_FRACTURE=none
ISMIP7_FRICTION=budd
ISMIP7_FSSA_REFERENCE=auto    # default (not exported)
ISMIP7_FSSA_THETA=1    # default (not exported)
ISMIP7_GEOMETRY_SPACE=dg0
ISMIP7_INVERSION=/N/project/ice_rheology/ISMIP7/antarctica/results/reinvert_2km/final/inversion_icepack2_budd_n3_dg0_logvel_2000_int5000_b20k_rho7500_floating_ef2_vf_relax2014.h5
ISMIP7_KEEP_CHECKPOINTS=3
ISMIP7_KSP_MAXIT=1000    # default (not exported)
ISMIP7_KSP_RTOL=1e-6    # default (not exported)
ISMIP7_LAKE_ICE_BASE=1
ISMIP7_LC=2000
ISMIP7_LC_COARSE=5000
ISMIP7_MAP_CLIP=0
ISMIP7_MASS_RESIDUAL_TOL_GT=5e-5    # default (not exported)
ISMIP7_MELT_SLOPE=ant    # default (not exported)
ISMIP7_MESH=checkpoint
ISMIP7_N_FLOW=3.0
ISMIP7_OBS_DATA_ROOT=/N/project/ice_rheology/ISMIP7/antarctica/data
ISMIP7_OCX_FORCING=protocol
ISMIP7_OCX_OCEAN=main
ISMIP7_OUTPUT=1
ISMIP7_RESCUE_ENABLED=1    # default (not exported)
ISMIP7_RESCUE_MAXIT=600    # default (not exported)
ISMIP7_RUN_TAG=r2kiurlcnoamb
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
2003.025000 | 56727.714464 | 24072705.96 | 2432.6632 | 1049.1455 | 4.5836 | 45.5004 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0
2006.300000 | 56732.927220 | 24073479.46 | 2580.1985 | 983.4486 | 1.5965 | 32.8015 | 0.0008 | 0.0000 | 0.0000 | 0 | 0 | 0
2009.600000 | 56738.072693 | 24074217.57 | 2557.9542 | 942.6112 | 1.6079 | 32.8624 | 0.0089 | 0.0000 | 0.0000 | 0 | 0 | 0
2012.875000 | 56743.019510 | 24075025.38 | 2391.1056 | 915.5973 | 1.4677 | 32.6490 | 0.0167 | -0.0000 | 0.0000 | 0 | 0 | 0
2016.150000 | 56747.996646 | 24075989.74 | 2753.3737 | 898.8570 | 1.4084 | 32.2625 | 0.0234 | 0.0000 | 0.0000 | 0 | 0 | 0
2019.425000 | 56753.622316 | 24077315.89 | 2494.7227 | 882.5107 | 1.3834 | 32.1039 | 0.0285 | 0.0000 | 0.0000 | 0 | 0 | 0
2022.725000 | 56761.067865 | 24079448.02 | 3271.6137 | 872.5322 | 1.2450 | 32.1942 | 0.0343 | 0.0000 | 0.0000 | 0 | 0 | 0
2026.000000 | 56767.755482 | 24081350.93 | 2725.6042 | 861.1629 | 1.2026 | 32.2008 | 0.0466 | 0.0000 | 0.0000 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ocx_r2kiurlcnoamb_2000_timeseries.csv
  920 steps, 2003.0->2026.0, dt=0.025 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2600.6   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt            919.8   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1301.8   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           392.6   [ -400.0,   200.0] Gt/yr     WARN
  dVAF/dt (post-2016)           1.8   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1824.6   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```
