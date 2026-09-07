# Core 9: ctrl2015_cesm2_waccm (32 km)

> **SUPERSEDED - the numbers below are INVALID.** This run predates the two
> Aug 2026 ice-front fixes: the calving terminus BC covered only ~5% of the
> ice front (the 2500 m boundary-id sidecar was used on the 32 km mesh), and
> the CG1 lumped-lift front-thickness bias, now fixed by the `dg0`
> `ISMIP7_GEOMETRY_SPACE` default. Both defects enter the MAP as well as the
> forward, so this core needs re-inversion and re-running. Kept for provenance
> only; see `MATRIX_STATUS.md`, the single owner of the invalidation detail.

- date: 2026-07-20
- git: 2f55628
- log: `antarctica/results/logs/core09_ctrl_cesm_20260719_233511.log`
- timeseries: `antarctica/results/ctrl2015_cesm2_waccm_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK

COMPLETE 2015-2300: the balanced control. RACMO climatology baseline + OI-climatology ocean; the ESM enters only via the acabf fallback, so cores 9 and 10 are physically near-identical (per preflight). This is the CTRL that the projections difference against (compare_ismip6 --ctrl-csv). It IS the control, so no ensemble overlay.

## Run environment

```
ISMIP7_APPARENT_MB=1
ISMIP7_DT=0.1
ISMIP7_FIXED_FRONT=1
ISMIP7_FRICTION=budd
ISMIP7_LC=32000
ISMIP7_OUTPUT_INTERVAL=100
ISMIP7_SUBCYCLES=1,4,16,64
OMP_NUM_THREADS=1
```

## Budget at marker years

year | vaf_mm_sle | mass_gt | smb_gtyr | melt_gtyr | outflux_gtyr | calv_gt | clamp_gt | resid_gt | amb_gtyr
--- | --- | --- | --- | --- | --- | --- | --- | --- | ---
2015.1 | 55931.177845 | 23897546.47 | 2526.4241 | 1549.0472 | 881.0612 | 0.4648 | 0.9056 | 0.0000 | -95.4160
2055.8 | 55934.070375 | 23901127.79 | 2526.4241 | 1713.7893 | 766.7838 | 0.5718 | 16.0051 | 0.0000 | -95.4160
2096.5 | 55938.923739 | 23902467.95 | 2526.4241 | 1742.6908 | 824.0601 | 0.7416 | 14.7972 | 0.0000 | -95.4160
2137.2 | 55945.379265 | 23904257.39 | 2526.4241 | 1668.6406 | 863.6733 | 0.1629 | 15.7543 | 0.0000 | -95.4160
2177.9 | 55967.637050 | 23911441.20 | 2526.4241 | 1630.6900 | 1031.1062 | 0.1638 | 152.7110 | 0.0000 | -95.4160
2218.6 | 55960.323796 | 23911081.05 | 2526.4241 | 1560.1449 | 1037.3634 | 0.6365 | 15.2570 | 0.0000 | -95.4160
2259.3 | 55946.354123 | 23906401.22 | 2526.4241 | 1791.0505 | 888.0263 | 0.1749 | 27.3108 | 0.0000 | -95.4160
2300.0 | 55951.935173 | 23910953.60 | 2526.4241 | 1576.4524 | 1619.3209 | 0.6074 | 18.5300 | 0.0000 | -95.4160

## Observational audit

```
ISMIP6-track audit: ctrl2015_cesm2_waccm_32000_timeseries.csv
  2850 steps, 2015.1->2300.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2526.4   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1668.0   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1062.4   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)            46.6   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           0.1   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway      51012.2   [peak < 6000, growth<1.5x/yr] FAIL

  OFF TRACK (1 FAIL row)
```
