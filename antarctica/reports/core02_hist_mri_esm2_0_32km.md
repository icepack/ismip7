# Core 2: hist_mri_esm2_0 (32 km)

> **SUPERSEDED - the numbers below are INVALID.** This run predates the
> annual-mean atmosphere-forcing fix (`fa7230c`): the ISMIP7 SDBN1 reader
> collapsed each year's 12 monthly slices with `isel(time=0)`, so JANUARY
> (peak austral summer, the maximum-ablation month) was applied as the whole
> year's forcing for every aSMB/ts field read through `get_field`. Ablation is
> overstated by a factor that grows with warming (MRI ssp585 continental aSMB
> integral: 2050 -229 vs +562 Gt/yr, sign flipped; 2108 -16583 vs -1276; 2300
> -82606 vs -15174), so the SMB, mass, VAF, sea-level and observational-audit
> numbers here are all wrong, as are the ISMIP6-envelope comparisons and the
> attribution of the late-century solver walls to the split scheme. It also
> predates the two Aug 2026 ice-front fixes - the calving terminus BC that
> covered only ~5% of the ice front, and the CG1 lumped-lift front-thickness
> bias now fixed by the `dg0` `ISMIP7_GEOMETRY_SPACE` default - which enter the
> MAP as well as the forward, so this core needs re-inversion, not just
> re-running with the corrected reader. Kept for provenance only; see
> `MATRIX_STATUS.md`, the single owner of the invalidation detail.

- date: 2026-07-19
- git: ec86c51
- log: `antarctica/results/logs/core02_v4_20260719_190422.log`
- timeseries: `antarctica/results/hist_mri_esm2_0_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK

Full 1850-2014 window, completed across two legs: cold to 1876.9 (x16 subcycles insufficient at the 1877 Amery emptying event) then a warm resume with the x64 rung that crossed it (one 553 s step) and cruised to 2014. Contrast with core 1: MRI-ESM2-0 GAINS ~+88 mm VAF over the window (stronger SMB) where CESM2-WACCM loses ~-87 mm. Final state hist_mri_esm2_0_32000_final.h5 (t_yr=2014.0) is the restart base for cores 4/6/8. Same open caveats as core 1 (aSMB units, 2014 handoff, peak-clause audit).

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
1850.1 | 55931.177845 | 23897546.47 | 2463.6371 | 1243.4912 | 881.0629 | 0.4562 | 0.9013 | 0.0000 | -338.2254
1873.5 | 55943.870354 | 23903705.39 | 2259.0525 | 1246.6241 | 910.3685 | 0.5220 | 12.7414 | -0.0000 | -338.2254
1896.9 | 55951.899844 | 23906406.28 | 2369.0454 | 1383.6448 | 841.0948 | 0.1212 | 15.9365 | 0.0000 | -338.2254
1920.3 | 55959.044099 | 23906486.57 | 2645.1251 | 1277.0399 | 921.3428 | 0.2657 | 13.1917 | 0.0000 | -338.2254
1943.8 | 55975.562462 | 23911639.82 | 2341.1089 | 1315.4594 | 875.8504 | 0.2384 | 13.6901 | -0.0000 | -338.2254
1967.2 | 55990.749337 | 23914342.59 | 2760.5285 | 1435.9177 | 869.5076 | 0.7074 | 17.5838 | 0.0000 | -338.2254
1990.6 | 56014.295116 | 23921742.93 | 2197.5016 | 1290.9527 | 1353.1465 | 0.2095 | 25.3632 | 0.0000 | -338.2254
2014.0 | 56019.025577 | 23918814.65 | 2092.4971 | 1622.0842 | 861.9645 | 0.2596 | 32.2344 | 0.0000 | -338.2254

## Observational audit

```
ISMIP6-track audit: hist_mri_esm2_0_32000_timeseries.csv
  1640 steps, 1850.1->2014.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2532.6   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1361.1   [  600.0,  1800.0] Gt/yr     PASS
  front discharge             927.1   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           130.0   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           0.5   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway      60446.4   [peak < 6000, growth<1.5x/yr] FAIL

  OFF TRACK (1 FAIL row)
```
