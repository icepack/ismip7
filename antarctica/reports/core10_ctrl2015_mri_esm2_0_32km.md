# Core 10: ctrl2015_mri_esm2_0 (32 km)

> **SUPERSEDED - the numbers below are INVALID.** This run predates the two
> Aug 2026 ice-front fixes: the calving terminus BC covered only ~5% of the
> ice front (the 2500 m boundary-id sidecar was used on the 32 km mesh), and
> the CG1 lumped-lift front-thickness bias, now fixed by the `dg0`
> `ISMIP7_GEOMETRY_SPACE` default. Both defects enter the MAP as well as the
> forward, so this core needs re-inversion and re-running. Kept for provenance
> only; see `MATRIX_STATUS.md`, the single owner of the invalidation detail.

- date: 2026-07-20
- git: 2f55628
- log: `antarctica/results/logs/core10_resume_20260720_013359.log`
- timeseries: `antarctica/results/ctrl2015_mri_esm2_0_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK

PARTIAL 2015-2040.5: walled at a front-emptying event that the split scheme + subcycle ladder cannot cross here (a warm resume re-walls at step 1). Its physical twin core 9 (same RACMO baseline; ESM enters only via the acabf fallback) reached 2300, so the CTRL behavior IS demonstrated - core 10 is the redundant ESM label and its 2040.5+ stretch is a split-scheme-edge artifact, not new physics. Monolithic coupling would close it.

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
2018.7 | 55931.463459 | 23898331.73 | 2526.4241 | 1552.0474 | 646.3359 | 0.6037 | 5.1375 | -0.0000 | -95.4160
2022.4 | 55931.696751 | 23899193.12 | 2526.4241 | 1694.5074 | 960.0517 | 0.6334 | 15.4496 | 0.0000 | -95.4160
2026.0 | 55931.991654 | 23899490.47 | 2526.4241 | 1654.7217 | 890.4382 | 0.5242 | 19.5074 | -0.0000 | -95.4160
2029.6 | 55932.176607 | 23899691.89 | 2526.4241 | 1678.7957 | 871.1498 | 0.1900 | 17.2068 | -0.0000 | -95.4160
2033.2 | 55932.637761 | 23900079.98 | 2526.4241 | 1731.1614 | 774.1029 | 0.2087 | 14.2390 | 0.0000 | -95.4160
2036.9 | 55932.884743 | 23900257.86 | 2526.4241 | 1677.4765 | 900.9694 | 0.2018 | 13.7516 | -0.0000 | -95.4160
2040.5 | 55932.822184 | 23900405.49 | 2526.4241 | 1694.5104 | 881.4232 | 0.5979 | 12.2411 | 0.0000 | -95.4160

## Observational audit

```
ISMIP6-track audit: ctrl2015_mri_esm2_0_32000_timeseries.csv
  255 steps, 2015.1->2040.5, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2526.4   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1649.3   [  600.0,  1800.0] Gt/yr     PASS
  front discharge             857.0   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           111.2   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           0.1   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1296.3   [peak < 6000, growth<1.5x/yr] PASS

  ON TRACK (0 FAIL rows)
```
