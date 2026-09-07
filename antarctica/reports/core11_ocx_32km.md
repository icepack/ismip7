# Core 11: ocx (32 km)

> **SUPERSEDED - the numbers below are INVALID.** This run predates the two
> Aug 2026 ice-front fixes: the calving terminus BC covered only ~5% of the
> ice front (the 2500 m boundary-id sidecar was used on the 32 km mesh), and
> the CG1 lumped-lift front-thickness bias, now fixed by the `dg0`
> `ISMIP7_GEOMETRY_SPACE` default. Both defects enter the MAP as well as the
> forward, so this core needs re-inversion and re-running. Kept for provenance
> only; see `MATRIX_STATUS.md`, the single owner of the invalidation detail.

- date: 2026-07-20
- git: 2f55628
- log: `antarctica/results/logs/core11_ocx_20260720_013359.log`
- timeseries: `antarctica/results/ocx_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK

COMPLETE 1990-2025: obs-constrained (RACMO actual-year SMB + OI-climatology ocean) from a ~2015 geometry, so 1990-2000 is relaxation and 2000-2025 the validation window. NOTE: the real OCX forcing product (RACMO2.3p2-ERA + cold/main/vary/warm ocean) is now on the share at /ISMIP7/AIS/OCX - switching this core off the RACMO+OI stopgap onto it is a follow-up.

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
1990.1 | 55931.177845 | 23897546.47 | 2329.5594 | 1549.0472 | 881.0608 | 0.4626 | 0.9053 | 0.0000 | 101.4300
1995.1 | 55933.515722 | 23899491.39 | 2452.7724 | 1617.8415 | 835.1687 | 0.1838 | 13.4450 | 0.0000 | 101.4300
2000.1 | 55935.084441 | 23900899.32 | 2431.5487 | 1647.2140 | 1020.6730 | 0.1680 | 46.6523 | 0.0000 | 101.4300
2005.1 | 55937.732648 | 23902095.34 | 2669.0142 | 1664.3017 | 902.4732 | 0.2336 | 38.3102 | 0.0000 | 101.4300
2010.0 | 55939.948059 | 23903335.36 | 2517.7362 | 1681.4336 | 863.8369 | 0.5927 | 34.8421 | -0.0000 | 101.4300
2015.0 | 55940.709411 | 23903966.29 | 2291.3376 | 1691.8999 | 815.9447 | 0.1829 | 17.5162 | -0.0000 | 101.4300
2020.0 | 55943.189717 | 23905238.93 | 2484.6944 | 1692.3668 | 918.1189 | 0.6006 | 15.4066 | -0.0000 | 101.4300
2025.0 | 55947.475389 | 23907196.48 | 2580.9458 | 1712.0333 | 850.9313 | 0.6457 | 15.6473 | 0.0000 | 101.4300

## Observational audit

```
ISMIP6-track audit: ocx_32000_timeseries.csv
  350 steps, 1990.1->2025.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2515.9   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1665.7   [  600.0,  1800.0] Gt/yr     PASS
  front discharge             856.8   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           279.4   [ -400.0,   200.0] Gt/yr     WARN
  dVAF/dt (post-2016)           0.5   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1058.4   [peak < 6000, growth<1.5x/yr] PASS

  ON TRACK (0 FAIL rows)
```
