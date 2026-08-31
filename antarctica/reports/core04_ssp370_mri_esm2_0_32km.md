# Core 4: ssp370_mri_esm2_0 (32 km)

- date: 2026-07-19
- git: 634f5ac
- log: `antarctica/results/logs/core04_ssp370_mri_esm2_0_20260719_212413.log`
- timeseries: `antarctica/results/ssp370_mri_esm2_0_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK

FIRST COMPLETE PROJECTION: full 2015-2100 ssp370 window from the core-2 historical final, zero rescue events (MRI forcing is mild enough for straight dt=0.1 cruising). ISMIP6-ensemble overlay follows once the CTRL pair (cores 9/10) lands - projections difference against the balanced control.

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
2014.1 | 56019.034152 | 23918806.87 | 2571.4159 | 1751.7189 | 821.4338 | 0.1184 | 26.3265 | -0.0000 | -338.2254
2026.4 | 56026.256479 | 23919068.64 | 2356.7635 | 1635.5444 | 841.5694 | 0.1202 | 29.1859 | 0.0000 | -338.2254
2038.6 | 56032.425374 | 23919785.18 | 2591.8527 | 1844.3048 | 815.2133 | 0.1365 | 37.1762 | 0.0000 | -338.2254
2050.9 | 56042.808584 | 23921810.97 | 3105.2911 | 1607.3195 | 813.1443 | 0.2323 | 31.5218 | 0.0000 | -338.2254
2063.2 | 56044.645516 | 23918893.99 | 2061.9408 | 1689.5353 | 810.1446 | 0.3151 | 34.4887 | 0.0000 | -338.2254
2075.5 | 56042.443057 | 23912407.74 | 2630.6985 | 1860.5219 | 806.2282 | 0.3362 | 37.5466 | 0.0000 | -338.2254
2087.7 | 56043.308308 | 23906183.69 | 2305.9393 | 2151.8215 | 741.1288 | 0.1193 | 48.2721 | 0.0000 | -338.2254
2100.0 | 56040.089817 | 23893929.27 | 2043.7651 | 2633.0672 | 769.0480 | 0.7215 | 85.0891 | 0.0000 | -338.2254

## Observational audit

```
ISMIP6-track audit: ssp370_mri_esm2_0_32000_timeseries.csv
  860 steps, 2014.1->2100.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2303.4   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1827.3   [  600.0,  1800.0] Gt/yr     WARN
  front discharge             807.9   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)          -292.9   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           0.2   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway        919.8   [peak < 6000, growth<1.5x/yr] PASS

  ON TRACK (0 FAIL rows)
```
