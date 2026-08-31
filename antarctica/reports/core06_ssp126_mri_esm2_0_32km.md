# Core 6: ssp126_mri_esm2_0 (32 km)

- date: 2026-07-20
- git: 223537b
- log: `antarctica/results/logs/core06_ssp126_mri_esm2_0_20260719_212415.log`
- timeseries: `antarctica/results/ssp126_mri_esm2_0_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- ISMIP6 ensemble: outside envelope (pool: exp01,exp02,exp03,exp04,exp05)

COMPLETE: full 2015-2300 in one leg, no rescues. Overlay against the (partial, to 2040.5) MRI CTRL clips to the common window; RCP8.5-class pool is a bounding comparison for ssp126.

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
2014.1 | 56019.036430 | 23918808.93 | 2571.4159 | 1738.6002 | 821.5763 | 0.1183 | 27.0933 | -0.0000 | -338.2254
2054.9 | 56045.548150 | 23922312.66 | 2798.8319 | 1651.7911 | 918.2689 | 0.3351 | 30.7930 | 0.0000 | -338.2254
2095.8 | 56072.626233 | 23928587.84 | 2624.7414 | 1397.4089 | 811.5173 | 0.1624 | 25.3504 | 0.0000 | -338.2254
2136.6 | 56102.880095 | 23937613.51 | 2424.4818 | 1724.9353 | 814.1727 | 0.2036 | 55.9157 | 0.0000 | -338.2254
2177.5 | 56139.909265 | 23949458.57 | 2341.0551 | 2007.4573 | 842.4044 | 0.3808 | 40.4288 | 0.0000 | -338.2254
2218.3 | 56180.221956 | 23960404.57 | 2648.9085 | 1739.8659 | 871.9367 | 0.1432 | 43.9441 | 0.0000 | -338.2254
2259.2 | 56219.378794 | 23971966.03 | 2828.3224 | 1654.5106 | 868.3694 | 0.1401 | 50.8607 | -0.0000 | -338.2254
2300.0 | 56267.362665 | 23987292.13 | 3179.1825 | 1746.6243 | 1322.8246 | 0.5873 | 153.2950 | 0.0000 | -338.2254

## Observational audit

```
ISMIP6-track audit: ssp126_mri_esm2_0_32000_timeseries.csv
  2860 steps, 2014.1->2300.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2687.2   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1651.7   [  600.0,  1800.0] Gt/yr     PASS
  front discharge             855.7   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           239.2   [ -400.0,   200.0] Gt/yr     WARN
  dVAF/dt (post-2016)           0.9   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       2624.0   [peak < 6000, growth<1.5x/yr] PASS

  ON TRACK (0 FAIL rows)
```

## ISMIP6 ensemble overlay

```
ISMIP6 ensemble comparison: ssp126_mri_esm2_0_32000_timeseries.csv - ctrl2015_mri_esm2_0_32000_timeseries.csv
  45 members from ['exp01', 'exp02', 'exp03', 'exp04', 'exp05'], overlap 2016.0-2040.5

     year     ours   ens p5   median      p95      min      max   [mm SLE vs ctrl]
   2016.0    -0.50    -0.90     0.09     1.06    -1.71     3.10
   2020.0    -1.56    -1.23     0.39     4.04    -1.60     5.94
   2025.0    -5.18    -1.88     0.39     4.66    -2.37     8.00
   2040.5   -15.05    -7.59     1.66     9.79    -9.95    12.42

  at 2040.5: ours -15.05 mm vs ensemble [-9.95, +12.42] (45 members) -> OUTSIDE ensemble envelope
```
