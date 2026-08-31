# Core 3: ssp370_cesm2_waccm (32 km)

- date: 2026-07-20
- git: 74f5bf8
- log: `antarctica/results/logs/core03_resume_20260719_230721.log`
- timeseries: `antarctica/results/ssp370_cesm2_waccm_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK

COMPLETE: full 2015-2100 ssp370 window in two legs (cold to 2069.7, warm resume through the wall to 2100; 17 rescued events total). End-state is strongly forced: net SMB -1987 Gt/yr (surface loss regime), melt 5064 Gt/yr, VAF -265 mm vs 2015 (~+0.26 m SLE contribution). Ensemble overlay lands with the CTRL pair.

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
2014.1 | 55843.737389 | 23846183.30 | 2522.4336 | 1718.2538 | 1341.5888 | 0.2041 | 33.7633 | -0.0000 | -776.7664
2026.4 | 55828.993218 | 23837286.37 | 2489.2076 | 1967.1982 | 1260.1075 | 0.1923 | 121.4826 | 0.0000 | -776.7664
2038.6 | 55817.713171 | 23826684.01 | 1974.3440 | 2334.5090 | 1183.7883 | 0.1889 | 124.1745 | 0.0000 | -776.7664
2050.9 | 55797.639114 | 23809746.55 | 2257.0193 | 2659.6811 | 1251.0101 | 0.2239 | 59.7118 | 0.0000 | -776.7664
2063.2 | 55774.587073 | 23788053.45 | 1563.4653 | 3252.3594 | 1350.4303 | 0.5442 | 100.0381 | 0.0000 | -776.7664
2075.5 | 55746.403605 | 23762367.19 | 1623.8723 | 3548.7191 | 901.1629 | 0.1592 | 162.1355 | 0.0000 | -776.7664
2087.7 | 55711.040816 | 23731271.11 | 679.4374 | 4473.8720 | 2586.3942 | 0.1924 | 223.4036 | -0.0000 | -776.7664
2100.0 | 55666.395398 | 23687648.97 | -1987.3181 | 5063.6682 | 489.4587 | 0.1513 | 291.2418 | 0.0000 | -776.7664

## Observational audit

```
ISMIP6-track audit: ssp370_cesm2_waccm_32000_timeseries.csv
  860 steps, 2014.1->2100.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        1932.3   [ 2000.0,  2900.0] Gt/yr     WARN
  shelf basal melt           3108.4   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge            1147.7   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)         -1858.2   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)          -2.1   [   -2.0,     2.0] mm SLE/yr WARN
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       4682.9   [peak < 6000, growth<1.5x/yr] FAIL

  OFF TRACK (3 FAIL rows)
```
