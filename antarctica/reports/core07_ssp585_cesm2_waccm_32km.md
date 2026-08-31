# Core 7: ssp585_cesm2_waccm (32 km)

- date: 2026-07-20
- git: 76f9a14
- log: `antarctica/results/logs/core07_resume_20260719_235733.log`
- timeseries: `antarctica/results/ssp585_cesm2_waccm_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK

PARTIAL 2015-2124.5 (of 2300): 108+ years of ssp585 - the configuration that died at year 9 before this session's initialization/rescue work. 19 rescued events; stops in the deep-collapse era (net SMB -6659 Gt/yr, melt 13,818 Gt/yr, VAF -584 mm vs 2015, clamp/limiter withholding ~1272 Gt/yr from emptying cells) where front eras become continuous and the split scheme saturates. The monolithic (u,M,tau,h) coupling (task: gia forward_monolithic port) owns 2124->2300.

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
2014.1 | 55843.743528 | 23846190.40 | 2522.4336 | 1648.2506 | 1341.7162 | 0.2040 | 33.8789 | -0.0000 | -776.7664
2029.9 | 55825.619810 | 23832282.86 | 2663.4601 | 2122.0560 | 1266.8845 | 0.2144 | 42.4222 | 0.0000 | -776.7664
2045.6 | 55807.394681 | 23814235.15 | 1834.5888 | 2202.3907 | 998.8225 | 0.2043 | 135.6049 | 0.0000 | -776.7664
2061.4 | 55777.800525 | 23788410.53 | 1485.0048 | 3310.4540 | 670.0142 | 0.2124 | 181.5846 | 0.0000 | -776.7664
2077.2 | 55721.435415 | 23740613.40 | 613.5372 | 4397.0203 | 820.8244 | 0.4907 | 188.1391 | 0.0000 | -776.7664
2093.0 | 55633.271398 | 23664380.15 | -1272.2931 | 7619.5752 | 617.7545 | 0.1910 | 453.7094 | -0.0000 | -776.7664
2108.7 | 55511.816968 | 23563802.65 | -3318.9635 | 9454.9702 | 613.3629 | 0.5644 | 708.8146 | -0.0000 | -776.7664
2124.5 | 55346.818200 | 23441769.69 | -6659.0158 | 13817.7404 | 506.1416 | 0.5796 | 1271.9444 | -0.0000 | -776.7664

## Observational audit

```
ISMIP6-track audit: ssp585_cesm2_waccm_32000_timeseries.csv
  1105 steps, 2014.1->2124.5, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                         114.1   [ 2000.0,  2900.0] Gt/yr     FAIL
  shelf basal melt           5280.0   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge             997.2   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)         -3688.2   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)          -4.5   [   -2.0,     2.0] mm SLE/yr WARN
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       8674.9   [peak < 6000, growth<1.5x/yr] FAIL

  OFF TRACK (4 FAIL rows)
```
