# Core 1: hist_cesm2_waccm (32 km)

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
- git: 235cdf6
- log: `antarctica/results/logs/core01_v3_20260719_180733.log`
- timeseries: `antarctica/results/hist_cesm2_waccm_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK

Full 1850-2014 window completed. 13 hard-era events (front-cell emptying, first at ~1873) crossed via the dt-subcycle rescue ladder (x4/x16 with trust-region + rescue speed limiter). Final state hist_cesm2_waccm_32000_final.h5 (t_yr=2014.0) is the projection restart base for cores 3/5/7. Known caveats carried in-report: the ~9 percent aSMB unit question and the 2014->2015 projection handoff (both open user decisions); observational-audit envelopes are present-day references applied to a historical-mean run - WARNs expected.

## Run environment

```
ISMIP7_APPARENT_MB=1
ISMIP7_DT=0.1
ISMIP7_FIXED_FRONT=1
ISMIP7_FRICTION=budd
ISMIP7_LC=32000
ISMIP7_OUTPUT_INTERVAL=100
ISMIP7_SUBCYCLES=1,4,16
OMP_NUM_THREADS=1
```

## Budget at marker years

year | vaf_mm_sle | mass_gt | smb_gtyr | melt_gtyr | outflux_gtyr | calv_gt | clamp_gt | resid_gt | amb_gtyr
--- | --- | --- | --- | --- | --- | --- | --- | --- | ---
1850.1 | 55931.177846 | 23897546.47 | 2783.1069 | 1124.1966 | 881.0670 | 0.4714 | 0.8948 | 0.0000 | -776.7664
1873.5 | 55919.353592 | 23893876.20 | 2863.9971 | 1309.0831 | 910.4889 | 0.1632 | 14.5754 | -0.0000 | -776.7664
1896.9 | 55903.353734 | 23886777.58 | 2703.3646 | 1418.4027 | 814.4768 | 0.6311 | 17.0863 | 0.0000 | -776.7664
1920.3 | 55893.279290 | 23880328.57 | 2436.7839 | 1336.5743 | 787.9429 | 0.5597 | 16.2098 | -0.0000 | -776.7664
1943.8 | 55880.383600 | 23872191.11 | 2170.1154 | 1387.7592 | 926.5330 | 0.1756 | 22.0024 | -0.0000 | -776.7664
1967.2 | 55876.696311 | 23868024.83 | 2302.3328 | 1369.4627 | 1189.8193 | 0.6122 | 18.1043 | 0.0000 | -776.7664
1990.6 | 55840.421160 | 23848005.70 | 2655.2884 | 1368.2127 | 808.1643 | 0.2915 | 30.7289 | 0.0000 | -776.7664
2014.0 | 55843.912492 | 23846281.16 | 3148.1438 | 1589.0361 | 1252.9472 | 0.1966 | 122.6249 | 0.0000 | -776.7664

## Observational audit

```
ISMIP6-track audit: hist_cesm2_waccm_32000_timeseries.csv
  1640 steps, 1850.1->2014.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2519.9   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1371.4   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1050.1   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)          -314.0   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)          -0.5   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.1   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway      89479.0   [peak < 6000, growth<1.5x/yr] FAIL

  OFF TRACK (1 FAIL row)
```
