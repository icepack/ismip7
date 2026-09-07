# Core 5: ssp126_cesm2_waccm (32 km)

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

- date: 2026-07-20
- git: 223537b
- log: `antarctica/results/logs/core05_resume_20260719_230542.log`
- timeseries: `antarctica/results/ssp126_cesm2_waccm_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK
- ISMIP6 ensemble: outside envelope (pool: exp01,exp02,exp03,exp04,exp05)

COMPLETE: full 2015-2300 in two legs (cold to the 2056.4 wall, warm resume through to 2300). Ensemble pool is the RCP8.5-class set (high-emission upper bound - treat the overlay as bounding, not analog, for ssp126).

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
2014.1 | 55843.737744 | 23846183.07 | 2522.4336 | 1724.9029 | 1341.4185 | 0.2017 | 34.1800 | 0.0000 | -776.7664
2054.9 | 55794.983695 | 23806509.57 | 2727.5606 | 2346.3413 | 1309.8575 | 0.5923 | 54.8746 | 0.0000 | -776.7664
2095.8 | 55723.352139 | 23743560.03 | 2571.4101 | 3221.3551 | 721.7129 | 0.2101 | 103.0170 | 0.0000 | -776.7664
2136.6 | 55614.053573 | 23669798.77 | 1991.2883 | 3852.1368 | 569.1959 | 0.1886 | 232.1471 | 0.0000 | -776.7664
2177.5 | 55556.564364 | 23624185.51 | 2383.7013 | 3313.3909 | 640.0487 | 1.4570 | 204.9473 | -0.0000 | -776.7664
2218.3 | 55494.691772 | 23576598.09 | 1999.4495 | 3164.5810 | 1687.6659 | 0.2918 | 179.0062 | -0.0000 | -776.7664
2259.2 | 55436.514623 | 23535627.18 | 2379.7392 | 3616.6360 | 1061.0458 | 0.1714 | 356.2701 | 0.0000 | -776.7664
2300.0 | 55383.516157 | 23493615.86 | 1862.2169 | 4458.8119 | 1482.4435 | 0.3159 | 277.5496 | 0.0000 | -776.7664

## Observational audit

```
ISMIP6-track audit: ssp126_cesm2_waccm_32000_timeseries.csv
  2860 steps, 2014.1->2300.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2162.4   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           3242.9   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge            1606.6   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)         -1235.3   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)          -1.6   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.5   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway      77453.3   [peak < 6000, growth<1.5x/yr] FAIL

  OFF TRACK (3 FAIL rows)
```

## ISMIP6 ensemble overlay

```
ISMIP6 ensemble comparison: ssp126_cesm2_waccm_32000_timeseries.csv - ctrl2015_cesm2_waccm_32000_timeseries.csv
  45 members from ['exp01', 'exp02', 'exp03', 'exp04', 'exp05'], overlap 2016.0-2101.0

     year     ours   ens p5   median      p95      min      max   [mm SLE vs ctrl]
   2016.0     0.81    -0.90     0.09     1.06    -1.71     3.10
   2020.0     8.28    -1.23     0.39     4.04    -1.60     5.94
   2025.0    12.03    -1.88     0.39     4.66    -2.37     8.00
   2050.0    42.77   -13.46     3.04    16.85   -16.97    28.58
   2075.0    85.67   -29.72     8.71    55.32   -43.32    91.83
   2100.0   132.11   -53.27    11.98   117.92   -84.83   161.22
   2101.0   133.74   -53.29     5.98    74.25   -59.37    80.55

  at 2101.0: ours +133.74 mm vs ensemble [-59.37, +80.55] (27 members) -> OUTSIDE ensemble envelope
```
