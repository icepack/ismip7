# Core 8: ssp585_mri_esm2_0 (32 km)

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
- log: `antarctica/results/logs/core08_resume_20260720_013359.log`
- timeseries: `antarctica/results/ssp585_mri_esm2_0_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK
- ISMIP6 ensemble: outside envelope (pool: exp01,exp02,exp03,exp04,exp05)

COMPLETE: full 2015-2300 in two legs (cold cruise to the 2129.0 wall - 114 clean years - then a warm resume through the deep era to 2300). Contrast: the CESM ssp585 twin saturates at 2124.5; MRI's milder ocean lets the split scheme finish.

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
2014.1 | 56019.035912 | 23918808.61 | 2571.4159 | 1733.8308 | 821.6118 | 0.1185 | 26.2979 | 0.0000 | -338.2254
2054.9 | 56050.551112 | 23923924.06 | 2889.5666 | 1679.6322 | 863.1947 | 0.3308 | 34.6587 | 0.0000 | -338.2254
2095.8 | 56017.640082 | 23865032.42 | -1847.8574 | 3372.5576 | 727.2144 | 0.1110 | 95.2704 | 0.0000 | -338.2254
2136.6 | 55648.671545 | 23473661.33 | -17047.5407 | 6525.4121 | 492.0922 | 0.7390 | 960.2526 | 0.0000 | -338.2254
2177.5 | 54927.013631 | 22939352.36 | -30671.9316 | 17350.8823 | 4357.4849 | 0.0817 | 3336.1917 | -0.0000 | -338.2254
2218.3 | 54053.701965 | 22426559.15 | -38719.5399 | 26149.5140 | 85.6790 | 0.0632 | 5281.9072 | 0.0000 | -338.2254
2259.2 | 53286.452785 | 22036157.78 | -38366.7454 | 27584.0303 | 3121.8797 | 0.0808 | 6145.4964 | 0.0000 | -338.2254
2300.0 | 52638.928580 | 21698288.99 | -56369.2442 | 38725.0869 | 82.6481 | 0.0503 | 8677.9074 | -0.0000 | -338.2254

## Observational audit

```
ISMIP6-track audit: ssp585_mri_esm2_0_32000_timeseries.csv
  2860 steps, 2014.1->2300.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                      -21387.1   [ 2000.0,  2900.0] Gt/yr     FAIL
  shelf basal melt          14602.5   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge             605.3   [  700.0,  2400.0] Gt/yr     WARN
  dM/dt (post-2016)         -7791.1   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)         -11.9   [   -2.0,     2.0] mm SLE/yr FAIL
  budget residual             113.9   [   -0.5,     0.5] Gt/yr     FAIL
  no discharge runaway      38372.7   [peak < 6000, growth<1.5x/yr] FAIL

  OFF TRACK (6 FAIL rows)
```

## ISMIP6 ensemble overlay

```
ISMIP6 ensemble comparison: ssp585_mri_esm2_0_32000_timeseries.csv - ctrl2015_mri_esm2_0_32000_timeseries.csv
  45 members from ['exp01', 'exp02', 'exp03', 'exp04', 'exp05'], overlap 2016.0-2040.5

     year     ours   ens p5   median      p95      min      max   [mm SLE vs ctrl]
   2016.0    -1.02    -0.90     0.09     1.06    -1.71     3.10
   2020.0    -4.64    -1.23     0.39     4.04    -1.60     5.94
   2025.0    -5.20    -1.88     0.39     4.66    -2.37     8.00
   2040.5   -16.73    -7.59     1.66     9.79    -9.95    12.42

  at 2040.5: ours -16.73 mm vs ensemble [-9.95, +12.42] (45 members) -> OUTSIDE ensemble envelope
```
