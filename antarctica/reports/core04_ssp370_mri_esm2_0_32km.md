# Core 4: ssp370_mri_esm2_0 (32 km)

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
