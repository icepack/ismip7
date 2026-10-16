# A 25 km rehearsal of IU's Budd forwards

IU Quartz, 9 October 2026. The recipe of the September rehearsal
(`rehearsal_25km.md`, issue 138) run again for IU's Budd line: a 25 km Budd
MAP re-inverted from IU's non-relaxed 2 km MAP refitted under the front-cell
rule (run record `inversion-2km-budd-b20k-ef2-vf`), on that MAP's own
objective and PR 163's code, then the melt refit, the MAP checks, all eleven
cores without and with the apparent mass-balance reference (issue #104), the
track audit, the core reports, the ISMIP7 writer, the compliance checker and
the scalar tool. The two rehearsals share the mesh, the recipe and the
post-processing and differ in more than the MAP. September's ran `b554238`
with scpc_gamg and no FSSA (two cores resumed under scpc_mumps), from Rice's
0241 snapshot under Rice's objective (sigma 0.3). This one ran `288382a` with
vertex_front sampling, exact_front 2, the vertex drag gate, the 2.5 m membrane
floor, FSSA theta 1 and scpc_mumps throughout, from IU's vf MAP under IU's
objective. A difference between IU's rehearsal and September's belongs to the
whole configuration.

Everything the rehearsal wrote is kept on Quartz project storage under
`antarctica/results/rehearsal_25km_iu/`, with `MANIFEST.md` and `SHA256SUMS`.
The CESM2-WACCM historical was also run at 2 km from the vf MAP itself, on its
own mesh, without and with the reference (The 2 km historicals, below).

## Configuration

| item | value |
|---|---|
| code | `288382a`, the head of PR 163 (main after PR 164, the front-cell rule of issue 167, exact_front 2, FSSA on by default) |
| mesh | `antarctica_250000_25000_buffered20000`, the September build (md5 `3e8b44b0`, 4,509 vertices, 7,615 cells) |
| warm start | IU's Budd MAP `inversion_icepack2_budd_n3_dg0_logvel_2000_int5000_b20k_rho7500_floating_ef2_vf.h5`, sha256 `d0ff592e` |
| MAP | 25 km Budd, 300 iterations on IU's objective; sha256 `551f0370` |
| melt | K 6.5e-5 with per-basin offsets refitted on this mesh under `ISMIP7_RASTER_SAMPLE=vertex_front`; sha256 `385c1f47` |
| forward | dt 0.025, scpc_mumps on 8 ranks, `ISMIP7_MAP_CLIP=0`, fixed front, no collapse, ISMIP7 output on, SMB-elevation feedback on, FSSA at its default (theta 1, with the `auto` reference resolving to `step` in attempt A and `start` in attempt B), 2003 start; exact_front 2, the vertex drag gate, the 2.5 m membrane floor and the vertex_front sampling from the MAP |
| attempt A | `ISMIP7_APPARENT_MB=0`, tag `rehiunoamb` |
| attempt B | `ISMIP7_APPARENT_MB=1`, tag `rehiuamb` |

Both attempts ran side by side. IU's objective is the one its 2 km chains
used: bilaplacian prior with sigma 30 and rho 7.5 km on both controls, no
dH/dt term, L-BFGS-B without a metric, fluidity on floating ice only, the
log-velocity weight 85380.4 held from the warm start (issue 68). Every
objective key the 25 km MAP records matches the warm start's.

## Stages and gates

| stage | job | verdict |
|---|---|---|
| forcing audit | login node | `audit_forcing_versions.py` exits 0: 481 mirror entries, none missing, behind, pinned or replaced; 93 files predate the download manifest |
| melt refit (G2) | 11888047, 11888663, 11888664 | every basin roots between -0.848 and +0.936 K for 1067.4 Gt/yr; the K50 alone melts 1320 Gt/yr on this mesh. The forward's own callback reproduces every basin at ratio 1.0000 (basin 15 at 1.0001). Against the OCX ocean main v2 no basin and no 256 km block is flagged in 2000, 2015 or 2025 |
| re-inversion (G3) | 11888050 | the startup ramp converged on its first rung. 300 iterations to the limit in 1 h 26 min on 8 ranks: misfit 5.28e5 to 3.74e3, total 4.23e3, \|grad\| still about 90. theta spans -23.3 to 17.9 and phi -27.0 to 9.9 on ice (September's MAP: -0.9 to 4.5 and -2.6 to 0.4); 61 vertices carry \|theta\| above 10 and 64 \|phi\| above 10 |
| discharge score (G4a) | 11888665 | 2058 Gt/yr across the grounding line against 1803 with the observed velocity, ratio 1.14 (September 1.52): 1.73 where the observed speed is under 100 m/yr, 1.12 at 100 to 500, 0.73 at 500 to 1500 and 0.57 above 1500 |
| census and self-consistency (G4b) | 11888666 | 0 of 1,919 floating cells carry friction under the HAF gate (the old sign test would have put it on 253); the forward re-solves the MAP's velocity to rel L2 2.5e-16 |
| front flux (G4d) | 11888667 | the t=0 front carries 1,231 Gt/yr, 861 out of the 12,713 km floating front and 370 out of grounded fronts, against 825 with the observed velocity (581 and 243); grounding-line flux 1,688 against 1,560 |
| probes (G4c) | 11888668, 11888669 | both stopped at step 1 by the runaway tripwire at 2e4 m/yr, September's setting: 2.2e4 m/yr without the reference and 3.2e4 with it, on one front vertex at (-2026, 485) km that carries no velocity observation. The step's solve had converged and the transport closed its mass. IU went to the matrix without a second probe; the matrix arms no tripwire |

The twelve fastest vertices of the 25 km MAP all lie on the ice front with no
velocity observation (`fast_cell_r25iu.py`, job 11888753). The 2 km warm start
has the same pattern at a larger speed: 1.4e5 m/yr on floating front vertices
near Shirase Glacier, and phi reaching -34 on ice (job 11888754).

## The matrix

Every core of both attempts reached its end year with the budget closed to
resid 0.00 on every row, on scpc_mumps throughout. Wall time is on 8 ranks.
Peak speed is the largest vertex speed at any step.

Attempt B, with the reference (`rehiuamb`):

| core | experiment | wall (min) | rescued steps | VAF change (mm SLE) | mass change (Gt) | peak speed (m/yr, year) | track audit | checker errors | scalars |
|---|---|---|---|---|---|---|---|---|---|
| C001 | historical CESM2-WACCM | 5 | 0 | +2.6 | -814 | 3.2e+04 (2005) | ON TRACK | 0 | exit 0 |
| C002 | historical MRI-ESM2-0 | 5 | 0 | -0.8 | -1,450 | 3.2e+04 (2006) | ON TRACK | 0 | exit 0 |
| C003 | ssp370 CESM2-WACCM | 19 | 0 | +53.8 | -140,087 | 3.3e+04 (2099) | OFF TRACK (dM/dt (post-2016) -1641.4) | 0 | exit 0 |
| C004 | ssp370 MRI-ESM2-0 | 18 | 0 | +46.4 | -15,801 | 3.2e+04 (2100) | ON TRACK | 0 | exit 0 |
| C005 | ssp126 CESM2-WACCM | 64 | 0 | +196.3 | -220,215 | 1.1e+05 (2161) | ON TRACK | 0 | exit 0 |
| C006 | ssp126 MRI-ESM2-0 | 63 | 0 | +139.1 | +6,356 | 2.2e+05 (2217) | ON TRACK | 0 | exit 0 |
| C007 | ssp585 CESM2-WACCM | 92 | 4 | -542.6 | -1,137,267 | 4.5e+06 (2127) | OFF TRACK (SMB -1394.5; shelf basal melt 51471.6; dM/dt (post-2016) -3989.6) | 0 | exit 0 |
| C008 | ssp585 MRI-ESM2-0 | 72 | 2 | +32.9 | -854,248 | 1.5e+06 (2149) | OFF TRACK (SMB 24.2; shelf basal melt 11243.8; dM/dt (post-2016) -2996.8) | 0 | exit 0 |
| C009 | control CESM2-WACCM | 56 | 1 | +82.4 | -11,510 | 2.3e+05 (2214) | ON TRACK | 0 | exit 0 |
| C010 | control MRI-ESM2-0 | 54 | 0 | +41.4 | +1,255 | 2.3e+05 (2227) | ON TRACK | 0 | exit 0 |
| C011 | OCX | 5 | 0 | +10.6 | +4,008 | 3.2e+04 (2005) | ON TRACK | 0 | exit 0 |

Attempt A, without the reference (`rehiunoamb`):

| core | experiment | wall (min) | rescued steps | VAF change (mm SLE) | mass change (Gt) | peak speed (m/yr, year) | track audit | checker errors | scalars |
|---|---|---|---|---|---|---|---|---|---|
| C001 | historical CESM2-WACCM | 7 | 0 | +30.3 | +6,699 | 2.2e+04 (2003) | OFF TRACK (dM/dt (post-2016) 577.9) | 0 | exit 0 |
| C002 | historical MRI-ESM2-0 | 6 | 0 | +29.6 | +7,423 | 2.2e+04 (2003) | OFF TRACK (dM/dt (post-2016) 633.7) | 0 | exit 0 |
| C003 | ssp370 CESM2-WACCM | 19 | 0 | +369.4 | -37,061 | 3.3e+03 (2089) | ON TRACK | 0 | exit 0 |
| C004 | ssp370 MRI-ESM2-0 | 18 | 0 | +375.0 | +90,485 | 2.3e+03 (2045) | OFF TRACK (dM/dt (post-2016) 1056.5) | 0 | exit 0 |
| C005 | ssp126 CESM2-WACCM | 62 | 0 | +1291.5 | +195,891 | 3.3e+03 (2145) | OFF TRACK (dM/dt (post-2016) 686.2) | 0 | exit 0 |
| C006 | ssp126 MRI-ESM2-0 | 61 | 0 | +1255.2 | +386,052 | 2.3e+03 (2047) | OFF TRACK (dM/dt (post-2016) 1352.1) | 0 | exit 0 |
| C007 | ssp585 CESM2-WACCM | 88 | 3 | +515.9 | -708,272 | 3.3e+07 (2100) | OFF TRACK (SMB -1028.3; shelf basal melt 46136.6; dM/dt (post-2016) -2487.0) | 0 | exit 0 |
| C008 | ssp585 MRI-ESM2-0 | 82 | 2 | +1244.3 | -360,061 | 8.4e+07 (2116) | OFF TRACK (SMB 267.4; shelf basal melt 9806.0; dM/dt (post-2016) -1265.8) | 0 | exit 0 |
| C009 | control CESM2-WACCM | 60 | 0 | +1138.2 | +357,135 | 2.3e+03 (2046) | OFF TRACK (dM/dt (post-2016) 1250.8) | 0 | exit 0 |
| C010 | control MRI-ESM2-0 | 58 | 0 | +1146.6 | +382,798 | 2.3e+03 (2047) | OFF TRACK (dM/dt (post-2016) 1340.4) | 0 | exit 0 |
| C011 | OCX | 5 | 0 | +78.0 | +21,169 | 2.2e+04 (2003) | OFF TRACK (dM/dt (post-2016) 946.5) | 0 | exit 0 |

The track audit judges against present-day envelopes, so a projection fails
its forced rows by design; the rows that gate a rehearsal, the budget residual
and the discharge runaway, pass on every core. Every core's 31 files pass
isschecker 0.5.1 with 0 errors (core 11 through `isschecker_ocx.py`), and
`compare_scalars.py` exits 0 on all of them.

## Forced response

Projection minus control in VAF, positive for sea-level rise (mm SLE), the
ISMIP6 ctrl_proj convention: the control and the projections of an ESM branch
from the same historical endpoint.

| ESM | scenario | IU A, without | IU B, with | September A | September B |
|---|---|---|---|---|---|
| CESM2-WACCM | ssp126 | -21 / -77 / -152 | -19 / -66 / -113 | -21 / -74 / -140 | -16 / -62 / -109 |
| CESM2-WACCM | ssp370 | -35 | -31 | -35 | -31 |
| CESM2-WACCM | ssp585 | -60 / -35 / +614 | -56 / +5 / +617 | -61 / -127 / +361 | -54 / -96 / +357 |
| MRI-ESM2-0 | ssp126 | -13 / -50 / -109 | -13 / -49 / -98 | -12 / -52 / -114 | -15 / -49 / -104 |
| MRI-ESM2-0 | ssp370 | -37 | -38 | -40 | -39 |
| MRI-ESM2-0 | ssp585 | -46 / -130 / -104 | -44 / -82 / +2 | -48 / -168 / -209 | -46 / -124 / -87 |

Values at 2100 / 2200 / 2300; the ssp370 runs end in 2100. A negative value
is a projection that gains VAF against its control: at 25 km the grounding
line cannot retreat through cells this coarse while the forced SMB rises
(September's finding 7).

Two figures go with this section, written by `plot_r25iu.py` (kept with the
rehearsal's scripts on Quartz) into the artifact tree's `figs/`:
`r25iu_slc.png`, projection minus control in VAF by ESM and attempt with
September's dashed beside IU's, and `r25iu_drift_calving.png`, the CESM2-WACCM
historical then control, VAF change since 2003 and the calving flux out of the
fixed front, for both rehearsals and both attempts.

## Findings

1. **The forced response hardly depends on the reference; the drift does.**
   Through 2100 the two attempts' projection minus control agrees to 4 mm SLE
   in every scenario, and by 2300 they part by 3 mm (CESM2-WACCM ssp585) to
   106 mm (MRI-ESM2-0 ssp585). The controls themselves differ by more than a
   metre: without the reference they gain 1,138 and 1,147 mm SLE of VAF by
   2300 and with it 82 and 41; the historicals gain 30.3 and 29.6 mm in 12
   years against 2.6 and -0.8.
2. **Without the reference the 25 km MAP's flow relaxes within a year, and the
   calving flux falls with it.** Under the fixed front the calving column is
   the ice the transport carries out of the t=0 extent, the flux across the
   t=0 front. Both attempts start at 1,251 Gt/yr. With the reference every
   cell's t=0 tendency is cancelled, the geometry holds (VAF constant to
   0.01 mm through 2004) and the MAP's velocity carries on: 1,245 Gt/yr in
   2015 and about 1,200 through 2300. Without it the calving falls 22 percent
   in the first year, to 777 Gt/yr in 2015 and about 470 by 2300. The front
   band's mass changes by under 1 percent; the speed drops. From
   `front_flux_check.py` on the CESM2-WACCM historical (job 11890115):

   | state | front flux | floating front | grounded front | grounded ice 0 to 5 km from the front | grounding-line flux | band mass (Gt) |
   |---|---|---|---|---|---|---|
   | MAP | 1,231 Gt/yr | 861 | 370 | 195 m/yr | 1,688 | 112,407 |
   | A, 2008 | 833 | 699 | 134 | 67 m/yr | 1,145 | 112,730 |
   | A, 2015 | 777 | 653 | 124 | 45 m/yr | 908 | 113,371 |
   | B, 2015 | 1,245 | 874 | 371 | 208 m/yr | 1,710 | 112,588 |
   | observed velocity on the MAP | 825 | 581 | 243 | 19 m/yr | 1,560 | |

   The MAP's t=0 velocity carries 1.5 times the observed-velocity flux out of
   the front (grounded ice within 5 km of it moves at 195 m/yr against 19
   observed); attempt B keeps that excess for three centuries, and attempt A
   loses it and undershoots, with the grounding-line flux at 58 percent of the
   observed-velocity flux by 2015. September's rehearsal shows the same split
   on a smaller scale (about 1,240 against 840 to 1,000 Gt/yr). Which cells
   drive the slowdown was not measured.
3. **IU against September.** Through 2100 the forced responses agree to 5 mm
   SLE in every scenario and attempt. By 2300 the CESM2-WACCM ssp585 response
   is +614 and +617 mm against September's +361 and +357, and the MRI-ESM2-0
   ssp585 response -104 and +2 against -209 and -87; ssp126 stays within
   13 mm. Without the reference IU's MAP drifts more than September's
   (controls +1,138 and +1,147 mm against +869 and +870); with it less (+82
   and +41 against +127 and +100). These differences belong to the
   configuration as a whole: code, front treatment, FSSA, solver and MAP.
4. **The 25 km MAP's controls are wide, and its fastest vertices are
   unobserved front vertices.** IU's prior (sigma 30, rho 7.5 km) allows about a
   hundred times the deviation of September's sigma 0.3 at the same rho, a
   penalty about 1e4 times weaker (it scales as 1/sigma^2), and at 25 km
   the bilaplacian's 7.5 km length is well under a cell. theta reaches -23 and
   phi -27 on ice, and the MAP runs 3.1e4 m/yr at a front vertex near
   (-2026, 485) km, which tripped both probes. In the matrix the same kind of
   single-vertex spike recurs and passes, as in September (finding 6 there):
   with the reference the controls and ssp126 runs peak at 1.1e5 to 2.3e5 m/yr
   near (-1888, 1023) km and the ssp585 runs at 4.5e6 and 1.5e6 m/yr near
   (-2140, 660) km; without it the ssp585 runs peak at 3.3e7 and 8.4e7 m/yr
   near (2690, -497) and (2607, -440) km. The budget closes through each.
   The friction control reaches the submission too: every core's `strbasemag`
   peaks at 7.8 MPa against the checker's accepted 1 MPa, on 0.05 percent of
   its values, and the long runs add velocity range warnings from the spikes
   (2 to 15 checker warnings a core, none of them errors).
5. **This configuration needed no rescue to speak of.** At most 4 rescued
   steps a core (September: up to 118), no solver switch, and 54 to 92 min a
   286-year core on 8 ranks (September: 77 to 421 min, with two cores resumed
   under scpc_mumps after scpc_gamg stalled). scpc_mumps ran together with
   FSSA, the front treatment and IU's MAP, so the drop in rescues belongs to
   all of them together.
6. **The late ssp585 melt is still mostly demand on cells that hold too little
   ice** (issue 136, September's finding 5): the CESM2-WACCM ssp585 run
   without the reference melts -39,132 Gt/yr at 2186 while the positivity
   clamp returns +39,595.
7. **Tooling.**
   - `core_report.py` imports `icepack2_tools.fssa`, which imports Firedrake
     at module level, so the report tool now runs only in the Firedrake
     environment; the rehearsal's `make_report_r25iu.sh` loads it.
   - `submit.sh` refuses one argument holding several `KEY=VALUE` pairs, the
     guard against zsh-joined submissions, and a `front_flux_check.py` list of
     `label=path` states is such an argument; the rehearsal's
     `ff_r25iu.script` takes the list comma-separated.
   - A non-interactive shell on Quartz reaches no GitHub remote, so the clone
     came from a git bundle of `288382a`.

## The 2 km historicals

Core 1 (CESM2-WACCM, 2003 to 2015) from IU's vf MAP itself, on its own 2 km
mesh (`ISMIP7_MESH=checkpoint`), without and with the reference, on the 25 km
forwards' settings at 32 ranks, with a melt refit on that mesh under the
front-cell rule (run record `calibration-melt-refit-2km-vertex-front`: every
basin roots between -0.709 and +1.220 K, and the forward's callback reproduces
every basin at 1.0000 or 1.0001). The start backdates the grounded ice 12 years
(+928 Gt on 1,287,008 cells). Records `core01-2km-hist-cesm2waccm-r2kiu{noamb,amb}`;
artifacts under `antarctica/results/rehearsal_2km_iu/`, with every yearly
checkpoint kept.

| run | wall | VAF change (mm SLE) | dM/dt in 2014 (Gt/yr) | calving at step 1, step 2 and in 2014 (Gt/yr) | track audit |
|---|---|---|---|---|---|
| 2 km, with the reference | 5 h 20 min, 38.5 s a solve | +3.4 | +118 | 1,704, 1,532, 1,496 | ON TRACK |
| 2 km, without | 8 h 18 min, 60.6 s a solve | +15.4 | +338 | 1,668, 1,428, 1,263 | ON TRACK |
| 25 km, with the reference | 5 min | +2.6 | +75 | 1,251, 1,247, 1,245 | ON TRACK |
| 25 km, without | 7 min | +30.3 | +835 | 1,251, 1,175, 775 | OFF TRACK (dM/dt) |

From `front_flux_check.py` on the 2 km states (job 11896233), fluxes in Gt/yr:

| state | front flux | floating front | floating u.n (m/yr) | grounded front | grounded ice 0 to 5 km from the front | floating ice 0 to 5 km in | grounding-line flux | band mass (Gt) |
|---|---|---|---|---|---|---|---|---|
| MAP | 1,493 | 1,313 | 219.5 | 180 | 53 m/yr | 432 m/yr | 2,256 | 6,845 |
| with the reference, 2004 | 1,531 | 1,340 | 224.0 | 191 | 55 m/yr | 446 m/yr | 2,451 | 6,845 |
| with the reference, 2015 | 1,496 | 1,300 | 219.5 | 196 | 54 m/yr | 435 m/yr | 2,478 | 6,874 |
| without, 2004 | 1,323 | 1,220 | 206.1 | 103 | 36 m/yr | 404 m/yr | 2,300 | 6,792 |
| without, 2015 | 1,264 | 1,150 | 200.6 | 114 | 34 m/yr | 385 m/yr | 2,355 | 6,784 |
| observed velocity on the MAP | 1,079 | 998 | 166.0 | 81 | 38 m/yr | 419 m/yr | 2,156 | |

The MAP row is the 2015 geometry; the runs start from the backdated 2003
geometry, which raises the grounding-line flux by about 200 Gt/yr.

1. **Every 2 km step converged on its first direct solve**: 480 of 480 in each
   run, 7.7 and 8.7 Newton iterations on average and 17 at most, no rescue or
   subcycle, scpc_mumps with no clip of the controls from the 12-year
   backdated start.
2. **At 2 km the MAP's velocity near the front is far closer to observed than
   at 25 km.** Grounded ice within 5 km of the front moves at 53 m/yr against
   38 observed (25 km: 195 against 19). The floating speed bands lie within 2
   percent of observed from 5 km in and 3 percent above it in the outer 5 km
   (432 against 419 m/yr). The floating front's own u.n is 219.5 against
   166.0 m/yr, 32 percent above, and the t=0 front carries 1.38 times the
   observed-velocity flux (25 km: 1.49).
3. **Without the reference the 2 km run drifts half as fast as the 25 km run**
   (+15.4 against +30.3 mm SLE in 12 years; dM/dt +338 against +835 Gt/yr in
   2014) and passes the track audit. Its calving falls 12 percent from step 2
   to 2014 (25 km: 34 percent) as the floating ice within 250 km of the front
   runs 7 to 16 percent slower than in the run with the reference by 2015,
   while the grounding-line flux holds (2,300 to 2,355 Gt/yr; 25 km: 1,145 to
   908).
4. **With the reference the calving loses 172 Gt/yr after step 1 and then
   holds.** The first step clears about 4 Gt of sub-1 m film from cells outside
   the t=0 extent, as the relaxation year from this MAP did
   (`inversion-2km-budd-b20k-ef2-vf-relax2014-year`); from step 2 the flux is
   1,532 Gt/yr, 1,496 by 2014.
5. **The Shirase front patch** that the MAP carries at 1.4e5 m/yr stays there
   for twelve years with the reference, its geometry held, and relaxes within
   three steps without it (3.0e4 to 1.6e4 m/yr).
6. **Cost.** 38.5 and 60.6 s a diagnostic solve on 32 ranks, 25.7 and 40.4 min
   a model year at dt 0.025: a 286-year core at those rates would take about 5
   and 8 days.

`plot_r2kiu.py` draws the two pairs side by side into the 2 km tree's `figs/`
(`r2kiu_hist.png`: VAF change, and the calving flux in quarter-year means).

## The relaxed line

On 10 October IU ran the relaxed initial state of issue #162 through the same
historicals: a relaxation year and a 250-iteration re-inversion at 25 km from
the rehearsal MAP (records `inversion-25km-budd-rehearsal-iuvf-relax2014-year`
and `inversion-25km-budd-rehearsal-iuvf-relax2014`), and forwards at 2 km and
1 km from IU's relaxed 2 km MAPs (`inversion-2km-{budd,rc}-b20k-ef2-vf-relax2014`).
Under PR 163 a forward on the relaxed MAP's own mesh took the relaxed
geometry, and on any other mesh the controls alone. Commit `c4c9eb8` (this
PR) makes the controls alone the rule on every mesh: BedMachine's 2015
geometry, rewound to 2003, with the relaxed MAP's friction and fluidity. The
runs from the relaxed geometry stay as evidence; the production suite starts
from the relaxed controls.

Core 1 (CESM2-WACCM), 2003 to 2015, fluxes in Gt/yr:

| run | wall | VAF change (mm SLE) | dM/dt in 2014 | calving at step 1, step 2 and in 2014 | record |
|---|---|---|---|---|---|
| 25 km, rehearsal MAP, without the reference | 7 min | +30.3 | +835 | 1,251, 1,175, 775 | `core01-25km-hist-cesm2waccm-rehiunoamb` |
| 25 km, relaxed geometry, without | 4 min | +30.0 | +843 | 1,002, 998, 769 | `core01-25km-hist-cesm2waccm-rehiurlxnoamb` |
| 25 km, relaxed geometry, with (reference net -401) | 3 min | +2.3 | +61 | 1,002, 1,002, 1,005 | `core01-25km-hist-cesm2waccm-rehiurlxamb` |
| 2 km, vf MAP, without | 8 h 18 min | +15.4 | +338 | 1,668, 1,428, 1,263 | `core01-2km-hist-cesm2waccm-r2kiunoamb` |
| 2 km, relaxed geometry, without | 7 h 07 min | +13.1 | +302 | 1,447, 1,432, 1,307 | `core01-2km-hist-cesm2waccm-r2kiurlxnoamb` |
| 2 km, relaxed geometry, with (reference net +25) | 5 h 55 min | +2.4 | +84 | 1,447, 1,436, 1,429 | `core01-2km-hist-cesm2waccm-r2kiurlxamb` |
| 2 km, relaxed controls, without (production core 1) | 7 h 16 min | +15.5 | +306 | 1,820, 1,514, 1,295 | `core01-2km-hist-cesm2waccm-r2kiurlcnoamb` |
| 1 km, relaxed controls, without (64 ranks) | 5 h 37 min | +15.7 | +219 | 1,889, 1,474, 1,363 | `core01-1km-hist-cesm2waccm-r1kiurlxnoamb` |
| 2 km RC, relaxed geometry, without | 34 h 49 min | +5.4 | +195 | 1,513, 1,466, 1,371 | `core01-2km-hist-cesm2waccm-r2kiurcrlxnoamb` |
| 2 km RC, relaxed geometry, with (reference net +99) | 5 h 13 min | +2.3 | +78 | 1,514, 1,503, 1,501 | `core01-2km-hist-cesm2waccm-r2kiurcrlxamb` |
| 2 km RC, relaxed controls, without | failed before step 1 | | | | `core01-2km-hist-cesm2waccm-r2kiurcrlcnoamb` |

Core 2 (MRI-ESM2-0) without the reference follows core 1 to within 0.7 mm:
+12.4 mm from the 2 km relaxed geometry, +14.8 (production core 2) and +15.1
at 1 km from the relaxed controls, +4.8 for RC from its relaxed geometry
(records `core02-*-r*rl*`). Every run that stepped closed its budget to
resid 0.00 on all 480 rows.

1. **Relaxation removes the first-step calving drop and leaves the drift.**
   Started from its relaxed geometry, the MAP sheds 15 (2 km) and 4 (25 km)
   Gt/yr of calving at step 2 where the unrelaxed MAP shed 240 and 76, and
   still gains 13.1 and 30.0 mm SLE in 12 years against 15.4 and 30.3.
2. **The relaxed controls on BedMachine geometry drift as the unrelaxed MAP
   does**, +15.5 and +14.8 mm at 2 km and +15.7 and +15.1 at 1 km. The
   relaxed geometry itself starts with +0.85 mm SLE of VAF and +72 Gt over
   BedMachine's. The first step on BedMachine geometry calves 1,820 (2 km)
   and 1,889 Gt/yr (1 km), against 1,447 from the relaxed geometry. Step 2's
   direct solve stopped at 200 Newton iterations in both 2 km runs and the
   rescue continuation completed it; every other step converged directly.
3. **Against Smith's dH/dt, 2005 to 2015** (`dh_maps_r2kiu.py`, jobs
   11903960, 11907101 and 11907102; `figs/w2005/` of the 2 km tree), in Gt
   of grounded mass change:

   | run | all grounded cells | Smith-observed cells | Smith | Amundsen | Peninsula | Amery to Wilkes | cell-scale RMS (m) |
   |---|---|---|---|---|---|---|---|
   | 2 km, vf MAP | +4,602 | +3,758 | -776 | -1,288 (Smith -1,320) | +1,309 (-215) | +1,407 (-332) | 5.7 |
   | 2 km, relaxed geometry | +3,858 | +3,024 | -772 | -1,411 (-1,320) | +1,299 (-215) | +1,238 (-332) | 5.4 |
   | 2 km, relaxed controls | +4,504 | +3,685 | -778 | -1,611 (-1,321) | +1,315 (-215) | +1,485 (-332) | 5.8 |
   | 1 km, relaxed controls | +4,294 | +3,192 | -536 | -1,772 (-1,209) | +857 (-150) | +1,525 (-308) | 5.8 |

   Basins are IMBIE2's (basin_numbers_ismip8km_v2.nc) by cell centroid: the
   Peninsula is 11 to 13 and Amery to Wilkes 2 to 4. The basin columns sum
   the model over the cells Smith observes, with Smith over the same cells in
   brackets; Smith covers fewer cells on the 1 km mesh. The
   excess sits in the Peninsula and East Antarctica from Amery to Wilkes at
   both resolutions, with Amundsen's loss reproduced; MRI-ESM2-0 gives the
   same pattern (+3,753 and +3,296 Gt on Smith-observed cells at 2 and 1 km).
4. **Without the reference, regularized Coulomb drifts least and costs
   most.** From its relaxed geometry it gains +5.4 and +4.8 mm, at 258 and
   250 s a solve, about five times Budd's 52 s, and 24 Newton iterations on
   average (181 at the most) against Budd's 9; each of these two cores took
   three links across the 11 to 12 October maintenance. With the reference,
   the RC core took one job (11899126) at 37.0 s a solve and drifted
   +2.3 mm, level with the other runs that use the reference. RC's controls
   on BedMachine geometry never started: the seed from the MAP's velocity
   diverged (||F|| 5.7e8 after 50 iterations) and the continuation from
   n = 1 stalled three times at 200 iterations with ||F|| 1.19e8, so both
   cores failed after 2 h 12 min. The production suite is Budd alone.
5. **1 km costs about what 2 km does per model year at twice the ranks**:
   5 h 37 min and 6 h 37 min on 64 ranks, 40 and 47 s a solve, against 7 h 16
   min and 7 h 44 min on 32 at 2 km. One step at 1 km (step 20) peaks at
   9.4e5 m/yr at (346.9, -1,797.1) km, with 2.0e4 at the steps either side.

## The 2 km production suite

IU's Budd suite at 2 km (tag `r2kiurlcnoamb`), started on 11 October from the
relaxed-controls historicals above: cores 1 and 2 are those runs, every control
and projection of an ESM branches from its historical's 2015 state, and OCX
starts from the MAP. It is production grade; the group chooses between this
line and Rice's buffer-0 line on 20 or 21 October.

| item | value |
|---|---|
| code | `c4c9eb8`, Quartz worktree `run2kp` of the rehearsal clone; wrapper `submit_r2kiu.sh proj rlcnoamb <experiment>` and `ocx rlcnoamb` |
| mesh | `antarctica_5000_2000_buffered20000`, 969,764 vertices and 1,928,285 cells, the MAP's own |
| MAP | IU's relaxed Budd MAP (md5 `0d9315d5`), its controls alone on BedMachine's 2015 geometry rewound to 2003 |
| melt | the 2 km vertex_front refit (sha256 `3d2e5348`) |
| forward | dt 0.025, scpc_mumps on 32 ranks and 96 GiB, `ISMIP7_MAP_CLIP=0`, no reference, fixed front, no collapse, SMB-elevation feedback, lake ice base, FSSA at its default, ISMIP7 output on; five-yearly checkpoints, three kept; 72 h links that resubmit themselves |
| output | `write_ismip7_output.py` at `c4c9eb8` with source_id IU and ism_id icepack2, isschecker 0.5.1, ismip7-scalars 0.1.0 and `compare_scalars.py`, one job per core |

The table describes the controls, the projections and OCX. Cores 1 and 2 ran
the same commit `c4c9eb8` from the worktree `run2kc`, with yearly checkpoints,
all kept, in one link each.

Before the suite IU tried scpc_gamg on the control's first step (record
`test-2km-control-scpc-gamg-probe-r2kiurlcnoamb`): the same 9 Newton
iterations and the same budget as scpc_mumps, but 1,533 s on 32 ranks and
1,804 s on 64 against 76 s. Every core runs scpc_mumps.

Every annual file stores its cells in the same order, the OCX run, both
historicals, the controls and both sides of a restart alike (centroids equal
entrywise, job 12249718), so one overlap operator, built once from OCX's first
year, regrids every core.

| core | experiment | jobs | status |
|---|---|---|---|
| C001, C002 | historicals | 11899766, 11899769 | done (records `core0{1,2}-2km-hist-*-r2kiurlcnoamb`) |
| C003 | ssp370 CESM2-WACCM | 11912184 | done 14 Oct, 53 h in one link |
| C004 | ssp370 MRI-ESM2-0 | 11912188 | done 14 Oct, 61 h in one link |
| C005 | ssp126 CESM2-WACCM | 11912183 | running |
| C006 | ssp126 MRI-ESM2-0 | 11912187 | running |
| C007 | ssp585 CESM2-WACCM | 11912185 | running |
| C008 | ssp585 MRI-ESM2-0 | 11912189, 11981754 | running; moved off a shared node at 2045.4, resumed from the 2045 checkpoint on the same 32 ranks |
| C009 | control CESM2-WACCM | 11912182 | running |
| C010 | control MRI-ESM2-0 | 11912186 | running |
| C011 | OCX | 11912190 | done 12 Oct, 13 h |

C008's first link averaged 89 s a solve over its last 500 steps on a node
shared with other users; IU cancelled it right after its 2045 checkpoint and
resubmitted with the same ranks, memory and solver, and the resumed link runs
52 to 55 s a solve. The resumed series holds no repeated or missing row.

Finished cores, through the writer, the checker and the scalar tool:

| core | VAF change (mm SLE) | resid | rescued steps | track audit | checker | scalars |
|---|---|---|---|---|---|---|
| C001 | +15.5 | 0.00 | 1 | ON TRACK | 0 errors | exit 0 |
| C002 | +14.8 | 0.00 | 1 | ON TRACK | 0 errors | exit 0 |
| C003 | +53.3 | 0.00 | 0 | OFF TRACK (dM/dt -1,462) | 0 errors | exit 0 |
| C004 | +120.3 | 0.00 | 0 | ON TRACK | 0 errors | exit 0 |
| C011 | +40.0 | 0.00 | 1 | ON TRACK | 0 errors | exit 0 |

The projections' VAF change runs from 2015, the historicals' and OCX's from
2003. Each checker log carries the same two warnings: the non-mandatory
variables the model does not write, and `strbasemag` above the accepted 1 MPa
(1.9 to 2.4 MPa on 0.001 to 0.004 percent of values). The scalar tool's slvaf
residual note reads 2.2 (C004) to 8.6 percent (C003) of its largest value.

**ssp585 by ESM.** By 2100 the CESM2-WACCM ssp585 run gains 54 mm SLE less VAF
than MRI-ESM2-0's (+69.8 against +124.1 since 2015), and 47.6 mm of that is
Filchner-Ronne (`dh_pair_r2kiu.py`, job 12250092; figure
`figs/dh585_2015_2100_v2.png` of the production tree). Filchner-Ronne loses
21.4 mm under CESM2-WACCM and gains 26.2 under MRI-ESM2-0; every other basin
differs by under 6 mm. The ice floating there in 2015 thins by 275 m on average
against 48 m, and the basin's grounded area shrinks by 41,000 against 6,000
km². The forcing makes the difference: the basin's mean thermal forcing at 200
to 600 m rises from 0.72 to 2.06 °C (2015 to 2024 against 2095 to 2100) under
CESM2-WACCM and from 0.55 to 0.94 °C under MRI-ESM2-0 (`tf_basins.py`, job
12250097), and total melt in the 2090s is 7,076 against 2,847 Gt/yr while
cumulative SMB differs by 2 percent.

## Cost

The 25 km rehearsal took 62 jobs and 141 core-hours on IU Quartz through the
scalar tool, both attempts included (September: 76 jobs and 358 core-hours):
123 for the matrix, 11 for the re-inversion and the rest in minutes-long
checks and post-processing. Each 286-year core took 54 to 92 min on 8 ranks at
dt 0.025, each historical 5 to 7 min. The output chain took 1 h 32 min an
attempt, most of it the writer on the 286-year cores. The 2 km pair took 5
jobs and 438 core-hours, 8 h 18 min and 5 h 20 min of wall time on 32 ranks
for its two historicals.

## Artifacts

On IU Quartz under `antarctica/results/rehearsal_25km_iu/`:

| path | holds |
|---|---|
| `MANIFEST.md`, `SHA256SUMS` | every job and verdict of the rehearsal, and a checksum of every file |
| `mesh/`, `calibration/`, `maps/`, `mapcheck/` | the September mesh copied in, the refit and its sidecar, the MAP and its timing, the MAP checks and front-flux JSON |
| `results/`, `logs/` | the forwards' timeseries, final and periodic checkpoints and annual files, and every log |
| `rehiunoamb/`, `rehiuamb/` | per attempt: `submission/AIS/RICE/icepack2/` (CORE/C001 to C011 and `params.nc`), `checker/`, `scalars/`, `reports/`, `movies/`, `summary.json`, `audits.json`, `output.json` |
| `figs/` | the two figures of the forced-response section |
| `scripts/` | the wrapper `submit_r25iu.sh` and every job script and tool the rehearsal ran, with `CODE_VERSIONS.txt` |

## Reproducing

Every submission went through `submit_r25iu.sh MODE [noamb|amb]` from a
scratch clone at `288382a`, which also ran the post-processing. The wrapper
spells out every knob, since the runner defaults are the 1 km mesh and
regularized Coulomb:

```
common   ISMIP7_LC=25000 ISMIP7_LC_COARSE=250000 ISMIP7_BUFFER_M=20000 ISMIP7_FRICTION=budd
         ISMIP7_N_FLOW=3.0 ISMIP7_GEOMETRY_SPACE=dg0 ISMIP7_BNDIDS=<the sidecar>
forward  ISMIP7_INVERSION=<the MAP> ISMIP7_MESH=checkpoint ISMIP7_DELTAT_PER_BASIN_NPZ=<the refit>
         ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_mumps ISMIP7_MAP_CLIP=0 ISMIP7_DT=0.025
         ISMIP7_APPARENT_MB=<0 or 1> ISMIP7_FIXED_FRONT=1 ISMIP7_FRACTURE=none
         ISMIP7_SMB_ELEVATION_FEEDBACK=1 ISMIP7_LAKE_ICE_BASE=1
refit    calibrate_deltaT.script, ISMIP7_RASTER_SAMPLE=vertex_front DELTAT_K=6.5e-5
inverse  submit.sh inversion, ISMIP7_WARM_START=<the vf MAP> ISMIP7_WARM_START_STRICT=0
         ISMIP7_RASTER_SAMPLE=vertex_front ISMIP7_PRIOR_FORM=bilaplacian
         ISMIP7_PRIOR_SIGMA_THETA=30 ISMIP7_PRIOR_SIGMA_PHI=30 ISMIP7_PRIOR_RHO=7500
         ISMIP7_GRAD_PRECOND=none ISMIP7_DHDT_WEIGHT=0 ISMIP7_LOG_VEL_WEIGHT=85380.44865839917
         ISMIP7_FLUIDITY_CONTROL=floating ISMIP7_EXACT_FRONT=2 ISMIP7_DRAG_GATE=vertex
         ISMIP7_RC_HVISC_FLOOR=2.5 ISMIP7_LAKE_ICE_BASE=1 ISMIP7_TRANSFER_FILL=extend
         ISMIP7_INVERSION_LINEAR_SOLVER=scpc_mumps ISMIP7_MAXITER=300
matrix   submit.sh projection, 8 ranks, 32 GiB, 12 h: hist_cesm_waccm with ISMIP7_CHAIN_THEN="control
         ssp126_cesm_waccm ssp370_cesm_waccm ssp585_cesm_waccm", the same for MRI-ESM2-0, and ocx
2 km     submit_r2kiu.sh from a worktree of the same clone: ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000,
         ISMIP7_INVERSION=<the vf MAP> ISMIP7_MESH=checkpoint, the 2 km refit, the forward knobs
         above, 32 ranks, 96 GiB, ISMIP7_CHECKPOINT_EVERY_YR=1 ISMIP7_KEEP_CHECKPOINTS=0
```
