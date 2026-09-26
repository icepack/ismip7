# ismip7-scalars against the model's own scalars, 32 km (issue #13)

IU Quartz, 23 September 2026. `scripts/batch_runners/scalar_processing.script`
at `b73d3ff` ran ismip7-scalars 0.1.0 (`3f36eb3`) and then
`scripts/compare_scalars.py` over the three p2 runs (runlog
`core07-32km-ssp585-cesm2waccm-p2-*`), jobs 10595618 to 10595620, about 3 min
each on one core. The trees were written by `write_ismip7_output.py` at
`91bf573` in the issue #12 job 10595165. Each run is its own reference, the
state stamped 2016 (`--hist ssp585 --refyear 2016`), so every sea-level series
starts at zero in 2015. The grids are the ones the downloader names, copied
from `/ISMIP7/Output-Processing/Data` with the Globus web app; their sha256
are in each summary below.

T is the tool, N the model. `compare_scalars.py` splits T - N into the area
factor, the maximum-extent mask, the writer's fill convention and a residual
against the mesh sums, and gates on identities only.

Causes 1 and 3 below were fixed on 24 September 2026, and `params.nc` moved
into the upload; "After the fixes" at the end has the runs that show it.

## What holds

Every gate, in all three runs:

- The tool's expressions, replayed on the same files, give its output exactly
  (worst T - R: 0).
- The grid sums of iareagr, iareafl, tendlicalvf and tendligroundf are the
  mesh sums, and so is tendacabf once weighted by pixel coverage. The
  writer's remap is conservative, and the calving booked by collapse (up to
  -14,173 Gt/yr) arrives whole.
- topg is the same in every year. sla20 equals slg20 (worst 2e-14 m), and
  slg20 - slvaf equals its closed form in the tool's own lim and limnsw.
- The yearly change of lim matches dlithkdt, which ties the state and flux
  time stamps together.

## What differs, at 2300

| | p2_none | p2_mask | p2_mask_front |
|---|---|---|---|
| lim, limnsw, iareagr: T - N | +2.3 to +2.6 % | +2.3 to +2.6 % | +2.3 to +2.6 % |
| iareafl: T - N | +2.6 % | +3.5 % | +3.5 % |
| slvaf N, T (mm) | -373.6, -350.8 | -484.9, -457.9 | -478.6, -451.7 |
| slvaf T - N: area factor, grid VAF (mm) | -16.6, +39.4 | -16.9, +43.9 | -16.9, +43.8 |
| slg20 = sla20 N, T (mm) | -324.9, -301.6 | -435.6, -408.2 | -429.3, -402.0 |
| tendacabf N, T (Gt/yr) | -9,492, -10,395 | -9,492, -10,395 | -9,492, -10,395 |
| tendlibmassbffl N, T (Gt/yr) | -123,427, -43,340 | -124,591, -4,647 | -124,561, -5,368 |
| tendligroundf T - N | +1.4 % | +1.1 % | +1.1 % |

Four causes account for all of it.

1. **The area factor.** The tool weights every pixel by af2 = (1/k)^2 of
   EPSG:3031, and the model sums map-plane cell areas. That is +2.2 to +2.6 %
   on the state scalars and -16 to -17 mm on slvaf. For the fluxes it depends on
   latitude, about +1 % on the grounding-line flux and up to +4 % on the calving
   of a collapse year. It is the whole of T - N for iareagr, iareafl,
   tendlicalvf and tendligroundf. (issue 97, closed on 24 September)
2. **Volume above flotation on the 8 km grid.** Before the area factor, the
   tool's sea level sits +39 to +44 mm above the mesh's by 2300, about 10 % of
   the signal. A diagnostic regridded the mesh integrand max(lithk - hf, 0) on
   grounded cells with the writer's own overlap operator. It set that against
   the same integrand on the submitted pixel means, and the drift sits in the
   partly grounded pixels. In a pixel that straddles the grounding line, the
   means of thickness and bed let the flotation deficit of its floating or
   ice-free part cancel the VAF of its grounded part. The deficit grows as the
   shelves thin away, and the number of such pixels grows with it. Grid minus
   mesh VAF, in Gt, for p2_none:

   | year | total | domain edge | grounding line | interior | GL pixels |
   |---|---|---|---|---|---|
   | 2015 | +1,595 | -1,933 | -5,801 | +9,329 | 4,455 |
   | 2100 | -1,305 | -2,172 | -8,423 | +9,290 | 4,510 |
   | 2200 | -8,080 | -2,116 | -14,680 | +8,716 | 5,257 |
   | 2300 | -12,686 | -1,531 | -18,730 | +7,575 | 6,364 |

   The change from 2015 to 2300, -14,281 Gt, is the +39.4 mm, and 90 % of it
   is in the grounding-line pixels (p2_mask: -15,898 Gt, 90 %). The 2300 total
   is the comparison's limnsw residual to the gigatonne. The interior term
   barely moves (+9,329 to +7,575 Gt). (issue #99)
3. **The writer's fill conventions.** acabf is a mean over the covered part
   of a pixel and libmassbffl over the part that floats at year end. The tool
   sums both over whole pixels, which adds -840 Gt/yr of SMB at the domain
   edge in 2300. For melt the partly floating pixels add -30,933 Gt/yr in
   p2_none and -3,300 to -3,800 Gt/yr in the mask modes. Forum thread 50 has
   the organisers reading a flux as a mass change per unit horizontal cell
   area. (issue 96, closed on 24 September)
4. **The model books the fluxes requested of every cell.** Melt and SMB are
   booked before the positivity limiter, on cells with no ice to take them.
   tendacabf is identical in all three modes, collapse or not, because it is
   the forcing over the whole mesh. From the yearly mesh files, in Gt/yr:

   | | year | SMB total | SMB on ice-free cells | melt total | melt off floating ice |
   |---|---|---|---|---|---|
   | p2_none | 2015 | +2,483 | +6 | -1,163 | -14 |
   | p2_none | 2200 | -2,715 | -2,916 | -47,318 | -38,508 |
   | p2_none | 2300 | -9,492 | -7,800 | -123,427 | -111,231 |
   | p2_mask | 2300 | -9,492 | -8,206 | -124,591 | -123,333 |

   The melt booked off floating ice is the comparison's tendlibmassbffl
   residual exactly. The scalar bounds in isschecker 0.5.1's table are
   ±1e9 kg/s for the tend* series, and 2300's -3.9e9 kg/s is four times that;
   0.5.1 does not range-check scalars.
   Pull request 95 (`55ab22b`) books the applied flux after the limiter,
   and the p4 pair below is the check.

The p2 runs predate `2a3a58e`, which books ligroundf in both directions, so
they exercise the chain and say nothing about the current booking.

## The p4 pair: the applied-flux booking

The issue 12 session ran a 32 km ctrl and ssp585 (tag p4) on its fixes
branch. That branch books the SMB and melt the transport applied (`55ab22b`)
and writes cells within 1 cm of flotation as grounded (`ae44194`); the same
branch's writer regridded them. Jobs 10597334 and 10597335 at `739cdd6` ran
the same comparison, and every gate holds. The near-flotation rule moved
floating area to grounded in 41 (ctrl) and 14 (ssp585) years, and the gate on
the two areas' sum absorbs it. Melt and SMB in Gt/yr, with "off the mask" the
melt the model books outside the writer's year-end floating mask:

| | year | tendlibmassbffl N | T | fill | off the mask | tendacabf N | T |
|---|---|---|---|---|---|---|---|
| p2_none | 2015 | -1,163 | -1,699 | -534 | -14 | +2,483 | +2,577 |
| p2_none | 2300 | -123,427 | -43,340 | -30,933 | -111,231 | -9,492 | -10,395 |
| p4 ssp585 | 2015 | -1,153 | -1,698 | -534 | -4 | +2,483 | +2,577 |
| p4 ssp585 | 2300 | -3,929 | -4,978 | -2,995 | -1,986 | -2,043 | -2,153 |
| p4 ctrl | 2300 | -1,047 | -1,592 | -545 | -9 | +2,458 | +2,547 |

Booking the applied flux takes the melt off the mask at 2300 from 111,231 to
1,986 Gt/yr, and what remains is melt the model applied to cells that are
not floating at year end. The pixel convention then carries most of T - N:
the tool's shelf melt is 47 % above the model's in 2015 (-1,698 against
-1,153 Gt/yr) and 50 % above it throughout the control. Sea level is the p2
run's to a tenth of a millimetre, since the booking is output only. The
control's grid VAF term stays small, +3.4 mm at 2300 against -378.9 mm, where
its shelves persist.

The sections below are the job summaries as written, with two Quartz paths
shortened: `<scratch>` is the scratch directory of the IU account that ran
them and `<checkout>` the shared checkout.

## Scalar comparison: p2_none_ssp585_C007

- submission: `<scratch>/ismip7_issue12/p2_none/AIS/RICE/icepack2/CORE/C007`
- tool output: `<scratch>/ismip7_issue13/p2_none/tool/nc/AIS/RICE/icepack2/CORE/C007`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2300 (286)
- pixel coverage: `<checkout>/antarctica/results/ssp585_cesm2_waccm_p2_none_32000_ismip7_annual.h5.overlap.npz`
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit b73d3ff
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -2.000e+00 against 4.036e+00 at tendacabf 2152 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | -4.889e+05 against 1.176e+06 at iareafl 2160 |
| tendacabf with coverage against N | pass | 5.310e+01 against 1.570e+02 at 2203 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | -2.159e-14 against 1.000e-06 at 2170 |
| slg20 - slvaf identity | pass | 4.968e-15 against 1.000e-09 at 2224 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | 1.039e+13 against 2.364e+13 at 2103 |
| lim residual sign | pass | 2.935e+12 against 2.652e+13 at 2073 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor (af2), largest share of max |N|: lim +2.57%, limnsw +2.58%, iareagr +2.32%, iareafl +2.34%, tendacabf -0.77%, tendlibmassbffl -0.19%, tendlicalvf -1.13%, tendligroundf +1.27%, slvaf -3.48%, slg20 -3.48%, sla20 -3.48%
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention, tendacabf: the tool's whole-pixel sum minus the coverage-weighted one is -8.63% of max |N| (2289)
- fill convention, tendlibmassbffl: the tool's whole-pixel sum minus the coverage-weighted one is -25.64% of max |N| (2289)
- tendlibmassbffl: the model's value carries -127167.9 Gt/yr booked off ice that floats at year end (2289), since it books the melt requested of every cell, ice-free ones included
- slvaf: residual +8.28% of max |N| (2297), above 2%
- slg20: residual +8.82% of max |N| (2297), above 2%
- sla20: residual +8.82% of max |N| (2297), above 2%

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.3691e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2036 (worst) | 2.3690e+19 | 2.4298e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2050 | 2.3690e+19 | 2.4298e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2100 | 2.3644e+19 | 2.4251e+19 | +2.561% | +2.561% | +0.000% | +0.000% | -0.000% |
| 2200 | 2.3315e+19 | 2.3912e+19 | +2.521% | +2.521% | +0.000% | +0.000% | -0.000% |
| 2300 | 2.3070e+19 | 2.3662e+19 | +2.501% | +2.502% | +0.000% | +0.000% | -0.000% |

#### limnsw (kg)

max |N| = 2.0835e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 (worst) | 2.0661e+19 | 2.1195e+19 | +2.560% | +2.552% | +0.000% | +0.000% | +0.008% |
| 2050 | 2.0683e+19 | 2.1216e+19 | +2.558% | +2.554% | +0.000% | +0.000% | +0.004% |
| 2100 | 2.0729e+19 | 2.1261e+19 | +2.552% | +2.559% | +0.000% | +0.000% | -0.006% |
| 2200 | 2.0831e+19 | 2.1359e+19 | +2.534% | +2.573% | +0.000% | +0.000% | -0.039% |
| 2300 | 2.0797e+19 | 2.1322e+19 | +2.520% | +2.581% | +0.000% | +0.000% | -0.061% |

#### iareagr (m2)

max |N| = 1.2078e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2077e+13 | 1.2357e+13 | +2.316% | +2.316% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.2076e+13 | 1.2356e+13 | +2.317% | +2.317% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.2070e+13 | 1.2350e+13 | +2.318% | +2.318% | +0.000% | +0.000% | -0.000% |
| 2107 (worst) | 1.2070e+13 | 1.2350e+13 | +2.319% | +2.319% | +0.000% | +0.000% | -0.000% |
| 2200 | 1.1854e+13 | 1.2130e+13 | +2.286% | +2.286% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.1528e+13 | 1.1798e+13 | +2.236% | +2.236% | +0.000% | +0.000% | -0.000% |

#### iareafl (m2)

max |N| = 1.4262e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4261e+12 | 1.4595e+12 | +2.338% | +2.338% | +0.000% | +0.000% | +0.000% |
| 2028 (worst) | 1.4262e+12 | 1.4595e+12 | +2.340% | +2.340% | +0.000% | +0.000% | -0.000% |
| 2050 | 1.4238e+12 | 1.4571e+12 | +2.335% | +2.335% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.3701e+12 | 1.4032e+12 | +2.327% | +2.327% | +0.000% | +0.000% | -0.000% |
| 2200 | 7.3472e+11 | 7.5905e+11 | +1.706% | +1.706% | +0.000% | +0.000% | +0.000% |
| 2300 | 1.6128e+11 | 1.6541e+11 | +0.289% | +0.289% | +0.000% | +0.000% | +0.000% |

#### tendacabf (kg s-1)

max |N| = 3.2158e+08

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 7.8686e+07 | 8.1663e+07 | +0.926% | +0.226% | +0.000% | +0.699% | -0.000% |
| 2050 | 8.1233e+07 | 8.4112e+07 | +0.895% | +0.176% | +0.000% | +0.719% | -0.000% |
| 2100 | 9.2626e+07 | 9.4499e+07 | +0.583% | +0.430% | +0.000% | +0.153% | -0.000% |
| 2200 | -8.6043e+07 | -1.0142e+08 | -4.781% | -0.044% | +0.000% | -4.737% | +0.000% |
| 2289 (worst) | -3.2158e+08 | -3.5179e+08 | -9.395% | -0.766% | +0.000% | -8.629% | -0.000% |
| 2300 | -3.0077e+08 | -3.2939e+08 | -8.900% | -0.622% | +0.000% | -8.279% | +0.000% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 286 years

#### tendlibmassbffl (kg s-1)

max |N| = 4.4761e+09

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.6856e+07 | -5.3847e+07 | -0.380% | -0.011% | +0.000% | -0.378% | +0.010% |
| 2050 | -5.3456e+07 | -7.4826e+07 | -0.477% | -0.018% | +0.000% | -0.489% | +0.030% |
| 2100 | -1.9834e+08 | -2.4773e+08 | -1.103% | -0.056% | +0.000% | -1.752% | +0.704% |
| 2200 | -1.4994e+09 | -8.1765e+08 | +15.231% | -0.078% | +0.000% | -11.952% | +27.261% |
| 2289 (worst) | -4.4761e+09 | -1.6007e+09 | +64.239% | -0.143% | +0.000% | -25.645% | +90.027% |
| 2300 | -3.9112e+09 | -1.3734e+09 | +56.697% | -0.149% | +0.000% | -21.899% | +78.745% |

#### tendlicalvf (kg s-1)

max |N| = 3.5241e+04

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -1.7758e+04 | -1.8032e+04 | -0.775% | -0.775% | +0.000% | +0.000% | +0.000% |
| 2050 | -2.7627e+04 | -2.7968e+04 | -0.968% | -0.968% | +0.000% | +0.000% | -0.000% |
| 2081 (worst) | -2.9600e+04 | -2.9999e+04 | -1.135% | -1.135% | +0.000% | +0.000% | -0.000% |
| 2100 | -2.4412e+04 | -2.4747e+04 | -0.952% | -0.952% | +0.000% | +0.000% | -0.000% |
| 2200 | -1.7086e+04 | -1.7381e+04 | -0.835% | -0.835% | +0.000% | +0.000% | -0.000% |
| 2300 | -2.9820e+03 | -3.0591e+03 | -0.219% | -0.219% | +0.000% | +0.000% | +0.000% |

#### tendlifmassbf (kg s-1)

zero on both sides, 286 years

#### tendligroundf (kg s-1)

max |N| = 6.9158e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.7015e+07 | 6.7769e+07 | +1.090% | +1.090% | +0.000% | +0.000% | -0.000% |
| 2050 | 6.7432e+07 | 6.8198e+07 | +1.107% | +1.107% | +0.000% | +0.000% | -0.000% |
| 2100 | 6.7569e+07 | 6.8407e+07 | +1.211% | +1.211% | +0.000% | +0.000% | -0.000% |
| 2176 (worst) | 5.6114e+07 | 5.6989e+07 | +1.266% | +1.266% | +0.000% | +0.000% | -0.000% |
| 2200 | 4.9400e+07 | 4.9981e+07 | +0.840% | +0.840% | +0.000% | +0.000% | +0.000% |
| 2300 | 4.8463e+07 | 4.9120e+07 | +0.950% | +0.950% | +0.000% | +0.000% | +0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -58.291 | -57.365 | +0.925 | -1.162 | +0.000 | +0.000 | +2.088 |
| 2100 | -186.180 | -181.905 | +4.275 | -3.715 | +0.000 | +0.000 | +7.989 |
| 2200 | -467.230 | -452.291 | +14.939 | -11.776 | +0.000 | +0.000 | +26.715 |
| 2293 (worst) | -389.327 | -366.209 | +23.119 | -16.358 | +0.000 | +0.000 | +39.477 |
| 2300 | -373.603 | -350.820 | +22.783 | -16.626 | +0.000 | +0.000 | +39.409 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -57.067 | -56.163 | +0.903 | -1.135 | +0.000 | +0.000 | +2.039 |
| 2100 | -179.003 | -174.758 | +4.245 | -3.558 | +0.000 | +0.000 | +7.803 |
| 2200 | -432.172 | -416.903 | +15.269 | -10.823 | +0.000 | +0.000 | +26.092 |
| 2293 (worst) | -341.037 | -317.480 | +23.557 | -14.999 | +0.000 | +0.000 | +38.556 |
| 2300 | -324.869 | -301.635 | +23.234 | -15.256 | +0.000 | +0.000 | +38.490 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -57.067 | -56.163 | +0.903 | -1.135 | +0.000 | +0.000 | +2.039 |
| 2100 | -179.003 | -174.758 | +4.245 | -3.558 | +0.000 | +0.000 | +7.803 |
| 2200 | -432.172 | -416.903 | +15.269 | -10.823 | +0.000 | +0.000 | +26.092 |
| 2293 (worst) | -341.037 | -317.480 | +23.557 | -14.999 | +0.000 | +0.000 | +38.556 |
| 2300 | -324.869 | -301.635 | +23.234 | -15.256 | +0.000 | +0.000 | +38.490 |

## Scalar comparison: p2_mask_ssp585_C007

- submission: `<scratch>/ismip7_issue12/p2_mask/AIS/RICE/icepack2/CORE/C007`
- tool output: `<scratch>/ismip7_issue13/p2_mask/tool/nc/AIS/RICE/icepack2/CORE/C007`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2300 (286)
- pixel coverage: `<checkout>/antarctica/results/ssp585_cesm2_waccm_p2_mask_32000_ismip7_annual.h5.overlap.npz`
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit b73d3ff
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -4.000e+00 against 8.068e+00 at tendligroundf 2076 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | 5.061e+01 against 1.164e+02 at tendlicalvf 2188 |
| tendacabf with coverage against N | pass | 5.310e+01 against 1.570e+02 at 2203 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | -2.542e-14 against 1.000e-06 at 2147 |
| slg20 - slvaf identity | pass | 5.440e-15 against 1.000e-09 at 2183 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | 1.024e+13 against 2.323e+13 at 2225 |
| lim residual sign | pass | 3.969e+12 against 2.653e+13 at 2030 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor (af2), largest share of max |N|: lim +2.57%, limnsw +2.58%, iareagr +2.32%, iareafl +2.48%, tendacabf -0.77%, tendlibmassbffl -0.20%, tendlicalvf -4.09%, tendligroundf +1.26%, slvaf -2.94%, slg20 -2.92%, sla20 -2.92%
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention, tendacabf: the tool's whole-pixel sum minus the coverage-weighted one is -8.63% of max |N| (2289)
- fill convention, tendlibmassbffl: the tool's whole-pixel sum minus the coverage-weighted one is -4.64% of max |N| (2202)
- tendlibmassbffl: the model's value carries -140550.2 Gt/yr booked off ice that floats at year end (2289), since it books the melt requested of every cell, ice-free ones included
- slvaf: residual +7.71% of max |N| (2282), above 2%
- slg20: residual +8.14% of max |N| (2282), above 2%
- sla20: residual +8.14% of max |N| (2282), above 2%

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.3690e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2049 (worst) | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | +0.000% |
| 2050 | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | +0.000% |
| 2100 | 2.3614e+19 | 2.4220e+19 | +2.562% | +2.562% | +0.000% | +0.000% | -0.000% |
| 2200 | 2.3300e+19 | 2.3896e+19 | +2.516% | +2.516% | +0.000% | +0.000% | -0.000% |
| 2300 | 2.3101e+19 | 2.3694e+19 | +2.501% | +2.502% | +0.000% | +0.000% | -0.000% |

#### limnsw (kg)

max |N| = 2.0870e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 (worst) | 2.0661e+19 | 2.1195e+19 | +2.555% | +2.548% | +0.000% | +0.000% | +0.008% |
| 2050 | 2.0683e+19 | 2.1216e+19 | +2.554% | +2.550% | +0.000% | +0.000% | +0.004% |
| 2100 | 2.0734e+19 | 2.1265e+19 | +2.542% | +2.554% | +0.000% | +0.000% | -0.012% |
| 2200 | 2.0862e+19 | 2.1388e+19 | +2.517% | +2.568% | +0.000% | +0.000% | -0.050% |
| 2300 | 2.0837e+19 | 2.1361e+19 | +2.508% | +2.577% | +0.000% | +0.000% | -0.069% |

#### iareagr (m2)

max |N| = 1.2078e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2077e+13 | 1.2357e+13 | +2.316% | +2.316% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.2076e+13 | 1.2356e+13 | +2.317% | +2.317% | +0.000% | +0.000% | -0.000% |
| 2084 (worst) | 1.2077e+13 | 1.2357e+13 | +2.318% | +2.318% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.2073e+13 | 1.2352e+13 | +2.318% | +2.318% | +0.000% | +0.000% | +0.000% |
| 2200 | 1.1862e+13 | 1.2139e+13 | +2.288% | +2.288% | +0.000% | +0.000% | +0.000% |
| 2300 | 1.1532e+13 | 1.1802e+13 | +2.237% | +2.237% | +0.000% | +0.000% | -0.000% |

#### iareafl (m2)

max |N| = 1.4262e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4261e+12 | 1.4595e+12 | +2.338% | +2.338% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.4106e+12 | 1.4441e+12 | +2.348% | +2.348% | +0.000% | +0.000% | -0.000% |
| 2100 | 1.2201e+12 | 1.2550e+12 | +2.450% | +2.450% | +0.000% | +0.000% | +0.000% |
| 2125 (worst) | 1.0578e+12 | 1.0932e+12 | +2.484% | +2.484% | +0.000% | +0.000% | -0.000% |
| 2200 | 4.6519e+11 | 4.8259e+11 | +1.221% | +1.220% | +0.000% | +0.000% | +0.000% |
| 2300 | 3.6185e+10 | 3.7451e+10 | +0.089% | +0.089% | +0.000% | +0.000% | +0.000% |

#### tendacabf (kg s-1)

max |N| = 3.2158e+08

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 7.8686e+07 | 8.1663e+07 | +0.926% | +0.226% | +0.000% | +0.699% | -0.000% |
| 2050 | 8.1233e+07 | 8.4112e+07 | +0.895% | +0.176% | +0.000% | +0.719% | -0.000% |
| 2100 | 9.2626e+07 | 9.4499e+07 | +0.583% | +0.430% | +0.000% | +0.153% | -0.000% |
| 2200 | -8.6043e+07 | -1.0142e+08 | -4.781% | -0.044% | +0.000% | -4.737% | +0.000% |
| 2289 (worst) | -3.2158e+08 | -3.5179e+08 | -9.395% | -0.766% | +0.000% | -8.629% | -0.000% |
| 2300 | -3.0077e+08 | -3.2939e+08 | -8.900% | -0.622% | +0.000% | -8.279% | +0.000% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 286 years

#### tendlibmassbffl (kg s-1)

max |N| = 4.4965e+09

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.6856e+07 | -5.3847e+07 | -0.378% | -0.011% | +0.000% | -0.377% | +0.010% |
| 2050 | -5.3483e+07 | -7.3928e+07 | -0.455% | -0.018% | +0.000% | -0.480% | +0.044% |
| 2100 | -1.8712e+08 | -2.0710e+08 | -0.444% | -0.068% | +0.000% | -1.452% | +1.076% |
| 2200 | -1.4841e+09 | -2.8222e+08 | +26.730% | -0.136% | +0.000% | -3.736% | +30.602% |
| 2289 (worst) | -4.4965e+09 | -1.4515e+08 | +96.772% | -0.064% | +0.000% | -2.214% | +99.049% |
| 2300 | -3.9481e+09 | -1.4726e+08 | +84.528% | -0.058% | +0.000% | -2.330% | +86.916% |

#### tendlicalvf (kg s-1)

max |N| = 4.1115e+08

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -1.7758e+04 | -1.8032e+04 | -0.000% | -0.000% | +0.000% | +0.000% | +0.000% |
| 2050 | -2.0384e+07 | -2.0269e+07 | +0.028% | +0.028% | +0.000% | +0.000% | -0.000% |
| 2100 | -3.0739e+07 | -3.0167e+07 | +0.139% | +0.139% | +0.000% | +0.000% | -0.000% |
| 2200 | -6.2767e+07 | -6.2979e+07 | -0.051% | -0.051% | +0.000% | +0.000% | -0.000% |
| 2237 (worst) | -4.1115e+08 | -4.2797e+08 | -4.092% | -4.092% | +0.000% | +0.000% | +0.000% |
| 2300 | -4.9485e+07 | -4.9734e+07 | -0.061% | -0.061% | +0.000% | +0.000% | +0.000% |

#### tendlifmassbf (kg s-1)

zero on both sides, 286 years

#### tendligroundf (kg s-1)

max |N| = 6.8061e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.7015e+07 | 6.7769e+07 | +1.108% | +1.108% | +0.000% | +0.000% | -0.000% |
| 2050 | 6.7232e+07 | 6.7998e+07 | +1.126% | +1.126% | +0.000% | +0.000% | +0.000% |
| 2088 (worst) | 6.6711e+07 | 6.7571e+07 | +1.263% | +1.263% | +0.000% | +0.000% | -0.000% |
| 2100 | 6.5420e+07 | 6.6175e+07 | +1.110% | +1.110% | +0.000% | +0.000% | -0.000% |
| 2200 | 4.7724e+07 | 4.8041e+07 | +0.466% | +0.466% | +0.000% | +0.000% | -0.000% |
| 2300 | 5.1191e+07 | 5.1743e+07 | +0.812% | +0.812% | +0.000% | +0.000% | -0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -58.345 | -57.357 | +0.989 | -1.161 | +0.000 | +0.000 | +2.150 |
| 2100 | -201.297 | -193.536 | +7.761 | -3.543 | +0.000 | +0.000 | +11.304 |
| 2200 | -554.014 | -532.224 | +21.790 | -11.495 | +0.000 | +0.000 | +33.285 |
| 2271 (worst) | -544.745 | -516.387 | +28.358 | -15.673 | +0.000 | +0.000 | +44.032 |
| 2300 | -484.859 | -457.891 | +26.968 | -16.913 | +0.000 | +0.000 | +43.881 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -56.968 | -56.005 | +0.963 | -1.137 | +0.000 | +0.000 | +2.100 |
| 2100 | -191.768 | -184.133 | +7.635 | -3.404 | +0.000 | +0.000 | +11.039 |
| 2200 | -515.956 | -493.919 | +22.036 | -10.470 | +0.000 | +0.000 | +32.506 |
| 2270 (worst) | -498.206 | -469.519 | +28.687 | -14.242 | +0.000 | +0.000 | +42.928 |
| 2300 | -435.559 | -408.239 | +27.320 | -15.535 | +0.000 | +0.000 | +42.854 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -56.968 | -56.005 | +0.963 | -1.137 | +0.000 | +0.000 | +2.100 |
| 2100 | -191.768 | -184.133 | +7.635 | -3.404 | +0.000 | +0.000 | +11.039 |
| 2200 | -515.956 | -493.919 | +22.036 | -10.470 | +0.000 | +0.000 | +32.506 |
| 2270 (worst) | -498.206 | -469.519 | +28.687 | -14.242 | +0.000 | +0.000 | +42.928 |
| 2300 | -435.559 | -408.239 | +27.320 | -15.535 | +0.000 | +0.000 | +42.854 |

## Scalar comparison: p2_mask_front_ssp585_C007

- submission: `<scratch>/ismip7_issue12/p2_mask_front/AIS/RICE/icepack2/CORE/C007`
- tool output: `<scratch>/ismip7_issue13/p2_mask_front/tool/nc/AIS/RICE/icepack2/CORE/C007`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2300 (286)
- pixel coverage: `<checkout>/antarctica/results/ssp585_cesm2_waccm_p2_mask_front_32000_ismip7_annual.h5.overlap.npz`
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit b73d3ff
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -4.000e+00 against 8.061e+00 at tendligroundf 2087 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | 5.241e+01 against 1.166e+02 at tendlicalvf 2165 |
| tendacabf with coverage against N | pass | 5.310e+01 against 1.570e+02 at 2203 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | -1.960e-14 against 1.000e-06 at 2256 |
| slg20 - slvaf identity | pass | 5.149e-15 against 1.000e-09 at 2290 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | -1.059e+13 against 2.327e+13 at 2210 |
| lim residual sign | pass | 3.969e+12 against 2.653e+13 at 2030 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor (af2), largest share of max |N|: lim +2.57%, limnsw +2.58%, iareagr +2.32%, iareafl +2.48%, tendacabf -0.77%, tendlibmassbffl -0.20%, tendlicalvf -2.48%, tendligroundf +1.25%, slvaf -2.96%, slg20 -2.94%, sla20 -2.94%
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention, tendacabf: the tool's whole-pixel sum minus the coverage-weighted one is -8.63% of max |N| (2289)
- fill convention, tendlibmassbffl: the tool's whole-pixel sum minus the coverage-weighted one is -4.72% of max |N| (2202)
- tendlibmassbffl: the model's value carries -140296.7 Gt/yr booked off ice that floats at year end (2289), since it books the melt requested of every cell, ice-free ones included
- slvaf: residual +7.79% of max |N| (2282), above 2%
- slg20: residual +8.23% of max |N| (2282), above 2%
- sla20: residual +8.23% of max |N| (2282), above 2%

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.3690e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2049 (worst) | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | +0.000% |
| 2050 | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2100 | 2.3625e+19 | 2.4232e+19 | +2.562% | +2.562% | +0.000% | +0.000% | -0.000% |
| 2200 | 2.3296e+19 | 2.3892e+19 | +2.516% | +2.516% | +0.000% | +0.000% | -0.000% |
| 2300 | 2.3099e+19 | 2.3691e+19 | +2.501% | +2.501% | +0.000% | +0.000% | -0.000% |

#### limnsw (kg)

max |N| = 2.0868e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 (worst) | 2.0661e+19 | 2.1195e+19 | +2.556% | +2.548% | +0.000% | +0.000% | +0.008% |
| 2050 | 2.0683e+19 | 2.1216e+19 | +2.554% | +2.550% | +0.000% | +0.000% | +0.004% |
| 2100 | 2.0733e+19 | 2.1264e+19 | +2.544% | +2.554% | +0.000% | +0.000% | -0.010% |
| 2200 | 2.0860e+19 | 2.1385e+19 | +2.518% | +2.568% | +0.000% | +0.000% | -0.050% |
| 2300 | 2.0835e+19 | 2.1359e+19 | +2.509% | +2.577% | +0.000% | +0.000% | -0.068% |

#### iareagr (m2)

max |N| = 1.2078e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2077e+13 | 1.2357e+13 | +2.316% | +2.316% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.2076e+13 | 1.2356e+13 | +2.317% | +2.317% | +0.000% | +0.000% | -0.000% |
| 2100 | 1.2073e+13 | 1.2353e+13 | +2.318% | +2.318% | +0.000% | +0.000% | -0.000% |
| 2106 (worst) | 1.2073e+13 | 1.2353e+13 | +2.318% | +2.318% | +0.000% | +0.000% | +0.000% |
| 2200 | 1.1860e+13 | 1.2136e+13 | +2.286% | +2.286% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.1529e+13 | 1.1800e+13 | +2.236% | +2.236% | +0.000% | +0.000% | -0.000% |

#### iareafl (m2)

max |N| = 1.4262e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4261e+12 | 1.4595e+12 | +2.338% | +2.338% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.4134e+12 | 1.4469e+12 | +2.348% | +2.348% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.2506e+12 | 1.2856e+12 | +2.453% | +2.453% | +0.000% | +0.000% | -0.000% |
| 2125 (worst) | 1.0601e+12 | 1.0955e+12 | +2.484% | +2.484% | +0.000% | +0.000% | -0.000% |
| 2200 | 4.6642e+11 | 4.8388e+11 | +1.225% | +1.225% | +0.000% | +0.000% | +0.000% |
| 2300 | 3.6937e+10 | 3.8230e+10 | +0.091% | +0.091% | +0.000% | +0.000% | -0.000% |

#### tendacabf (kg s-1)

max |N| = 3.2158e+08

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 7.8686e+07 | 8.1663e+07 | +0.926% | +0.226% | +0.000% | +0.699% | -0.000% |
| 2050 | 8.1233e+07 | 8.4112e+07 | +0.895% | +0.176% | +0.000% | +0.719% | -0.000% |
| 2100 | 9.2626e+07 | 9.4499e+07 | +0.583% | +0.430% | +0.000% | +0.153% | -0.000% |
| 2200 | -8.6043e+07 | -1.0142e+08 | -4.781% | -0.044% | +0.000% | -4.737% | +0.000% |
| 2289 (worst) | -3.2158e+08 | -3.5179e+08 | -9.395% | -0.766% | +0.000% | -8.629% | -0.000% |
| 2300 | -3.0077e+08 | -3.2939e+08 | -8.900% | -0.622% | +0.000% | -8.279% | +0.000% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 286 years

#### tendlibmassbffl (kg s-1)

max |N| = 4.4882e+09

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.6856e+07 | -5.3847e+07 | -0.379% | -0.011% | +0.000% | -0.377% | +0.010% |
| 2050 | -5.3560e+07 | -7.4565e+07 | -0.468% | -0.018% | +0.000% | -0.485% | +0.035% |
| 2100 | -1.8939e+08 | -2.1409e+08 | -0.550% | -0.069% | +0.000% | -1.473% | +0.991% |
| 2200 | -1.4839e+09 | -2.8506e+08 | +26.710% | -0.136% | +0.000% | -3.800% | +30.646% |
| 2289 (worst) | -4.4882e+09 | -1.5341e+08 | +96.582% | -0.058% | +0.000% | -2.413% | +99.053% |
| 2300 | -3.9471e+09 | -1.7010e+08 | +84.153% | -0.065% | +0.000% | -2.679% | +86.897% |

#### tendlicalvf (kg s-1)

max |N| = 4.4912e+08

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -1.7758e+04 | -1.8032e+04 | -0.000% | -0.000% | +0.000% | +0.000% | +0.000% |
| 2050 | -1.1888e+07 | -1.1764e+07 | +0.028% | +0.028% | +0.000% | +0.000% | -0.000% |
| 2100 | -4.4700e+07 | -4.3777e+07 | +0.206% | +0.206% | +0.000% | +0.000% | +0.000% |
| 2200 | -5.1568e+07 | -5.1935e+07 | -0.082% | -0.082% | +0.000% | +0.000% | -0.000% |
| 2268 (worst) | -2.6941e+08 | -2.8054e+08 | -2.479% | -2.479% | +0.000% | +0.000% | -0.000% |
| 2300 | -4.8258e+07 | -4.8477e+07 | -0.049% | -0.049% | +0.000% | +0.000% | +0.000% |

#### tendlifmassbf (kg s-1)

zero on both sides, 286 years

#### tendligroundf (kg s-1)

max |N| = 6.8612e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.7015e+07 | 6.7769e+07 | +1.099% | +1.099% | +0.000% | +0.000% | -0.000% |
| 2050 | 6.7325e+07 | 6.8092e+07 | +1.117% | +1.117% | +0.000% | +0.000% | +0.000% |
| 2088 (worst) | 6.7260e+07 | 6.8118e+07 | +1.251% | +1.251% | +0.000% | +0.000% | +0.000% |
| 2100 | 6.7960e+07 | 6.8729e+07 | +1.122% | +1.122% | +0.000% | +0.000% | -0.000% |
| 2200 | 4.7256e+07 | 4.7566e+07 | +0.452% | +0.452% | +0.000% | +0.000% | +0.000% |
| 2300 | 5.0869e+07 | 5.1415e+07 | +0.795% | +0.795% | +0.000% | +0.000% | +0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -58.321 | -57.387 | +0.934 | -1.162 | +0.000 | +0.000 | +2.096 |
| 2100 | -198.015 | -191.383 | +6.632 | -3.540 | +0.000 | +0.000 | +10.172 |
| 2200 | -546.929 | -525.020 | +21.909 | -11.395 | +0.000 | +0.000 | +33.303 |
| 2282 (worst) | -520.802 | -492.569 | +28.233 | -16.157 | +0.000 | +0.000 | +44.390 |
| 2300 | -478.623 | -451.705 | +26.918 | -16.871 | +0.000 | +0.000 | +43.788 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -56.966 | -56.057 | +0.910 | -1.137 | +0.000 | +0.000 | +2.046 |
| 2100 | -189.303 | -182.773 | +6.530 | -3.404 | +0.000 | +0.000 | +9.934 |
| 2200 | -508.768 | -486.609 | +22.159 | -10.365 | +0.000 | +0.000 | +32.524 |
| 2282 (worst) | -471.918 | -443.353 | +28.565 | -14.787 | +0.000 | +0.000 | +43.351 |
| 2300 | -429.294 | -402.023 | +27.271 | -15.492 | +0.000 | +0.000 | +42.764 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -56.966 | -56.057 | +0.910 | -1.137 | +0.000 | +0.000 | +2.046 |
| 2100 | -189.303 | -182.773 | +6.530 | -3.404 | +0.000 | +0.000 | +9.934 |
| 2200 | -508.768 | -486.609 | +22.159 | -10.365 | +0.000 | +0.000 | +32.524 |
| 2282 (worst) | -471.918 | -443.353 | +28.565 | -14.787 | +0.000 | +0.000 | +43.351 |
| 2300 | -429.294 | -402.023 | +27.271 | -15.492 | +0.000 | +0.000 | +42.764 |

## Scalar comparison: p4_ctrl_C009

- submission: `<scratch>/ismip7_issue12/fixed/p4/AIS/RICE/icepack2/CORE/C009`
- tool output: `<scratch>/ismip7_issue13/p4_ctrl/tool/nc/AIS/RICE/icepack2/CORE/C009`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2300 (286)
- pixel coverage: `<checkout>/antarctica/results/ctrl2015_cesm2_waccm_p4_32000_ismip7_annual.h5.overlap.npz`
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit 739cdd6
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -4.000e+00 against 8.075e+00 at tendligroundf 2241 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | -5.319e+06 against 1.511e+07 at grounded area gained 2139 |
| tendacabf with coverage against N | pass | -3.718e+00 against 8.747e+01 at 2215 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | 2.037e-14 against 1.000e-06 at 2292 |
| slg20 - slvaf identity | pass | -4.897e-15 against 1.000e-09 at 2057 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | -9.805e+12 against 2.371e+13 at 2072 |
| lim residual sign | pass | 3.799e+12 against 2.655e+13 at 2057 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor (af2), largest share of max |N|: lim +2.56%, limnsw +2.57%, iareagr +2.32%, iareafl +2.35%, tendacabf +0.97%, tendlibmassbffl -0.91%, tendlicalvf -1.23%, tendligroundf +1.17%, slvaf -1.99%, slg20 -1.99%, sla20 -1.99%
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention, tendacabf: the tool's whole-pixel sum minus the coverage-weighted one is +2.66% of max |N| (2300)
- fill convention, tendlibmassbffl: the tool's whole-pixel sum minus the coverage-weighted one is -52.11% of max |N| (2299)
- near flotation: 41 years write floating area as grounded, at most 1.152e+09 m2 (2151)
- tendlibmassbffl: the model's value carries -13.8 Gt/yr outside the writer's year-end floating mask (2213)

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.3808e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.3688e+19 | 2.4296e+19 | +2.553% | +2.553% | +0.000% | +0.000% | +0.000% |
| 2050 | 2.3702e+19 | 2.4310e+19 | +2.554% | +2.554% | +0.000% | +0.000% | -0.000% |
| 2100 | 2.3723e+19 | 2.4331e+19 | +2.556% | +2.556% | +0.000% | +0.000% | -0.000% |
| 2200 | 2.3765e+19 | 2.4375e+19 | +2.559% | +2.559% | +0.000% | +0.000% | -0.000% |
| 2300 (worst) | 2.3808e+19 | 2.4418e+19 | +2.563% | +2.563% | +0.000% | +0.000% | -0.000% |

#### limnsw (kg)

max |N| = 2.0799e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.0661e+19 | 2.1195e+19 | +2.564% | +2.556% | +0.000% | +0.000% | +0.008% |
| 2050 | 2.0680e+19 | 2.1213e+19 | +2.565% | +2.558% | +0.000% | +0.000% | +0.006% |
| 2100 | 2.0705e+19 | 2.1239e+19 | +2.566% | +2.561% | +0.000% | +0.000% | +0.006% |
| 2200 | 2.0752e+19 | 2.1286e+19 | +2.569% | +2.565% | +0.000% | +0.000% | +0.004% |
| 2300 (worst) | 2.0799e+19 | 2.1334e+19 | +2.571% | +2.570% | +0.000% | +0.000% | +0.002% |

#### iareagr (m2)

max |N| = 1.2086e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2077e+13 | 1.2357e+13 | +2.315% | +2.315% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.2078e+13 | 1.2358e+13 | +2.316% | +2.316% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.2085e+13 | 1.2366e+13 | +2.325% | +2.319% | +0.000% | +0.000% | +0.006% |
| 2151 (worst) | 1.2075e+13 | 1.2356e+13 | +2.326% | +2.316% | +0.000% | +0.000% | +0.010% |
| 2200 | 1.2074e+13 | 1.2354e+13 | +2.315% | +2.315% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.2074e+13 | 1.2354e+13 | +2.314% | +2.314% | +0.000% | +0.000% | -0.000% |

#### iareafl (m2)

max |N| = 1.4261e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4261e+12 | 1.4594e+12 | +2.338% | +2.338% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.4223e+12 | 1.4556e+12 | +2.329% | +2.329% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.4117e+12 | 1.4437e+12 | +2.248% | +2.297% | +0.000% | +0.000% | -0.049% |
| 2200 | 1.4151e+12 | 1.4484e+12 | +2.337% | +2.337% | +0.000% | +0.000% | +0.000% |
| 2268 (worst) | 1.4137e+12 | 1.4472e+12 | +2.348% | +2.348% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.4109e+12 | 1.4443e+12 | +2.341% | +2.341% | +0.000% | +0.000% | -0.000% |

#### tendacabf (kg s-1)

max |N| = 7.7877e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 7.7877e+07 | 8.0705e+07 | +3.631% | +0.970% | +0.000% | +2.662% | -0.000% |
| 2050 | 7.7877e+07 | 8.0705e+07 | +3.631% | +0.970% | +0.000% | +2.662% | -0.000% |
| 2100 | 7.7877e+07 | 8.0705e+07 | +3.631% | +0.970% | +0.000% | +2.662% | -0.000% |
| 2200 | 7.7877e+07 | 8.0705e+07 | +3.631% | +0.970% | +0.000% | +2.662% | -0.000% |
| 2300 (worst) | 7.7877e+07 | 8.0705e+07 | +3.631% | +0.970% | +0.000% | +2.662% | -0.000% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 286 years

#### tendlibmassbffl (kg s-1)

max |N| = 3.3263e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.3069e+07 | -4.9455e+07 | -49.261% | -0.885% | +0.000% | -48.772% | +0.396% |
| 2050 | -3.2987e+07 | -4.9410e+07 | -49.375% | -0.884% | +0.000% | -48.949% | +0.458% |
| 2100 | -3.2872e+07 | -4.9446e+07 | -49.827% | -0.861% | +0.000% | -49.539% | +0.574% |
| 2200 | -3.3150e+07 | -5.0076e+07 | -50.887% | -0.899% | +0.000% | -50.822% | +0.835% |
| 2299 (worst) | -3.3201e+07 | -5.0548e+07 | -52.149% | -0.897% | +0.000% | -52.106% | +0.854% |
| 2300 | -3.3189e+07 | -5.0459e+07 | -51.919% | -0.888% | +0.000% | -51.910% | +0.879% |

#### tendlicalvf (kg s-1)

max |N| = 3.2724e+04

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -1.9721e+04 | -2.0003e+04 | -0.861% | -0.861% | +0.000% | +0.000% | +0.000% |
| 2050 | -2.4211e+04 | -2.4553e+04 | -1.046% | -1.046% | +0.000% | +0.000% | -0.000% |
| 2100 | -3.1519e+04 | -3.1865e+04 | -1.057% | -1.057% | +0.000% | +0.000% | -0.000% |
| 2200 | -2.9860e+04 | -3.0259e+04 | -1.219% | -1.219% | +0.000% | +0.000% | -0.000% |
| 2205 (worst) | -3.1547e+04 | -3.1949e+04 | -1.228% | -1.228% | +0.000% | +0.000% | -0.000% |
| 2300 | -2.8915e+04 | -2.9311e+04 | -1.209% | -1.209% | +0.000% | +0.000% | +0.000% |

#### tendlifmassbf (kg s-1)

zero on both sides, 286 years

#### tendligroundf (kg s-1)

max |N| = 6.7960e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.2242e+07 | 6.2889e+07 | +0.951% | +0.951% | +0.000% | +0.000% | -0.000% |
| 2050 | 6.2915e+07 | 6.3569e+07 | +0.961% | +0.961% | +0.000% | +0.000% | +0.000% |
| 2100 | 6.4576e+07 | 6.5284e+07 | +1.041% | +1.041% | +0.000% | +0.000% | +0.000% |
| 2170 (worst) | 6.7372e+07 | 6.8167e+07 | +1.171% | +1.171% | +0.000% | +0.000% | +0.000% |
| 2200 | 6.6658e+07 | 6.7421e+07 | +1.123% | +1.123% | +0.000% | +0.000% | -0.000% |
| 2300 | 6.7775e+07 | 6.8566e+07 | +1.165% | +1.165% | +0.000% | +0.000% | +0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -50.702 | -51.021 | -0.319 | -1.021 | +0.000 | +0.000 | +0.702 |
| 2100 | -120.082 | -121.335 | -1.254 | -2.454 | +0.000 | +0.000 | +1.200 |
| 2200 | -249.433 | -252.080 | -2.647 | -4.959 | +0.000 | +0.000 | +2.312 |
| 2300 (worst) | -378.899 | -383.065 | -4.165 | -7.557 | +0.000 | +0.000 | +3.392 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -50.410 | -50.740 | -0.330 | -1.015 | +0.000 | +0.000 | +0.686 |
| 2100 | -119.523 | -120.792 | -1.269 | -2.442 | +0.000 | +0.000 | +1.173 |
| 2200 | -248.595 | -251.277 | -2.683 | -4.940 | +0.000 | +0.000 | +2.258 |
| 2300 (worst) | -377.745 | -381.960 | -4.216 | -7.528 | +0.000 | +0.000 | +3.312 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -50.410 | -50.740 | -0.330 | -1.015 | +0.000 | +0.000 | +0.686 |
| 2100 | -119.523 | -120.792 | -1.269 | -2.442 | +0.000 | +0.000 | +1.173 |
| 2200 | -248.595 | -251.277 | -2.683 | -4.940 | +0.000 | +0.000 | +2.258 |
| 2300 (worst) | -377.745 | -381.960 | -4.216 | -7.528 | +0.000 | +0.000 | +3.312 |

## Scalar comparison: p4_ssp585_C007

- submission: `<scratch>/ismip7_issue12/fixed/p4/AIS/RICE/icepack2/CORE/C007`
- tool output: `<scratch>/ismip7_issue13/p4_ssp585/tool/nc/AIS/RICE/icepack2/CORE/C007`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2300 (286)
- pixel coverage: `<checkout>/antarctica/results/ssp585_cesm2_waccm_p4_32000_ismip7_annual.h5.overlap.npz`
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit 739cdd6
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -2.000e+00 against 4.029e+00 at tendacabf 2240 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | -5.174e+06 against 1.337e+07 at grounded area gained 2259 |
| tendacabf with coverage against N | pass | -5.112e+00 against 3.812e+01 at 2194 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | -1.898e-14 against 1.000e-06 at 2134 |
| slg20 - slvaf identity | pass | 5.173e-15 against 1.000e-09 at 2071 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | 1.040e+13 against 2.364e+13 at 2103 |
| lim residual sign | pass | 2.028e+12 against 2.652e+13 at 2074 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor (af2), largest share of max |N|: lim +2.57%, limnsw +2.58%, iareagr +2.32%, iareafl +2.34%, tendacabf +2.14%, tendlibmassbffl -2.67%, tendlicalvf -1.13%, tendligroundf +1.30%, slvaf -3.48%, slg20 -3.48%, sla20 -3.48%
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention, tendacabf: the tool's whole-pixel sum minus the coverage-weighted one is -4.00% of max |N| (2259)
- fill convention, tendlibmassbffl: the tool's whole-pixel sum minus the coverage-weighted one is -51.00% of max |N| (2203)
- near flotation: 14 years write floating area as grounded, at most 1.076e+09 m2 (2144)
- tendlibmassbffl: the model's value carries -2441.8 Gt/yr outside the writer's year-end floating mask (2283)
- slvaf: residual +8.27% of max |N| (2297), above 2%
- slg20: residual +8.82% of max |N| (2297), above 2%
- sla20: residual +8.82% of max |N| (2297), above 2%

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.3691e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2036 (worst) | 2.3690e+19 | 2.4298e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2050 | 2.3690e+19 | 2.4298e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2100 | 2.3644e+19 | 2.4251e+19 | +2.561% | +2.561% | +0.000% | +0.000% | -0.000% |
| 2200 | 2.3315e+19 | 2.3912e+19 | +2.521% | +2.521% | +0.000% | +0.000% | -0.000% |
| 2300 | 2.3070e+19 | 2.3662e+19 | +2.501% | +2.502% | +0.000% | +0.000% | -0.000% |

#### limnsw (kg)

max |N| = 2.0835e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 (worst) | 2.0661e+19 | 2.1195e+19 | +2.560% | +2.552% | +0.000% | +0.000% | +0.008% |
| 2050 | 2.0683e+19 | 2.1216e+19 | +2.558% | +2.554% | +0.000% | +0.000% | +0.004% |
| 2100 | 2.0729e+19 | 2.1261e+19 | +2.552% | +2.559% | +0.000% | +0.000% | -0.006% |
| 2200 | 2.0831e+19 | 2.1359e+19 | +2.534% | +2.573% | +0.000% | +0.000% | -0.039% |
| 2300 | 2.0797e+19 | 2.1322e+19 | +2.520% | +2.581% | +0.000% | +0.000% | -0.061% |

#### iareagr (m2)

max |N| = 1.2078e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2077e+13 | 1.2357e+13 | +2.316% | +2.316% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.2076e+13 | 1.2356e+13 | +2.317% | +2.317% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.2070e+13 | 1.2350e+13 | +2.318% | +2.318% | +0.000% | +0.000% | -0.000% |
| 2144 (worst) | 1.2029e+13 | 1.2310e+13 | +2.324% | +2.315% | +0.000% | +0.000% | +0.009% |
| 2200 | 1.1854e+13 | 1.2130e+13 | +2.286% | +2.286% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.1528e+13 | 1.1798e+13 | +2.236% | +2.236% | +0.000% | +0.000% | -0.000% |

#### iareafl (m2)

max |N| = 1.4262e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4261e+12 | 1.4595e+12 | +2.338% | +2.338% | +0.000% | +0.000% | +0.000% |
| 2016 (worst) | 1.4260e+12 | 1.4593e+12 | +2.338% | +2.338% | +0.000% | +0.000% | -0.000% |
| 2050 | 1.4238e+12 | 1.4571e+12 | +2.335% | +2.335% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.3701e+12 | 1.4032e+12 | +2.327% | +2.327% | +0.000% | +0.000% | -0.000% |
| 2200 | 7.3468e+11 | 7.5900e+11 | +1.706% | +1.706% | +0.000% | +0.000% | +0.000% |
| 2300 | 1.6128e+11 | 1.6541e+11 | +0.289% | +0.289% | +0.000% | +0.000% | +0.000% |

#### tendacabf (kg s-1)

max |N| = 9.8191e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 7.8686e+07 | 8.1664e+07 | +3.032% | +0.742% | +0.000% | +2.291% | +0.000% |
| 2050 | 8.1236e+07 | 8.4117e+07 | +2.935% | +0.577% | +0.000% | +2.358% | -0.000% |
| 2100 | 9.3556e+07 | 9.5785e+07 | +2.270% | +1.389% | +0.000% | +0.882% | -0.000% |
| 2200 | 2.3996e+06 | -1.2703e+06 | -3.738% | -0.203% | +0.000% | -3.535% | +0.000% |
| 2273 (worst) | -7.2732e+07 | -7.6963e+07 | -4.309% | -0.438% | +0.000% | -3.872% | +0.000% |
| 2300 | -6.4723e+07 | -6.8223e+07 | -3.565% | -0.067% | +0.000% | -3.498% | -0.000% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 286 years

#### tendlibmassbffl (kg s-1)

max |N| = 2.5420e+08

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.6534e+07 | -5.3816e+07 | -6.799% | -0.194% | +0.000% | -6.652% | +0.047% |
| 2050 | -5.2311e+07 | -7.4606e+07 | -8.771% | -0.325% | +0.000% | -8.538% | +0.093% |
| 2100 | -1.6483e+08 | -2.2567e+08 | -23.934% | -1.133% | +0.000% | -24.024% | +1.223% |
| 2150 (worst) | -2.4137e+08 | -3.3581e+08 | -37.151% | -2.313% | +0.000% | -41.644% | +6.807% |
| 2200 | -1.9415e+08 | -2.7460e+08 | -31.651% | -1.839% | +0.000% | -45.335% | +15.523% |
| 2300 | -1.2451e+08 | -1.5774e+08 | -13.070% | -0.495% | +0.000% | -37.336% | +24.761% |

#### tendlicalvf (kg s-1)

max |N| = 3.5241e+04

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -1.7758e+04 | -1.8032e+04 | -0.775% | -0.775% | +0.000% | +0.000% | +0.000% |
| 2050 | -2.7626e+04 | -2.7967e+04 | -0.968% | -0.968% | +0.000% | +0.000% | +0.000% |
| 2081 (worst) | -2.9600e+04 | -2.9999e+04 | -1.135% | -1.135% | +0.000% | +0.000% | +0.000% |
| 2100 | -2.4412e+04 | -2.4747e+04 | -0.952% | -0.952% | +0.000% | +0.000% | -0.000% |
| 2200 | -1.7086e+04 | -1.7381e+04 | -0.835% | -0.835% | +0.000% | +0.000% | +0.000% |
| 2300 | -2.9820e+03 | -3.0591e+03 | -0.219% | -0.219% | +0.000% | +0.000% | +0.000% |

#### tendlifmassbf (kg s-1)

zero on both sides, 286 years

#### tendligroundf (kg s-1)

max |N| = 6.5987e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.2253e+07 | 6.2899e+07 | +0.979% | +0.979% | +0.000% | +0.000% | -0.000% |
| 2050 | 6.3325e+07 | 6.4009e+07 | +1.037% | +1.037% | +0.000% | +0.000% | +0.000% |
| 2100 | 6.4722e+07 | 6.5491e+07 | +1.164% | +1.164% | +0.000% | +0.000% | +0.000% |
| 2176 (worst) | 5.5413e+07 | 5.6271e+07 | +1.301% | +1.301% | +0.000% | +0.000% | +0.000% |
| 2200 | 4.9180e+07 | 4.9759e+07 | +0.877% | +0.877% | +0.000% | +0.000% | +0.000% |
| 2300 | 4.8429e+07 | 4.9086e+07 | +0.996% | +0.996% | +0.000% | +0.000% | +0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -58.291 | -57.365 | +0.925 | -1.162 | +0.000 | +0.000 | +2.088 |
| 2100 | -186.180 | -181.905 | +4.275 | -3.715 | +0.000 | +0.000 | +7.990 |
| 2200 | -467.200 | -452.279 | +14.921 | -11.776 | +0.000 | +0.000 | +26.697 |
| 2293 (worst) | -389.297 | -366.172 | +23.125 | -16.357 | +0.000 | +0.000 | +39.482 |
| 2300 | -373.573 | -350.785 | +22.789 | -16.625 | +0.000 | +0.000 | +39.414 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -57.067 | -56.164 | +0.903 | -1.135 | +0.000 | +0.000 | +2.039 |
| 2100 | -179.003 | -174.758 | +4.245 | -3.558 | +0.000 | +0.000 | +7.803 |
| 2200 | -432.142 | -416.890 | +15.252 | -10.823 | +0.000 | +0.000 | +26.074 |
| 2293 (worst) | -341.006 | -317.443 | +23.563 | -14.999 | +0.000 | +0.000 | +38.561 |
| 2300 | -324.839 | -301.599 | +23.240 | -15.255 | +0.000 | +0.000 | +38.495 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -57.067 | -56.164 | +0.903 | -1.135 | +0.000 | +0.000 | +2.039 |
| 2100 | -179.003 | -174.758 | +4.245 | -3.558 | +0.000 | +0.000 | +7.803 |
| 2200 | -432.142 | -416.890 | +15.252 | -10.823 | +0.000 | +0.000 | +26.074 |
| 2293 (worst) | -341.006 | -317.443 | +23.563 | -14.999 | +0.000 | +0.000 | +38.561 |
| 2300 | -324.839 | -301.599 | +23.240 | -15.255 | +0.000 | +0.000 | +38.495 |

## After the fixes, 24 September 2026

Three commits on `claude/scalar-output-fixes` answer causes 1 and 3 and ship
`params.nc`. The sections above stay as they were measured.

- `677b469` (issue 96). Every gridded flux is a mean over the whole pixel of
  what the cells under it booked, and the fill policy says only where it is
  fill. `libmassbffl` keeps the melt of every cell in a pixel that still
  floats at year end. The flux files carry `flux_pixel_mean = whole_pixel`.
- `7d03212` (issue 97). The forward integrates the ten scalars over true
  area, each cell's map-plane area times af2 at its centroid. The writer
  checks every year of the scalars CSV against the annual files and stamps
  the scalar files `scalar_area`, and `compare_scalars.py --native-af2` has
  to agree with the stamp.
- `df1a2c0` (issue 98). `ismip7-scalars-set-params` writes `params.nc` into
  the upload, beside `CORE/`, and the tool reads it there with no
  `--params-path`.

IU Quartz jobs at `df1a2c0`, from a scratch clone of the branch:

| job | what | result |
|---|---|---|
| 10604339 | the p4 pair regridded again, `params.nc` into the tree, isschecker 0.5.1, af2 against the organisers' grid | 0 errors in every test group for C009 and C007; af2 within 5.95e-8 at all 579,121 pixel centres |
| 10604341 | the new comparison over the issue 12 trees, with and without `--native-af2` | `comparison.csv` byte-identical to jobs 10597334 and 10597335; the switch on map-plane scalars fails, exit 1 |
| 10604342, 10604343 | scalar processing over the regridded p4 pair | exit 0, every gate holds, fill term zero |
| 10604612 | the p4 control for five years (tag p5), scalars over true area | 5 annual files; `resid` 0.0000; timeseries equal to p4's |
| 10604344 | p5 regridded, isschecker 0.5.1 | scalars stamped `true_area`; errors only in the 93 length checks a five-year run fails |
| 10604347, 10604348 | scalar processing over p5, with the switch and without | exit 0 with the area term zero; exit 2 |

### Whole-pixel means

Melt and SMB in Gt/yr. T before is jobs 10597334 and 10597335, T after jobs
10604342 and 10604343. "Left out" is the melt the model books in pixels with
no floating ice at year end, which the `no_floating_ice` fill policy leaves as
fill:

| | year | tendlibmassbffl N | T before | T after | left out | tendacabf N | T before | T after |
|---|---|---|---|---|---|---|---|---|
| p4 ctrl | 2015 | -1,044 | -1,561 | -1,048 | -2 | +2,458 | +2,547 | +2,482 |
| p4 ctrl | 2300 | -1,047 | -1,592 | -1,048 | -6 | +2,458 | +2,547 | +2,482 |
| p4 ssp585 | 2015 | -1,153 | -1,698 | -1,163 | -2 | +2,483 | +2,577 | +2,507 |
| p4 ssp585 | 2150 | -7,617 | -10,597 | -7,399 | -381 | +1,667 | +1,645 | +1,719 |
| p4 ssp585 | 2300 | -3,929 | -4,978 | -2,181 | -1,771 | -2,043 | -2,153 | -2,045 |

The fill term is zero and every gate holds, the tendacabf gate included,
which a whole-pixel tree passes without the overlap cache; what T - N keeps
for tendacabf is the area factor. The control's shelf melt as the tool reads
it is the model's to 0.4 %, where it read 50 % high. In the ssp585 run T
falls below N from about 2150, by the melt booked in pixels with no floating
ice at year end: 2 Gt/yr in 2015, 1,771 in 2300 and 2,223 at most, in 2283.
The request fills such a pixel under any pixel convention, and isschecker
0.5.1 reports a value there as an error. Keeping the melt of every cell in a
pixel that still floats holds 216 Gt/yr more at 2300 than the year-end
floating cells alone would (1,986, the "off the mask" of the p4 section). The
writer's summary line reports the same numbers as the comparison's residual,
and the near-flotation rule accounts for at most 0.6 Gt/yr of them (2.0 in
the control). What that melt is follows below.

isschecker 0.5.1 finds 0 errors in every test group on the regridded pair,
with `AIS/RICE/icepack2/params.nc` in the tree. Its warnings are the
non-mandatory variables and, for the ssp585, `libmassbffl` below
-0.008 kg m-2 s-1 in 0.00545 % of values (0.0382 % before) and `strbasemag`
above 1e6 Pa in 0.000813 %. The control's `libmassbffl` warning (0.00767 %
before) is gone.

### What the left-out melt is

Jobs 10604946 and 10605000 read the p4 ssp585 annual files with the writer's
overlap cache and split each year's left-out melt by the state, at the start
of that year, of the cells that booked it. The script, `leftout_melt.py`, ran
from the scratch clone at `df1a2c0` and is attached to issue 105. Its totals
match the writer's to the digit. Melt in Gt over 2016 to 2300, and in Gt/yr
at 2300:

| cells that booked it | 2016 to 2300 | 2300 |
|---|---|---|
| floating at the start of the year, no floating ice at its end | 18,196 | 87 |
| grounded at the start, afloat and emptied within the year | 711 | 0 |
| no ice (1 m or less) at both ends of the year | 180,082 | 1,684 |
| all left out | 198,990 | 1,771 |

Only the first row is shelf ice that melted away during the year. It carries
the peak: in 2283 a single 15,384 km² cell at 84.3°S, 142.9°W melted through
46 m of ice and holds 27 % of that year's left-out melt.

The cells of the third row cover 1.30 million km² by 2300. The positivity
limiter lets a cell lose only the ice it held when a step began, so a cell
with no ice books as melt its share of whatever arrives in the step. Their
thickness budgets say what arrives:

| arriving in those cells | 2016 to 2300 | 2300 |
|---|---|---|
| the frozen apparent-MB reference | 134,363 | 1,373 |
| ice flowing in | 72,394 | 640 |
| of which across the grounding line | 66,198 | 627 |
| taken by their negative SMB | -26,671 | -329 |

Under a pinned front the reference is never cleared, so a shelf cell that melt
has emptied keeps receiving it and books it as melt in the same step. Shared
in proportion with the ice flowing in, the reference accounts for about
117,000 Gt of the left-out melt over the run (59 %) and 1,150 Gt/yr at 2300,
29 % of the native `tendlibmassbffl`. That is melt of ice the model's state
never holds. The fill drops it from the tool's sum except in pixels that
still float, where the whole-pixel mean keeps it as part of the 216 Gt/yr
above, and the native scalar keeps all of it. The rest, about 620 Gt/yr at
2300, is real melt the gridded field cannot carry: ice that crosses the
grounding line into an empty cell and melts on arrival, and the shelf ice of
the first row. Whether the production runs keep the reference is a group
decision (issue 104), the melt it books is issue 105, and where the
submission reports the real melt is issue 109.

Note of 25 September 2026: the group chose option 3 of issue 109, and the
forward now writes the melt of ice flowing into marine cells holding no ice
at either end of the year as `lifmassbf`, the inflow's share of what those
cells melt; the reference's share stays in `libmassbffl`. The same split
applied to the p4 annual files moves 131, 451 and 571 Gt/yr out of the
pixels the fill blanks in 2150, 2250 and 2300 (against the class-wide 535 at
2300 above) and leaves 250, 978 and 1,200 Gt/yr there: 216, 937 and 1,094 of
them the reference's share on emptied marine cells, and 0.4, 8.7 and 20 the
reference melted on emptied land cells, which the melt law counts as afloat
at zero thickness. The runs in this report predate the change, and
their `lifmassbf` is zero.

### True area

`regrid.area_factor` matches `af2_AIS_08000m_v1.nc` at all 579,121 pixel
centres to 5.95e-8, the file's float32 rounding (mean 2.06e-8). The p5 run
logged af2 from 0.9526 to 1.0567 on the 32 km mesh and wrote five years. Its
timeseries is p4's to the last digit over the same 50 rows, with `resid`
0.0000 throughout, so the mass budget stays map-plane and the physics is
untouched. The writer found every year of its scalars CSV on true area: the
forward's `iareagr` and `iareafl` equal the writer's own true-area mesh sums
to the CSV's seven digits.

The p5 native scalars over p4's in 2015, against the tool's area factor on
p4 (R over the replay with af2 = 1):

| | lim | limnsw | iareagr | iareafl | tendacabf | tendlibmassbffl | tendlicalvf | tendligroundf |
|---|---|---|---|---|---|---|---|---|
| p5 over p4 | +2.5743 % | +2.5815 % | +2.3234 % | +2.3391 % | +1.0022 % | +0.6194 % | +1.4290 % | +1.0385 % |
| tool's af2 on p4 | +2.5659 % | +2.5733 % | +2.3164 % | +2.3381 % | +0.9991 % | +0.6214 % | +1.4284 % | +1.0382 % |

With `--native-af2` the area term is zero and every gate holds. The exact sums
agree to 6.26e-5 of their L1, against an allowance of 5.93e-4, af2's largest
change between neighbouring pixels. The difference is where each side samples
af2. The model takes it at the centroid of each cell, and this mesh's interior
cells reach 320 km, where the curvature of af2 puts the centroid value up to
about 1e-4 above the cell mean. The tool takes it at 8 km pixel centres. The
difference barely moves over the years (lim +2.5743 % in 2015, +2.5742 % in
2019), so it cancels in sea-level changes. Without the switch the stamp
refuses the comparison (exit 2), and the switch on the map-plane p4 scalars
fails the exact sums by 2.27e-2 of their L1 (job 10604341).

The first p5 attempt, job 10604433, lost rank 0 to the job's 16 GiB limit
while loading the RACMO climatology, before any output code ran. The rerun
asked for 32 GiB, and its step peaked at 1.2 GiB per rank, as the p4
control's did.

### params.nc in the upload

`ismip7-scalars-set-params --rhoi 917 --rhow 1024 --rhof 1000 --modelpath
<tree>/AIS` put `params.nc` into both trees. The tool read it there with no
`--params-path` in all four scalar jobs, and the density gate holds.
isschecker reads one set-counter directory, so the file beside `CORE/` never
reaches it.

The sections below are the three new job summaries as written, with the same
two paths shortened and their headings dated.

## Scalar comparison: p4_ctrl_C009, 24 September

- submission: `<scratch>/ismip7_issue96_98/trees/p4/AIS/RICE/icepack2/CORE/C009`
- tool output: `<scratch>/ismip7_issue96_98/scalars/p4_ctrl_C009/tool/nc/AIS/RICE/icepack2/CORE/C009`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2300 (286)
- flux means: whole pixel (`flux_pixel_mean`)
- pixel coverage: none given
- native scalars: over map-plane area; the writer's stamp: map_plane (`<scratch>/ismip7_issue96_98/trees/p4/AIS/RICE/icepack2/CORE/C009`)
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit df1a2c0
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -4.000e+00 against 8.075e+00 at tendligroundf 2241 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | -5.319e+06 against 1.511e+07 at grounded area gained 2139 |
| tendacabf against N | pass | -3.722e+00 against 8.723e+01 at 2215 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | 2.037e-14 against 1.000e-06 at 2292 |
| slg20 - slvaf identity | pass | -4.897e-15 against 1.000e-09 at 2057 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | -9.805e+12 against 2.371e+13 at 2072 |
| lim residual sign | pass | 3.799e+12 against 2.655e+13 at 2057 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor (af2), largest share of max |N|: lim +2.56%, limnsw +2.57%, iareagr +2.32%, iareafl +2.35%, tendacabf +1.00%, tendlibmassbffl -0.62%, tendlicalvf -1.23%, tendligroundf +1.17%, slvaf -1.99%, slg20 -1.99%, sla20 -1.99%
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention: acabf and libmassbffl are whole-pixel means (flux_pixel_mean), the tool's own, so nothing is undone and the fill term is zero
- near flotation: 41 years write floating area as grounded, at most 1.152e+09 m2 (2151)
- tendlibmassbffl: the model's value carries -8.2 Gt/yr in pixels with no floating ice at year end, which the fill leaves out (2213)

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.3808e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.3688e+19 | 2.4296e+19 | +2.553% | +2.553% | +0.000% | +0.000% | +0.000% |
| 2050 | 2.3702e+19 | 2.4310e+19 | +2.554% | +2.554% | +0.000% | +0.000% | -0.000% |
| 2100 | 2.3723e+19 | 2.4331e+19 | +2.556% | +2.556% | +0.000% | +0.000% | -0.000% |
| 2200 | 2.3765e+19 | 2.4375e+19 | +2.559% | +2.559% | +0.000% | +0.000% | -0.000% |
| 2300 (worst) | 2.3808e+19 | 2.4418e+19 | +2.563% | +2.563% | +0.000% | +0.000% | -0.000% |

#### limnsw (kg)

max |N| = 2.0799e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.0661e+19 | 2.1195e+19 | +2.564% | +2.556% | +0.000% | +0.000% | +0.008% |
| 2050 | 2.0680e+19 | 2.1213e+19 | +2.565% | +2.558% | +0.000% | +0.000% | +0.006% |
| 2100 | 2.0705e+19 | 2.1239e+19 | +2.566% | +2.561% | +0.000% | +0.000% | +0.006% |
| 2200 | 2.0752e+19 | 2.1286e+19 | +2.569% | +2.565% | +0.000% | +0.000% | +0.004% |
| 2300 (worst) | 2.0799e+19 | 2.1334e+19 | +2.571% | +2.570% | +0.000% | +0.000% | +0.002% |

#### iareagr (m2)

max |N| = 1.2086e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2077e+13 | 1.2357e+13 | +2.315% | +2.315% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.2078e+13 | 1.2358e+13 | +2.316% | +2.316% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.2085e+13 | 1.2366e+13 | +2.325% | +2.319% | +0.000% | +0.000% | +0.006% |
| 2151 (worst) | 1.2075e+13 | 1.2356e+13 | +2.326% | +2.316% | +0.000% | +0.000% | +0.010% |
| 2200 | 1.2074e+13 | 1.2354e+13 | +2.315% | +2.315% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.2074e+13 | 1.2354e+13 | +2.314% | +2.314% | +0.000% | +0.000% | -0.000% |

#### iareafl (m2)

max |N| = 1.4261e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4261e+12 | 1.4594e+12 | +2.338% | +2.338% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.4223e+12 | 1.4556e+12 | +2.329% | +2.329% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.4117e+12 | 1.4437e+12 | +2.248% | +2.297% | +0.000% | +0.000% | -0.049% |
| 2200 | 1.4151e+12 | 1.4484e+12 | +2.337% | +2.337% | +0.000% | +0.000% | +0.000% |
| 2268 (worst) | 1.4137e+12 | 1.4472e+12 | +2.348% | +2.348% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.4109e+12 | 1.4443e+12 | +2.341% | +2.341% | +0.000% | +0.000% | -0.000% |

#### tendacabf (kg s-1)

max |N| = 7.7877e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 7.7877e+07 | 7.8655e+07 | +0.999% | +0.999% | +0.000% | +0.000% | -0.000% |
| 2050 | 7.7877e+07 | 7.8655e+07 | +0.999% | +0.999% | +0.000% | +0.000% | -0.000% |
| 2100 | 7.7877e+07 | 7.8655e+07 | +0.999% | +0.999% | +0.000% | +0.000% | -0.000% |
| 2126 (worst) | 7.7877e+07 | 7.8655e+07 | +0.999% | +0.999% | +0.000% | +0.000% | -0.000% |
| 2200 | 7.7877e+07 | 7.8655e+07 | +0.999% | +0.999% | +0.000% | +0.000% | -0.000% |
| 2300 | 7.7877e+07 | 7.8655e+07 | +0.999% | +0.999% | +0.000% | +0.000% | -0.000% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 286 years

#### tendlibmassbffl (kg s-1)

max |N| = 3.3263e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.3069e+07 | -3.3202e+07 | -0.397% | -0.616% | +0.000% | +0.000% | +0.219% |
| 2024 (worst) | -3.2855e+07 | -3.2989e+07 | -0.403% | -0.620% | +0.000% | +0.000% | +0.217% |
| 2050 | -3.2987e+07 | -3.3095e+07 | -0.325% | -0.609% | +0.000% | +0.000% | +0.284% |
| 2100 | -3.2872e+07 | -3.2946e+07 | -0.223% | -0.586% | +0.000% | +0.000% | +0.363% |
| 2200 | -3.3150e+07 | -3.3182e+07 | -0.098% | -0.611% | +0.000% | +0.000% | +0.513% |
| 2300 | -3.3189e+07 | -3.3215e+07 | -0.079% | -0.608% | +0.000% | +0.000% | +0.529% |

#### tendlicalvf (kg s-1)

max |N| = 3.2724e+04

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -1.9721e+04 | -2.0003e+04 | -0.861% | -0.861% | +0.000% | +0.000% | +0.000% |
| 2050 | -2.4211e+04 | -2.4553e+04 | -1.046% | -1.046% | +0.000% | +0.000% | -0.000% |
| 2100 | -3.1519e+04 | -3.1865e+04 | -1.057% | -1.057% | +0.000% | +0.000% | -0.000% |
| 2200 | -2.9860e+04 | -3.0259e+04 | -1.219% | -1.219% | +0.000% | +0.000% | -0.000% |
| 2205 (worst) | -3.1547e+04 | -3.1949e+04 | -1.228% | -1.228% | +0.000% | +0.000% | -0.000% |
| 2300 | -2.8915e+04 | -2.9311e+04 | -1.209% | -1.209% | +0.000% | +0.000% | +0.000% |

#### tendlifmassbf (kg s-1)

zero on both sides, 286 years

#### tendligroundf (kg s-1)

max |N| = 6.7960e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.2242e+07 | 6.2889e+07 | +0.951% | +0.951% | +0.000% | +0.000% | -0.000% |
| 2050 | 6.2915e+07 | 6.3569e+07 | +0.961% | +0.961% | +0.000% | +0.000% | +0.000% |
| 2100 | 6.4576e+07 | 6.5284e+07 | +1.041% | +1.041% | +0.000% | +0.000% | +0.000% |
| 2170 (worst) | 6.7372e+07 | 6.8167e+07 | +1.171% | +1.171% | +0.000% | +0.000% | +0.000% |
| 2200 | 6.6658e+07 | 6.7421e+07 | +1.123% | +1.123% | +0.000% | +0.000% | -0.000% |
| 2300 | 6.7775e+07 | 6.8566e+07 | +1.165% | +1.165% | +0.000% | +0.000% | +0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -50.702 | -51.021 | -0.319 | -1.021 | +0.000 | +0.000 | +0.702 |
| 2100 | -120.082 | -121.335 | -1.254 | -2.454 | +0.000 | +0.000 | +1.200 |
| 2200 | -249.433 | -252.080 | -2.647 | -4.959 | +0.000 | +0.000 | +2.312 |
| 2300 (worst) | -378.899 | -383.065 | -4.165 | -7.557 | +0.000 | +0.000 | +3.392 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -50.410 | -50.740 | -0.330 | -1.015 | +0.000 | +0.000 | +0.686 |
| 2100 | -119.523 | -120.792 | -1.269 | -2.442 | +0.000 | +0.000 | +1.173 |
| 2200 | -248.595 | -251.277 | -2.683 | -4.940 | +0.000 | +0.000 | +2.258 |
| 2300 (worst) | -377.745 | -381.960 | -4.216 | -7.528 | +0.000 | +0.000 | +3.312 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -50.410 | -50.740 | -0.330 | -1.015 | +0.000 | +0.000 | +0.686 |
| 2100 | -119.523 | -120.792 | -1.269 | -2.442 | +0.000 | +0.000 | +1.173 |
| 2200 | -248.595 | -251.277 | -2.683 | -4.940 | +0.000 | +0.000 | +2.258 |
| 2300 (worst) | -377.745 | -381.960 | -4.216 | -7.528 | +0.000 | +0.000 | +3.312 |

## Scalar comparison: p4_ssp585_C007, 24 September

- submission: `<scratch>/ismip7_issue96_98/trees/p4/AIS/RICE/icepack2/CORE/C007`
- tool output: `<scratch>/ismip7_issue96_98/scalars/p4_ssp585_C007/tool/nc/AIS/RICE/icepack2/CORE/C007`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2300 (286)
- flux means: whole pixel (`flux_pixel_mean`)
- pixel coverage: none given
- native scalars: over map-plane area; the writer's stamp: map_plane (`<scratch>/ismip7_issue96_98/trees/p4/AIS/RICE/icepack2/CORE/C007`)
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit df1a2c0
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -2.000e+00 against 4.029e+00 at tendacabf 2240 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | -5.174e+06 against 1.337e+07 at grounded area gained 2259 |
| tendacabf against N | pass | -5.108e+00 against 3.769e+01 at 2194 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | -1.898e-14 against 1.000e-06 at 2134 |
| slg20 - slvaf identity | pass | 5.173e-15 against 1.000e-09 at 2071 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | 1.040e+13 against 2.364e+13 at 2103 |
| lim residual sign | pass | 2.028e+12 against 2.652e+13 at 2074 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor (af2), largest share of max |N|: lim +2.57%, limnsw +2.58%, iareagr +2.32%, iareafl +2.34%, tendacabf +2.11%, tendlibmassbffl -2.34%, tendlicalvf -1.13%, tendligroundf +1.30%, slvaf -3.48%, slg20 -3.48%, sla20 -3.48%
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention: acabf and libmassbffl are whole-pixel means (flux_pixel_mean), the tool's own, so nothing is undone and the fill term is zero
- near flotation: 14 years write floating area as grounded, at most 1.076e+09 m2 (2144)
- tendlibmassbffl: the model's value carries -2223.3 Gt/yr in pixels with no floating ice at year end, which the fill leaves out (2283)
- slvaf: residual +8.27% of max |N| (2297), above 2%
- slg20: residual +8.82% of max |N| (2297), above 2%
- sla20: residual +8.82% of max |N| (2297), above 2%

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.3691e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.3688e+19 | 2.4296e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2036 (worst) | 2.3690e+19 | 2.4298e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2050 | 2.3690e+19 | 2.4298e+19 | +2.566% | +2.566% | +0.000% | +0.000% | -0.000% |
| 2100 | 2.3644e+19 | 2.4251e+19 | +2.561% | +2.561% | +0.000% | +0.000% | -0.000% |
| 2200 | 2.3315e+19 | 2.3912e+19 | +2.521% | +2.521% | +0.000% | +0.000% | -0.000% |
| 2300 | 2.3070e+19 | 2.3662e+19 | +2.501% | +2.502% | +0.000% | +0.000% | -0.000% |

#### limnsw (kg)

max |N| = 2.0835e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 (worst) | 2.0661e+19 | 2.1195e+19 | +2.560% | +2.552% | +0.000% | +0.000% | +0.008% |
| 2050 | 2.0683e+19 | 2.1216e+19 | +2.558% | +2.554% | +0.000% | +0.000% | +0.004% |
| 2100 | 2.0729e+19 | 2.1261e+19 | +2.552% | +2.559% | +0.000% | +0.000% | -0.006% |
| 2200 | 2.0831e+19 | 2.1359e+19 | +2.534% | +2.573% | +0.000% | +0.000% | -0.039% |
| 2300 | 2.0797e+19 | 2.1322e+19 | +2.520% | +2.581% | +0.000% | +0.000% | -0.061% |

#### iareagr (m2)

max |N| = 1.2078e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2077e+13 | 1.2357e+13 | +2.316% | +2.316% | +0.000% | +0.000% | +0.000% |
| 2050 | 1.2076e+13 | 1.2356e+13 | +2.317% | +2.317% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.2070e+13 | 1.2350e+13 | +2.318% | +2.318% | +0.000% | +0.000% | -0.000% |
| 2144 (worst) | 1.2029e+13 | 1.2310e+13 | +2.324% | +2.315% | +0.000% | +0.000% | +0.009% |
| 2200 | 1.1854e+13 | 1.2130e+13 | +2.286% | +2.286% | +0.000% | +0.000% | -0.000% |
| 2300 | 1.1528e+13 | 1.1798e+13 | +2.236% | +2.236% | +0.000% | +0.000% | -0.000% |

#### iareafl (m2)

max |N| = 1.4262e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4261e+12 | 1.4595e+12 | +2.338% | +2.338% | +0.000% | +0.000% | +0.000% |
| 2016 (worst) | 1.4260e+12 | 1.4593e+12 | +2.338% | +2.338% | +0.000% | +0.000% | -0.000% |
| 2050 | 1.4238e+12 | 1.4571e+12 | +2.335% | +2.335% | +0.000% | +0.000% | +0.000% |
| 2100 | 1.3701e+12 | 1.4032e+12 | +2.327% | +2.327% | +0.000% | +0.000% | -0.000% |
| 2200 | 7.3468e+11 | 7.5900e+11 | +1.706% | +1.706% | +0.000% | +0.000% | +0.000% |
| 2300 | 1.6128e+11 | 1.6541e+11 | +0.289% | +0.289% | +0.000% | +0.000% | +0.000% |

#### tendacabf (kg s-1)

max |N| = 9.8191e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 7.8686e+07 | 7.9439e+07 | +0.766% | +0.766% | +0.000% | +0.000% | +0.000% |
| 2050 | 8.1236e+07 | 8.1833e+07 | +0.608% | +0.608% | +0.000% | +0.000% | -0.000% |
| 2100 | 9.3556e+07 | 9.4916e+07 | +1.384% | +1.384% | +0.000% | +0.000% | -0.000% |
| 2116 (worst) | 8.7757e+07 | 8.9827e+07 | +2.109% | +2.109% | +0.000% | +0.000% | +0.000% |
| 2200 | 2.3996e+06 | 2.1682e+06 | -0.236% | -0.236% | +0.000% | +0.000% | +0.000% |
| 2300 | -6.4723e+07 | -6.4818e+07 | -0.096% | -0.096% | +0.000% | +0.000% | -0.000% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 286 years

#### tendlibmassbffl (kg s-1)

max |N| = 2.5420e+08

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.6534e+07 | -3.6839e+07 | -0.120% | -0.147% | +0.000% | +0.000% | +0.027% |
| 2050 | -5.2311e+07 | -5.2842e+07 | -0.209% | -0.259% | +0.000% | +0.000% | +0.050% |
| 2100 | -1.6483e+08 | -1.6607e+08 | -0.487% | -0.998% | +0.000% | +0.000% | +0.511% |
| 2200 | -1.9415e+08 | -1.6635e+08 | +10.935% | -1.528% | +0.000% | +0.000% | +12.464% |
| 2283 (worst) | -1.4930e+08 | -7.9575e+07 | +27.431% | -0.285% | +0.000% | +0.000% | +27.716% |
| 2300 | -1.2451e+08 | -6.9111e+07 | +21.795% | -0.278% | +0.000% | +0.000% | +22.073% |

#### tendlicalvf (kg s-1)

max |N| = 3.5241e+04

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -1.7758e+04 | -1.8032e+04 | -0.775% | -0.775% | +0.000% | +0.000% | +0.000% |
| 2050 | -2.7626e+04 | -2.7967e+04 | -0.968% | -0.968% | +0.000% | +0.000% | +0.000% |
| 2081 (worst) | -2.9600e+04 | -2.9999e+04 | -1.135% | -1.135% | +0.000% | +0.000% | +0.000% |
| 2100 | -2.4412e+04 | -2.4747e+04 | -0.952% | -0.952% | +0.000% | +0.000% | -0.000% |
| 2200 | -1.7086e+04 | -1.7381e+04 | -0.835% | -0.835% | +0.000% | +0.000% | +0.000% |
| 2300 | -2.9820e+03 | -3.0591e+03 | -0.219% | -0.219% | +0.000% | +0.000% | +0.000% |

#### tendlifmassbf (kg s-1)

zero on both sides, 286 years

#### tendligroundf (kg s-1)

max |N| = 6.5987e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.2253e+07 | 6.2899e+07 | +0.979% | +0.979% | +0.000% | +0.000% | -0.000% |
| 2050 | 6.3325e+07 | 6.4009e+07 | +1.037% | +1.037% | +0.000% | +0.000% | +0.000% |
| 2100 | 6.4722e+07 | 6.5491e+07 | +1.164% | +1.164% | +0.000% | +0.000% | +0.000% |
| 2176 (worst) | 5.5413e+07 | 5.6271e+07 | +1.301% | +1.301% | +0.000% | +0.000% | +0.000% |
| 2200 | 4.9180e+07 | 4.9759e+07 | +0.877% | +0.877% | +0.000% | +0.000% | +0.000% |
| 2300 | 4.8429e+07 | 4.9086e+07 | +0.996% | +0.996% | +0.000% | +0.000% | +0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -58.291 | -57.365 | +0.925 | -1.162 | +0.000 | +0.000 | +2.088 |
| 2100 | -186.180 | -181.905 | +4.275 | -3.715 | +0.000 | +0.000 | +7.990 |
| 2200 | -467.200 | -452.279 | +14.921 | -11.776 | +0.000 | +0.000 | +26.697 |
| 2293 (worst) | -389.297 | -366.172 | +23.125 | -16.357 | +0.000 | +0.000 | +39.482 |
| 2300 | -373.573 | -350.785 | +22.789 | -16.625 | +0.000 | +0.000 | +39.414 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -57.067 | -56.164 | +0.903 | -1.135 | +0.000 | +0.000 | +2.039 |
| 2100 | -179.003 | -174.758 | +4.245 | -3.558 | +0.000 | +0.000 | +7.803 |
| 2200 | -432.142 | -416.890 | +15.252 | -10.823 | +0.000 | +0.000 | +26.074 |
| 2293 (worst) | -341.006 | -317.443 | +23.563 | -14.999 | +0.000 | +0.000 | +38.561 |
| 2300 | -324.839 | -301.599 | +23.240 | -15.255 | +0.000 | +0.000 | +38.495 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2050 | -57.067 | -56.164 | +0.903 | -1.135 | +0.000 | +0.000 | +2.039 |
| 2100 | -179.003 | -174.758 | +4.245 | -3.558 | +0.000 | +0.000 | +7.803 |
| 2200 | -432.142 | -416.890 | +15.252 | -10.823 | +0.000 | +0.000 | +26.074 |
| 2293 (worst) | -341.006 | -317.443 | +23.563 | -14.999 | +0.000 | +0.000 | +38.561 |
| 2300 | -324.839 | -301.599 | +23.240 | -15.255 | +0.000 | +0.000 | +38.495 |

## Scalar comparison: p5_ctrl_C009, 24 September

- submission: `<scratch>/ismip7_issue96_98/trees/p5/AIS/RICE/icepack2/CORE/C009`
- tool output: `<scratch>/ismip7_issue96_98/scalars/p5_ctrl_C009/tool/nc/AIS/RICE/icepack2/CORE/C009`
- grids: `<checkout>/ISMIP7/Output-Processing/Data/AIS/af2_AIS_08000m_v1.nc` (float32), `<checkout>/ISMIP7/Output-Processing/Data/AIS/maxmask1_AIS_08000m_v0.nc` (int32)
- densities in params.nc: 917 / 1024 / 1000
- reference: the state stamped 2016 (nominal 2015), the run's own
- years: 2015 to 2019 (5)
- flux means: whole pixel (`flux_pixel_mean`)
- pixel coverage: none given
- native scalars: over true area; the writer's stamp: true_area (`<scratch>/ismip7_issue96_98/trees/p5/AIS/RICE/icepack2/CORE/C009`)
- comparison: numpy 2.5.2, netCDF4 1.7.4
- commit df1a2c0
- ismip7-scalars 0.1.0, isschecker 0.5.1
- af2_AIS_08000m_v1.nc sha256 e62c8d274cae4c262b495211ad6e8ba6fc870786e45ee31021317e4e9c3a6b0f
- iaf2_GIC_AIS_08000m_v0.nc sha256 a133c411b899b84082625d8b2e78350294a7a76dbc8722adb7a5a46d79393999
- maxmask1_AIS_08000m_v0.nc sha256 5a13d364ba3cbc7fdbe28bc7d54164c157ccbd33075c178efc77f7031fb7d310

Exit status 0.

### Gates

| gate | result | worst |
|---|---|---|
| densities | pass | 0.000e+00 against 0.000e+00 at rhoi |
| scalar files against the CSV | pass | -6.509e+04 against 1.751e+05 at iareafl 2019 |
| T against the replay R | pass | 0.000e+00 against 2.430e+10 at lim 2015 |
| forbidden-policy sums against N | pass | -8.645e+08 against 8.202e+09 at ice area 2019 |
| tendacabf against N | pass | -2.383e+03 against 4.671e+04 at 2015 |
| zero fluxes | pass | 0.000e+00 against 0.000e+00 at tendlibmassbfgr 2015 |
| sla20 = slg20 (fixed bed) | pass | 1.215e-14 against 1.000e-06 at 2019 |
| slg20 - slvaf identity | pass | -4.764e-15 against 1.000e-09 at 2016 |
| topg constant in time | pass | 0.000e+00 against 0.000e+00 at 2015 |
| lim change against dlithkdt | pass | -3.502e+12 against 2.488e+13 at 2019 |
| lim residual sign | pass | 0.000e+00 against 1.442e+16 at 2015 |
| T differs from N | pass | 0.000e+00 against 0.000e+00 at lim 2015 |
| af2 > 0 under ice and flux | pass | 0.000e+00 against 0.000e+00 at all years |

### Named differences

- area factor: the native scalars carry it too (--native-af2), so the area term is zero; the exact sums' worst |C - N| is 6.26e-05 of their L1 (ice area 2019), against an allowance of 5.93e-04
- maximum-extent mask: 0 pixels carry ice or flux outside maxmask1; lim +0.000%, limnsw +0.000%, slvaf +0.000%, slg20 +0.000%
- fill convention: acabf and libmassbffl are whole-pixel means (flux_pixel_mean), the tool's own, so nothing is undone and the fill term is zero
- tendlibmassbffl: the model's value carries -2.3 Gt/yr in pixels with no floating ice at year end, which the fill leaves out (2019)

### By scalar

T - N and its parts, as shares of max |N| over the run; sea level in mm.

#### lim (kg)

max |N| = 2.4299e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 (worst) | 2.4298e+19 | 2.4296e+19 | -0.008% | +0.000% | +0.000% | +0.000% | -0.008% |
| 2019 | 2.4299e+19 | 2.4297e+19 | -0.008% | +0.000% | +0.000% | +0.000% | -0.008% |

#### limnsw (kg)

max |N| = 2.1197e+19

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 2.1195e+19 | 2.1195e+19 | -0.000% | +0.000% | +0.000% | +0.000% | -0.000% |
| 2019 (worst) | 2.1197e+19 | 2.1197e+19 | -0.001% | +0.000% | +0.000% | +0.000% | -0.001% |

#### iareagr (m2)

max |N| = 1.2358e+13

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.2358e+13 | 1.2357e+13 | -0.007% | +0.000% | +0.000% | +0.000% | -0.007% |
| 2019 (worst) | 1.2357e+13 | 1.2357e+13 | -0.007% | +0.000% | +0.000% | +0.000% | -0.007% |

#### iareafl (m2)

max |N| = 1.4594e+12

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 1.4594e+12 | 1.4594e+12 | -0.001% | +0.000% | +0.000% | +0.000% | -0.001% |
| 2017 (worst) | 1.4593e+12 | 1.4593e+12 | -0.001% | +0.000% | +0.000% | +0.000% | -0.001% |
| 2019 | 1.4591e+12 | 1.4590e+12 | -0.001% | +0.000% | +0.000% | +0.000% | -0.001% |

#### tendacabf (kg s-1)

max |N| = 7.8657e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 (worst) | 7.8657e+07 | 7.8655e+07 | -0.003% | +0.000% | +0.000% | +0.000% | -0.003% |
| 2019 | 7.8657e+07 | 7.8655e+07 | -0.003% | +0.000% | +0.000% | +0.000% | -0.003% |

#### tendlibmassbfgr (kg s-1)

zero on both sides, 5 years

#### tendlibmassbffl (kg s-1)

max |N| = 3.3274e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -3.3274e+07 | -3.3202e+07 | +0.218% | +0.000% | +0.000% | +0.000% | +0.218% |
| 2019 (worst) | -3.3101e+07 | -3.3027e+07 | +0.223% | +0.000% | +0.000% | +0.000% | +0.223% |

#### tendlicalvf (kg s-1)

max |N| = 2.3476e+04

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -2.0003e+04 | -2.0003e+04 | +0.001% | +0.000% | +0.000% | +0.000% | +0.001% |
| 2019 (worst) | -2.3476e+04 | -2.3476e+04 | +0.001% | +0.000% | +0.000% | +0.000% | +0.001% |

#### tendlifmassbf (kg s-1)

zero on both sides, 5 years

#### tendligroundf (kg s-1)

max |N| = 6.3013e+07

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | 6.2889e+07 | 6.2889e+07 | -0.000% | +0.000% | +0.000% | +0.000% | -0.000% |
| 2016 (worst) | 6.2952e+07 | 6.2952e+07 | -0.000% | +0.000% | +0.000% | +0.000% | -0.000% |
| 2019 | 6.3013e+07 | 6.3012e+07 | -0.000% | +0.000% | +0.000% | +0.000% | -0.000% |

#### slvaf (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2019 (worst) | -6.018 | -5.906 | +0.112 | +0.000 | +0.000 | +0.000 | +0.112 |

#### slg20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | -0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2019 (worst) | -5.974 | -5.865 | +0.109 | +0.000 | +0.000 | +0.000 | +0.109 |

#### sla20 (mm)

| year | N | T | T - N | area | max mask | fill | residual |
|---|---|---|---|---|---|---|---|
| 2015 | -0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| 2019 (worst) | -5.974 | -5.865 | +0.109 | +0.000 | +0.000 | +0.000 | +0.109 |
