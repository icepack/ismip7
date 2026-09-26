# Geometry discretization and the calving front

`ISMIP7_GEOMETRY_SPACE` selects the finite-element space for the geometry
(`h`, `s`, `b`, and everything derived from them: `phi_eff`, `N_ref`, `C_w0`,
`H_init`, and the forcing fields `accum` / `ocean_melt`). The apparent-MB
reference `a_ref_mb` is a per-cell transport source and stays DG0 either way.

- `dg0` (default): cell-wise geometry. The momentum solve and the mass
  transport share one thickness field.
- `cg1`: the pre-August-2026 behaviour, kept for A/B comparison.

The control fields (`theta`, `phi`, `fluidity_prior`, `obs_mask`) and the
velocity stay CG1 in both modes. The controls are regularized by a
Whittle-Matern prior on `|grad theta|`, which needs a differentiable field;
the geometry is not differentiated in the residual (see below).

## Why this changed

Under CG1 geometry there were two thicknesses. The transport evolved a DG0
field; the momentum solve wanted CG1; a lumped-mass lift bridged them:

```
h_CG1[i] = sum_K |K∩supp(phi_i)| h_DG0[K] / sum_K |K∩supp(phi_i)|
```

The lift conserves volume exactly and reconstructs well in the interior. On
the domain boundary its stencil is one-sided: a boundary node averages only
interior cells, so it is pulled toward interior values. On an unbuffered mesh
that boundary is the calving front, where the terminus traction

```
f_I - f_W = 1/2 rho_I g h^2 - 1/2 rho_W g d^2
```

goes as `h^2`. Measured at 32 km with the n=3 Budd MAP:

| | front `<h>` | terminus force | front `<u.n>` | outflux |
|---|---|---|---|---|
| before the lift | 105 m | 1.0x | +238 m/yr | 848 Gt/yr |
| after the lift | 200 m | **3.3x** | +698 m/yr | 3980 Gt/yr |

BedMachine's ice-front thickness is 145 m over all fronts and 167 m for
floating fronts (median 152, p90 306), so the unlifted CG1 field ran 30% thin
and the lifted one 25% thick. That is a factor 3.6 spread in the front driving
force from discretization alone, larger than any physics knob in the model.
Force and mass disagreed too: with one velocity field, the outflux from the
momentum-solve thickness and from the transport carrier differed by 22% while
the volume integrals agreed to 0.000%.

Under DG0 geometry `h` and `h_dg` are the *same Function*. Re-measured:

```
front <h>                      158.9 m      (BedMachine: 145 all / 167 floating)
outflux from momentum-solve h  identical
outflux from transport carrier identical
force/mass disagreement        0.000%
```

## What it costs

The driving stress is represented as broken cell gradient + facet jump:

```
cell : -rho_I g h grad(s) dx
facet:  rho_I g avg(h) <jump(s, nu), avg(v)> dS
```

Together these are the distributional gradient of a possibly discontinuous
surface. With CG1 `s` the jump vanishes and the cell term carries everything.
With DG0 `s` the cell gradient is identically zero (UFL folds it away) and the
facet term carries everything. Both are consistent, and they differ in
accuracy. Against the exact continuum driving force on a smooth manufactured
problem:

| n | CG1 error | rate | DG0 error | rate |
|---|---|---|---|---|
| 8 | 1.20e-02 | - | 5.98e-02 | - |
| 32 | 7.55e-04 | 2.00 | 1.45e-02 | 1.02 |
| 128 | 4.72e-05 | 2.00 | 3.59e-03 | 1.00 |

CG1 is second order and DG0 first order on the interior driving stress, so the
trade is 1.4% against 0.08% error at 32x32 in the interior, bought for an
unbiased calving front. At Antarctic resolutions the bed and thickness data
error dominates that gap, and the CG1 path never reached its asymptotic rate at
the front: the lift smooths `h` by up to 1.4 km at the steep PIG grounding-zone
gradient.

Two further consequences to watch, neither yet a demonstrated problem:

1. **The grounding line is a staircase.** `He = smooth_heaviside(haf(H, b))` is
   cell-constant, so GL migration proceeds cell-by-cell rather than
   sub-element. Coarse at 32 km; much less so at 2500 m.
2. **The lift was incidentally damping grid-scale noise.** Removing it removes
   that filter. First-order upwind transport is strongly diffusive at the grid
   scale and should cover it. Worth measuring.

## Sampling the raster: cell average, not centroid

A DG0 dof sits at the cell centroid, so the obvious
`icepack.interpolate(raster, Q_g)` takes a one-point sample of a 500 m
BedMachine raster per cell, so at 32 km one pixel stands in for a 32 km cell.
The DG0 driving stress is entirely the facet jump in `s`, so that sampling
noise reads as slope. Measured at 32 km:

| construction | rms \|jump s\| | peakedness (L4/L2) | front `<h>` | \|driving force\| vs CG1 |
|---|---|---|---|---|
| cg1 | 0 (continuous) | - | 106.9 m | - |
| dg0, centroid sample | 366 m | 2.04 | 209.5 m | +7.7% |
| **dg0, cell average** | **257 m** | **1.49** | **152.9 m** | **+1.0%** |

The centroid version failed to converge in 200 Newton iterations. The cell
average, meaning the L2 projection of the CG1 interpolant, is 42% smoother and
much less peaked, lands on BedMachine's front thickness (median 152 m), and
reproduces the CG1 driving force to 1%. `geometry.sample_to_geometry` does
this; do not replace it with a direct interpolate onto the DG0 space. Cell
average here means that L2 projection (`ISMIP7_RASTER_SAMPLE=vertex`, the
default). The raster's true cell mean is a separate option (`cell_mean`) and
measures rougher across exactly these facet jumps, so it is kept for the
record. See its row in `antarctica/README.md`.

The same rule applies to the RACMO SMB climatology, which sets the mass budget
and the `a_ref` balance: `forcing.load_racmo_smb_climatology` cell-averages onto
a DG0 space (`forcing._sample_raster`) and stays the nodal interpolant on CG1.

It also applies when an inversion is used on a different timing mesh. A direct
DG0-to-DG0 cross-mesh interpolation samples the discontinuous source field at
each target-cell centroid; non-nested mesh boundaries then alias into artificial
target-cell jumps. This escaped the first cached timing campaign because the
controls and velocity are continuous CG1 fields, where the same transfer is
appropriate, while `h`, `b`, and `s` are not. On the 2500/25000 timing mesh the
aliased cache had an RMS surface jump of 67.2 m versus 51.5 m for target-native
BedMachine cell averages, and its first transient steps ran from 19,092 m/yr to
4.05e6 and then 5.10e9 m/yr. The later transport-budget failure was a symptom of
that diagnostic runaway. Exact-mesh caches therefore construct `h` and `b` from
BedMachine on the target mesh, recompute hydrostatic `s`, transfer only the
continuous inversion products, and recompute geometry-dependent anchors before
the cold diagnostic continuation.

## Where a CG1 reconstruction is still used, and why that is legitimate

A DG0 field has no pointwise gradient, so three places reconstruct one via
`geometry.cg1_lift` (the same lumped lift, used deliberately and locally):

- `geometry.surface_slope`, used by `dual_friction.weertman_anchor`. `C_w0` is
  a fixed reference scaling that defines what `theta = 0` means, and enters no
  force.
- `forcing.compute_sin_alpha`. Feeds the Burgard melt *parameterization*, not
  the momentum residual.
- The thermomechanical fluidity prior in `inversion_icepack2.py`, which runs on
  CG1-lifted geometry throughout.

None of these appear in the momentum residual or in a mass flux, which is the
property that matters: the terminus traction and the transport flux must
integrate the same field, and they do.

The fluidity prior must be lifted rather than L2-projected. An L2 DG0 to CG1
projection overshoots at the front and produced a negative fluidity (`A_prior`
min -9.88 against a DG0 range of `[1.0, 446.7]`), which poisons
`log(A/A_prior)`. `cg1_lift` is a convex combination of cell values, so it
cannot overshoot.

## MAPs are not interchangeable

The inversion absorbs the front treatment into `theta` and `phi`: whatever the
momentum balance gets wrong at the terminus, the optimizer compensates for by
adjusting friction and fluidity until the modelled velocity matches
observations. The t=0 velocity misfit therefore looks fine either way, and only
prognostic behaviour separates them: drift, calving flux against the observed
1300 Gt/yr, and front-region velocity against MEaSUReS.

Because of that, MAPs carry a geometry tag (`map_geom_tag()`):

- CG1: `inversion_icepack2_budd_n3_<lc>.h5` (legacy, untagged)
- DG0: `inversion_icepack2_budd_n3_dg0_<lc>.h5`

A forward prefers the tagged MAP for its own geometry space. With only the
legacy one present it loads, projects and warns loudly, which is a smoke test
rather than a result. Driving a DG0 forward with the CG1 MAP at 32 km raises
the initial misfit from 8.6e3 to 1.5e5, the size of the inconsistency.

Run the DG0 inversion with:

```
OMP_NUM_THREADS=1 ISMIP7_FRICTION=budd ISMIP7_LC=32000 ISMIP7_LC_COARSE=320000 \
  ISMIP7_N_FLOW=3.0 \
  ISMIP7_GEOMETRY_SPACE=dg0 \
  ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_320000_32000.msh ISMIP7_MAXITER=500 \
  mpiexec -n 8 python antarctica/scripts/inversion_icepack2.py
```

## The ice front of a buffered mesh (issue #115)

On a buffered mesh the ice front sits inside the domain, beside ice-free buffer
cells. The floor-cell ocean drag (`ISMIP7_OCEAN_DRAG`, 1e-2 MPa yr/m, ramping to
zero at `ISMIP7_H_OCEAN` = 10 m) acts in those cells, and the CG1 front nodes
they share with the ice take it. The inversion and the legacy fixed-front mask
(`ISMIP7_FIXED_FRONT`) keep the drag everywhere. A level-set front
(`ISMIP7_CALVING`) keeps it only where the signed distance to the current ice
extent exceeds one cell diameter, so every ice cell and the first ring or two
of water lose it.

Measured on the v4 2500/25000 timing cache, the 2.5 km Budd MAP's own state,
with `antarctica/scripts/front_flux_check.py`. The front is the set of facets
between t=0 ice (`H_init >= 1 m`) and the buffer, and fluxes are upwind:

| state | floating front u.n (m/yr) | floating front flux (Gt/yr) | grounding-line flux (Gt/yr) |
|---|---|---|---|
| `velocity_obs` on the same facets | +82 | 168 | 2,139 |
| t=0, drag on (the MAP) | +3.5 | 10.5 | 2,266 |
| first solve with the level-set gate | +768 | 1,793 | 2,599 |

At t=0 the deficit sits at the front nodes: floating ice 5 to 100 km behind
the front moves at 419 to 496 m/yr, against 383 to 428 m/yr in
`velocity_obs`. The gated solve speeds floating ice up at every distance, to
1,362 m/yr within 5 km of the front and 1,519 to 1,950 m/yr from 25 to 250 km.
The MAP's front band is slightly stiffer than the rest of the ice (mean log
fluidity -0.13 against -0.04).

Ten years from the same cache through `run_timing.py`, with no forcing and
`ISMIP7_APPARENT_MB=div` (equal to balance without forcing), the two runs
differing only in `ISMIP7_CALVING`; records `test-2500m-budd-legacy-front` and
`test-2500m-budd-levelset-fixed` in `antarctica/runlog/`:

| | legacy mask | level set, `ISMIP7_CALVING=fixed` |
|---|---|---|
| calving, 2016 / 2020 / 2024 means (Gt/yr) | 12.4 / 12.5 / 12.4 | 3,534 / 3,135 / 2,839 |
| grounding-line flux, 2025 (Gt/yr) | 2,266 | 4,028 |
| mass change, 2015 to 2025 (Gt) | -2 | -30,908 |
| VAF change, 2015 to 2025 (mm SLE) | +0.004 | -18.8 |
| fastest node (m/yr) | 4,611 | 81,018 at step 1, 14,161 in 2025 |
| Newton iterations a step, mean / max | 5.0 / 9 | 8.8 / 16 |

Every step of both runs converged on its first direct solve. The timing
harness's speed tripwire (2e4 m/yr) stops the level-set run at step 1, so the
ten-year run leaves it unset, as the production runners do. The apparent-MB
reference is built from the drag-on velocity before the level set exists, so
step 1 of both runs moves ice with the same velocity and the gated velocity
first moves ice in step 2.

A restart checks the loaded state against the residual with the drag on
everywhere. From the level-set run's 2024.0 checkpoint that residual was 2.19e10
against an acceptance of 7.84, so the restart re-solved the state with the
drag on, and its first step calved 21.7 Gt/yr where the uninterrupted run
calved 2,865.

## The melt calibration follows the forward's melt path

The calibrations melt on the same `ISMIP7_GEOMETRY_SPACE` as the forward.
Under `dg0` they evaluate the forward's own cell by cell melt path: bed and
thickness sampled onto the cells, the surface from flotation, the slope of
`forcing.compute_sin_alpha` (the constant under the default
`ISMIP7_MELT_SLOPE=ant`, the uncapped cell slope under `local`), forcing at
each centroid and its own draft, the forward's melt set
`forcing.melt_receiving` (the seawater floating test `forcing.is_floating` on
cells holding ice, `h > 0`) and cell areas. A K or an offset fitted there is
the one the forward applies, by construction. The runs melt with the tracked
calibration of item 4.

The earlier K files were fitted under `cg1`: BedMachine on CG1 nodes with its
raster mask, the nodal slope capped at 5e-3, lumped-mass areas. The K file
records `geometry_space`, `melt_slope`, `sin_alpha_ant` and `sin_alpha_cap`,
and `load_K_per_basin` warns once when a run melts on a geometry other than the
one its K was fitted on, once when its slope convention or constant differs
from the file's (an untagged file reads as `local`), and, under `local` only,
once when the file was fitted against a capped slope while the run applies none.
An offsets file records the same, and `load_deltaT_per_basin` refuses any of
those differences instead.

Two things were found by fitting through the forward's path, both measured on
the 2500 m mesh (`inversion_icepack2_budd_2500.h5`, BedMachine v4.1
vertex-sampled, the 865.0 Gt/yr table, 21 September 2026):

1. **The melt callbacks floated ice with fresh-water density**, issue #66,
   closed on 22 September 2026 once the shared test landed.
   Their `haf <= 0` test used the SMB conversion's 1000 kg/m3 where the forward
   floats ice at 917/1024, which read every floating cell thicker than 78
   percent of its flotation thickness as grounded: 364 038 km2, 24 percent of
   BedMachine's shelf area, got no melt (the seawater test grounds 5 160 km2 of
   it on nodes). The callbacks, the calibration and `check_melt_bound.py` now
   share `forcing.is_floating`, with seawater. The rows `check_melt_bound.py`
   records for the 2 km mesh (1732 and 646 Gt/yr) were measured through the
   fresh-water test with the legacy K files.

2. **With a consistent floating test, the geometry space barely matters; the
   slope convention does.** What the forward's uncapped DG0 path applies with
   each K file, seawater flotation, ice-present cells:

   | K fitted on | slope cap | K* | K total-match | forward applies today (uncapped) | basins in Burgard K5..K95 |
   |---|---|---|---|---|---|
   | cg1 nodes (the production file) | 5e-3 | 4.26e-5 | 5.34e-5 | 3555 Gt/yr | 4 of 16 |
   | dg0 cells | 5e-3 | 4.37e-5 | 5.68e-5 | 3679 Gt/yr, 865 if the forward caps too | 4 of 16 |
   | dg0 cells | none | 6.67e-6 | 1.45e-5 | 865 Gt/yr, by construction | 0 of 16 |

   The capped DG0 fit reproduces the nodal per-basin K within about 10 percent
   per basin, 22 percent in basin 7 (basin 9, Amundsen, 1.47e-4 on both), over
   a floating area of 1 577 840 km2 against 1 509 122 on nodes. The uncapped cell slope has a median of 1.5e-2
   over floating cells and integrates 3.7 times the capped nodal melt at K = 1,
   so a K fitted to it absorbs mesh slope noise and lands every basin below
   Burgard's range, while the forward as it runs today would apply 4.1 times
   the target with the production K once the flotation test is right. The
   slope convention was settled in issue 26. The ISMIP7 reference example is one
   constant mean-Antarctic slope, "no slope dependency", and the toolbox's K
   percentiles (July 2026: K05 4.75e-5, K50 8.5e-5, K95 1.375e-4) were sampled
   with sin(alpha) = 5.115e-3, back-computed from the notebook's own gamma_T
   conversion. That is now the forward's and the calibration's default
   (`ISMIP7_MELT_SLOPE=ant`, `ISMIP7_SIN_ALPHA_ANT`); the local slope stays as
   `local`, capped in the calibration and uncapped in the forward as before.

3. **Under the constant slope the geometry space is immaterial and K lands
   below the toolbox's K50 on the July table.** Same mesh, same table,
   seawater flotation:

| slope | geometry | K* | K total-match | melt at K* | basins in K05..K95 |
|---|---|---|---|---|---|
| constant 5.115e-3 | dg0 cells | 4.06e-5 | 4.44e-5 | 791 Gt/yr | 7 of 16 |
| constant 5.115e-3 | cg1 nodes | 4.04e-5 | 4.44e-5 | 787 Gt/yr | 7 of 16 |

   The two fits agree to 0.5 percent. Fitted to the July 2026 table (1067.4
   Gt/yr) on the 2 km MAP mesh, DG0 cells (24 September 2026, run record
   `calibration-melt-2km-1067`), K* is 4.46e-5, just under the toolbox's K05,
   and the total-match K is 5.59e-5, between K05 and K50; K* melts 851 Gt/yr
   and 7 of 16 basin K fall in K05..K95. Both sit below K50, as a
   term-1-only fit should: the notebook's own K50 applies 1571 Gt/yr against
   the 1067 observed, because terms 2 to 4 pull K up. The protocol's
   per-basin adjustment is a temperature offset at fixed K
   (`calibrate_deltaT.py`, antarctica/README.md section 5).

4. **The tracked calibration.** The group chose on 25 September 2026 (issue
   26) the K50 of the toolbox objective run through this path on the
   1000 m / 10 km production mesh, with the offsets fitted for every K first
   and the objective restricted to the K whose offsets keep present-day
   thermal forcing plausible (`select_melt_parameters.py`): K = 6.5e-5, with
   offsets from -0.68 to +1.20 K, tracked as
   `antarctica/calibration/deltaT_per_basin_1000_K6.500e-05.npz`. Every run
   reads it unless another file is named. The forward melts
   `forcing.melt_receiving`, the set the fit summed over. Measured through
   the forward's own callback on the production mesh (run record
   `calibration-melt-forward-1km-k50`): 1067.390 Gt/yr against the 1067.386
   the offsets were fitted to, every basin within 0.006 Gt/yr. The
   callbacks' earlier melt set, every `haf <= 0` cell, also covered 257 687
   ice-free open-ocean cells at draft 0 there, booking 156.3 Gt/yr of melt
   and 23.1 Gt/yr of refreezing on them, and 6 666 cells of bare land, where
   the climatology melts nothing. At 32 km the same offsets put the
   basins at 0.33 to 1.74 times their totals, 1069.5 Gt/yr in all.

## Incompatibilities

- `ISMIP7_LEGACY_TRANSPORT=1` requires `ISMIP7_GEOMETRY_SPACE=cg1` (it
  re-projects CG1 `h` <-> DG0 `h_dg` every step; under DG0 they are the same
  Function and the projection would be self-referential). Enforced.
- A run must not change geometry space mid-trajectory. Restart checkpoints
  record `geometry_space`; a mismatch projects and warns.
