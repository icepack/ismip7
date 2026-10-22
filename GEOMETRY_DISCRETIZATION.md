# Geometry discretization and the calving front

`ISMIP7_GEOMETRY_SPACE` selects the finite-element space for the geometry
(`h`, `s`, `b`, and everything derived from them: `phi_eff`, `N_ref`, `C_w0`,
`H_init`, and the forcing fields `accum` / `ocean_melt`). The apparent-MB
reference `a_ref_mb` is a per-cell transport source and stays DG0 either way.

- `dg0` (default): cell-wise geometry. The momentum solve and the mass
  transport share one thickness field.
- `cg1`: the pre-August-2026 behaviour, kept for A/B comparison.
- `dg1`: cell-wise linear geometry, discontinuous across cells, on vertex
  nodes. One thickness field as under `dg0`, at second order. See "DG1
  geometry" below for what it carries and what it refuses.

The control fields (`theta`, `phi`, `fluidity_prior`, `obs_mask`) and the
velocity stay CG1 in every mode. The controls are regularized by a
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
much less peaked, lands on BedMachine's front thickness (median 152 m) on this
buffer-0 mesh (a buffered mesh's front is thinned; see "Front cells on a
buffered mesh"), and
reproduces the CG1 driving force to 1%. `geometry.sample_to_geometry` does
this; do not replace it with a direct interpolate onto the DG0 space. Cell
average here means that L2 projection (`ISMIP7_RASTER_SAMPLE=vertex`, which
the default `vertex_front` keeps everywhere but the marine front). The
raster's true cell mean is a separate option (`cell_mean`) and measures
rougher across exactly these facet jumps, so it is kept for the record. See its row in `antarctica/README.md`.

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

## Front cells on a buffered mesh (issue #167)

On a buffered mesh BedMachine's ice front crosses cells. Vertex sampling gives
each of them the mean of its three vertex samples, and a vertex over open
water samples zero, so the model's front is a rim of partly covered cells at
a fraction of the front's thickness. Every cell of 1 m or more counts as ice
(`ISMIP7_FRONT_HMIN`), and the transport is upwind, so the flux out of the
front is the rim's thickness times a normal speed the inversion fits to the
observations. IU's final 2 km MAPs carried 343 Gt/yr out of a front band of
39.5 m, where BedMachine's own 500 m front is 163 m thick and carries 1,138 to
1,297 Gt/yr under MEaSUREs v2 (`bm_front_flux.py`). Without apparent MB the
shelves then gained about 850 Gt/yr (run record
`inversion-2km-budd-b20k-ef2-relax2014-year`). A buffer-0 mesh puts the front
on the mesh boundary and keeps 139 m at 2 km.

`ISMIP7_RASTER_SAMPLE=vertex_front` (`geometry.front_cells`) rebuilds the
marine front from BedMachine's mask, sampled on the `cell_samples` lattice
and again four times finer in mixed cells:

- a cell holds ice when at least half its samples are BedMachine ice;
- a cell without ice that holds an ocean sample is water, `H = 0`;
- a cell with ice that shares a vertex with a water cell, or holds an ocean
  sample itself, takes BedMachine's mean thickness and bed over its ice
  samples;
- every other cell keeps its vertex samples, so interior ice, land margins
  and the open buffer are unchanged.

It is the default. A MAP records its sampling (`raster_sample`, an objective
key), a forward on the MAP's own mesh follows that record, a transfer to
another mesh rebuilds under the run's sampling, and the melt calibration is
the file fitted under the run's sampling (`runconfig.MELT_CALIBRATIONS`).

The t = 0 floating front under `velocity_obs` (run record
`test-i167-front-probes`):

| mesh | sampling | front thickness | front band mean | front flux |
|---|---|---|---|---|
| BedMachine 500 m | raster | 163 m | | 1,138 to 1,297 Gt/yr |
| 2 km, 20 km buffer | `vertex` | 47.7 m | 39.5 m | 179 Gt/yr |
| 2 km, 20 km buffer | `vertex_front` | 190.9 m | 129.5 m | 998 Gt/yr |
| 1 km, 20 km buffer | `vertex` | 36.3 m | 25.6 m | 179 Gt/yr |
| 1 km, 20 km buffer | `vertex_front` | 178.4 m | 90.3 m | 1,037 Gt/yr |

On the 2 km mesh the rule rebuilds 27,854 cells and empties 16,637; the
floating area falls 3.4 % and the floating mass 0.26 %, and the
melt-receiving area falls 3.6 % (2.4 % at 1 km), which is why the rule has
its own melt calibration. Meshes whose nodes follow BedMachine's marine front
(`ISMIP7_MESH_FRONT=bm`, `mesh_antarctica.py --front bm`) carry the same front
flux within 3 % under the rule at 6 to 8 % more cells, and the final Budd
controls started their refit at misfit 4,250 there against 1,882 on the
current mesh (`test-i167-frontbm-meshes`,
`test-2km-budd-ef2-vertex-front-eval1`), so the production meshes stay the
current ones.

The controls absorb the front they were fitted against, so a MAP is refitted
under the rule before it drives a forward on a buffered mesh. IU's Budd MAP
refitted under it (`inversion-2km-budd-b20k-ef2-vf`) carries 1,313 Gt/yr out
of a 191 m front at t = 0, its relaxation year without apparent MB holds the
front at 1,420 to 1,317 Gt/yr with dM/dt near zero
(`inversion-2km-budd-b20k-ef2-vf-relax2014-year`), and it transfers onto the
1 km mesh with an apparent MB of +161 Gt/yr. Rice's Budd snapshot of 24
September, fitted on the buffer-0 mesh, transfers onto the 1 km mesh near
balance under `vertex` (apparent MB -473 Gt/yr) and with its shelves at three
times the observed speed under the rule (+9,004 Gt/yr;
`test-1km-transfer-t0-vertex-front`), so a Rice MAP drives a buffered-mesh
forward only after a refit under the rule (IU, 9 October). The two lines run
side by side until the submission is chosen on 20 and 21 October: Rice's
MAPs and forwards on buffer-0 meshes, IU's on the buffered meshes under
`vertex_front`. On a buffer-0 mesh the rule is unmeasured, and
`ISMIP7_RASTER_SAMPLE=vertex` keeps a buffer-0 run as it was before it.

## The drag gate and the membrane floor (issue #153)

On the 20 km buffered 2 km mesh the ocean drag reached the ice front by two
routes. A water cell touching the ice at one vertex drags that front node,
since the velocity is CG1, and `front.ocean_drag_cells` skipped only the cells
sharing an edge with ice (`ISMIP7_DRAG_GATE=facet`). The membrane floor
`h_visc_floor` (`ISMIP7_RC_HVISC_FLOOR`) gives the first water row that much
ghost ice, sheared between the moving front and the dragged water one cell
out. `ISMIP7_DRAG_GATE=vertex`, the default, closes the first route; a lower
floor narrows the second.

RC's stage-1 controls on that mesh at evaluation 1, model minus observed speed
(mean; front band: floating ice within 25 km of open water):

| gate, floor | floating ice (median ratio) | front band |
|---|---|---|
| facet, 10 m | -119 m/yr (-27 %) | -138 m/yr |
| vertex, 10 m | -39.2 (-11 %) | -52.2 |
| vertex, 5 m | -22.9 (-7 %) | -34.0 |
| vertex, 2.5 m | -12.7 (-4 %) | -22.2 |
| vertex, 1 m | -5.6 (-2 %) | -13.6 |
| drag off, 10 m | -4.7 (-1 %) | -17 |

Over 300 stage-2 iterations the inversion absorbs most of this into phi. The
front band ends at -13.1, -14.0, -15.5 and -18.5 m/yr at 1, 2.5, 5 and 10 m,
and the total misfit falls with the floor, 2,833 to 2,791. Against the 1 m MAP
the floating phi is softer by +0.14, +0.34 and +0.64 at its 99th percentile
(2.5, 5, 10 m), a dipole on Ross and Filchner-Ronne: the shelf interior runs
faster and the front band slower, by 22, 51 and 95 m/yr at its 1st
percentile.

Forward cost against 10 m: scpc_gamg, five unforced years drifting by the flux
divergence (one year at 1 km), every solve converged with no rescue:

| floor | 2 km, 32 ranks, from the MAP state | 2 km, cold start | 1 km, 64 ranks, cold start | 1 km start-up |
|---|---|---|---|---|
| 1 m | 2.5x | 1.6x, start-up 4.8 h | 1.53x | 19 min |
| 2.5 m | 1.9x | 1.4x | 1.14x | 17 min |
| 5 m | 1.3x | 1.0x | 1.10x | 15 min |
| 10 m | 14.4 s a solve | 11.4 s | 18.2 s | 11 min |

The inversion runs 2.5 m (`runconfig.HVISC_FLOOR_DEFAULT`, decided 6 October
2026 for the submission schedule): it keeps about a fifth of the 10 m
artifact for 14 % more forward time at 1 km, where 1 m takes 53 %. The MAP
records `drag_gate` and `h_visc_floor`, and a forward runs both
(`runconfig.forward_drag_gate`, `forward_hvisc_floor`); a MAP older than the
record was inverted at 10 m and keeps it. At rho 75 km scpc_gamg stalls at
every floor, so the floor leaves the inversion's solver choice there
unchanged. Records: `test-2km-rc-b20k-eval1-*`, `test-2km-rc-b20k-vgate-floor*`,
`test-2km-rc-b20k-forward-floor*`, `test-1km-rc-forward-floor*-drift`.

## The push at an ice edge (issues #153, #166)

On a facet between a cell holding ice and one holding none, the DG0 driving
stress is the facet term `rho_I g avg(H) (s - s_other)`. It equals the
depth-integrated front push for floating ice and on flat land, falls short at
a grounded marine cliff by `g D (rho_I H - rho_W D) / 2`, and on land follows
the step in bed height. `ISMIP7_EXACT_FRONT` replaces it on every such facet
with the push of its version (`dual_friction.front_cliff_correction`, the
push itself in `cliff_push`). The MAP records the version as `exact_front`,
and the forward runs the MAP's version and refuses a knob that differs.

| ice-free neighbour's bed `B` | version 1 | version 2 |
|---|---|---|
| at or below the ice base (ocean, lower ground) | free-cliff push `g (rho_I H^2 - rho_W d^2) / 2` | the same |
| seabed between the ice base and sea level | free-cliff push | the face above `B`, water down to `B` |
| land between the ice base and surface | free-cliff push | the face above `B`, `h_e = s - B` |
| rock at or above the ice surface | free-cliff push, into the rock | 0 |

For 500 m of ice on a 100 m bed, rock at 0, 300 and 700 m: the facet term
alone pushes 1.2, 0.6 and -0.2 times the free-cliff push, version 1 pushes
1.0 in each case, and version 2 pushes 1.0, 0.36 and 0.

Census of IU's final Budd MAP (2 km buffered mesh,
`antarctica/scripts/check_cliff_facets.py`), push in MN/m:

| class | facets | length | edge `H` p50 / p90 / p99 | version 1 push p90 / p99 | version 2 minus version 1, summed |
|---|---|---|---|---|---|
| ocean, bed below the ice base | 11,128 | 30,362 km | 40 / 98 / 162 m | 5.0 / 15.5 | 0 |
| ocean, shoal | 123 | 300 km | 63 / 147 / 231 m | 25.1 / 88.4 | -926 GN |
| land below the ice base | 2,113 | 6,167 km | 3 / 30 / 87 m | 4.1 / 33.9 | 0 |
| land, partial wall | 239 | 697 km | 47 / 129 / 273 m | 61.5 / 334.3 | -14,382 GN |
| rock above the ice surface | 2,756 | 7,597 km | 5 / 83 / 271 m | 30.0 / 315.1 | -122,276 GN |

Most rock walls carry a thin edge, and 288 changed facets with edge ice
thicker than 100 m (849 km) carry 83 % of the change; the largest lie in the
Transantarctic and Prince Charles Mountains. In 124 grounded edge cells
(424 km^2) the change exceeds the cell's own basal drag. On the 1 km
production mesh the rock walls are longer (16,871 km) and their edges
thinner (p90 41 m), with 108 km of changed facets thicker than 200 m against
195 km at 2 km.

At the final Budd MAP's own controls, version 2 against version 1:

| | version 1 | version 2 |
|---|---|---|
| inversion: misfit, total, \|grad\| | 1,108.24, 1,388.36, 0.66 | 1,113.05, 1,393.17, 1.70 |
| inversion: Newton iterations, first forward | 6 | 12 |
| forward t = 0: Newton iterations (no FSSA) | 15 | 22 |
| forward t = 0: grounded discharge | 2,761.3 Gt/yr | 2,758.7 Gt/yr |

Restricting version 1 to ocean facets instead (2441818) gave misfit 1,146.8
and \|grad\| 3.28 at the same controls. By IMBIE basin the discharge changes
at most 0.99 % (basin 6) and 0.63 % (basin 7). The speed change stays at the
walls: grounded ice within 2 km of a changed facet moves p99 22.9 m/yr,
2 to 6 km 12.6, 6 to 20 km 2.6, and beyond 20 km 0.00; floating ice beyond
20 km 1.9. IU's final Budd MAP was refitted under version 2 from these
controls (`inversion-2km-budd-b20k-rho7500-floating-ef2`). IU's final RC MAP
was refitted under version 2 from its version 1 chain
(`inversion-2km-rc-b20k-rho7500-floating-ef2`); at the version 1 controls
version 2 lowered RC's misfit from 1,318.2 to 1,316.8. The sub-element path
(icepack_tools) has version 1 only and refuses 2. Records:
`test-2km-budd-b20k-exact-front-*`.

## The melt calibration follows the forward's melt path

The calibrations melt on the same `ISMIP7_GEOMETRY_SPACE` as the forward.
Under `dg0` they evaluate the forward's own cell by cell melt path: bed and
thickness sampled onto the cells, the surface from flotation, the slope of
`forcing.compute_sin_alpha` (the constant under the default
`ISMIP7_MELT_SLOPE=ant`, the uncapped cell slope under `local`), forcing at
each centroid and its own draft, the forward's melt set
`forcing.melt_receiving` (the seawater floating test `forcing.is_floating` on
cells holding ice, `h > 0`, on a bed below sea level) and cell areas. A K or
an offset fitted there is the one the forward applies, by construction. The
runs melt with the tracked calibration of item 4.

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
   thermal forcing plausible (`select_melt_parameters.py`): K = 6.5e-5. The
   offsets at that K, from -0.68 to +1.20 K, are refitted on Rice's build of
   the mesh, the submission mesh (issue 20), and tracked as
   `antarctica/calibration/deltaT_per_basin_1000_K6.500e-05.npz`, which a run
   sampled with `vertex` reads unless another file is named; the default
   sampling reads the `_vertex_front` file refitted at the same K ("Front
   cells on a buffered mesh"). The forward melts
   `forcing.melt_receiving`, the set the fit summed over. Measured through
   the forward's own callback on Rice's build (run record
   `calibration-melt-refit-1km-rice-k50`): 1067.389 Gt/yr against the
   1067.386 the offsets were fitted to, every basin within 0.006 Gt/yr. The
   callbacks' earlier melt set, every `haf <= 0` cell, also covered 257 836
   ice-free open-ocean cells at draft 0 there, booking 156.2 Gt/yr of melt
   and 23.1 Gt/yr of refreezing on them, and 6 805 cells of bare land, where
   the climatology melts nothing. The fit carries between builds only to a
   few tenths of a percent: the one on IU's build, where the K was selected,
   missed basins 1 and 7 on Rice's by 0.18 and 0.17 percent, with offsets at
   most 0.0008 K apart. At 32 km the offsets put the basins at 0.33 to 1.74
   times their totals, 1069.5 Gt/yr in all.

## DG1 geometry

`ISMIP7_GEOMETRY_SPACE=dg1` holds `h`, `s` and `b` cell-wise linear, on the
three vertices of each cell, with jumps across facets. Built for the question
whether a second-order geometry with sharp fronts removes the cell-scale dH
streaks of long DG0 forwards.

**Vertex nodes.** The space is `FunctionSpace(mesh, "DG", 1,
variant="equispaced")` (`geometry.geometry_function_space`). Firedrake's
default DG1 puts the nodes inside the cell, where a nodal floor does not bound
the field at the vertices and a vertex sample is not the CG1 interpolant.
`cellview.CellView` refuses any other DG1.

**Sampling.** `vertex` holds the CG1 vertex interpolant exactly, so a cold
start's geometry is continuous and its cell means are the DG0 sample: on the
buffer-0 meshes the means of `h` and `b` differ from DG0's by 9.1e-13 m at
32 km and 1.4e-12 m at 8 km. `vertex_front` classifies cells as under DG0 and
writes a rebuilt or water cell at all three nodes, so it is constant and its
mean is again the DG0 cell's value; every other cell keeps its linear sample.
`cell_mean` has no DG1 rule and is refused.

**Driving stress.** Both terms of `build_rc_residual` carry force: the cell
gradient the smooth part, the facet jump the rest. Against a smooth test field
on an L2-projected surface (`tests/test_dg1_geometry.py`), relative error of
the driving force:

| n | DG0 | DG1 | CG1 |
|---|---|---|---|
| 16 | 1.99e-2 | 7.08e-3 | 7.25e-3 |
| 32 | 1.21e-2 | 1.80e-3 | 1.79e-3 |
| 64 | 6.55e-3 | 4.54e-4 | 4.47e-4 |
| rate, 32 to 64 | 0.88 | 1.99 | 2.00 |

For floating DG1 ice the two terms are exactly the distributional gradient of
`rho_I g (1 - rho_I/rho_W) H^2 / 2` (cell: `H grad H = grad(H^2) / 2`; facet:
`avg(H) jump(H) = jump(H^2) / 2`), the balance the terminus traction uses. The
exact-front push decides which cells hold ice on their cell means
(`H_cell`), and on piecewise-constant DG1 fields assembles DG0's residual to
roundoff. The Budd HAF gate is pointwise: zero at every point of a floating
cell, and a grounding line crossing a cell gates it point by point. The Budd
reference pressure is pointwise as well: `N_ref` is
`dual_friction.budd_reference_pressure` evaluated on the frozen t = 0
thickness and surface, which a restart and a full-state MAP carry as
`N_ref_thickness` and `N_ref_surface`, so `N_hat` starts at the inversion's
own `N/N`. A nodal DG1 interpolant of `N` misses the kinks `N` has inside a
cell (the bed at sea level, the 1 m overburden floor, flotation); on the
32 km probe's DG1 MAP it put `N_hat` at a median 1.04, up to 1.7 in grounded
cells whose bed crosses sea level and up to the cap of 3 near the grounding
line, where clean grounded cells sit at exactly 1. DG0 keeps its per-cell
`N_ref` Function built from the same expression.

**Transport.** The same backward-Euler upwind solve with the DG volume term
`-h u.grad(phi)` (`transport.transport_form`); the L2 error of its steady
solution on a manufactured problem converges at rate 1.96 against DG0's 0.99
(32 to 64 cells a side). Testing with 1 still leaves only the boundary
outflow, so the transport residual and `resid` close exactly as under DG0. The
only limiter is a positivity floor (`cellview.clamp_cells`): a cell whose mean
falls below the floor is filled to it and booked as clamp, and a cell whose
lowest vertex does has its slope scaled until that vertex sits on the floor,
keeping the mean. Every advance holds the floor's assembled mass change to
the fill its pre-floor cell means call for (`cellview.floor_fill`), within the
transport tolerance, so mass the slope scaling moved stops the run. Firedrake's `VertexBasedLimiter` is kept
out: thickness peaks wherever ice fills a bed trough under a smooth surface,
and clipping those peaks rebuilds the surface roughness DG1 removes. On
BedMachine's `vertex` sample:

| buffer-0 mesh | rms jump of `s` over ice facets: DG0 / DG1 / DG1 after `VertexBasedLimiter` | ice cells the limiter moves by more than 1 m |
|---|---|---|
| 32 km | 238.4 / 0 / 77.4 m | 6,310 of 9,035 (p99 508 m) |
| 8 km | 106.1 / 0 / 67.1 m | 38,360 of 64,156 (p99 558 m) |

**Per-cell rules** read the cell means of `h` and `b`, with the surface the
DG0 formula builds from them (`cellview.CellView`, `forcing.forcing_geometry`),
so each is the DG0 rule applied to those means:

| rule | under `dg1` |
|---|---|
| SMB, melt, collapse mask | the forcing lives on the DG0 cells (`forcing.forcing_space`) and reads the cell means; the draft slope of `ISMIP7_MELT_SLOPE=local` is the DG0 slope of the cell-mean draft |
| melt calibration | a `dg0` file applies as fitted (`runconfig.melt_geometry_space`); the calibration tools run under `dg1` fit and stamp `dg0` |
| SMB-elevation feedback | the surface change of the cell means since `H_init`'s |
| grounded cells, forcing mask, drag gate, tripwire | the flotation test and the thresholds on the cell means |
| level-set front and its law | the level set and the law read the cell means (`h_cell`, the bed's means) |
| removals (fixed front, level set, slivers, collapse) | whole cells, booked at mean times area; the sub-cell shed scales the whole cell |
| dH/dt term of the inversion (`ISMIP7_DHDT_WEIGHT`) | one DG1 upwind step (`transport.transport_form`) forced on the DG0 cells, with the grounded test on the cell means (`obs_dhdt.grounded_cells`); the term scores the cell means of the tendency against Smith binned onto the cells (`obs_dhdt.cell_mean_form`), leaving the slope inside a cell free |
| backdate | the DG0 rule on the cell means; a changed cell moves by its mean's change, keeping its slope, floored by the positivity operator (`obs_dhdt.backdate_cells`) |
| ISMIP7 output | the per-cell fields are the cell means and their surface; the grounding-line and boundary fluxes are the transport's own, on the DG1 traces |

A t = 0 DG1 geometry has the DG0 cell values as its means to 1e-12 m, so at
t = 0 every one of these rules sees what it sees under DG0. A forward on
another mesh than its MAP's rebuilds the DG1 geometry from BedMachine with the
run's sampling, as DG0 does.

**Apparent mass balance.** The reference lives on the transport space. Under
DG1 it is the field whose integral against every DG1 basis function is the
t = 0 upwind divergence, the transport's own volume term included:
`a_ref = M^-1 R(h0, u0)` with the cell-local P1 mass matrix
(`cellview.local_mass_inverse`). A transport step from the t = 0 state with
it as the only source returns that state at every node; a cell-constant
reference, the DG0 rule on the cell means, cancels each cell's net flux alone
and lets the slopes move (`tests/test_dg1_cell_rules.py`). The fixed-front and
free-front masks and `balance`'s forcing act on whole cells.

The reference holds a DG1 state less tightly than a DG0 one at a steep front
cell. On the 32 km probe MAPs (unforced, fixed front, `balance`; run record
`test-32km-dg1-stage2-amb`) `resid` and `clamp` are 0 at every step, and the
neighbour-scale rms of the cell-mean dH is:

| years | DG0, apparent MB | DG1, apparent MB | DG1, without |
|---|---|---|---|
| 1 | 0.0007 m | 0.0034 m | 3.0 m |
| 2 | 0.0012 m | 0.0265 m | |
| 5 | 0.0023 m | 0.0410 m | 9.3 m |
| 10 | 0.0205 m | 0.2196 m | |

The DG1 excess starts at one floating front cell (2188, -1296 km), whose
vertices hold 583, 96 and 0 m and which moves at 1,306 m/yr. Without the
reference that cell refills to 304 m and slows to 115 m/yr within a year, in
both spaces; with it, DG0 holds the cell exactly for 15 years, and under DG1
it thins by 35 m over three years, speeds up to 1,600 to 2,000 m/yr, and
passes the change to its floating neighbours (dh -38 to +51 m by year 10)
while VAF moves by -0.04 mm. With the slope scaling of the positivity floor
turned off (an experiment, `test-32km-dg1-stage2-amb-noscale`) the cell holds
through the first year and then drifts the same way, its zero vertex falling
to -49 m: the floor brings the onset forward, and the held state itself is
unstable at that cell under DG1.

**Provenance.** A DG1 MAP is named `_dg1` and records `geometry_space=dg1`; a
DG1 state seeds only a DG1 run and a DG1 run starts only from one, on every
path (`runconfig.geometry_space_mismatch`). An inversion is the one exception:
it may take a warm start's controls and fluidity prior from the other space
and refit them on its own geometry, never taking the warm start's geometry
(`runconfig.warm_start_across_spaces`). At 32 km a DG1 refit from the DG0
probe MAP's controls started at misfit 6,938, against 35,793 cold.
`geometry_space` is an objective key of the handoff, so the refit's first
link runs with `ISMIP7_WARM_START_STRICT=0` (`inversion.sbatch` exports 1),
as any first link that changes the objective does.

**What `dg1` refuses** (`runconfig.dg1_refusals` and the drivers):
sub-element friction and `ISMIP7_LEGACY_TRANSPORT`, which integrate a
continuous surface; an inversion taking a warm start's geometry from another
mesh, where a discontinuous linear field interpolated node by node samples one
source cell per target node; adaptive transfer;
`ISMIP7_SUBSTEP_TRACE`, which prints one value per cell; and timing caches:
one is prepared on another mesh and its contract
(`timing_campaign.cache_required_fields`) names the DG0 `N_ref`, so
`prepare_timing_cache.py`, `redistribute_checkpoint.py --manifest` and a
forward reading a cache all stop under `dg1`.

**The t = 0 imbalance is the data's.** The cell-scale part of the initial flux
divergence, which the apparent mass balance cancels, is nearly the same in
both spaces. Neighbour-scale rms: the area-weighted rms over interior ice of
a cell's upwind divergence minus the mean of its three facet neighbours',
BedMachine `vertex` sample:

| velocity | mesh | DG0 | DG1 |
|---|---|---|---|
| MEaSUREs v2 | 32 km | 3.81 m/yr | 3.45 m/yr |
| MEaSUREs v2 | 8 km | 3.73 m/yr | 3.98 m/yr |
| the 32 km DG0 Budd MAP's, in both | 32 km | 9.95 m/yr | 10.43 m/yr |

DG0's own truncation error accounts for 2.31 m/yr rms of the last row; the rest
is BedMachine's thickness against the inverted velocity at the cell scale.

**The 32 km probe.** Budd MAPs inverted in each space under one objective
(`vertex` sampling on the buffer-0 32 km mesh, `exact_front` 2), then 40
unforced years from each with the fixed front and no apparent mass balance,
and DG0 again with it as the smooth reference
(`antarctica/scripts/geometry_noise.py`; run records `inversion-32km-budd-*-dg1probe`
and `test-32km-dg1probe-*`). Neighbour-scale rms of the cell-mean dH over
interior ice:

| years | DG0 | DG1 | DG0, apparent MB |
|---|---|---|---|
| 1 | 4.1 m | 3.0 m | 0.0007 m |
| 5 | 12.4 m | 9.3 m | 0.002 m |
| 15 | 22.2 m | 17.4 m | 0.15 m |
| 40 | 33.7 m | 27.4 m | 0.36 m |

| | DG0 | DG1 |
|---|---|---|
| MAP velocity misfit | 1,919 | 1,368 |
| correlation with year 1, at 15 / 40 years | 0.81 / 0.67 | 0.75 / 0.62 |
| outflux 2015 / 2035 / 2054 (Gt/yr) | 1,503 / 1,133 / 1,049 | 1,306 / 1,500 / 1,301 |
| Newton iterations a solve; rescues | 8.2; 1 | 4.3; 0 |
| diagnostic solve | 0.53 s | 0.27 s |

`resid` was 0.0000 Gt at every one of the 1,600 steps of each run. The
apparent mass balance removes 99 % of the cell-scale noise, and in both spaces
the rest grows where the first year put it (correlation 0.75 to 0.81 at 15
years): at this resolution the noise is the initial imbalance. DG1 starts with
27 % less of it, from a MAP that fits the velocity 29 % better, and the gap
narrows to 19 % by year 40, with the same pattern at the margins. The
production noise is at cells 16 to 32 times finer; this probe measures the
32 km cell scale.

**Forced, at 25 km.** IU's 25 km buffered mesh, `vertex_front`, the CESM2-WACCM
historical from a 12-year backdated start and then ssp585 to 2100 with the
collapse mask (`mask_front`), the SMB-elevation feedback on, the 25 km melt
calibration fitted on `dg0` cells, fixed front, from point-velocity MAPs
refit in each space (records `test-25km-dg1-stage2-*`):

| | DG0 | DG1 | DG1, apparent MB |
|---|---|---|---|
| backdate (cells, mass) | 2,992, +566 Gt | 2,992, +566 Gt | |
| t = 0 melt | 997.8 Gt/yr | 997.8 Gt/yr | |
| historical VAF, 2003 to 2015 | +30.8 mm | +28.1 mm | +2.0 mm |
| ssp585 VAF, 2015 to 2100 | +387.6 mm | +334.8 mm | |
| collapse at 2100 (flagged, removed) | 417, 328 cells | 433, 377 cells | |
| isschecker 0.5.1, historical | 0 errors | 0 errors | |

`resid` is 0 at every step of every run, with no rescue. The backdate and the
t = 0 melt are DG0's to the reported digit, as the cell means make them; the
later differences carry the two MAPs' (velocity misfit 3,339 under DG0, 4,475
under DG1). The ssp585 output passes the checker's naming, numerical and
spatial tests and fails only its time axis, which stops at 2100.

## Incompatibilities

- `ISMIP7_LEGACY_TRANSPORT=1` requires `ISMIP7_GEOMETRY_SPACE=cg1` (it
  re-projects CG1 `h` <-> DG0 `h_dg` every step; under DG0 they are the same
  Function and the projection would be self-referential). Enforced.
- A run must not change geometry space mid-trajectory. Restart checkpoints
  record `geometry_space`; a mismatch between `dg0` and `cg1` projects and
  warns, and one involving `dg1` stops the run
  (`runconfig.geometry_space_mismatch`).
