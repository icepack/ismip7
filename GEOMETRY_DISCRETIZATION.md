# Geometry discretization and the calving front

`ISMIP7_GEOMETRY_SPACE` selects the finite-element space for the geometry
(`h`, `s`, `b`, and everything derived from them: `phi_eff`, `N_ref`, `C_w0`,
`H_init`, and the forcing fields `accum` / `ocean_melt`). The apparent-MB
reference `a_ref_mb` is a per-cell transport source and stays DG0 either way.

- `dg0` (**default**) - cell-wise geometry. The momentum solve and the mass
  transport use **one** thickness field.
- `cg1` - the pre-August-2026 behaviour, kept for A/B comparison.

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

That lift conserves volume exactly, and in the interior it is a reasonable
reconstruction. On the **domain boundary** its stencil is one-sided: a
boundary node averages only interior cells, so it is pulled toward interior
values. On an unbuffered mesh the domain boundary *is* the calving front, and
the terminus traction

```
f_I - f_W = 1/2 rho_I g h^2 - 1/2 rho_W g d^2
```

goes as `h^2`. Measured at 32 km with the n=3 Budd MAP:

| | front `<h>` | terminus force | front `<u.n>` | outflux |
|---|---|---|---|---|
| before the lift | 105 m | 1.0x | +238 m/yr | 848 Gt/yr |
| after the lift | 200 m | **3.3x** | +698 m/yr | 3980 Gt/yr |

BedMachine's actual ice-front thickness is **145 m** over all fronts and
**167 m** for floating fronts (median 152, p90 306), so the unlifted CG1 field
was ~30% thin and the lifted one ~25% thick - a factor 3.6 spread in the
front driving force from discretization alone, larger than any physics knob in
the model. Force and mass also disagreed: with one velocity field, the outflux
computed from the momentum-solve thickness and from the transport carrier
differed by **22%**, while the volume integrals agreed to 0.000%.

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
surface. With CG1 `s` the jump vanishes and the cell term carries everything;
with DG0 `s` the cell gradient is identically zero (UFL folds it away) and the
facet term carries everything. Both are consistent, but not equally accurate.
Against the exact continuum driving force on a smooth manufactured problem:

| n | CG1 error | rate | DG0 error | rate |
|---|---|---|---|---|
| 8 | 1.20e-02 | - | 5.98e-02 | - |
| 32 | 7.55e-04 | 2.00 | 1.45e-02 | 1.02 |
| 128 | 4.72e-05 | 2.00 | 3.59e-03 | 1.00 |

**CG1 is 2nd-order, DG0 is 1st-order** on the interior driving stress. This is
the real trade: ~1.4% vs ~0.08% error at 32x32 in the interior, bought against
an unbiased calving front. At Antarctic resolutions the bed and thickness data
error dominates that gap, and the CG1 path never achieved its asymptotic rate
at the front anyway (the lift smooths `h` by up to ~1.4 km at the steep PIG
grounding-zone gradient).

Two further consequences to watch, neither yet a demonstrated problem:

1. **The grounding line is a staircase.** `He = smooth_heaviside(haf(H, b))` is
   cell-constant, so GL migration proceeds cell-by-cell rather than
   sub-element. Coarse at 32 km; much less so at 2500 m.
2. **The lift was incidentally damping grid-scale noise.** Removing it removes
   that filter. First-order upwind transport is strongly diffusive at the grid
   scale and should cover it, but this needs measuring rather than assuming.

## Sampling the raster: cell average, not centroid

A DG0 dof sits at the cell centroid, so the obvious
`icepack.interpolate(raster, Q_g)` takes a **one-point sample** of a 500 m
BedMachine raster per cell - at 32 km, one pixel standing in for a 32 km cell.
Because the DG0 driving stress is *entirely* the facet jump in `s`, that
sampling noise is read as slope. Measured at 32 km:

| construction | rms \|jump s\| | peakedness (L4/L2) | front `<h>` | \|driving force\| vs CG1 |
|---|---|---|---|---|
| cg1 | 0 (continuous) | - | 106.9 m | - |
| dg0, centroid sample | 366 m | 2.04 | 209.5 m | +7.7% |
| **dg0, cell average** | **257 m** | **1.49** | **152.9 m** | **+1.0%** |

The centroid version failed to converge in 200 Newton iterations. The cell
average - the L2 projection of the CG1 interpolant, which is what a DG0 field
*means* - is 42% smoother, much less peaked, lands on BedMachine's true front
thickness (median 152 m), and reproduces the CG1 driving force to 1%.
`geometry.sample_to_geometry` does this; do not replace it with a direct
interpolate onto the DG0 space.

The same rule applies to the RACMO SMB climatology, which sets the mass budget
and the `a_ref` balance: `forcing.load_racmo_smb_climatology` cell-averages onto
a DG0 space (`forcing._sample_raster`) and stays the nodal interpolant on CG1.

## Where a CG1 reconstruction is still used, and why that is legitimate

A DG0 field has no pointwise gradient, so three places reconstruct one via
`geometry.cg1_lift` (the same lumped lift, used deliberately and locally):

- `geometry.surface_slope`, used by `dual_friction.weertman_anchor`. `C_w0` is a fixed
  reference *scaling* that defines what `theta = 0` means. It is not a force.
- `forcing.compute_sin_alpha`. Feeds the Burgard melt *parameterization*, not
  the momentum residual.
- The thermomechanical fluidity prior in `inversion_icepack2.py`, which runs on
  CG1-lifted geometry throughout.

None of these appear in the momentum residual or in a mass flux, which is the
property that matters: the terminus traction and the transport flux must
integrate the same field, and they do.

Note the fluidity prior must be lifted, **not** L2-projected. An L2 DG0->CG1
projection overshoots at the front and produced a *negative* fluidity
(`A_prior` min -9.88 against a DG0 range of `[1.0, 446.7]`), which is
unphysical and poisons `log(A/A_prior)`. `cg1_lift` is a convex combination of
cell values, so it cannot overshoot.

## MAPs are not interchangeable

The inversion absorbs the front treatment into `theta`/`phi`: whatever the
momentum balance gets wrong at the terminus, the optimizer compensates for by
adjusting friction and fluidity until the modelled velocity matches
observations. **The t=0 velocity misfit therefore cannot validate the front
treatment** - it will look fine either way. Only prognostic behaviour can:
drift, calving flux against the observed ~1300 Gt/yr, and front-region velocity
against MEaSUREs.

Because of that, MAPs carry a geometry tag (`map_geom_tag()`):

- CG1: `inversion_icepack2_budd_n3_<lc>.h5` (legacy, untagged)
- DG0: `inversion_icepack2_budd_n3_dg0_<lc>.h5`

A forward prefers the tagged MAP for its own geometry space. If only the legacy
one exists it will load, project, and warn loudly - runnable as a smoke test,
but not a result. Driving a DG0 forward with the CG1 MAP at 32 km raises the
initial misfit from 8.6e3 to 1.5e5, which is the size of the inconsistency.

Run the DG0 inversion with:

```
OMP_NUM_THREADS=1 ISMIP7_FRICTION=budd ISMIP7_LC=32000 ISMIP7_LC_COARSE=320000 \
  ISMIP7_N_FLOW=3.0 \
  ISMIP7_GEOMETRY_SPACE=dg0 \
  ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_320000_32000.msh ISMIP7_MAXITER=500 \
  mpiexec -n 8 python antarctica/scripts/inversion_icepack2.py
```

## Open: the per-basin melt K was calibrated under CG1

`antarctica/scripts/calibrate_melt.py` is CG1 throughout - it builds its own
`FunctionSpace(mesh, "CG", 1)` and vertex-samples BedMachine, `sin_alpha` and
the floating mask - and it is not affected by `ISMIP7_GEOMETRY_SPACE`. The DG0
forward evaluates that same per-basin `K` with a cell-wise draft, a cell-wise
`sin_alpha` (CG1 lift, then cell-sampled) and a cell-wise `haf <= 0` floating
mask, so the melt-receiving area shifts by roughly a one-cell band at the
grounding line and at the ice front. At 32 km, where shelves are only a few
cells wide, that can move the integrated shelf melt by a non-trivial fraction
of the 865 Gt/yr the `K_SCALE = 1.26` calibration targets.

Nothing in this change compensates for that, deliberately: the integrated DG0
melt total needs to be checked against the observational target and `K`
recalibrated. The recalibration is tracked separately and belongs against the
2026-07-31 re-release of the ISMIP7 AIS ocean-melt parameterisation toolbox
(new observational constraint datasets, new cold/warm term-3 targets; the
guidance is to re-run the calibration notebook rather than adjust the model),
and the `sin(alpha)` treatment in the quadratic parameterisation is itself
still unsettled upstream. Until then, treat DG0 melt totals as uncalibrated.

## Incompatibilities

- `ISMIP7_LEGACY_TRANSPORT=1` requires `ISMIP7_GEOMETRY_SPACE=cg1` (it
  re-projects CG1 `h` <-> DG0 `h_dg` every step; under DG0 they are the same
  Function and the projection would be self-referential). Enforced.
- A run must not change geometry space mid-trajectory. Restart checkpoints
  record `geometry_space`; a mismatch projects and warns.
