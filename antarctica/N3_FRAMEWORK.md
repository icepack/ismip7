# Branch `antarctica-n3`: standard Glen n=3

This branch is a parallel line to `antarctica`, identical in every respect
**except the ice flow-law exponent**: it runs the composite viscous rheology
at **n = 3 (standard Glen)** instead of n = 4 (Goldsby-Kohlstedt dislocation
creep). n = 3 is the ISMIP6/ISMIP7 default and the more directly comparable
configuration; n = 4 remains on `antarctica` untouched.

## What changed vs `antarctica`

The n = 4 assumption lived entirely in two environment-variable defaults, so
the branch difference is small and self-documenting:

| knob | `antarctica` (n=4) | `antarctica-n3` (n=3) | why |
|---|---|---|---|
| `ISMIP7_N_FLOW` | `4.0` | `3.0` | the flow exponent |
| `ISMIP7_A4_FACTOR` | `10.0` | `1.0` | `A0 = rate_factor(260 K)` **is** the n=3 fluidity, so the composite main term is plain Glen and needs no rescale. At n=4 the factor lifts `A_3` to `A_4` so `A_4·τc⁴ ≈ A_3·τc³` at `τc`. |

Both defaults are module constants (`N_FLOW_DEFAULT`, `A4_FACTOR_DEFAULT` /
`A4_FACTOR_N4`) in `simulation.py` and `inversion_icepack2.py`; a run still
overrides them with the env vars. The `a4_factor` default is **derived from**
`ISMIP7_N_FLOW` (`a4_factor_default()`: 10.0 at n=4, 1.0 otherwise), so
`ISMIP7_N_FLOW=4` alone is enough to reproduce the `antarctica` rheology and
cannot silently pair the legacy untagged n=4 MAP with a 1.0 prefactor; an
explicit `ISMIP7_A4_FACTOR` still wins. The inversion and the forward that
loads its MAP **must** use the same pair.

## MAPs

The inversion is n-specific (the fluidity, the `τc^{n-1}` linearization, and
the n-continuation target all depend on n), so n = 3 needs its **own** MAPs.
They carry an `_n3` filename tag so they coexist on disk with the n = 4 MAPs:

- n = 4: `inversion_icepack2_budd<geom>_<lc>.h5` (flow exponent untagged, legacy)
- n = 3: `inversion_icepack2_budd_n3<geom>_<lc>.h5`

`<geom>` is the geometry-space tag that follows the n tag - `_dg0` under the
`ISMIP7_GEOMETRY_SPACE` default, empty for `cg1` - so the n = 3 Budd MAP is
`inversion_icepack2_budd_n3_dg0_<lc>.h5` unless the geometry space is
overridden; see `../GEOMETRY_DISCRETIZATION.md` for why a MAP is only valid for
the geometry space it was inverted under. Both tags are produced by
`icepack2_tools/naming.py` (`map_n_tag()`, empty at n = 4 for backward
compatibility, and `map_basename()`), which the inversion, the forward, the
preflight and the launch gates all go through. MAP h5 files are gitignored and
regenerated per machine.

## Physical fluidity prior (the n=3 inversion method)

At n = 3 a constant fluidity baseline forces the control `phi = log(A/A0)` to
carry all the spatial fluidity structure, which blows up (`phi` to ±36) since
n = 3 lacks the n = 4 `a4×10` boost. So the n = 3 inversion follows the
`mismip_time-dependent-da` (Recinos/fenics_ice) method: it regularizes the
**deviation from a physical prior mean**, not amplitude.

- **Fluidity**: `phi = log(A / A_prior)` where `A_prior(x)` is a
  thermomechanical field from a fixed-velocity depth-averaged enthalpy solve
  (`icepack2_tools/thermo_model.py`), driven by the observed geometry/velocity
  and the ISMIP7 `tas` mean-annual surface temperature. Inspect it with
  `antarctica/scripts/thermo_prior.py`.
- **Friction**: `theta = log(C / C_w0)` on the balance anchor (already physical).
- **Regularization**: Whittle-Matérn `gamma·(theta² + L²|∇theta|²)`
  (`icepack2_tools/prior.py`) at a physical `gamma` whose default is coupled to
  the misfit normalization `ISMIP7_MISFIT_NORM`, because normalizing the misfit
  by the observational `sigma²` would otherwise weaken the prior by the same
  factor. Env: `ISMIP7_GAMMA_THETA`, `ISMIP7_GAMMA_PHI`, `ISMIP7_L_REG`
  (7.5 km), `ISMIP7_FLUIDITY_PRIOR` (thermo|legacy); the defaults are tabulated
  in `README.md` §4 under "Environment knobs (inversion)".

The MAP stores `A_prior`; the forward loads it and rebuilds `A = A_prior·exp(phi)`.
Result (32 km): controls physical (theta/phi p99 ≈ 3-4, was 14-17), misfit
n=4-comparable, forward reproduces the inversion. Those misfit figures are on
the legacy dimensional scale (`ISMIP7_MISFIT_NORM=none`); the default chi^2
misfit is dimensionless and its values do not compare to them.

**Run it** (small mesh → FEW ranks; MUMPS is fragile at ~90 vertices/rank, and
a failed first solve now raises rather than writing a garbage theta=phi=0 MAP):

```
# 32 km n=3 Budd MAP (-n 8; runs the thermo prior + Whittle-Matern by default)
OMP_NUM_THREADS=1 ISMIP7_FRICTION=budd ISMIP7_LC=32000 ISMIP7_LC_COARSE=320000 \
  ISMIP7_N_FLOW=3.0 \
  ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_320000_32000.msh ISMIP7_MAXITER=500 \
  mpiexec -n 8 python antarctica/scripts/inversion_icepack2.py
```

The forward path (CTRL, projections, matrix) then runs as on `antarctica`,
with `ISMIP7_RUN_TAG=n3` set so the n=3 outputs (`hist_<esm>_n3_...`,
`ctrl2015_<esm>_n3_...`, ...) coexist with the n=4 results and the
historical -> projection / CTRL restart chain stays within the n=3 line
(see the run-management flags in `antarctica/README.md` §6).
The MAP clip default is n-aware: 10 when `ISMIP7_N_FLOW != 4` (the physical
n = 3 controls reach ~8) and 6 at n = 4 (tuned for garbage outliers); set
`ISMIP7_MAP_CLIP` explicitly to override.

See `COMPOSITE_RHEOLOGY.md` for the full formulation.
