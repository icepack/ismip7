# DRAFT ISSUE — not yet posted

**STATUS UPDATE (Aug 30): the net term was implemented, validated, and then
EXCLUDED from the production objective by decision.** Validation held (net
+604 → −54 Gt/yr vs observed −73 at 2500 m, ~10% pointwise-RMS cost, stable
10-yr forward), but three arguments removed it from the default: assimilating
the integrated trend forfeits it as independent validation; proj − CTRL
differencing under the apparent-MB scheme cancels the shared frozen drift in
the reported signal anyway; and the integral is bought with regionally
structured adjustments to what is demonstrably div(h·u) data noise (interior
response times leave no room for real decadal interior signal). The term
remains in-tree behind ISMIP7_DHDT_NET_SIGMA (default 0) as an experiment
knob, and the per-iteration net= diagnostic always prints. If posted, this
issue should be framed as a record of the investigation and the decision, not
as a proposal.

Target repo: `icepack/ismip7` (or `hoffmaao/ismip7` if we want it fork-local
first). Delete this header when posting.

---

**Title: Initialized state gains mass where observations show loss; proposal:
net mass-balance term in the inversion objective**

## Problem

Our forward runs gain mass where Antarctica should lose it. We have now
isolated this to the **initial state itself**, before any forcing or transport
enters: taking one prognostic step of the DG0 upwind transport from the
converged 2500 m transient MAP and integrating the thickness tendency over
grounded, altimetry-observed ice gives

| | model | observed (Smith et al. 2020) |
|---|---|---|
| net dH/dt | **+604 Gt/yr** | **−73 Gt/yr** |
| grounding-line band (HAF 0–50 m) | +127 | −5.9 |
| interior (HAF > 200 m) | +447 | −52 |

(225,717 cells, 11.2 M km²; the observed integral is consistent with IMBIE-3
for the 2003–2019 era, so the target field is sound.)

Two earlier links in this chain are already fixed and are not the residual
cause: the calving-terminus BC that covered only ~5 % of the front
(`0cd72c9`), and the initialization tendency error that `ISMIP7_APPARENT_MB`
freezes into a constant per-cell source. What remains is that **nothing in the
objective ever asks the initial state to lose mass at the observed rate**.

## Why the existing pointwise dH/dt term is not enough

The transient inversion (`ISMIP7_DHDT_WEIGHT`) penalises
`(dh/dt_model − dh/dt_obs)²` cell-by-cell, and at 32 km it did fix the net
(−518 → −26 Gt/yr). At 2500 m it converged properly — final pointwise RMS
2.2 m/yr, consistent with its χ² — **and still left +604 Gt/yr of net bias**.
The arithmetic explains it: +604 Gt/yr over 11.2 M km² is a mean of
**+0.06 m/yr per cell, i.e. ~3 % of the 2.2 m/yr local residual**. A pointwise
L2 misfit reduces variance; its gradient barely feels a mean that small
relative to the noise. At 2500 m the local residual is dominated by
grid-scale flux-divergence structure at the grounding line that we have
measured to be non-convergent under refinement (pointwise max |div(h u)|
grows 293 → 1769 m/yr from 32 km to 2500 m while area-integrals stay stable),
so the optimizer spends the controls on unfixable local noise and the
systematic bias survives.

## Proposal

Add an **integral (net mass balance) term** to the inversion objective,
alongside the pointwise term:

```
J_net = 1/2 · [ ( ∫_grounded,observed (dh/dt_model − dh/dt_obs) dA ) / σ_net ]²
```

with `σ_net` of order 25 Gt/yr (IMBIE-scale uncertainty on the integrated
balance). This penalises exactly the quantity that is wrong — the net — with
a gradient that does not vanish into the pointwise noise.

Implementation notes (already prototyped on `antarctica-n3`):

- `dh/dt_model` is the same one-step implicit-Euler DG0 upwind tendency the
  pointwise term uses (the model's own transport operator, so the inversion
  is charged for the divergence its own scheme will produce), so `J_net`
  costs **no additional solve**.
- The squared integral is taped by projecting the masked residual onto the
  `R` (Real) function space — a one-dof solve whose solution is the domain
  mean — and integrating a spatially-constant quadratic of it. Everything
  stays in UFL and differentiates through tlm_adjoint's existing
  `EquationSolver` machinery; no functional algebra is required.
- Knob: `ISMIP7_DHDT_NET_SIGMA` (Gt/yr; `0` disables). Recorded in the MAP
  checkpoint attributes alongside the other objective provenance
  (`misfit_norm`, `gamma_*`, `dhdt_weight`), since MAPs produced with and
  without it are not interchangeable.
- Restricted to grounded, observation-covered ice for the same reasons as the
  pointwise term: shelf altimetry is firn/tide/ocean-confounded, floating
  thickness change does not move VAF, and ocean melt is identically zero on
  grounded ice so no melt model enters the constraint.

## Natural refinement (follow-up, not this change)

A per-IMBIE-basin version — 16 scalars instead of 1 — would prevent
compensating errors between basins (e.g. a spurious East Antarctic gain
cancelling a real Amundsen loss inside a single net). The single-integral
version is the minimal change that addresses the observed artifact and is
already testable.

## Validation status

- 32 km adjoint test (cold start, sigma_net = 25 Gt/yr, 7 L-BFGS
  iterations): the adjoint descends through the R-space projection at
  unchanged cost (~1.3 s vs an ~80 s forward), and the net moves an order of
  magnitude toward the observed value while the velocity and pointwise terms
  also fall -- no term is traded away:

  | iter | vel chi^2 | dhdt chi^2 | net (Gt/yr) |
  |---|---|---|---|
  | 1 | 8154 | 3977 | -3289 |
  | 3 | 4137 | 2415 | -1430 |
  | 4 | 3338 | 1981 | +304 |
  | 7 | 1675 | 1120 | -269 |

  One implementation note: tlm_adjoint's linear-solver cache cannot copy the
  R-space python-type PETSc Mat, so the one-dof solve opts out of caching
  (`cache_jacobian=False, cache_adjoint_jacobian=False`) -- free, since a 1x1
  solve gains nothing from a cache.
- The protocol note on the Smith 2003–2019 observation window overlapping the
  post-2015 projection era (declared in `icepack2_tools/obs_dhdt.py`) applies
  to this term exactly as to the pointwise one.

## Relation to other work

- Builds on the transient (dH/dt) inversion merged in hoffmaao/ismip7#6.
- Complements — does not replace — `ISMIP7_APPARENT_MB`: with a correctly
  signed initial tendency, the frozen `a_ref` correction becomes genuinely
  small instead of masking a bias of hundreds of Gt/yr.
- Orthogonal to the mesh/sidecar convention merge from
  `upstream/integration` (Dan/David).
