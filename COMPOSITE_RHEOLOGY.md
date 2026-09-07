# Composite rheology for icepack2 dual-form SSA

> **Branch `antarctica-n3`: standard Glen n=3.** This branch runs the flow
> law at `n = 3` with `a4_factor = 1` (the composite "main" term is just
> Glen, since `A0 = rate_factor(260 K)` is already the n=3 fluidity). The
> n=4 Goldsby-Kohlstedt composite (`a4_factor ≈ 10`) is on the `antarctica`
> branch. Everything below is written for the general exponent `n_flow`;
> read `n_flow = 3`, `a4_factor = 1` here. MAPs inverted at n=3 carry an
> `_n3` filename tag so they coexist with the untagged n=4 MAPs (naming rule:
> `antarctica/N3_FRAMEWORK.md`); the forward and its inversion must use the
> same `ISMIP7_N_FLOW` / `ISMIP7_A4_FACTOR`. The `a4_factor` default is
> derived from `ISMIP7_N_FLOW` (10.0 at n=4, 1.0 otherwise), so setting
> `ISMIP7_N_FLOW=4` alone reproduces the `antarctica` rheology; an explicit
> `ISMIP7_A4_FACTOR` still wins.


The standard formulation lets the membrane stress and basal shear stress
become unconstrained where the local ice thickness vanishes, so the SNES
Jacobian goes singular at the calving front. To allow `h → 0` cleanly we
add a small linear (n=1) regularization to both the viscous and friction
contributions of the action functional, with the regularization using a
**constant reference thickness** `H_ref > 0` for the viscous term.

This mirrors `icepack2/test/dome_test.py` and is used in
`antarctica/scripts/simulation.py` (forward runs) and
`antarctica/scripts/inversion_icepack2.py` (MAP estimation).

## Notation

| symbol | meaning | units |
|---|---|---|
| `u`     | depth-averaged velocity                                  | m/yr      |
| `M`     | depth-integrated membrane stress tensor                  | MPa       |
| `τ`     | basal shear stress vector                                | MPa       |
| `h`     | ice thickness (geometry space: DG0 by default)            | m         |
| `s`     | upper surface elevation                                  | m         |
| `b`     | bed elevation                                            | m         |
| `A`     | depth-averaged ice fluidity = `A_prior · exp(φ)`         | MPa⁻ⁿ·yr⁻¹|
| `K`     | sliding coefficient = `K_base · exp(−n·θ)`               | (yr/m)·MPa⁻ⁿ |
| `K_base`| baseline sliding coefficient = `u_c / (φ_eff · τ_c)ⁿ`    | (yr/m)·MPa⁻ⁿ |
| `φ_eff` | effective-pressure fraction in `[0.01, 1]`               | —         |
| `θ`     | log-friction control field                                | —         |
| `φ`     | log-fluidity control field                                | —         |
| `n_flow`| dislocation-creep flow exponent (= 4, Goldsby-Kohlstedt) | —         |
| `m_slide`| Weertman sliding exponent (= 3)                          | —         |
| `τ_c`   | reference stress for normalization (= 0.1 MPa)            | MPa       |
| `u_c`   | reference speed = mean of `|u_obs|`                       | m/yr      |
| `α`     | composite regularization weight                           | —         |
| `H_ref` | constant reference thickness for the linear viscous term  | m         |

The Glen rate factor `A₀ = A(260 K)` follows `icepack.rate_factor`.

## Action functional

The SSA momentum-balance problem is recovered as the stationary point of

```
L(z) =   ψ_visc^{(n_flow)}(M, h; A)            # Goldsby-Kohlstedt dislocation creep (n=4)
       + α · ψ_visc^{(1)}(M, H_ref; A_lin)     # diffusion-like regularizer (n=1)
       + ψ_fric^{(m_slide)}(τ; K)              # Weertman sliding (m=3)
       + α · ψ_fric^{(1)}(τ; K_lin)            # linear sliding regularizer (m=1)
       + Π_mom(u, M, τ; h, s)
       + Π_term(u, h, s)                       # only on calving boundaries
```

The flow exponent and the sliding exponent are decoupled (we use `n_flow=4`
for the depth-averaged GK dislocation creep but keep the basal-slip law at
Weertman `m_slide=3`).

with `z = (u, M, τ)` the mixed unknown.

### Glen power dissipation terms

For exponent `p ∈ {1, n}` the (dual) viscous and friction powers used by
`icepack2.model.minimization` are

```
ψ_visc^{(p)}(M, h; A) = ∫  2 · h · A / (p+1) · |M_dev|^{p+1}  dx

ψ_fric^{(p)}(τ; K)   = ∫    K / (p+1) · |τ|^{p+1}            dx
```

where `|M_dev|² = ½(M:M − tr(M)²/(d+1))` is the squared deviatoric
membrane stress norm. The case `p = 1` collapses to a quadratic form in
`M` and `τ` (positive definite for any `K, A > 0`), which is what makes
the regularization useful.

### Momentum balance

```
Π_mom = ∫ [ −h M : ∇_sym u  +  ( τ − ρ_I g h ∇s ) · u ] dx
      + ∫ ρ_I g ⟨h⟩ ⟦s · ν⟧ ⟨u⟩ dS
```

`⟨·⟩` is the facet average, `⟦·⟧` is the jump. The boundary contribution
is automatic from the action.

### Calving terminus (optional)

```
Π_term = ∫_{calving} ½ [ ρ_I g h² − ρ_W g d² ] (u · ν) ds
```
with `d = min(0, s − h)` the draft below sea level.

## Composite specification

The main (nonlinear) rheology uses the spatially-varying inverted controls:

```
A      = A_prior · exp(φ)                      A_prior(x) = fluidity prior mean
K      = K_base · exp(−m_slide · θ)
K_base = u_c / (φ_eff · τ_c)^{m_slide}
φ_eff  = max(0.01, 1 − ρ_W g max(0, −b) / (ρ_I g max(H, 1)))
```

`A_prior(x)` is the thermomechanical fluidity prior mean the inversion
computes and stores in the MAP checkpoint, so the control
`φ = log(A / A_prior)` is a deviation from a physical field (see
`antarctica/N3_FRAMEWORK.md` for the prior method). MAPs that predate the
physical prior carry no `fluidity_prior` and the forward falls back to the
legacy constant baseline `a4_factor · A₀` (`a4_factor ≈ 10` on `antarctica`,
so `A_4·τ_c⁴ ≈ A_3·τ_c³` at `τ = τ_c`; `1` on this branch).

The linear rheology is obtained by linearizing the main forms about the
reference stress `τ_c`:

```
A_lin = A · τ_c^{n_flow − 1}                   # so A_lin · τ_c = A · τ_c^{n_flow} at τ_c
K_lin = u_c / (φ_eff · τ_c) · exp(−θ)          # so K_lin · τ_c = u_c at τ_c
```

so that at `|τ| = τ_c` and `|M_dev| = τ_c` the linear and main
contributions agree up to the `α` factor. With `α ≪ 1` the dislocation
creep part dominates everywhere `h > 0`; the linear part is what holds
`M` and `τ` in check where `h → 0`.

**Key trick:** the linear viscous term uses `H_ref = 100 m` (a
constant) instead of the evolving `h`. So at every node the linear
viscous action contributes

```
α · 2 · H_ref · A_lin / 2 · |M_dev|² = α · H_ref · A_lin · |M_dev|²
```

which is positive definite in `M` regardless of the local `h`. The Glen
viscous term is multiplied by `h`, so it vanishes at the calving front;
the linear regularization is the only thing pinning `M` there.

The linear friction term uses the same `K_lin` everywhere (no thickness
dependence — `ψ_fric` doesn't carry an `h` factor in either form).

## Defaults

| env var                  | default (this branch) | used in                              |
|--------------------------|-----------------------|--------------------------------------|
| `ISMIP7_N_FLOW`          | `3.0` (n=4 on `antarctica`) | both                           |
| `ISMIP7_M_SLIDE`         | `3.0`                 | both                                 |
| `ISMIP7_A4_FACTOR`       | derived from `ISMIP7_N_FLOW`: `10.0` at n=4, else `1.0` | both |
| `ISMIP7_COMPOSITE_ALPHA` | `1e-2` (forward, `budd`/`regularized_coulomb`; `1e-4` legacy) | `simulation.py` / `control/run.py` |
| `ISMIP7_COMPOSITE_ALPHA` | `1e-2` (inversion) | `inversion_icepack2.py`                 |
| `ISMIP7_H_REF`           | `100.0` m          | both                                    |
| `ISMIP7_H_CLAMP_INIT`    | `0.0` m (`budd`/`regularized_coulomb`; `10.0` legacy) | `simulation.py` initial diagnostic only |
| `ISMIP7_H_CLAMP`         | `0.0` m            | both, advection floor + inversion       |

With the residual friction laws (`ISMIP7_FRICTION=budd`, the default, or
`regularized_coulomb`) the forward runs share the inversion's `α = 1e-2`
and start from the true h=0 BedMachine geometry (`h_clamp_init = 0`),
because the MAP was inverted against that geometry. The legacy action-form
path keeps the original split: a 100× stronger `α` in the inversion (every
L-BFGS-B step does a full diagnostic solve, so SNES robustness is worth a
small bias in the recovered `θ` and `φ`) and a 10 m initial clamp in the
forward run.

## Sanity checks

- **Smoke test (synthetic disk, 50 km, h=500 max, h=0 at edge):** SNES
  converges through `n=1 → 3` continuation, `h_min = 0` exactly,
  velocities finite. See `antarctica/scripts/composite_smoke_test.py`
  (in commit history) or the inline test in the change log.

- **Forward run with composite + `h_clamp_init = 10` + `h_clamp = 0`
  (advection):** runs to completion at dt=0.1 over 10 yr on the 2500 m
  mesh. Cells reach 0 cleanly during melt but buffer cells with
  `h_clamp_init = 10` feed mass back through the diagnostic velocity, so
  this is *not* a clean drift test. See
  `antarctica/results/ctrl2015_cesm2_waccm_2500_timeseries.csv` (2026-05-15).

- **Inversion with composite + `h_clamp = 0`:** done - the `_budd` / `_rc`
  MAP checkpoints are inverted against the true BedMachine geometry (h=0
  over the buffered ocean region), and the forward runs that load them
  start with `h_clamp_init = 0` (no initial-thickness clamp).
