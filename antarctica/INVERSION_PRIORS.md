# The inversion's priors, controls and likelihood

Where the friction and rheology inversion stands after the 26 and 27
September 2026 work on branch `feat/recinos-friction-prior` (merged into
`feat/recinos-upstream-sync` on 27 September). This is the topic doc for the
prior and control choices; the chains themselves are tracked on the board
(issue #24, issue #31). Numbers below were measured in that session unless
marked as read from a log.

## 1. What the inversion minimises

The objective is a velocity misfit plus two Whittle-Matern priors, all inside
the tape the optimiser differentiates (`inversion_icepack2.py`,
`_prior_energy_form`, `forward_total`):

* misfit: the per-node chi-squared against MEaSUREs with the observed sigma
  (floor `ISMIP7_SIGMA_U_FLOOR`, 3 m/yr), plus the log-speed term with its
  frozen weight, multiplied by `ISMIP7_MISFIT_SCALE`. `1` (the old default)
  is the area MEAN chi-squared; `nodes` is the SUM over observed nodes, a
  Gaussian likelihood, which is the balance Recinos et al. (2023) use;
* priors: bi-Laplacian `A M^-1 A` with `A = delta M + gamma K`, delta and
  gamma from a marginal standard deviation sigma and a correlation length rho
  by the closed forms of Villa et al. (2021), `icepack2_tools/prior.py`.

The prior operator does what it says. At the 32 km stage-1 MAP of 26
September a 47-sigma excursion of phi over one shelf cost `reg_phi` 4.9e5,
which is what sigma_phi 0.3 over the roughly 570 correlation patches of that
shelf predicts.

## 2. Controls

| `ISMIP7_FRICTION_CONTROL` | friction | prior mean | notes |
|---|---|---|---|
| `log` | `C = C_w0 exp(theta)` on the balance anchor | the anchor | the pre-September control |
| `sqrt` | `C = alpha^2` | zero | Recinos et al. (2023); a step is additive in sqrt(C), small where C is small |
| `exp` | `C = C_ref exp(alpha)`, one scalar `C_ref` (`ISMIP7_C_REF`, auto = grounded median of the start) | zero, sigma in log units (`ISMIP7_PRIOR_SIGMA_ALPHA`, auto = 1) | 27 September; a step multiplies C by one factor everywhere, so weak beds are as reachable as stiff ones |

The MAP records the control (`friction_control`), and under `exp` the
reference (`friction_c_ref`) with alpha saved as `log_friction` and `C_ref` as
a constant `C_w0`; the forward and `plot_map.py` take `C = C_ref exp(alpha)`
pointwise with no anchor, checkpoints carry it as `friction_exp`, and warm
starts rebase between the three controls (an exp-control chain link resumes
its alpha unchanged). Every forward state checkpoint copies the MAP's
`friction_control`, `friction_c_ref`, `subelement_friction`,
`subelement_scheme`, `subelement_scheme_version`, `exact_front` and
`fluidity_control`, so a restarted link rebuilds the same residual.

`ISMIP7_SUBELEMENT_SCHEME` picks the sub-element grounding scheme under
`ISMIP7_SUBELEMENT_FRICTION=1`: `sep1` (the default, ISSM's
SubelementFriction1: whole-cell quadrature, Weertman drag times the grounded
fraction) or `sep2` (quadrature over the grounded part). The MAP records it as
`subelement_scheme` and the forward follows it; a sub-element MAP without the
record runs SEP2.

The MAP also records which form of its scheme it was inverted under, as
`subelement_scheme_version` (`icepack2_tools/handoff.py`,
`SUBELEMENT_SCHEME_VERSIONS`), a decision of 2 October:

| scheme | version | form |
|---|---|---|
| `sep2` | 1 | the only form |
| `sep1` | 1 | `exp(He θ)` through the shared Weertman closure (PR 158 before `1c13b38`) |
| `sep1` | 2 | `exp(θ)` ungated, as under SEP2 (this code) |

The two SEP1 forms agree where θ is zero, under the sqrt and exp controls. A
record without the version stands for version 1 under SEP2 and for the current
version under SEP1 with the sqrt or exp control. A log-control SEP1 record
without it cannot say which form it was. A forward refuses a MAP of another
version, and that unversioned log-control case; a chain link counts a changed
version as a changed objective. Rice reports no log-control SEP1 MAP written
before the record (PR 158). A change to a scheme's residual at the same
control fields takes a new version.

The fluidity control is `A = A_prior exp(phi)` under every friction control.

`ISMIP7_FLUIDITY_PRIOR`: `pattyn` (default since 27 September) is the rate
factor of the Pattyn (2013) depth-averaged temperature; `thermo` is the
enthalpy model. `ISMIP7_FLUIDITY_CONTROL=floating` lets phi act on floating
ice only, the split the data-assimilating ISMIP6 groups initialised with
(Seroussi et al. 2020, Appendix C); the forward masks phi with the grounded
indicator of the live thickness, in the residual and in the calving law's
`A_map`, so the shelf rheology follows the grounding line. The inversion holds
phi at exactly zero on every node whose cells all have `1 - He` below 1e-3
(`dual_friction.floating_control_nodes`, about 35 m above flotation): it zeroes
the warm start's values there and projects them out of the gradient, and the
MAP records `phi_grounded=zero`. Its `log_fluidity` therefore has no grounded
variation to pair with the friction by mistake, and ice that goes afloat in a
forward takes the prior fluidity (issue #153). Only the L-BFGS-B path
(`ISMIP7_GRAD_PRECOND=none` or `mass`) runs it.
`ISMIP7_INVERT=phi|theta|both` moves one control with the other held as a
fixed coefficient, the objective unchanged, so a staged inversion
warm-starts each stage through the handoff check.

## 3. What was measured

**The prior operator is right; its weight against the data is uncalibrated.** Under
`misfit_scale=1` the prior outweighs the data by about the number of nodes
and every MAP sat on its prior. Under `nodes` with the 3 m/yr sigma a 400 m/yr
error costs 2e4 per node; chi-squared per node at every 32 km MAP is about
3e4, the likelihood is a hundred times the prior, and phi runs to
`exp(phi)` of 1e6 on the shelves unpunished (32 km stage 1; the 2 km sqrt
thermal MAP reached 5e6). Neither weight is calibrated. The discrepancy
principle (chi-squared per node of order one at the MAP, i.e. a sigma that
carries representation error) is the principled middle. The likelihood
weight is uncalibrated, and the error model is Rice's choice.

**The shelves are the misfit.** On the 2 km Recinos chains the grounded ice
is fitted to a median speed error of 3 to 4 m/yr against an observed median
of 8, and 96 to 99 percent of the chi-squared is on floating ice: under the
thermal prior Ross, Filchner-Ronne and Amery run more than 1000 m/yr too
fast. Under the Pattyn prior the same mesh and friction start at a misfit ten
times lower (2 km, iteration 0: 1.67e10 against 1.70e11 node-summed). This
is why the default moved to Pattyn.

**Joint descents die by line search at 32 km.** Every joint (theta, phi)
run reached 8 or 9 iterations and then failed the line search (TAO reason
-6) at trial points where the forward diverges and the continuation rescue
cannot recover; the thermal prior failed at 5. Cutting the iteration cap or
raising it changes nothing; step control on the line search is the lever.

**Staging as run overshoots.** phi-only on the shelves with friction frozen
(20 iterations, no failures) cut the shelf chi-squared from 8.3e8 to 1.4e8
but raised the grounded chi-squared from 4.1e8 to 1.0e9: the shelves went to
`exp(phi)` of 1e6, buttressing vanished and the streams sped up. The
theta-only stage that followed died by line search at 13 iterations. The
mode is sound; the likelihood weight is what let phi run.

**exp against sqrt, early.** On the thermal prior at 2 km the two controls
descend at the same rate per iteration (exp 15 percent in 7, sqrt 20 percent
in 9). At 32 km on the Pattyn prior, exp had the misfit down 51 percent at
iteration 9 against 63 percent for sqrt, with phi moving 100 times more
under exp. No evidence yet that either converges faster; both are limited by
the line search. The Pattyn-prior twins at Rice (section 5) are the clean
comparison.

**Published state.** Two 32 km MAPs saved a velocity with residual 3e35 and
4e30: the state backup tracked the last EVALUATED trial, so after a failed
line search the publishing solve started from a state of other controls and
"converged" in 0 iterations at an atol scaled from its own residual. Fixed:
the evaluation ring keeps the mixed state, the monitor restores the accepted
evaluation's controls and state (pinned), and a forward that reported
convergence at more than `ISMIP7_FNORM_CEILING` (1e4) times the last
accepted residual, floored at the warm start's recorded residual, is a
failed trial. The fix was verified on a 12-iteration
joint run (published residual 76). The stage-2 run showed the second half:
its accepted iterate had a forward that converged only relatively from a
starting residual of 1e17, hence the ceiling.

**Restart handoff.** Every checkpoint records the objective's settings and
the objective at the accepted iterate; a warm start freezes the auto weights
(log-speed weight, and sigma_alpha and C_ref within one friction control,
section 5) from the record, refuses a changed
objective under `ISMIP7_WARM_START_STRICT=1`, and checks that its first
evaluation reproduces the recorded objective. Verified: the stage-2 handoff
reproduced 5.2008401e7 to all printed digits.

**Cost.** The sub-element scheme costs nothing per iteration (790 to 860 s
per iteration at 2 km on 32 ranks, the same as the cell-wise law) but its
kernels take about 4.4 hours to compile, and Firedrake's kernel cache
defaults to the node-local `/tmp`, cold on every allocation. A single-node job
seeds its node-local cache from the shared store `ISMIP7_SHARED_JIT_CACHE`
and copies its new kernels back (`site_core.sh`; the Rice site file names
`/projects/ah301/sw/pyop2-cache`); with the cache warm the first iteration
arrives in 20 minutes.

**L-BFGS-B against TAO at one objective (issue #157).** Three arms
minimised one objective (bi-Laplacian prior, the log friction control, the
`legacy` fluidity prior, SEP2 sub-element friction with the exact front push,
`full_mumps`, `ISMIP7_EVAL_CONTINUATION=0`) on PR 155's `4cf7f0e`: scipy
L-BFGS-B under `ISMIP7_GRAD_PRECOND=none` and `mass`, and TAO lmvm under
`mass_consistent`, the path Rice's chains take. At 32 km they ran
40 iterations from a cold start on 4 ranks of the IU workstation; at 2 km, 15
iterations from Rice's snapshot 0948 on 32 Quartz ranks (records
`test-32km-inversion-opt-*`, `test-2km-inversion-opt-*`). The objective
starts at 7.8499e4 (32 km) and 5.6211e4 (2 km).

| mesh | arm | evaluations | evaluation time | best objective | best at the `none` arm's time | forward failures |
|---|---|---|---|---|---|---|
| 32 km | `none` | 43 | 375 s | 1.6342e4 | 1.6342e4 | 0 |
| 32 km | `mass` | 52 | 450 s | 1.6148e4 | 1.6155e4 | 0 |
| 32 km | `mass_consistent` | 42 | 730 s | 1.7489e4 | 2.2946e4 | 7 |
| 2 km | `none` | 18 | 3356 s | 5.0444e4 | 5.0444e4 | 0 |
| 2 km | `mass` | 21 | 3820 s | 5.1462e4 | 5.2370e4 | 0 |
| 2 km | `mass_consistent` | 21 | 5801 s | 5.3800e4 | 5.5115e4 | 1 |

* TAO ends highest at both resolutions, after equal iterations, after equal
  evaluations and at equal time. On `4cf7f0e` its median iteration cost 1.7
  times the L-BFGS-B arms' median evaluation at 32 km (14.5 s against 8.4 s)
  and 1.3 times at 2 km (245 s against 184 s). Most of that was a residual
  check TAO recorded on its tape, which PR 155 runs untaped since `71a6809`:
  there the same 32 km iterations take 7.9 s a TAO iteration against 7.2 s
  an L-BFGS-B evaluation (`test-32km-inversion-merged-*-check`), so the time
  columns above overstate TAO's cost and its deficit an evaluation stands.
  Its forward failed at eight trial points over the two meshes, each
  recovered by the re-ramp rescue; no L-BFGS-B arm had a forward failure.
* `mass` starts slowly: each early decrease is about four times the last, the
  line search growing a short first step, so the first eight evaluations at
  32 km moved the objective by 4.9 percent and the first six at 2 km by 0.06
  percent. It then descends fastest: over the last six evaluations at 2 km it
  gained 1,877 against 887 for `none`, and at 32 km it finished lowest. Which
  L-BFGS-B arm ends lower over a whole chain is open.
* The TAO arm's iterations 0 to 5 reproduce the issue #156 probe of the same
  configuration (job 10824069) to 7e-15 relative.

The 32 km arms again on PR 155 at `432c831`, under SEP1 (version 2) and the
direct forward, PR 158's defaults, with the log control and with the sqrt and
exp controls at their `auto` prior scales; 40 iterations from a cold start on
4 ranks of the IU workstation, the three arms of a control run together
(records `test-32km-inversion-sep1-opt-*`). Time includes the checkpoint the
scipy path now writes before each evaluation (`ISMIP7_CHECKPOINT_EVERY_IT=1`,
1.5 to 1.7 s here); the last column is the lowest objective an arm had
reached by the `none` arm's total time.

| control | arm | evaluations | time | best objective | best at the `none` arm's time |
|---|---|---|---|---|---|
| log, from 4.9239e4 | `none` | 44 | 311 s | 1.3944e4 | 1.3944e4 |
| log | `mass` | 48 | 336 s | 1.3788e4 | 1.3796e4 |
| log | `mass_consistent` | 41 | 294 s | 1.4022e4 | 1.4022e4 |
| sqrt, from 1.4221e5 | `none` | 43 | 312 s | 4.9213e4 | 4.9213e4 |
| sqrt | `mass` | 57 | 411 s | 4.8482e4 | 4.9457e4 |
| sqrt | `mass_consistent` | 41 | 313 s | 5.5072e4 | 5.5287e4 |
| exp, from 3.8615e4 | `none` | 43 | 297 s | 7.0787e3 | 7.0787e3 |
| exp | `mass` | 51 | 354 s | 6.4786e3 | 6.8559e3 |
| exp | `mass_consistent` | 41 | 307 s | 1.3284e4 | 1.3735e4 |

* No forward failed and no arm took a rescue; TAO's seven rescued failures at
  32 km were under SEP2 on `4cf7f0e`.
* The order holds under every control, `mass` lowest and TAO highest. TAO
  ends 0.6 % above `none` under the log control (7.0 % under SEP2), 11.9 %
  under the sqrt control and 88 % under the exp control, the control new
  chains take.
* A TAO iteration and an L-BFGS-B evaluation with its checkpoint cost the
  same here: 6.75 to 7.20 s against 6.60 to 7.01 s (medians).
* The exp control's cold start climbed the startup ramp on its first rung;
  under SEP2 it diverged near n = 2.1 on every rung (README, "Inversion
  solver").

At 2 km on Quartz under the production settings (`72ac7a1`: the exp control
from Rice's snapshot 0948 with its log-velocity weight, SEP1, the direct
forward, 32 ranks), 60 iterations of each arm, under `full_mumps` and under
`scpc_gamg`, which took the same iterates (README, "Inversion solver"; records
`test-2km-inversion-final-*`). Best objective from 1.2928e5, and the time and
evaluations each arm needed under `scpc_gamg` to reach TAO's 60-iteration
objective:

| arm | evaluations | after 1 h | after 60 iterations | to 8.31e4 |
|---|---|---|---|---|
| `none` | 65 | 6.530e4 | 5.820e4 | 0.47 h, 14 evaluations |
| `mass` | 71 | 8.363e4 | 7.083e4 | 1.05 h, 34 evaluations |
| `mass_consistent` | 79 | 9.485e4 | 8.308e4 | 2.38 h, 79 evaluations |

* No forward failed and no arm took a rescue.
* `none` is lowest at every time and every evaluation count, `mass` second;
  at 32 km `mass` led under every control. An evaluation costs the same on
  all three paths (106 to 115 s under `scpc_gamg`, 155 to 163 s under
  `full_mumps`), so the arms differ in progress an evaluation.
* Every arm was still descending at the 60-iteration cap.
* At 1 km under `scpc_gamg` (from 6.2094e5, records
  `test-1km-inversion-final-{tao,none}-scpc-gamg`), an evaluation costs the
  same on both paths (69 and 70 s); in the 7 evaluations and 499 s of TAO's
  whole run `none` reached 2.554e5 against TAO's 3.156e5, passing TAO's final
  objective at its sixth evaluation, and 1.497e5 after 13.
* 3 October: the chains take `none` (section 4).

## 4. Decisions taken

* 26 September: no dH/dt term in the inversion (the inversion carries no
  forcing); the chain default is `ISMIP7_DHDT_WEIGHT=0`.
* 26 September: the sub-element grounding scheme and the exact front push
  are on for every new inversion and the MAP records them.
* 27 September: the exp friction control replaces the square for new
  chains; the square chains already running were left to finish.
* 27 September: the Pattyn temperature is the fluidity prior; the thermal
  model stays available.
* 30 September: each objective evaluation solves directly at n=3 from the last
  converged state (`ISMIP7_DIRECT_FORWARD`), a failed trial gets one rescue
  rung (`ISMIP7_TRIAL_RESCUE_RUNGS`), and the inversion checkpoints every
  accepted iterate with its mixed state, so a chain link skips the ramp.
* 1 October: SEP1 is the inversion's sub-element scheme, and the SEP2 chains
  were stopped. On the same 2 km start SEP1 reached misfit 3.42e9 in 21
  accepted iterations with no failed line-search trial; SEP2 reached 4.21e9
  with nine (read from the NOTS logs).
* 3 October: the inversion's linear solver defaults to `scpc_gamg`
  (`ISMIP7_INVERSION_LINEAR_SOLVER`), and the chains run L-BFGS-B without a
  metric (`ISMIP7_GRAD_PRECOND=none`, the code default) in place of TAO with
  the mass-consistent metric. On the final Quartz round at 2 km and 1 km
  (section 3, issue #157; README, "Inversion solver") `scpc_gamg` took the
  same iterates 1.35 to 1.5 times faster at 2 km and 4.6 times at 1 km, and
  L-BFGS-B reached TAO's 60-iteration objective in 14 evaluations against 79.
  A chain switching optimizer mid-way runs one link with
  `ISMIP7_WARM_START_STRICT=0`: `grad_precond` is an objective key.

## 5. The forward this hands to: the 2003 chain

The forward protocol moved while this work was under way, and the MAPs here
feed it. Indiana's work since 22 September, merged into `upstream/main` and
into this branch: the historicals and OCX start in 2003 from the 2015
geometry backdated by the Smith mean dH/dt (issue #117, `ISMIP7_GEOMETRY_BACKDATE`);
a control or projection branches only from a historical that reached 2015,
and the runner queues them once it has (`simulation.historical_endpoint`);
the controls are forced with each ESM's ctrl ocean (issue #107); the 32 km
cores and the controls with the surface-elevation feedback on and off from
the 2003 start are recorded in `antarctica/reports`, and the 25 km rehearsal
of the whole matrix (issue #138) is on a branch of `icepack/ismip7` that is
not in this merge. What it means for a MAP from this branch: it is inverted
on the 2015 geometry, the forward backdates that geometry on a cold start,
and every MAP attribute the residual needs (section 2) travels through the
chain's restart checkpoints. The 32 km probes in section
3 start in 1850 and therefore ran with the backdating off.

Two measurements from those probes for whoever restarts a sub-element MAP:
on a restart the saved mixed state's residual under the rebuilt form is 1e8
before the runner re-solves it (2.7e8 for the exp MAP, 1.5e8 for the sqrt one,
recorded residuals 0.2 to 4), because the sub-element quadrature is rebuilt
from the loaded geometry only when the diagnostic solve starts; the re-solve
converges and the budget continues without a break, so the cost is one extra
Newton solve per link. An auto `prior_sigma_alpha` (and `friction_c_ref`) is
reused across a warm start only within one friction control; a warm start
from another control derives its own auto value.

## 6. The relaxed re-inversion

An initial-state option beside each production MAP (issue #162), settled
at IU on 6 October 2026, in `icepack2_tools/relaxation.py`,
`antarctica/scripts/relaxation/run.py` and the inversion's warm start; the
commands are in the README, "The relaxed initial state".

| step | what | settings |
|---|---|---|
| relaxation year | cold start from the MAP on its own mesh, 2014.0 to 2015.0 | one year of Smith dH/dt undone on grounded ice; half the production step (0.0125 yr at 0.025); OCX's 2014 forcing; apparent MB off; front pinned |
| re-inversion | warm start from the year's end state on the same mesh | 250 L-BFGS-B iterations; the MAP's objective, carried by the end state and held by the strict handoff; the relaxed geometry; the MAP's `θ` and `φ` |
| forwards | `ISMIP7_INVERSION` names the relaxed MAP | 2003 start backdated 12 years, as from any MAP |

Choices, IU's on 6 October:

* **OCX's 2014 forcing** for the year: observational, the same for every ESM,
  and what OCX itself runs through 2014.
* **The front stays pinned**, so the re-inversion keeps the MAP's extent.
* **`θ` stays the MAP's.** On grounded ice `C = C_w0 exp(θ)`, and the
  anchor `C_w0 = τ_d / |u_obs|^(1/m)` follows the geometry through the
  driving stress `τ_d`. The re-inversion builds its anchor on the relaxed
  geometry, so the friction it starts from is the MAP's times
  `R = τ_d(relaxed) / τ_d(2015)`. Where basal drag carries the driving
  stress the Weertman speed is `|u_obs| exp(-m θ)`, independent of the
  geometry, so the starting velocity stays the MAP's. Rebasing `θ` by
  `ln(1/R)` would hold the MAP's friction instead; at fixed friction that
  speed scales as `R^m`, and with `m = 3` a 10 % change in driving stress
  moves it by about 30 %, and the rebased `θ` would carry grid-scale
  structure the bi-Laplacian prior penalises. Under the exp and sqrt
  friction controls the friction has no anchor and the two are the same.
  The re-inversion logs the distribution of `ln R` over grounded dofs.
* **The re-inversion runs on the MAP's mesh** (2 km), and every forward takes
  only its controls, onto its own mesh's BedMachine geometry (`init_state`
  `relaxed-controls`): a 1 km forward through the usual transfer, a 2 km one
  directly. The relaxed thickness serves the re-inversion alone, so every
  resolution starts from the 2015 geometry the unrelaxed MAP's forwards start
  from.

The record: the end state carries `relaxation_end_state`, the source MAP and
its sha256, the year, the step, the forcing and the free-surface
stabilization the year ran under (`relax_fssa_theta`, 0 when off, and the
resolved `relax_fssa_reference`, `none` when off) as `relax_*`, and the MAP's
`handoff.OBJECTIVE_KEYS`, leaving the objective value behind, so the
re-inversion skips the handoff gap check and logs its first objective
against the source MAP's. The relaxed MAP records
`geometry_source_method = relaxed-forward-v1` and the end state's name and
sha256 on every checkpoint, and a chain link that resumes it inherits them.

## 7. Where the chains are (27 September, 22:00 EDT)

| job | mesh | control, prior | state |
|---|---|---|---|
| 1638543 | 2 km | exp, Pattyn | running, iteration 6 with 4 rescues (read from the log) |
| 1643735 | 2 km | sqrt, Pattyn | running, iteration 2 (read from the log) |
| 1638542 | 2 km | exp, thermal | cancelled at iteration 7 when the prior moved to Pattyn |
| 1635853 | 2 km | sqrt, thermal | ended at iteration 9, MAP saved |
| 1632739 | 2 km | sqrt, thermal, with dH/dt | ended at iteration 23, MAP saved |
| 1631512-14 | 2 km | log, thermal, prior sweep sigma 0.3, 1, 3 | hit the 24 h wall at iterations 81, 91, 95; checkpoints at 80 |
| 1623790, 1624269 | 1 km | log, thermal | link 1 out of memory at iteration 41; link 2 running from the iteration 40 checkpoint |

Local 32 km probes live under the session scratchpad; the analysis scripts
that produced the numbers above (`phi_stats.py`, `misfit_split.py`,
`plot_conv.py`) are there as well and are small enough to recreate from
this description.
