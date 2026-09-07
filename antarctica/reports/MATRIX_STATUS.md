# ISMIP7 Antarctica core-experiment matrix - status

> **Cores 1-8 are SUPERSEDED - their results are INVALID and are being
> re-run.** They were produced before the annual-mean atmosphere-forcing fix
> (`fa7230c`): the ISMIP7 SDBN1 reader collapsed each year's 12 monthly slices
> with `isel(time=0)`, applying JANUARY (peak austral summer, the
> maximum-ablation month) as the whole year's forcing for every aSMB/ts field.
> Ablation is overstated by a factor that grows with warming (MRI ssp585
> continental aSMB integral: 2050 -229 vs +562 Gt/yr, sign flipped; 2108
> -16583 vs -1276; 2300 -82606 vs -15174). Everything below for cores 1-8 -
> SMB, mass, VAF, sea-level contribution, audit verdicts, envelope
> comparisons, and where the runs stopped - reflects that bug. Cores 9-11
> (CTRL2015 on the RACMO climatology, OCX) do not read the aSMB path and are
> unaffected by *that* bug, but the CTRL differencing they support cannot be
> redone until the projections are re-run.
>
> **All 11 cores additionally predate the two ice-front fixes (Aug 2026), so
> cores 9-11 are superseded as well.** Every core ran with the 2500 m
> boundary-id sidecar on the 32 km mesh, which tagged ~5% of the ice front, so
> the calving-terminus back-pressure was missing over the other ~95%; and with
> CG1 geometry, whose lumped-mass lift biased the front thickness ~2x high and
> the terminus force ~3.3x. Both defects enter the MAP as well as the forward -
> the inversion absorbs a wrong front treatment into `θ`/`φ`, where the t=0
> misfit cannot reveal it - so the matrix needs **re-inversion and re-running**,
> not just re-running. Only the 32 km DG0 Budd MAP has been rebuilt so far. See
> `../../GEOMETRY_DISCRETIZATION.md`.

**Resolution:** 32 km (`lc=32000`), the validated demonstration resolution.
**Configuration (all cores):** Budd N_hat friction (exact-zero shelf), balanced
apparent-MB init (`ISMIP7_APPARENT_MB=1`), fixed calving front, `dt=0.1`,
transport-first stepping with the dt-subcycle rescue ladder
(`ISMIP7_SUBCYCLES=1,4,16,64`). Budget residual closes to ~0 on every run.
**Flow exponent n = 4** on the untagged legacy MAPs: every core here ran before
this branch flipped the `ISMIP7_N_FLOW` default to 3, so the per-core run-env
blocks record no `ISMIP7_N_FLOW` and the same env today would mean n = 3 (see
`../N3_FRAMEWORK.md`).

Results h5/CSV are gitignored; each core's `coreNN_<name>_32km.md` is the tracked
record (run env, budget at marker years, observational audit, ensemble overlay).

## Matrix

| Core | Experiment | Window | Reached | Status |
|---|---|---|---|---|
| 1 | historical CESM2-WACCM | 1850-2014 | 2014.0 | **superseded** (ran full window) |
| 2 | historical MRI-ESM2-0 | 1850-2014 | 2014.0 | **superseded** (ran full window) |
| 3 | ssp370 CESM2-WACCM | 2015-2100 | 2100.0 | **superseded** (ran full window) |
| 4 | ssp370 MRI-ESM2-0 | 2015-2100 | 2100.0 | **superseded** (ran full window) |
| 5 | ssp126 CESM2-WACCM | 2015-2300 | 2300.0 | **superseded** (ran full window) |
| 6 | ssp126 MRI-ESM2-0 | 2015-2300 | 2300.0 | **superseded** (ran full window) |
| 7 | ssp585 CESM2-WACCM | 2015-2300 | 2124.5 | **superseded** (partial, saturation) |
| 8 | ssp585 MRI-ESM2-0 | 2015-2300 | 2300.0 | **superseded** (ran full window) |
| 9 | CTRL2015 CESM2-WACCM | 2015-2300 | 2300.0 | **superseded** (ran full window; the control) |
| 10 | CTRL2015 MRI-ESM2-0 | 2015-2300 | 2040.5 | **superseded** (partial, edge; twin of core 9) |
| 11 | OCX obs-constrained | 1990-2025 | 2025.0 | **superseded** (ran full window) |

**Every core is invalidated by the defects in the banner above (cores 1-8 by
the January-forcing bug as well) and the matrix is being re-inverted and
re-run; the "Reached" column records only how far the superseded run got.**
Of the runs as they stood, 9 of 11 covered their full windows. The two that
did not:

- **Core 7 (ssp585 CESM):** 108+ years - the exact configuration that died at
  *year 9* before this session's initialization work. Stopped at 2124.5 in
  what the buggy forcing made a deep-collapse era (net SMB -6659 Gt/yr, melt
  ~13,800 Gt/yr) where front-emptying events become continuous and the
  operator-split diagnostic saturates. With ablation overstated several-fold,
  that wall is very likely forcing-induced rather than the structural boundary
  of the split scheme it was read as; the re-run settles it.
- **Core 10 (CTRL MRI):** stops at 2040.5 at a front-emptying event a warm
  resume also can't cross. Not affected by the forcing bug (RACMO baseline).
  Its physical twin **core 9** (same RACMO baseline; the ESM enters CTRL only
  via the acabf fallback) reached 2300, so CTRL behavior *is* demonstrated -
  core 10's tail is a split-scheme-edge artifact, not missing physics.

Core 10 points at the known fix for continuous front-emptying eras: the
monolithic implicit `(u, M, tau, h)` coupling (gia `forward_monolithic` port).
Whether core 7 still needs it is an open question for the re-run.

## What this run demonstrates

The forcing-magnitude claims below do NOT survive the January-forcing bug; the
numerical/conservation claims do, since they are independent of the forcing
values fed in.

- The **initialization + transport machinery is sound end to end**: every core
  starts from the balanced control, applies the full ISMIP7 protocol forcing
  path (re-referenced aSMB + per-year Burgard ocean melt, both ESMs, all
  scenarios), and conserves mass to machine precision through
  century-to-multi-century integrations. The *reader* on that path was wrong
  (January, not the annual mean); the budget closure is not evidence it was
  right.
- ~~**Inter-ESM and inter-scenario structure is physical**~~: withdrawn.
  The historical VAF contrast and the ssp585-vs-ssp126 spread were computed
  from January-only ablation, which biases both ESMs and grows with the
  scenario's warming, so the structure has to be re-established from the
  re-run.
- The **rescue ladder + dt-subcycling** carries the split scheme across
  discrete front-emptying events (dozens per run) that previously stopped it.
  How much of the window that buys under correct forcing is a re-run
  question - many of those events were driven by the spurious ablation.

## Open items (carried in the per-core reports)

1. **~9% aSMB unit inflation** (`forcing.py:smb_kgm2s_to_myr`) - a real
   pending decision; inflates every ESM SMB anomaly. (The matrix sitting
   above the ISMIP6 envelope is a real forced-response bias - a
   melt/dynamics parametrization target - not a differencing artifact.
   `control/run.py` branches the CTRL from the historical endpoint because
   that is the protocol-correct same-state differencing (the control and
   projections must share the historical endpoint so their common
   relaxation drift cancels in proj-CTRL), but an isolation test at fixed
   n=4 shows this correct control makes proj-CTRL LARGER, not smaller
   (ssp126 +161 vs +132 mm): a_ref is a t=0 balancing correction that
   zeroes the initial tendency, not a net sink, so the earlier "spurious
   ~160 mm a_ref trend" attribution was wrong. n=3 rheology roughly halves
   the overshoot (ssp126 +85 mm, still above the [-14,+50] mm envelope).
   Note the archived n=4 cores 9/10 CTRL records predate the hist-branch
   fix: they are cold-start controls branched from the 2015 inversion, a
   different initial state than the projections, so proj-CTRL differences
   computed against them do not cleanly isolate the forced response.)
   **Every number in this item, including the isolation test and the
   above-envelope overshoot it explains, was produced with the January-only
   forcing.** That bug overstates late-century ablation several-fold and is on
   its own a sufficient cause of the overshoot, so the magnitude and the
   attribution both have to be re-derived from the re-run before the ~9% unit
   question can be judged against them.
2. **2014→2015 projection handoff** starts projections at the historical
   final's 2014.0 with a one-year zero-anomaly gap.
3. **Runaway-detector peak clause** flags isolated one-step discharge spikes
   during emptying events as FAIL even though the budget closes; the audit
   verdict is otherwise ON TRACK. Worth refining to sustained-growth only.
4. **Monolithic forward** for cores 7 (and the 10 tail) beyond saturation.
5. **500 m / 2500 m production resolution**: the 2500 m `_budd` MAP on disk
   (`inversion_icepack2_budd_2500.h5`) is the untagged n=4, CG1-geometry one;
   this matrix is the 32 km demonstration. An n=3 production line needs its own
   `inversion_icepack2_budd_n3_dg0_2500.h5` (see `../N3_FRAMEWORK.md` for the
   naming rule), inverted with a 2500 m boundary-id sidecar.

For the pipeline and its knobs see `antarctica/README.md`; for the rheology see
`COMPOSITE_RHEOLOGY.md` and `antarctica/N3_FRAMEWORK.md`. (The deeper
solver/data history lives in the local, untracked `CLAUDE.md` working notes -
gitignored, so it is not part of this repo.)
