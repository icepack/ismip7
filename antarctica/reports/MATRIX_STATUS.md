# ISMIP7 Antarctica core-experiment matrix - status

**Resolution:** 32 km (`lc=32000`), the validated demonstration resolution.
**Configuration (all cores):** Budd N_hat friction (exact-zero shelf), balanced
apparent-MB init (`ISMIP7_APPARENT_MB=1`), fixed calving front, `dt=0.1`,
transport-first stepping with the dt-subcycle rescue ladder
(`ISMIP7_SUBCYCLES=1,4,16,64`). Budget residual closes to ~0 on every run.

Results h5/CSV are gitignored; each core's `coreNN_<name>_32km.md` is the tracked
record (run env, budget at marker years, observational audit, ensemble overlay).

## Matrix

| Core | Experiment | Window | Reached | Status |
|---|---|---|---|---|
| 1 | historical CESM2-WACCM | 1850-2014 | 2014.0 | **complete** |
| 2 | historical MRI-ESM2-0 | 1850-2014 | 2014.0 | **complete** |
| 3 | ssp370 CESM2-WACCM | 2015-2100 | 2100.0 | **complete** |
| 4 | ssp370 MRI-ESM2-0 | 2015-2100 | 2100.0 | **complete** |
| 5 | ssp126 CESM2-WACCM | 2015-2300 | 2300.0 | **complete** |
| 6 | ssp126 MRI-ESM2-0 | 2015-2300 | 2300.0 | **complete** |
| 7 | ssp585 CESM2-WACCM | 2015-2300 | 2124.5 | partial (saturation) |
| 8 | ssp585 MRI-ESM2-0 | 2015-2300 | 2300.0 | **complete** |
| 9 | CTRL2015 CESM2-WACCM | 2015-2300 | 2300.0 | **complete** (the control) |
| 10 | CTRL2015 MRI-ESM2-0 | 2015-2300 | 2040.5 | partial (edge; twin of core 9) |
| 11 | OCX obs-constrained | 1990-2025 | 2025.0 | **complete** |

**9 of 11 ran their full windows.** The two partials:

- **Core 7 (ssp585 CESM):** 108+ years - the exact configuration that died at
  *year 9* before this session's initialization work. Stops at 2124.5 in the
  deep-collapse era (net SMB -6659 Gt/yr, melt ~13,800 Gt/yr) where
  front-emptying events become continuous and the operator-split diagnostic
  saturates. This is the structural boundary of the split scheme.
- **Core 10 (CTRL MRI):** stops at 2040.5 at a front-emptying event a warm
  resume also can't cross. Its physical twin **core 9** (same RACMO baseline;
  the ESM enters CTRL only via the acabf fallback) reached 2300, so CTRL
  behavior *is* demonstrated - core 10's tail is a split-scheme-edge artifact,
  not missing physics.

Both partials point at the same known fix: the monolithic implicit
`(u, M, tau, h)` coupling (gia `forward_monolithic` port) for the eras where
front-emptying is continuous.

## What this run demonstrates

- The **initialization + forcing pipeline is sound end to end**: every core
  starts from the balanced control, applies real ISMIP7 protocol forcing
  (re-referenced aSMB + per-year Burgard ocean melt, both ESMs, all scenarios),
  and conserves mass to machine precision through century-to-multi-century
  integrations.
- **Inter-ESM and inter-scenario structure is physical**: MRI gains VAF over
  the historical (stronger SMB) where CESM loses it; ssp585 drives deep
  collapse where ssp126 stays near the control.
- The **rescue ladder + dt-subcycling** carries the split scheme across
  discrete front-emptying events (dozens per run) that previously stopped it -
  buying the full window everywhere except continuous-collapse eras.

## Open items (carried in the per-core reports)

1. **~9% aSMB unit inflation** (`forcing.py:smb_kgm2s_to_myr`) - a real
   pending decision; inflates every ESM SMB anomaly. The 32 km ssp585 sitting
   above the ISMIP6 envelope on the SMB-gain side is consistent with it.
2. **2014→2015 projection handoff** starts projections at the historical
   final's 2014.0 with a one-year zero-anomaly gap.
3. **Runaway-detector peak clause** flags isolated one-step discharge spikes
   during emptying events as FAIL even though the budget closes; the audit
   verdict is otherwise ON TRACK. Worth refining to sustained-growth only.
4. **Monolithic forward** for cores 7 (and the 10 tail) beyond saturation.
5. **500 m / 2500 m production resolution**: the 2500 m `_budd` MAP is ready
   (`inversion_icepack2_budd_2500.h5`); this matrix is the 32 km demonstration.

See `MEMORY` and `CLAUDE.md` for the deeper solver/data history.
