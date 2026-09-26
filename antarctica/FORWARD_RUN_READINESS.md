# Forward-run readiness for the ISMIP7 submission

Swept from the ISMIP discussion board (https://github.com/orgs/ismip/discussions,
45 threads) on 13 September 2026, and again on 19 September (48 threads, plus
the two in `ismip7-antarctic-ocean-forcing`). This file is the tracked home of
both sweeps. Sections 1 to 5 are the first, corrected where the second found
them out of date; section 6 is the second, with the thread-by-thread
checklist. Section 7 is a third, read-only pass on 21 September (49 threads)
over the threads that had moved; it updates that checklist in place. Ordered by
what blocks a submission first.

## 1. Protocol changes since 12 August

**Data freeze, 3 September (#37).** The forcing is frozen unless something is
badly wrong. The freeze copy is on Source Cooperative; Globus stays the archive
of record. Only the current version of each product is kept, and version
numbers are independent per product, so an atmosphere v2 with an ocean v3 is
fine. Use the latest of each and write the versions into the README. AIS
fracture for CESM2-WACCM ssp585 may still change, having been built from the
replaced SDBN1 v2.

**Source Cooperative mirror (#40).** Products `ismip7-ais-forcing`,
`ismip7-ais-observations`, `ismip7-ais-melt-calibration`. Re-synced with Globus
on 11 and 21 September, by hand every week or two. The forcing layout is
`data/<ESM>/<scenario>/<product>/<variable>/<file>`, with no `AIS/` level and
no version directories; the version sits in the filename. Anonymous HTTPS
listing needs a browser-like User-Agent. This is the route for NOTS, which has
no Globus endpoint; `antarctica/scripts/download_mirror.py` takes it with no
AWS tooling.

**Renamed MRI atmosphere (#37).** `SDBN1-*` became `GEMB-SDBN1-*` for
MRI-ESM2-0, a name change with the data unchanged. The core experiment uses
`SDBN1` for CESM2-WACCM and `GEMB-SDBN1` for MRI-ESM2-0; `dEBM2` is the second
downscaling for perturbed ensembles. Handled: `atmosphere_product` accepts
whichever exists.

**Control experiment (#28, #15).** ctrlclim is the 2000-2029 climatology of the
last 15 years of `historical` and the first 15 of `ssp126`, per ESM. Files are
time varying so the setup matches a projection. Fracture and excess meltwater
hold at 2015 conditions; SMB-height feedback, calving and GIA stay free.

**Historical start is free (#34).** A steady initial state may be assigned to
1960 or 1975 rather than 1850, with the 1960-1989 anomaly reference driving
changes from there. This closes the 1 K cooling question at an 1850 start.

**OCX (#32, #33, #45, #48).** Independent of CMIP by design, spin-up method is
the group's choice, no fracture forcing for historical or OCX (use the observed
front positions in `obs/`), ice mask ends 2021. The spatially shifted AIS OCX
`dacabfdz` of #45 was replaced IN PLACE around 8 September, confirmed on the
17th to be OCX only; the OCX gradients are v2 on the mirror. Nothing here reads
them (no SMB-height feedback), but it is the case that showed a re-sync could
not see a replaced file, see section 6. Open since 17 September: the OCX `main`
thermal forcing departs from the Zhou climatology around Mertz (#48). Confirmed
upstream on 20 September: OCX was built on an earlier extrapolation and is
being regenerated, with no date set, section 7.

**Fracture masks (#29, #30, #33).** Both ESMs' SSPs (CESM ssp585 at v2.1), none
for historical or OCX. Floating ice only. The masks light up near the grounding
line on Ross and FRIS, and the focus group suggests pairing them with a stress
criterion (Lai et al. 2020). Available as `ISMIP7_FRACTURE=mask` and, since the
second sweep, `mask_front`; the default and the core matrix are `none`. Which
the submission uses is a decision, section 6.

**NaN forcing outside the downscaled mask (#39)** is intended. Zero, nearest or
a large melt are all acceptable, stated in the README. Ours fills with zero.

**CESM2-WACCM ends in 2299 (#8, #49).** The empty 2300 atmosphere files were
removed, and by 22 September 2300 was back as the 2290-2299 mean, which the
reader opens as given. The ocean stops at 2299, while a 2015-2300 run needs
2300. Handled: the reader persists the last year on disk exactly one year past
the end and reports it once per variable; a gap inside the series stays an
error. The collapse masks end in 2299 too, and 2300 reads the 2299 mask
(section 6). On 23 September the organisers answered that one forcing year at
the end makes no significant difference and that duplicating 2299 into 2300 is
acceptable, so runs end in 2300 (section 8).

**Melt toolbox re-release (#25).** Rerun the calibration notebook with the new
constraint datasets. No ice-model rerun.

## 2. Local forcing, audited against the mirror (13 September)

`antarctica/scripts/audit_forcing_versions.py` compares every product the
mirror publishes for the two core ESMs against `ISMIP7/AIS` (313 GB on disk).

| what | mirror | local | verdict |
|---|---|---|---|
| CESM2-WACCM atmosphere `SDBN1-8000m` acabf, acabf-anomaly (hist, ssp126/370/585) | v2 | v2 | current |
| CESM2-WACCM ocean so, tf | v3 | v3 (ssp585 also v1 leftovers) | current, delete the leftovers |
| CESM2-WACCM fracture ssp126/370/585 | v2.1 | v2 | behind, re-download |
| MRI-ESM2-0 atmosphere | `GEMB-SDBN1-8000m` v1 | `SDBN1-8000m` v1 | data current, directory renamed |
| MRI-ESM2-0 ocean so, tf | v3 | v3 | current |
| MRI-ESM2-0 fracture | v1 | v1 | current |
| `ctrl` trees, both ESMs (atmosphere and ocean, 2015-2300) | v2, v3 | absent | download for cores 9 and 10 |
| OCX atmosphere `OCX/RACMO2.3p2-ERA/SDBN1-8000m` (acabf, tas, ts, pr, mrro, gradients, 1979-2025) and `OCX/ocean/{main,cold,warm,vary}` (1950-2025) | v1 | absent | download for core 11; `main` is the core variant |
| ocean `extras` (ct/sa bias and climatology), `grid/ocean/ISMIP7/8km-60m` | v3 | present | fine |
| `tas`, `ts`, `pr`, `mrro`, gradients for scenarios other than ssp585; 2 km atmosphere; `dEBM2` | v2, v1 | mostly absent | needed only for SMB-height feedback or a dEBM2 member |

318 of the mirror's 358 entries are absent locally and six fracture rows are
behind. Everything the current forwards read is at the freeze version.

**Staging since that audit** (the table stays the point-in-time record): the
CESM2-WACCM fracture v2.1 files and the OCX set were downloaded on 13
September, and the ssp585 CESM2-WACCM set was staged on NOTS on 14 September.
On 14 September the two `ctrl` trees came down as well, 8 km atmosphere and
ocean for both ESMs, 77 GB, every transfer size-checked. A re-run of the audit
that afternoon reads 440 mirror entries, 64 of them present locally at the
mirror's own version and none behind. Every `ctrl` row for
`SDBN1-8000m`, `GEMB-SDBN1-8000m` and `ocean` (so, tf, thetao) is current, so
cores 9 and 10 have their forcing. The 376 absent entries are the 2 km
atmospheres, `dEBM2`, and the per-scenario fields that only an SMB-height
feedback or a perturbed member reads.

## 3. Output and submission

**Status.** `ISMIP7_OUTPUT=1` makes the forward accumulate the yearly flux
means and snapshot the state each year (`icepack2_tools/ismip7_output.py`, one
Firedrake checkpoint per year, written atomically).
`antarctica/scripts/write_ismip7_output.py` regrids conservatively to the 8 km
grid through a cached supermesh overlap operator, applies the request's fill
policies and units, encodes time, and writes the 21 gridded and 10 scalar files
under `AIS/<source_id>/<ism_id>/CORE/<exp>/`. At full length, 2015 to 2300,
a 32 km CESM2-WACCM control and ssp585 pass `ismip7-compliance-checker` 0.5.1
with zero errors in every test group, the experiment-length checks included
(run records `core09-32km-ctrl2015-cesm2waccm-p4` and
`core07-32km-ssp585-cesm2waccm-p4`, 23 September 2026). Two output rules made
that possible, both found by the first full-length check (action 10).

Conventions chosen: the fluxes are what the transport applied, after the
positivity limiter (action 10), `acabf` is that SMB with the apparent-MB
correction travelling separately as `acabf_correction`, `ligroundf` is booked
into the first floating cell and signed positive from grounded to floating
(section 9), a floating cell within 1 cm of the bed is written as grounded,
`lithk` is zero where the ice mask is zero, and `base = orog - lithk` on the
grid.

What a submission needs (#5, #16, #17, #18, #19, #20, #22, #23):

- **Variable request:** `isschecker/data/ISMIP7_variable_request.csv` in
  `ismip/ISM_SimulationChecker`; the mandatory set is ISMIP6's. Fluxes
  (`acabf`, `libmassbf`, `licalvf`, `ligroundf`) need conservative regridding
  from the DG0 cells.
- **Grids:** the description files in `ISM_SimulationChecker/gdfs` define x and
  y for every supported resolution and feed CDO directly.
- **Time encoding:** `days since 1850-01-01`. State variables are snapshots at
  1 January of the following year (2015 stamped 2016-01-01), fluxes are yearly
  means stamped 1 July, the initial state is not requested, and the filename
  year range is the years actually run. A run's `t_end` is 1 January of the
  year after its last, so projections run to 2301 and the historical to 2015,
  putting the handoff checkpoint at 2015.0.
- **Names:** `<var>_AIS_<source_id>_<ism_id>_m001_<ESM>_f001_<scenario>_C0NN_<years>.nc`
  under `AIS/<source_id>/<ism_id>/<set_id>/<set_counter>/`, `ism_id` without
  dots or underscores.
- **Checker:** bounds were relaxed in July (topg to 5500 m). `licalvf` negative
  means loss. `ligroundf` is a specific mass flux booked into the last grounded
  or first floating cell, bounded to [-10, 10] kg m-2 s-1 at 0.5.1. A wholly
  floating pixel must sit more than 1 cm above `topg`. `topg` and `lithk` must not
  be masked to the evolving ice sheet, or the scalar tool's sea-level numbers
  break.
- **Scalars:** the model integrates the ten scalars itself on the native mesh.
  `ismip/ismip7-scalar-processing` also produces them, and the sea-level
  estimates (`sla20`, `slg20`, `slvaf`) which nothing here computes, so it
  still has to be run on the gridded files; 2D fields and scalars share the
  experiment folder. `scripts/compare_scalars.py` sets its scalars against
  ours and names every part of the difference.
- **README:** the template is the Google document from discussion #6, recording
  forcing versions, the NaN rule, the historical start year and the spin-up.
  Drafted in `ISMIP7_README_AIS_RICE_icepack2.md`, which owns those answers;
  open items are marked `[confirm]`.
- **Upload:** email ismip6 at gmail.com with the Globus id, `AIS`, group name
  and `ism_id` to receive a folder. (issue #39)
- **Experiment ids:** C001 CESM2-WACCM historical, C002 MRI-ESM2-0 historical,
  C003/C004 ssp370, C005/C006 ssp126, C007/C008 ssp585 (CESM first), C009/C010
  ctrl, C011 OCX. Output goes on the standard ISMIP7 grid closest to the native
  grid, 8 km here. 3D fields are requested at a few times only and we have none.
- **Deadline:** the site still describes the 30 June round; the board refers to
  the end of September 2026 for the next one, scope to confirm. (issue #40)

## 4. Model-side state

**Inversions.** RC and Budd MAPs exist on the adaptive-preset mesh
(`inversion_icepack2_{rc,budd}_n3_dg0_logvelnet_ua2000.h5`). Every Budd MAP
older than 13 September carries the shelf-friction defect and is unusable.
The 2 km RC and Budd MAPs inverting at Rice replace them (issue #24); their
22 September snapshots go through `make map-check`, native and transferred
onto the 1 km / 10 km mesh, and the numbers land in `MAP_CHECK.md`.

The 13 September Budd MAP was inverted while the shelf gate still multiplied
through by the grounded indicator `He`, and the shipped gate is height above
flotation alone, so it was re-inverted under the shipped law as NOTS 1390416
(200 iterations, final masked misfit 1.041e4, 14 September). The census
justifies that on its own: on the adaptive mesh the old sign gate puts 13 647 of
103 233 floating cells at the friction cap, the `He` form 8 781, and the
shipped gate 0. The superseded file is kept as
`inversion_icepack2_budd_n3_dg0_logvelnet_ua2000_hegate.h5`. The RC MAP never
carried the gate.

**Resolved 21 September: the MAPs reproduce their own velocity once the
forward assembles the residual they were inverted under.**
`check_budd_map.py --forward` re-solves the diagnostic at a MAP's controls and
compares against the velocity the MAP saved; `ISMIP7_CHECK_FRICTION` (default
`budd`) selects the law. Measured at Rice, serially, under the `ismip7-pr6`
clone:

| MAP | law | forward stabilizers | relative L2 | solved mean speed | saved mean speed |
|---|---|---|---|---|---|
| `inversion_icepack2_budd_n3_dg0_logvelnet_ua2000.h5`, 14 September | Budd | as shipped | 0.665 | 56.3 m/yr | 105.3 m/yr |
| same (NOTS 1569253) | Budd | `ISMIP7_OCEAN_DRAG=0 ISMIP7_U_LIM=0` | 1.858e-8 | 105.3 m/yr | 105.3 m/yr |
| `inversion_icepack2_rc_n3_dg0_logvelnet_ua2000.h5`, 14 September (NOTS 1568627) | regularized Coulomb | as shipped | 0.6854 | 68.4 m/yr | 137.1 m/yr |
| same (NOTS 1569253) | regularized Coulomb | `ISMIP7_OCEAN_DRAG=0 ISMIP7_U_LIM=0` | 9.049e-8 | 137.1 m/yr | 137.1 m/yr |
| `inversion_icepack2_budd_n3_dg0_logvelnet_ua2000_pr6.h5`, 18 September (NOTS 1563053) | Budd | as shipped | 1.104e-7 | 84.3 m/yr | 84.3 m/yr |
| fresh 32 km Budd, 10 iterations, current code | Budd | as shipped | 1.04e-7 | 131.76 m/yr | 131.76 m/yr |

The RC and `_pr6` rows are the same on both mesh routes, the `.msh` named and
the mesh taken from the checkpoint. The two 14 September MAPs came from
`ismip7-next`, whose `inversion_icepack2.py` passes no `residual_stabilizers`
while its forward applies `ocean_drag` and `u_lim`. The commit that shares the
stabilizers with the inversion (571d1c9) was authored at 19:04 CDT on 14
September, after the Budd re-inversion NOTS 1390416 started at 15:29. Inversion
and forward therefore assembled different residuals, and the 15 September mean
speeds were real: the drags alone account for them.

A forward from a MAP inverted before 571d1c9 starts from the inverted state
only with `ISMIP7_OCEAN_DRAG=0 ISMIP7_U_LIM=0`. MAPs inverted after it, the
`_pr6` Budd and the 2 km RC and Budd now inverting under the prior metric,
reproduce as shipped. The 14 September RC MAP needs re-inverting only to match
the default forward. The 10-year run of job 1368723 used the default
stabilizers, so it stays a pipeline exercise.

A separate defect affected the 32 km measurement alone. When the forward
builds its mesh from the `.msh`, the checkpoint carries its own copy numbered
differently, the direct load of the saved field fails, and the fallback copied
raw `dat` arrays between the two, so it compared permuted fields and reported
0.570 for the fresh 32 km MAP. The two mean speeds were equal to every printed
digit, and the nodal speeds differed by up to 7,023 m/yr in the given order and
by 0.012 m/yr once both were sorted, with identical sums. The comparison now
aligns nodes by coordinate and refuses two meshes that are not the same vertex
set (`check_budd_map.node_permutation`, `tests/test_forward_check_ordering.py`).
The adaptive-mesh measurements above used the checkpoint mesh, so the defect did not touch
them. On the 32 km MAP `probe_forward_consistency.py` also found `N_ref=None`
(the inversion's own call, via `ISMIP7_BUDD_NREF=none`) and the inversion's
composite alpha inert to every digit.

The check needs a MAP's final save: the every-20-iterate checkpoints carry the
controls but not the velocity.

**Forward.** The RC control on the adaptive mesh runs and holds (1 yr, resid 0). A
10-year CESM2-WACCM ssp585 on that mesh (NOTS job 1368723) took 10.5 minutes on
32 Sapphire Rapids ranks, 6 s per 0.1-year step, so a 2015-2300 projection is
about 5 node-hours and eleven cores about 2.5 node-days.

**Melt calibration, rerun 14 September.** The re-released calibration product
combines Paolo (2023), Davison (2023) and Adusumilli (2020), and its integrated
target is 1067.4 Gt/yr against the 865.0 Gt/yr of the Paolo plus Adusumilli
table the old calibration used. Both tables went through `calibrate_melt.py` on
the same adaptive mesh, so the comparison isolates the observations:

| observations | integrated target | K* | melt at K* |
|---|---|---|---|
| Paolo, Adusumilli | 865.0 Gt/yr | 4.347e-5 | 677.0 Gt/yr |
| Paolo, Davison, Adusumilli | 1067.4 Gt/yr | 4.700e-5 | 732.0 Gt/yr |

As summary statistics of the fit, K* rises 8% and the total-match K rises with
the target by 23%, 5.553e-5 to 6.853e-5. The old table's K* on this mesh sits
within 2% of the 2500 m result that predates it, so the mesh is not what moved.
`calibrate_melt.py` takes the newer table by default, `ISMIP7_MELT_OBS_CSV`
names either, and the saved npz records which one produced it.

A forward reads the per-basin `K_basin`, each basin fitted to its own
observation, so the melt it applies moves basin by basin. The ratio of new to
old K_b spans 0.86 to 3.58 with a median near 1.15:

| basin | old K_b | new K_b |
|---|---|---|
| 1, Antarctic Peninsula fringe | 2.100e-5 | 7.515e-5 |
| 7 | 3.724e-5 | 6.314e-5 |
| 9, Amundsen | 1.464e-4 | 1.935e-4 |
| 13, the one that falls | 7.922e-5 | 6.850e-5 |

**Still open:** nothing from the August list. Delivered here: the `ssp126`
ctrlclim default (`icepack2_tools/climatology.py`), the collapse mask in the
thickness update, the 2300 forcing year, the `GEMB-SDBN1` path, the
forcing-version audit, the output writer, and the melt calibration above.

## 5. Actions, in order

1. Drive the full-length ssp585 (NOTS 1390452, 2015 to 2301 on the adaptive mesh with
   `ISMIP7_OUTPUT=1`) through the writer and the compliance checker. It is the
   first run at experiment length, so it is what clears the checker's remaining
   length checks. Record which K calibration it read: a run picks up whichever
   `calibrated_K_per_basin_*.npz` is staged when it starts. The job is held in
   the queue until the slope convention shared by calibration and forward is
   settled, taken up in action 5, since it decides which K the run should read.
   It no longer waits on MAP self-consistency (section 4, 21 September);
   its MAP must pass action 3 in the configuration it runs. Once the
   convention is settled,
   `scontrol release 1390452` starts it. (issue #27)
2. Settle the `[confirm]` items in the submission README draft with the group. (issue #38)
3. Run the corrected `check_budd_map.py --forward` on every MAP the matrix
   will use, on its final save, and record the number beside the MAP. A MAP
   inverted before 571d1c9 must run with `ISMIP7_OCEAN_DRAG=0 ISMIP7_U_LIM=0`
   or be re-inverted (section 4). The 2 km RC and Budd MAPs now inverting under
   the prior metric are next when they finish. (issue #24)
4. Re-run `audit_forcing_versions.py` immediately before the production matrix
   and cite it in the README. The `ctrl` pull for cores 9 and 10 is done, and
   the mirror is re-synced with Globus by hand every week or two, so the freeze
   versions can still move under a long campaign. (issue #41)
5. Settle where the `libmassbffl` bound violation comes from. The request's AIS
   minimum is -0.008 kg m-2 s-1, which is 275.3 m/yr of ice, and the 10-year adaptive-mesh
   ssp585 of job 1368723 reached -0.0117, or 402.6 m/yr.

   `check_melt_bound.py` measures the slope side of it. The CG1 calibration
   behind the existing K files caps the draft slope `sin(alpha)` at 5e-3 and
   the forward applies no cap, so the melt the forward applies is a different
   field from the melt the per-basin K was fitted against. The script evaluates two halves at the reference geometry
   with `calibrated_K_per_basin_2000.npz` on the adaptive 2 km mesh, each capped and
   uncapped. The calibration half reproduces `calibrate_melt.py` on CG1 nodes,
   with BedMachine's raster surface and mask and the cap on the nodal slope.
   The forward half reproduces the forward on DG0 cells, with the surface from
   flotation, `forcing.compute_sin_alpha`'s cell slope, forcing at each cell's
   own draft, the callback's `haf <= 0` floating test and the cap on the cell
   slope. The calibration half floats 1 512 899 km2 over 47 288 nodes and the
   forward half 1 631 466 km2 over 85 820 cells, and they differ in mask and
   quadrature as well, so their totals compare in magnitude.

   The rows below predate the seawater flotation test, which landed with
   issue #66, and the
   `h > 0` test in the forward half, and are to be re-measured with the DG0
   melt total (issue #30); the current numbers are in
   `GEOMETRY_DISCRETIZATION.md`.

   | slope | max, m/yr | p99, m/yr | area mean, m/yr | integrated, Gt/yr | past the bound |
   |---|---|---|---|---|---|
   | calibration half, capped, what K was fitted to | 71.1 | 22.2 | 0.77 | 1067 | 0 |
   | calibration half, uncapped | 1804.9 | 256.3 | 4.18 | 5803 | 421 nodes, 2920.6 km2 |
   | forward half, uncapped, as the forward runs today | 1144.1 | 59.0 | 1.16 | 1732 | 97 cells, 388.2 km2 |
   | forward half, capped | 57.8 | 13.0 | 0.43 | 646 | 0 |

   The same script with `calibrated_K_per_basin_2500.npz`, the coefficient file
   the 10-year run of job 1368723 read, everything else unchanged, gives the
   like-for-like comparison against that run:

   | slope, 2500 file | max, m/yr | area mean, m/yr | integrated, Gt/yr | past the bound |
   |---|---|---|---|---|
   | calibration half, capped, what K was fitted to | 54.0 | 0.61 | 841 | 0 |
   | calibration half, uncapped | 1522.7 | 3.34 | 4634 | 320 nodes, 2193.7 km2 |
   | forward half, uncapped, as the forward runs today | 869.8 | 0.92 | 1380 | 63 cells, 248.6 km2 |
   | forward half, capped | 43.9 | 0.34 | 510 | 0 |

   The calibration half reproduces the target its own K was fitted to: 1067
   against 1067.4 Gt/yr for the 2000 file, and 841 against 865 Gt/yr for the
   2500 file, within 3%. That is the internal consistency check. With the same
   K, the forward applies 1.62 times the target with the 2000 file, 1732
   against 1067, and 1.64 times with the 2500 file, 1380 against 841. The ratio
   is the durable result, stable across both coefficient files.

   The 10-year run booked 1860 Gt/yr with the 2500 file, above the 1380 Gt/yr
   its forward half applies at the reference state, since that run carries
   warmer ssp585 thermal forcing over evolving geometry where the script holds
   the OI climatology at the initial state. That gap is a consistent residual.
   An earlier reading of this section put the uncapped forward half within 7%
   of the run; it compared the 2000 file's 1732 Gt/yr with a run on the 2500
   file and does not hold.

   Capping the forward's own slope gives 646 Gt/yr with the 2000 file, 39%
   under the target. Neither convention on its own reconciles the halves,
   which also differ in floating area, mask and quadrature.

   An earlier form of the script lifted the forward's slope onto CG1 nodes and
   melted it with CG1 forcing and the raster mask. Its forward rows, 4293 Gt/yr
   and 298 nodes past the bound uncapped and 1028 Gt/yr capped, reproduced
   neither half and are superseded, and with them the reading that capping the
   forward lands within 4% of the target.

   Both mechanisms for the bound violation are visible at the reference state.
   Uncapped, 97 forward cells melt past the bound over 388.2 km2, so the
   parameterisation itself exceeds it there. Their median area is 3.61 km2
   against 64 km2 for an 8 km pixel, so the gridded value comes from a small
   hot cell filling its pixel under the request's `no_floating_ice` fill
   policy. The full-length 32 km runs of 23 September measured the two further
   contributions (action 10). The bookkeeping was the larger: `book_advance`
   booked the melt requested of a step while a nearly ice-free floating cell
   can only lose what it holds, and it now books what the transport applied.
   With the evolved geometry and its warmer projected thermal forcing the
   ssp585 then falls below the bound in 0.0382 % of its values, a warning.
   Since 24 September every flux is a whole-pixel mean, which weights a hot
   cell by its share of the pixel. The p4 pair regridded that way falls below
   the bound in 0.00545 % of the ssp585's values, and the control no longer
   does.

   `calibrate_melt.py` now fits K through the forward's own melt path under
   `ISMIP7_GEOMETRY_SPACE=dg0` (issue #30). The forward and the calibration
   default to the ISMIP7 reference slope, one constant `sin(alpha)` =
   5.115e-3 (`ISMIP7_MELT_SLOPE=ant`); the local slope, capped or not, stays
   as `local` and is tied to the unsettled upstream local-slope question.
   `load_K_per_basin` warns once per run when the K file it reads was fitted
   under another convention, and under `local` when it records a cap the
   forward does not apply (issue #26, `GEOMETRY_DISCRETIZATION.md`).
6. Optional: read the provided `ctrl` trees in place of the `ssp126`
   reference-climate pool. Closed as icepack/ismip7#43, not planned for
   September 2026.
7. Optional: the stress criterion (Lai et al. 2020) alongside the collapse
   mask. Closed as icepack/ismip7#44, not planned for September 2026.

## 6. Second sweep, 19 September

### What moved on the board since the 13th

- **#30, 16 to 18 September.** Groups compared how they apply the collapse
  mask. Emptying every flagged floating cell opens holes behind the Ross and
  Filchner-Ronne fronts, the holes act as open ocean, and the inflowing
  glaciers lose their buttressing: 3.5 m to nearly 5 m of sea level by 2300 in
  one model, a doubling in another. The other end-member removes flagged ice
  only at the front. Ours was the first, with no connectivity test.
- **#46, 16 and 17 September.** isschecker 0.3.0 to 0.5.0: a range finding is a
  warning below 1 % of the values and an error above, units go through UDUNITS,
  `hfgeoubed` is soft. The record of section 3 ("pass every content check")
  predates all three and names no version.
- **#48, 17 September, open.** The OCX `main` thermal forcing against the Zhou
  climatology around Mertz: 50 % less melt there than the calibration, possibly
  an older extrapolated climatology underneath OCX.
- **#22, 18 September, open.** Whether `ligroundf` should carry the sign of the
  other fluxes. Ours is positive for ice leaving the grounded sheet.
- **#45, 17 September.** Closed as OCX only. **#36, 18 September.** Either
  elevation gradient is accepted if the README says which.

### Every Antarctic thread

`[x]` handled or not applicable, `[ ]` still owed (see the actions),
`[~]` open upstream.

Forcing data:

- [x] #2 size, 8 km product: only `*-8000m` is read.
- [x] #3 Globus authentication: the mirror route needs none.
- [x] #7 `pr-anomaly` withdrawn: precipitation is never read.
- [x] #8 CESM2-WACCM ends in 2299, empty 2300 gradient files: one year is held
      past the end of a series, atmosphere and now ocean alike, logged, and
      anything further raises. A file with no time slices is named, with the
      cause.
- [x] #49 (the forum thread) CESM2-WACCM `thetao` ends in 2299 for ssp126 and
      ssp585: the #8 rule, the ocean reader holds 2299 for the single year 2300.
      Upstream will not add 2300, and on 23 September the organisers accepted
      a 2299 duplicate for 2300, section 8.
- [x] #9, #24 time stamps and calendars differ between products: the year comes
      from the filename and only the year of a time value is ever read.
- [x] #10, #39 NaN fill and NaN outside the downscaled mask: zero-filled, in the
      README.
- [x] #37 MRI `GEMB-SDBN1` rename: the reader, the audit and now the Globus
      downloader, which used to skip all of MRI-ESM2-0.
- [x] #37, #40 data freeze, version directories dropped on the mirror: handled;
      each run now logs the product and version it opened and
      `core_report.py` carries that into the committed report.
- [x] #45, #41 files replaced in place under the same name and version: the
      mirror downloader keeps ETags in `.mirror_manifest.json` and refetches a
      moved one, the audit reports `REPLACED` and `PINNED` and fails on them
      (a version newer than the mirror's, fetched from Globus first, reads
      `AHEAD` and passes, and a collapse mask below
      `forcing.FRACTURE_MIN_VERSION` reads `OUTDATED` and fails), the Globus
      route has `--resync`.
- [x] #41 item 9, dotted fracture versions: the Globus downloader no longer
      passes over `v2.1` for `v2`.
- [x] #41 items 1, 3, 4, 5, 7, 8, 10 to 14: none in a tree that is read. The
      misfiled-scenario case (item 3) is what made `available_years` strict.
- [~] #37 CESM2-WACCM ssp585 fracture may be re-cut: re-audit before production. (issue #16)
- [x] #29, #33 no fracture forcing for historical, control or OCX: refused.
- [x] #15, #28 control definition: window and pool centralised. The control
      reads its ESM's `ctrl` ocean (icepack/ismip7#107); reading the `ctrl`
      atmosphere in place of the `ssp126` pool stays optional (action 6).
- [x] #34 1960-1989 anomaly reference: the anomaly is re-referenced to the
      2000-2029 pool and there is no temperature forcing, so no jump at the
      start of a historical. The README now says the historical starts in 1850
      from the 2015 state.
- [x] #35, #36 runoff-gradient sign, which gradient: no SMB-height feedback, and
      the README names both gradients as unused.
- [x] #32, #33, #41 item 6 OCX: the readers open the real tree and core 11 runs
      on it by default (`ISMIP7_OCX_FORCING`).
- [~] #48 Mertz, and Cook by the product author's account: OCX confirmed built
      on an earlier extrapolation and being regenerated, no date set (section 7).
      `check_melt_bound.py --ocx` is the tripwire. [ ] Run it. (issue #11)
- [~] #11 Ross warm stripe: a feature of the climatology. Thermal forcing is
      used unsmoothed, in the README. A perturbed member, not a fix.
- [x] #25 melt toolbox re-release: recalibrated on 14 September, section 4.
      The `06_nov` climatology's packaging fault is fixed upstream (`so` and
      `thetao` at `v4`) and the reader now finds it; `calibrate_melt.py` still
      reads `30_sep` whatever `ISMIP7_OI_VERSION` says, so `30_sep` stays the
      default.
- [x] ocean-forcing #61, #84, issue #124: `tf` and `so` are read one chunk file
      at a time; nothing to do.

Ice-shelf collapse:

- [x] #30 `ISMIP7_FRACTURE=mask_front` is the second end-member.
      [ ] Which mode the submission uses is a decision; the README's item 9 (issue #10)
      had said the mask is applied while the core matrix runs `none`.
- [x] #41 item 2 wrong MRI ssp126 `lake_properties`: never opened.

Output and submission:

- [x] #14, #16, #20 time encoding, no initial state, filename years.
- [x] #16 `licalvf` negative for loss. [x] #22 `ligroundf` sign: settled by
      the group on 22 September with the grounded sheet as the reference,
      positive for grounded ice going afloat, section 9.
- [x] #23 bounds. [x] #46 the bundled request is 0.5.1's and records its tag;
      `audit_variable_request.py` finds drift. [x] The checker at 0.5.1 passes
      a full-length 32 km control and ssp585, action 10. At 0.5.1 the `tend*`
      totals are bounded to [-1e9, 1e9] kg s-1, and scalars are still not
      range-checked upstream.
- [x] #19 scalar-tool pitfalls. [x] The tool ran on the three 32 km p2 runs and
      the p4 pair, and its scalars are compared with ours, `reports/scalar_comparison_32km.md`.
      [ ] Run it on the submitted files and clear the README confirm. (issue #13)
- [x] #17 names: the core counter follows from the forcing and the ids are
      validated. [~] What goes in the forcing field of an OCX filename is
      unsettled: isschecker checks it against CMIP model names and has no `ocx`
      experiment row, so it rejects the organisers' own GrIS example (`ERA5`,
      `ocx`), and no AIS example exists, section 7. (issue #18)
- [x] #5, #13, #21, #6, #18, #1, #12, #38.

### Actions added, continuing section 5

8. **Decide the collapse mode for the submission** (`none`, `mask`,
   `mask_front`) and settle README item 9. A 32 km ssp585 under each of the two
   mask modes is the evidence: the log prints the mode under all three, and
   `ctx["collapse_held_cells"]` counts the flagged floating cells `mask_front`
   is holding back, written with the flagged and removed counts to the
   timeseries (`collapse_*_cells`), the budget lines and the core report.
   `mask_front` has unit tests for the rule and a 1-against-3
   rank check of the facet sweep, and has not yet run inside a forward. (issue #10)
9. **Run `check_melt_bound.py --ocx` on the production mesh** before core 11
   runs on the OCX product, and hold that run until #48 is answered if the
   Mertz block is flagged. `ISMIP7_OCX_FORCING=stopgap` reproduces the old
   core 11 meanwhile. The 20 September answer names the cause and promises
   regenerated files with no date, and places the difference at Cook, so read
   the Cook block as well as Mertz, section 7. (issue #11)
10. **Re-run isschecker** on full-length 32 km control and ssp585 outputs and
    record the version in the README. Done on 23 September at 0.5.1, on IU
    Quartz, and closed as icepack/ismip7#12. The first pass, over the three p2
    ssp585 runs, found two errors the two-year rehearsal could not show:

    | collapse mode | `libmassbffl` below -0.008 kg m-2 s-1 | `base` within 1 cm of `topg`, wholly floating |
    |---|---|---|
    | none | 1.43 % of values, an error | 13 pixel-years, an error |
    | mask | 0.341 %, a warning | 24 pixel-years |
    | mask_front | 0.35 %, a warning | 17 pixel-years |

    The first came from booking the melt requested of the transport: by 2200
    the positivity limiter held back 99 % of it. `book_advance` now takes the
    withheld sink off the SMB, melt and reference in proportion, so the fluxes
    are what the transport applied. The second is a cell 1.9 to 9.2 mm afloat,
    inside the checker's `ELEVATION_TOLERANCE`; `write_ismip7_output.py` writes
    such a cell as grounded. The p4 control and ssp585, run with both, report
    zero errors. Their warnings are the non-mandatory variables, `strbasemag`
    in 0.000813 % of the ssp585's values, and the `libmassbffl` excursion of
    action 5, now 0.00767 % of the control's values and 0.0382 % of the
    ssp585's. The booked basal melt of the ssp585 peaks near 7600 Gt/yr in
    2150 and is 3929 Gt/yr in 2300, where the requested booking reported
    123427.
11. **Run `ismip7-scalar-processing`** on the submitted files, for `sla20`,
    `slg20` and `slvaf`, compare its scalars with the native ones, and clear
    the submission README's confirm. `scripts/batch_runners/scalar_processing.script`
    runs the tool and `scripts/compare_scalars.py` as one job. Done on
    23 September on the three 32 km p2 runs and the p4 pair, each its own
    reference: every identity holds, and T - N is the tool's area factor
    (+2.2 to +2.6 % on the state scalars), the volume above flotation of
    partly grounded 8 km pixels (+39 to +44 mm of sea level by 2300 in the p2
    runs, about 10 % of the signal), the writer's fill conventions for acabf
    and libmassbffl, and melt and SMB booked where no ice takes them (111,231
    Gt/yr of melt at 2300 in p2 none, 1,986 with the applied-flux booking of
    action 10), `reports/scalar_comparison_32km.md`. Each has its issue: the
    pixel means (issue 96), the area factor (issue 97) and `params.nc` in the
    upload (issue 98), all three closed on 24 September when PR 101 merged,
    and the grid VAF (issue #99).
    The first three were fixed on 24 September on `claude/scalar-output-fixes`,
    measured on the p4 pair regridded again and a five-year p5 control. Every
    flux is a whole-pixel mean: the fill term is zero, and the control's shelf
    melt as the tool reads it is the model's to 0.4 %. The native scalars
    integrate over true area, the writer stamps the scalar files with it, and
    with `--native-af2` the area term is zero. `params.nc` sits in the upload
    beside `CORE/`, where the tool reads it. isschecker 0.5.1 still reports zero
    errors on the pair. Late in the ssp585 run the tool reads less melt than
    the model books: 1,771 Gt/yr at 2300 sits in pixels with no floating ice
    at year end, which the request fills. Nine tenths of it over the run is
    booked by cells with no ice at either end of the year, and the frozen
    apparent-MB reference supplies about two thirds of what those cells melt:
    117,000 of the 198,990 Gt left out over 2016 to 2300, and 1,150 of the
    model's 3,929 Gt/yr at 2300. Whether the production runs keep the
    reference is a group decision (issue #104); the melt it books under a
    pinned front is tracked on its own (issue #105). The rest, about 82,000 Gt
    over the run and 620 Gt/yr at 2300, is melt of real ice gone by year end,
    most of it grounded ice that goes afloat into an empty cell and melts on
    arrival. The group chose on 25 September to report that as front melt
    (issue #109): the forward books the inflow's share of the melt in marine
    cells holding no ice at either end of the year as `lifmassbf`, which the
    request never fills, and leaves the reference's share in `libmassbffl`.
    Applied to the p4 annual files, the split moves 571 Gt/yr at 2300 out of
    the pixels the fill blanks and leaves 1,200 there, 1,094 of them the
    reference's share; the submission README states the booking. A main and
    branch pair restarted from p4 at 2294.0 checked it end to end on Quartz
    (runlog `test-32km-ssp585-front-melt-main` and `-branch`): the forward's
    split matches the function applied offline to 1.5e-14 m/yr, the tool's
    `tendlifmassbf` matches the model's to 0.000 %, and isschecker 0.5.1
    finds no error outside the length checks of a five-year series.
    What is left is the submitted files, paired with their historical. (issue #13)
12. **Adopt or refetch the forcing that predates the manifest.** Done on 21
    September, and the premise above was wrong. The first
    `audit_forcing_versions.py` run counts none of the Globus-era tree as
    `older`. `plan()` reaches `OLDER` only for a file that already matches the
    mirror object's byte length, and then compares mtimes; every file on Quartz
    carries an mtime at or after its object's, so the whole tree comes back
    `adopt`. `--older refetch` and `--older adopt` therefore name the same run
    here, and neither re-downloads anything. Checked twice: the audit prints an
    `older` line only when the count is nonzero and printed none, and replaying
    `plan()` offline over all 89,767 mirror objects against a 109,931 file
    inventory of the tree gave 87,012 `adopt`, 2,755 `fetch`, 0 `OLDER`.
    The manifest is what the run is for. `download_mirror.py --older refetch`
    over `data/CESM2-WACCM/`, `data/MRI-ESM2-0/` and `data/OCX/` wrote
    `ISMIP7/AIS/.mirror_manifest.json` across those 89,767 objects, which gives
    `REPLACED` (same name, new content, #45 and #41) something to compare
    against from here on. Recording an entry and fetching are one knob: the
    script records only the keys under the prefixes it is given, and it fetches
    whatever is missing under them, so full coverage also pulled every file
    behind the 16 `MISSING` rows the audit had been reporting, 2,755 files and
    84 GB, mostly the `ctrl` ocean `thetao`, `tf` and `so` and the `ctrl`
    `mrro` and `mrro-anomaly`. The run ended 87,012 `adopt`, 2,755 fetched, 0
    failed. Eight `MISSING` rows survived it with their files on disk, because
    `local_versions` stopped at the top of a row while the mirror nests `extra`
    and `extras` one level deeper; that is fixed here, and those rows read `ok`
    against the Quartz tree.

    **The mirror is not frozen.** Between 10:22 and 10:48 UTC on 21 September,
    while the pass above was running, it gained 4,576 objects and 226 GB with
    none removed: `pr`, `pr-anomaly`, `tas` and `tas-anomaly` for `ctrl`, at
    `2000m` and `8000m`, for both core ESMs. Two listings two hours apart
    settle it, and the manifest covers the 89,767 objects that existed at the
    time of the run. The `ctrl` atmosphere was the gap this opens, and those
    4,576 objects and 226 GB are fetched, nothing failed; all sixteen `ctrl`
    `pr`, `pr-anomaly`, `tas` and `tas-anomaly` rows read `ok` over 478 mirror
    entries with 0 behind, 0 pinned, 0 replaced and 0 older. Re-list before
    trusting any earlier listing, and see issue #41 for the audit immediately
    before the production matrix.

    **The audit sees three of the mirror's ten prefixes.** Its default is
    `CESM2-WACCM`, `MRI-ESM2-0` and `OCX`, so `ACCESS-CM2`, `CanESM5`,
    `GFDL-ESM4`, `IPSL-CM6A-LR` and `MPI-ESM1-2-HR` have never appeared in an
    audit here, and all five are absent from the Quartz tree in full: 14,298
    objects and 155 GB. `grid` and `parameterisations` are unaudited and
    present. The `ismip7-ais-melt-calibration` product is 17.2 GB absent of
    25.2 GB, `meltMIP` being the part that is local. None of this reads as
    `BEHIND` or `MISSING`, because a prefix the audit never lists cannot.

    A second reason the audit cannot answer "is every file here": `MISSING` is
    per row, and a row counts as present when any version of it is on disk. The
    run that settled action 12 reported eight rows missing while thousands of
    objects were absent under rows reading `ok`. Only a `download_mirror.py`
    pass over every prefix settles file-level completeness, and
    `download_mirror.py --dry-run` does it read-only. The audit's own docstring
    now says so, and issues #41 and #16 carry the same note where their exit
    criteria lean on it. Which prefixes the submission needs, and which a
    routine re-sync covers, was icepack/ismip7#49, closed on 22 September as
    not required for the 30 September submission and to be reopened for an
    October ESM submission; the five ESMs and the melt-calibration product
    stay unfetched.

    Closed as icepack/ismip7#14.
13. **Bring the Quartz forcing tree up to the mirror.** Done on 21 September.
    The 19 September reading (listing and NetCDF headers only, nothing run)
    undercounted the rows: the audit reports seven `BEHIND`. OCX
    `dacabfdz`, `dmrrodz` and `dtsdz` each stood at `v1` against `v2` on the
    mirror, at both `SDBN1-2000m` and `SDBN1-8000m`, the `dacabfdz` being the
    spatially shifted file of #45; and the CESM2-WACCM ssp585 fracture stood at
    `v2` against `v2.1`. Nothing reads the gradients and the core matrix runs
    without fracture, so no result was affected. Those 286 files, 3.94 GB, are
    fetched, and `audit_forcing_versions.py` over its 462 mirror entries now
    reports 0 behind, 0 pinned, 0 replaced and 0 older, exiting 0. The OCX
    `acabf` (47 years, 1979-2025) and the four OCX oceans were already there in
    the layout the readers expect, so core 11 passes its coverage gate on
    Quartz. Action 4 and issue #41 still call for a fresh audit immediately
    before the production matrix. Closed as icepack/ismip7#15.

### Read off the real files on 19 September

Three assumptions nothing here could check without the forcing tree, all
settled from NetCDF headers on Quartz:

- **Collapse mask time axis.** `time` is `int32` with `units = "year"`, 1950 to
  2299, so xarray leaves it as plain years: the reader handled it before and
  handles it now, and 2300 reads the 2299 slice. The variable is `mask`
  (`standard_name = ice_shelf_collapse_mask`), beside a scalar `mapping` and
  2-D `lon`/`lat` that are coordinates only because `mask` names them. The
  reader took the first data variable, which worked by position; it now goes
  by name. The files say the shelves "should collapse on January 1st", which
  is when the first advance of a year applies the year's mask.
- **Depth axis of the ocean climatology.** `z` is height relative to the sea
  surface, positive up, -30 to -1770 m in 30 levels, in the OI climatology and
  in the Zhou re-release alike. `forcing.py` clips the draft into that range
  assuming exactly this.
- **The Zhou `06_nov` re-release.** `tf` is at `v3`, `so` and `thetao` at `v4`,
  and the `v4` files hold their own variables: the July fault (the tf field
  shipped inside the so file) is fixed upstream. The path builder said `v3` for
  all three and could not have opened the fix; it now takes each variable's
  highest version. `30_sep` stays the default, since every K is fitted to it
  and `calibrate_melt.py` reads it whatever `ISMIP7_OI_VERSION` says. Moving
  the forward to `06_nov` without recalibrating shifts the melt.

## 7. Third sweep, 21 September

Read-only, two days after section 6. The board holds 49 threads. Read in full
through the GitHub API: #17, #22, #30, #37, #40, #48 and the new #49.
`ismip7-antarctic-ocean-forcing` has no item updated since the 18th. Nothing
was posted upstream, isschecker was read at its tag and never run, and no
mirror listing was taken. Thread and issue numbers now collide: `#49` below is
the forum thread about `thetao`, and icepack/ismip7#49 in section 6 is the
board item about mirror prefixes.

### What moved on the board since the 19th

- **#48, 20 September, open.** The author of the OCX ocean product confirmed
  that OCX was built from an earlier version of the extrapolated climatology,
  v1 or v2 and unlabelled in the files, while `so` and `thetao` of the Zhou
  climatology stand at v4. The author places the difference beneath the Cook
  ice shelf near 152.5°E, in the Wilkes Land basin, about 1.5 °C colder at
  500 m in OCX, and expects zero melt or refreezing there from any
  parameterisation calibrated to the climatology. Cook is near x +1090 km,
  y -2090 km in EPSG:3031 and Mertz near x +1440 km, y -2030 km, both
  converted here from approximate positions. The OCX files are being
  regenerated from the current climatology. A maintainer of the forcing mirror
  called the update necessary; its timing is the steering committee's call and
  was undecided on the 21st. Both question whether the v4 extrapolation is the
  more realistic one at Cook, 1.5 °C above the surface freezing point under a
  shelf observed to melt at 1.3 m/yr, and no change to the climatology is
  planned. The regenerated OCX will follow the current re-release, while every
  K here is fitted to `30_sep`; whether those two differ at Cook was not
  checked. (issue #11)
- **#30, 21 September.** UFEMISM (IMAU/KNMI) posted both end-members under
  ssp585 for both core ESMs: removing flagged ice only at a front that touches
  open ocean loses about 8 % less mass by 2300 than removing it everywhere.
  The 3.5 m to nearly 5 m of section 6 is that model with the mask everywhere
  against no fracture forcing, and the doubling is PISM's everywhere against
  margin-only sensitivity runs from July. The two published end-member spreads
  are therefore about 8 % and about 100 %. The organisers replied twice the
  same day: the ISMIP6-style application can behave unrealistically on the
  large shelves, a front-connected rule is a good way to limit that, no
  front-connected mask can be supplied since every model's front differs, and
  they welcome groups choosing different approaches so that the projections
  carry the spread. Upstream prescribes no mode. (issue #10)
- **#22, 21 September, open.** The sign question of the 18th drew one reply,
  from a respondent who speaks for the organisers in #30: gain positive and
  loss negative holds in general, and for `ligroundf` the sign depends on the
  reference, since the grounded ice loses what the shelves gain. No sign is
  prescribed. The organiser the board named on the 18th as best placed to
  answer has not replied. The accepted answer is a 30 June reply, marked on
  18 August, about converting the flux to a per-area quantity. The bundled
  request table, whose `ligroundf` row is identical to isschecker 0.5.0's,
  bounds `ligroundf` to [-1e9, 1e11] kg m-2 s-1 and `licalvf` to [-1e11, 0].
  The thread opened because the checker then required `ligroundf` to be
  nonnegative, and the lower bound was relaxed to 1 % of the upper for ice
  rumples. The table therefore expects `ligroundf` positive for ice crossing
  from grounded to floating, which is how the writer books it. Settled on
  the 22nd, section 9.
- **#40, 21 September.** The mirror was re-synced with Globus, after 29 August
  and 11 September. New: `ctrl` `pr`, `tas` and their anomalies for both core
  ESMs at all resolutions, which matches the 4,576 objects action 12 saw
  arrive that morning; dEBM2 atmosphere for the five additional ESMs,
  `historical` and `ssp370`, AIS at 8 km, their SDBN1 downscaling being absent
  from Globus so far; and the melt-calibration READMEs and licence files on
  Globus. The additional ESMs are announced as incomplete in both places, with
  the ISMIP7 protocol overview sheet as the place to check what is ready.
  Nothing in the announcement touches the fracture product. (issue #41)
- **#49, 21 September, open, new.** A group reports CESM2-WACCM `thetao` for
  ssp126 and ssp585 ending in 2299 where MRI-ESM2-0 runs to 2300, and asks for
  the year on Globus. No reply yet. This is the #8 case, already in section 1:
  `ISMIP7Ocean._chunk_for` holds 2299 for the single year 2300 and says so once
  per variable. If the year is added under the same version it shows as
  `fetch` in `download_mirror.py --dry-run` while the audit table stays `ok`;
  if the version is bumped, the audit reads `BEHIND`. Answered on the 21st
  and 23rd, section 8. (issue #41)
- **Unmoved.** #17 since 29 June, #37 since 4 September. isschecker is still
  0.5.0, tagged on 17 September at the head of its default branch, with no open
  pull request and no open issue, so nothing upstream is adding an `ocx` row.

### Read this pass in #17 and in the checker source

The organisers' list of core filenames of 19 June gives core 11, for GrIS, as
`iareafl_GrIS_NORCE_CISM3_m001_ERA5_f001_ocr_C011_2015-2300.nc`, and the ISMIP7
web page carries the same list with `ocx` in place of `ocr` (read there through
an automated page summary). No AIS example exists, and the conventions PDF the
thread links was unreadable without a sign-in. At isschecker 0.5.0: `ERA5` is
absent from `VALID_ESM_NAMES`, so field 5 is an error, and
`experiments_ismip7.csv` has five rows (`historical`, `ssp370`, `ssp126`,
`ssp585`, `ctrl`), so an `ocx` file set draws a naming error and its compliance
check is skipped. The OCX atmosphere here is `RACMO2.3p2-ERA`
(`forcing.OCX_ATMOSPHERE_SOURCE`), which the list lacks as well. Core 11 cannot
pass 0.5.0 under any forcing name, the experiment row being what is missing.
(issue #18)

## 8. Fourth sweep, 22 September

Read-only, one day after section 7. The board holds 50 threads. Read in full
through the GitHub API: #22, #30, #48, #49 and the new #50. The Greenland
thread #47 was skipped. `ismip7-antarctic-ocean-forcing`,
`ismip7-scalar-processing` and the documentation repository itself have no
commit or issue since the 20th. isschecker moved to 0.5.1. Nothing was posted
upstream. On Quartz one targeted listing of the CESM2-WACCM tree was taken,
and six 2300 files were read with h5py.

### What moved since the 21st

- **isschecker 0.5.1, 22 September.** Only the variable request changed,
  from checker issue 35 (opened and closed the same morning) and #50. In
  kg m-2 s-1, `licalvf` and `lifmassbf` go from [-1e11, 0] to [-10, 0], and
  `ligroundf` from [-1e9, 1e11] to [-10, 10]. The six `tend*` totals go from
  [0, 1e25] to [-1e9, 1e9] kg s-1. The experiment table is unchanged, so
  core 11 still has no `ocx` row (issue #18). Our fluxes are annual means in
  kg m-2 s-1 on the native cells before the conservative remap. Removing
  1000 m of ice in one year is -0.029, and `ligroundf` for u = 4 km/yr,
  H = 1000 m into a 500 m cell is 0.23, so the gridded fields sit one to two
  orders of magnitude inside the new bounds. A few cells outside are a
  warning, since the grading of 0.5.0 still applies. At 0.5.1 the scalars
  are still not range-checked (`_check_numerical` calls `_check_range` only
  when the file is not a scalar). If a later release does check them, a year
  in which the mask removes a Ross-sized shelf (about 1.5e17 kg) would put
  `tendlicalvf` near -5e9 kg s-1, past the new bound. The bundled
  `icepack2_tools/ismip7_variable_request.csv` was still 0.5.0's; it was
  re-vendored at 0.5.1 on 23 September.
- **#50, 22 September, open, new.** Are the front and grounding-line fluxes
  given over the face or averaged over the cell? An organiser answered: a
  mass change per unit horizontal cell area, pointing to checker issue 35.
  That is how the writer books all three (`icepack2_tools/ismip7_output.py`:
  facet flux or removed thickness divided by the native cell area, then the
  conservative remap). Nothing to change.
- **#49, 21 September, after section 7 was written.** The ocean-forcing
  maintainer does not plan to add 2300. The CMIP archives of CESM2-WACCM
  ssp126 and ssp585 end in December 2299, and a forcing author confirmed
  that from the raw files. MRI-ESM2-0 runs to December 2300. The same
  maintainer reported that the atmosphere group padded 2300 with the
  2290-2299 mean, and would not personally endorse that. The reporting group
  leans to ending its runs in 2299 and has asked the organisers whether that
  meets the protocol. No organiser has replied. Checked on Quartz: the
  CESM2-WACCM SDBN1-8000m v2 directories hold 286 years including 2300.
  `acabf`, `acabf-anomaly` and `tas` for 2300 match the 2290-2299 mean to
  2e-7 relative for both scenarios, against 3e-2 to 5e-1 for 2299 alone. The
  ocean v3 directories end in `2291-2299`. So a CESM2-WACCM run to 2300 here
  reads the padded decadal mean for SMB and holds 2299 for the ocean. Section
  1's "the 2300 atmosphere files were removed" and the reader's comment in
  `forcing.py` are stale, and the SMB sentence of the submission README
  claimed the 2299 hold for both. That sentence is corrected in this
  change. Answered on the 23rd, below. (issue #41)
- **#49, 23 September, answered.** Read that day through the GitHub API. At
  00:03 UTC a member of the ISMIP7 Antarctica team at Dartmouth College gave
  the team's answer to the reporting group's question: a single forcing
  year at the very end should make no significant difference, and
  duplicating the 2299 data into 2300 is acceptable. The question had been
  put to the ISMIP7 lead. Runs therefore end in 2300, and the runners, the
  time axis and the filenames of cores 5 and 7 stand. For 2300 a CESM2-WACCM
  run reads the distributed atmosphere file, the 2290-2299 mean, for SMB;
  holds 2299 for the ocean `tf` and `so`, which is the duplicate the reply
  accepts; and under a mask mode reads the 2299 collapse mask. The v2.1
  masks for ssp126 and ssp585 end in 2299; their time axes were read on
  Quartz that day with h5py. No reader changes. Closed as
  icepack/ismip7#78.
- **#30, 22 September.** PISM posted its margin-only runs. Under CESM2-WACCM
  ssp585, the mask applied to all floating ice doubles the sea-level
  contribution, 1.25 to 2.48 m of ice above flotation. Applied only at the
  shelf margins it has hardly any effect. Under MRI-ESM2-0 the mask has almost
  no effect either way. So the two published end-member spreads are now
  about 8 % (UFEMISM) and about 100 %, with PISM's margin-only runs close to
  no fracture forcing. The same author followed up on how a model should
  treat fractured ice inside the shelf, as melange or as a reduced-viscosity
  region rather than open ocean. No one has answered that, and the
  organisers prescribe no mode. (issue #10)
- **#30, 23 September: the MRI-ESM2-0 ssp585 mask was faulty.** At 00:00 UTC,
  under the PISM post above, a member of the ISMIP7 Antarctica team at
  Dartmouth College reported an issue with the MRI-ESM2-0 ssp585 fracture
  mask and an update on Globus. The update is v2, whose README credits
  improved wind forcing in the excess meltwater, which changes all three
  fracture files. IU fetched the MRI-ESM2-0 fracture tree from Globus to
  Quartz that day, while the mirror still served the 29 August v1 (checked
  19:30 UTC): v2 arrived for ssp585 and ssp534-over, and the v1 files of
  ssp126, ssp370 and ssp534-over came back byte-identical, their collapse
  masks matching the manifest's ETags. Flagged 8 km cells in the ssp585
  mask, v1 against v2: 428 against 3,449 in 2100, 583 against 20,158 in 2200
  and 593 against 25,149 in 2250, with every v1 cell flagged in v2 at those
  years. So PISM's MRI-ESM2-0 result above ran on the faulty mask. Both
  versions end in 2299, as the CESM2-WACCM masks do. The reader opens v2, and
  the audit reads the row `AHEAD` and passes. `FRACTURE_MIN_VERSION` in
  `icepack2_tools/forcing.py` makes v2 the floor for MRI-ESM2-0 ssp585, and
  v2.1 the floor for CESM2-WACCM ssp126, ssp370 and ssp585: a run under a mask
  mode refuses an older mask at startup, and the audit reads it `OUTDATED`
  and fails, whatever the mirror serves. NOTS and Midway therefore need v2
  from Globus before a mask-mode core 8 can run there. (issue #16)
- **Unmoved.** #22 since 13:32 UTC on the 21st (the new 0.5.1 bounds stop
  taking a side on the `ligroundf` sign, but the thread prescribes none).
  #48 since the 20th: no date for the regenerated OCX. #17 and #37 as in
  section 7.

## 9. The `ligroundf` sign, settled 22 September

Thread #22 was read again through the GitHub API on 22 September. Its last
reply is still the 13:32 UTC one of the 21st: gain positive and loss negative
in general, and for `ligroundf` "it depends which part you are considering as
your reference". isschecker 0.5.1 bounds the field to [-10, 10] kg m-2 s-1.
The group settled the reference on 22 September: the grounded ice sheet, so
`ligroundf` is positive for grounded ice going afloat and negative where
floating ice flows onto grounded ice.

The writer already booked the grounded-to-floating direction positive. It
dropped the other direction: `book_advance` took the upwind outflow of the
grounded cell alone, so ice flowing from a shelf onto a pinning point was
never booked, and an ice rumple contributed its whole throughput to
`tendligroundf` as discharge. The form now books the full upwind facet flux,
the same flux the DG0 transport moved, signed from grounded to floating and
still landed in the floating cell. `tests/test_ismip7_annual.py` checks both
directions on a strip of cells, and the rumple case, whose net is zero.
Series banked by earlier runs carry the one-sided booking; they differ from
the new one only in cells with a facet across which ice flowed from floating
to grounded. The README's `[confirm]` on this convention is cleared and the
board item for it is closed.
