# MAP checks: one released MAP, native and transferred

`make -C antarctica map-check MAP_CHECK_FRICTION=regularized_coulomb|budd`
takes one released MAP through a fixed ladder of checks and prints where each
stands. It was built for two group decisions. The submission mesh and time
step were decided on 25 September 2026 (issue 20):
`antarctica_10000_1000_buffered20000` at dt 0.025 yr. The basal friction law,
regularized Coulomb or Budd, which the inversions at Rice carry, is still open
(issue #24). The ladder measures what
a forward from the MAP does on the MAP's own mesh and after transfer onto the
production mesh, so the group can read cost, stability, t = 0 fidelity and
ten-year drift for each (mesh, law) row side by side. It records evidence; it
recommends nothing.

## The MAPs it was built for

Rice released mid-descent snapshots of both 2 km MAPs on 22 September 2026
(release `maps-2km-snap-2026-09-22`, one file per law, 581 MB each, md5 in the
release notes):

| file | law | iterations kept |
|---|---|---|
| `inversion_icepack2_rc_n3_dg0_logvelnet_2000_int5000_bilap_snap20260922_0948.h5` | regularized Coulomb | 60 + 80 across two links |
| `inversion_icepack2_budd_n3_dg0_logvelnet_2000_int5000_bilap_snap20260922_0948.h5` | Budd | 120 + 120 across two links |

Both were inverted on `antarctica_5000_2000_buffered0.msh` (2 km / 5 km gmsh,
925,183 vertices, no ocean buffer), DG0 geometry, n = 3, bilaplacian prior
with sigma 0.3 and rho 7500 m, gamma 1e5 on both controls, code `be5d685`.
Each carries `log_friction`, `log_fluidity`, `fluidity_prior`, bed, thickness,
surface, `velocity_obs` and `obs_mask`, and no `velocity`: a forward re-solves
the diagnostic on start. The final MAPs, with velocity, follow when the chains
end (issue #24); the ladder then runs again on them and the snapshot rows stay
in the table marked superseded.

## The transfer, and the hazard it closes

With `ISMIP7_INVERSION` naming a MAP and `ISMIP7_MESH` naming another mesh,
`simulation.setup_model` keeps the MAP's mesh as the interpolation source and
solves on `ISMIP7_MESH`: the continuous fields (`log_friction`,
`log_fluidity`, `fluidity_prior`, `velocity_obs`) are interpolated across, the
DG0 geometry is rebuilt from BedMachine cell averages on the target, and the
frozen anchors (`N_ref`, `C_w0`, `phi_eff`, `H_init`) are rebuilt on the target
geometry. That is the path the timing matrix's lanes take from the 2.5 km MAP.

The production mesh carries a 20 km ocean buffer and the 2 km MAPs carry none,
so every target dof in that ring lies outside the source mesh. Firedrake's
cross-mesh interpolate wrote `0.0` there. Zero is the prior for theta and phi
and a singular block for the fluidity prior: `A_eff = A_prior * exp(phi)`
multiplies the dislocation term, the `lin_reg` regularizer and the `alpha_gl`
collar in `dual_friction.build_rc_residual`, all at once. The loader now
fills each field with a stated value (`icepack2_tools/transfer.py`,
`interpolate_with_fill`):

| field | fill outside the source mesh |
|---|---|
| `log_friction`, `log_fluidity` | 0, the prior |
| `fluidity_prior` | the constant baseline `A0 * a4_factor`, the value the code uses when a MAP carries no prior at all |
| `velocity_obs` | the raster sample the forward makes on its own mesh |

A second artefact comes with the first. Firedrake locates a target point in a
source cell up to half a reference cell outside it (`mesh.tolerance`, 0.5 by
default) and evaluates that cell's linear basis there, so ring points within
about a kilometre of the 2 km front are extrapolated. Measured on the first
Quartz pass (22 September, jobs 10569251 and 10569254): the transferred prior
spanned [-218.69, 1028.05] from a source spanning [1.00, 783.69]; the same
snapshot onto the buffered 5 km production mesh gave [-164.7, 921.1] (issue
#81). The transfer therefore locates strictly: the source mesh's tolerance is
1e-8 for the interpolation and restored afterwards, so a dof the source does
not contain is a fill whatever its distance from the front. Behind that,
inside a cell linear interpolation stays within the cell's vertex values, so
a value beyond the source field's range can only be an extrapolation;
located dofs are still clamped to the source's range, component by
component, and counted, as a guard on the tolerance. The snapshot rows below
were measured before strict location, with the clamp alone: the 996 clamped
prior dofs they report are the band that now takes the fill.

Every filled or clamped field prints one `Transfer fill:` line with both
counts, the counts go into the context (`ctx["transfer_fill"]`), the cache
manifest and the score JSON, and a fluidity prior whose minimum is not
positive after loading aborts the run with the reason (which is what caught
the extrapolation). A same-mesh load misses nothing and prints nothing.

`ISMIP7_MESH=checkpoint` names the mesh embedded in the MAP or restart file.
`site_env.sh` always exports a derived `.msh` path, so this sentinel is the
only way a job submitted through `submit.sh projection` runs MAP-native, and
Quartz holds no `antarctica_5000_2000_buffered0.msh` at all: the native rows
below take their mesh from the checkpoint.

## The ladder

One MAP, one law, both meshes. Every stage stamps
`results/map_check/<map stem>/status_<stage>.txt`; the table reads each stage
off its artifacts first and its stamp second, and `make map-check` submits the
first stage whose dependencies passed (`MAP_CHECK_STAGES=all` submits every
runnable one, `MAP_CHECK_STAGES=<stage>` one, `FORCE_TIMING=1` resubmits).
`make map-check-dry-run` prints the table and every stage's composed sbatch
line. The manager is `scripts/manage_map_check.py`; it composes its requests
through `submit.sh` like the timing campaign and shares that code
(`SlurmStageRunner` in `scripts/manage_timing_campaign.py`).

Native means `ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_BUFFER_M=0
ISMIP7_MESH=checkpoint`, dt 0.1 (the matrix rule `0.125 * LC / 2500`), the
tracked `boundary_ids_antarctica_5000_2000_buffered0.json`. Transferred means
`ISMIP7_LC=1000 ISMIP7_LC_COARSE=10000 ISMIP7_BUFFER_M=20000`, the production
mesh, dt 0.05 for the strict lane (the matrix rule, so the lane compares with
the matrix row) and the production step, 0.025, for the control. The 22
September pass recorded below ran its 1 km control at 0.05, before the step
was chosen. Every job carries `ISMIP7_FRICTION=<law>`, DG0, n = 3, and the
lanes and controls carry the runaway tripwire defaults of
`timing_campaign.TRIPWIRE_DEFAULTS`.

| # | stage | job | passed means |
|---|---|---|---|
| 0 | `fetch` | curl of the release asset on the login node, md5 against the release notes (`MAP_CHECK_MD5=auto` reads them; a hex value pins it), the h5 root attributes read where h5py is importable | md5 equal; `friction`, `geometry_space`, `n_flow`, `mesh_basename`, `lc`, `lc_coarse`, `buffer_m` as expected |
| 1 | `repack` | `timing_redistribute.script`, one rank, 32G: the one-rank rewrite every later job loads at its own rank count; the repack keeps the release basename, so every record names the release file | the repack exists |
| 2 | `score_native` | `map_check_score.script`: `score_map.py --json` on the MAP's mesh, 16 ranks; under Budd also the `check_budd_map.py` shelf-gate census | the continuation converges, the discharge ratio is finite, no dof was filled, the fluidity prior minimum is positive |
| 3 | `prepare_transfer` | `timing_prepare.script` with the law and the target mesh, cache role `map-check-initial-state`, under `<stem>/cache/`, 16 ranks 96G, 1 to 3 h | `validate_map_check_manifest` passes with the MAP's and the mesh's sha256 and the `scpc_mumps` fingerprint; the fill counts are in the manifest |
| 4 | `audit_cache` | `timing_cache_audit.script`: the t = 0 finite-volume tendency of the transferred state | the audit JSON exists (record only) |
| 5 | `score_transfer` | `map_check_score.script --restart <cache>`: the same discharge score on the transferred state | finite ratio, positive prior minimum |
| 6 | `lane_transfer` | `timing_transient.script`, kind `map_check`, restart from the cache, 64 ranks under `MAP_CHECK_SOLVER` (`scpc_gamg`), the strict contract, 10 steps of dt 0.05 | the `make qualify` rule (whole interval, no diverged solve, mass residual at or under 5e-5 Gt), no tripwire, rescue off, `initial_state_source` is the release file |
| 7 | `lane_native` | the same lane on the MAP's mesh, cold start from the MAP inside the lane (setup is timed apart from the steps), 10 steps of dt 0.1 | the same rule |
| 8 | `control_transfer` | `submit.sh projection ISMIP7_EXPERIMENT=control ISMIP7_T_END=2025 ISMIP7_OUTPUT=1`, the production defaults (`ISMIP7_APPARENT_MB=1`, `ISMIP7_FIXED_FRONT=1`, `scpc_gamg`, dt 0.025, self-chaining), cold start from the MAP through the transfer inside the job, the K file `MAP_CHECK_K_NPZ` names, `ISMIP7_RUN_TAG=mapcheck_<law>_<snap>_<lc>` | the timeseries reaches 2025 with `resid` at 0.00 on every row and the final state written |
| 9 | `control_native` | the same control on the MAP's mesh | the same |
| 10 | `audit_controls` | `map_check_audit.script`: `check_ismip6_track.py` on both series (exit codes kept), `compare_runs.py` overlay, `region_budget.py` at each final state under its own mesh triple | the audit JSON and the figure exist |
| 11 | `summary` | the manager writes `<stem>/summary.md` from whatever JSON exists | always |

Dependencies: `fetch`, then `repack`, then `score_native`, `prepare_transfer`
and `lane_native` at once; `prepare_transfer` gates `audit_cache`,
`score_transfer` and `lane_transfer`; each lane gates its control; both
controls gate `audit_controls`. The controls cold-start from the MAP because
`setup_model` refuses `ISMIP7_APPARENT_MB=1` on a restart that carries no
`a_ref_mb`, and because that is the path production takes. The control driver
prints its "no historical endpoint" warning on a cold start: these controls
measure drift and are not a projection baseline. A lane that fails on the
state itself (the snapshot's own runaway cells, on both meshes) would leave
the controls blocked; `MAP_CHECK_CONTROLS_AFTER_FAILED_LANE=1` runs them
anyway, because the production configuration cancels the t = 0 tendency the
strict contract exposes and its cost per simulated year and drift are what
the mesh question needs. `summary.md` then says which lane failed, and those
control numbers are production numbers and no stability verdict.

Stage 12, by hand when the final MAPs land: `ISMIP7_CHECK_FRICTION=<law>
check_budd_map.py <final> --forward` on the MAP's mesh
(`budd_map_census.script`). A MAP inverted under the forward's residual
reproduces its own velocity to about 1e-7 (issue #24 carries the exit
criterion); the snapshots have no velocity to compare against.

## How to read the table

- **Does the transfer work.** `prepare_transfer` passing with a positive
  prior minimum and a continuation that converged, `score_transfer` within
  0.05 of `score_native` overall and 0.10 per band (0.05 is the observed
  discharge uncertainty, 100 of 2050 Gt/yr, Rignot 2019), and
  `lane_transfer` passing. A filled dof inside ice (a cell thicker than 1 m
  beyond the source outline) means the buffer-0 source does not cover the
  buffered target's ice; the fallback is an `antarctica_10000_1000_buffered0`
  mesh, to be built and re-timed.
- **Cost.** `seconds_per_step` from each lane record, on 64 ranks, against the
  matrix row for the 1 km mesh from the 2.5 km Budd MAP: 30.0 s per step,
  10.0 minutes per simulated year, 47.5 h per 285 years, 1.0 GiB per rank
  (`TIMING_MATRIX_QUARTZ_SCPC_GAMG.md`). Minutes per year is
  `seconds_per_step / dt / 60`; the eleven-experiment matrix is about 2260
  simulated years of transient loop.
- **Stability.** The lane verdicts, the tripwire, the Newton iterations per
  step (12.9 on the 1 km reference row) and the condensed iterations.
- **t = 0 fidelity.** The discharge ratio overall and per speed band on each
  mesh, the initial misfit against the observations, the filled dof counts,
  and the `Apparent MB: a_ref in [lo, hi] m/yr, net X Gt/yr` line of each
  control's log. The misfit line is the mean squared velocity error over the
  whole mesh, open ocean cells and the buffer ring included, where the
  thin-ice velocity is unconstrained and large; it is context for one mesh
  and no measure across the two (the Budd snapshot gives 1.29e5 on its own
  mesh and 2.02e4 transferred). The discharge ratio is the fidelity number.
- **Ten-year drift.** From each control's timeseries over 2016 to 2025:
  dVAF/dt, dM/dt, the discharge in 2016 and 2025 and its block growth, the
  grounded area from `iareagr`, melt and SMB, and the track verdict per row
  (`check_ismip6_track.py`: dVAF/dt warns outside 2 and fails outside 5 mm
  SLE/yr; the peak clause raises a known false positive, issue #33, written
  as such).

Confounders that go with every table: the snapshots are unconverged and at
different iteration counts under log-velocity weights re-derived per chain
link (issue #68); Budd carries the `ISMIP7_ALPHA_GL=0.5` grounding-line collar
and a frozen `N_ref` that regularized Coulomb has no counterpart to; the K
file was fitted under the local slope on a 2500 m mesh and the forward now
defaults to the constant Antarctic slope (issues #26, #30); apparent mass
balance zeroes the t = 0 tendency, so only the later drift and the size of
the correction separate rows; the native mesh has no buffer and the
production mesh 20 km of it, so front bookkeeping differs and the native
against transferred comparison carries both the mesh and the transfer.

## Results

### Snapshot pass, 22 September 2026 (issue #20, issue #24)

Release `maps-2km-snap-2026-09-22`: RC md5 `30d64e08ec9dd651d8a67eb00da3e30d`
(140 iterations kept, source sha256 `6f2c8694f360dd6f8017837913604d0412c05d9618034437a1ec89e57d850bbf`),
Budd md5 `1c5d1031651f873e386c66aea34a69fe` (240 iterations, source sha256
`fdc66d42ad8357a0cdb66ede45c8ce3b266ec32739c5e075e520c082bcdb3310`), both on
`antarctica_5000_2000_buffered0.msh` (925,183 vertices), MAP code `be5d685`.
Forward code: this branch at 7a0c7e0 for the physics (the later commits on
the branch touch the site setup, the cache publisher and the manager); the
Budd native score ran before the clamp commit, on the MAP's own mesh, where
no transfer happens. IU Quartz, `general` partition, one 128-core node per
job; `scpc_mumps` for the scores and the prepares on 16 ranks, `scpc_gamg`
for the lanes and controls on 64 ranks. Melt slope the `ant` default, K file
`K_issue11_mesh2500.npz`, `ISMIP7_N_FLOW=3.0`, `ISMIP7_A4_FACTOR=1.0`, DG0
geometry, tripwires `u_max 2e4 m/yr, h_max 5000 m, (dh/h)/dt 20/yr over
cells thicker than 100 m`. Banners: RC `regularized Coulomb (c0=0.5,
h_visc_floor=10m, cw0_floor=0.0e+00, eps_tauc=0.0e+00 MPa, alpha=1.0e-02)`;
Budd `Budd N_hat (N_ref=reference; exact-zero shelf; delta=0.020,
N_hat_cap=3.0, alpha_gl=0.50, h_visc_floor=10m, ocean_drag=1e-02@h<10m,
u_lim=2e+04)`. The rows below are for the snapshots and are superseded when
the final MAPs are run.

**Transfer (both laws, `prepare_transfer`).** 164,735 of 1,869,088 vertex
dofs per continuous field lie in the 20 km ring outside the source mesh and
took the stated fill; the clamp caught 996 located dofs on the fluidity prior
and 2 on log fluidity, none on log friction or the observed velocity. The
prior on the target is [1.00, 783.69], the source range. The 8-step
continuation converged in about 15 minutes on 16 ranks for both laws (Budd
job 10569349, RC 10569531). Every transferred dof inside ice came from the
source: no cell thicker than 1 m lies beyond the source outline (the ring is
open ocean).

**t = 0 (`score_native`, `score_transfer`).** Q(u_model)/Q(u_obs) across the
grounding line, overall and per speed band of the observations; `misfit0` is
the whole-mesh mean squared velocity error and is context for one mesh only.

| law | mesh | Q ratio | bands (<100, 100-500, 500-1500, >1500 m/yr) | Q(u_obs) Gt/yr | misfit0 | job |
|---|---|---|---|---|---|---|
| Budd | 2 km / 5 km, native | 1.06 | 1.38, 0.99, 1.08, 0.92 | 2151 | 1.29e5 | 10569253 |
| Budd | 1 km / 10 km, transferred | 0.67 | 0.90, 0.66, 0.64, 0.48 | 2378 | 2.02e4 | 10569518 |
| RC | 2 km / 5 km, native | 2.75 | 6.62, 2.62, 1.78, 0.68 | 2151 | 1.19e5 | 10569455 |
| RC | 1 km / 10 km, transferred | 2.67 | 6.41, 2.54, 1.74, 0.52 | 2378 | 2.97e4 | 10569708 |

The Budd snapshot reproduces the observed grounding-line discharge on its
own mesh to within the observational uncertainty overall, with the slow band
a third high. Transferred, it discharges a third less in every band. Two
things the transfer rebuilds explain the direction: the friction anchor
`C_w0 = tau_d / |u_obs|^(1/m)` (`weertman_anchor`) is recomputed from the
1 km BedMachine cell averages and the interpolated observations while the
inverted log adjustment `theta` was fitted against the 2 km anchor, so the
transferred friction is `C_w0(1 km) exp(theta(2 km))` and not the inverted
field; and the grounding line and its thickness are rebuilt at 1 km, which
moves Q(u_obs) itself by a tenth (2151 to 2378 Gt/yr). Budd alone also has
`N_ref`, the shelf gate and the `alpha_gl` collar that follow the rebuilt
grounding line. The RC snapshot discharges 2.7 times the observations on its
own mesh, 6.4 times in the slow band, and the transfer leaves that ratio
where it is: at 140 iterations the RC descent is far from the observations,
or the forward's RC form and the inversion's differ; the final MAP's
self-consistency check (stage 12, issue #24) separates the two.

**Strict 10-step lanes (`lane_transfer`, `lane_native`).** All four lanes
failed at their first step on the runaway tripwire, on hotspots the states
carry at t = 0 (the cache audit's `no_forcing_dhdt`, the thickness tendency
with no forcing and no apparent-mass-balance correction, names the same
cells before the step is taken). The failures are the snapshots' own: the
Budd hotspot is the same Lambert Glacier confluence cell on both meshes.

| law | mesh | step-1 diagnostic | tripwire | job |
|---|---|---|---|---|
| Budd | 2 km / 5 km, native | 91 Newton, 21,871 condensed iterations, 596 s; 25.6 s for the step | speed 1.00e5 m/yr at (1698255, 700000), Lambert confluence; dh +522 m in 0.1 yr on a 1365 m grounded cell there; (dh/h)/dt 11/yr on a 110 m floating cell at (-1608931, -328070), Pine Island Bay | 10569571 |
| Budd | 1 km / 10 km, transferred | 19 Newton, 566 condensed iterations, 42 s; 43.6 s for the step | (dh/h)/dt 25/yr on a 126 m grounded cell at (-1244404, 136859), Rutford Ice Stream area; dh +420 m in 0.05 yr on a 1186 m cell at (1692481, 701685), Lambert confluence; speed max 1.65e4 m/yr there | 10569519 |
| RC | 2 km / 5 km, native | no transient step: the gamg continuation converged its first step after 64 Newton and 54,867 condensed iterations (1374 s) and diverged on the second (200 Newton iterations, residual 1.2e6, 1587 s), then restarted with 16 steps; the lane was stopped there after 53 minutes. The mumps continuation of the native score converged in 8 steps (42, 10, 16, 8, 8, 9, 8, 9 Newton iterations) | none reached | 10569456 |
| RC | 1 km / 10 km, transferred | 24 Newton, 1,136 condensed iterations, 76 s; 78.2 s for the step | speed 1.09e5 m/yr at (-2412787, 1261392), northern Antarctic Peninsula; (dh/h)/dt 57/yr on a 209 m grounded cell at (-1551062, 836317); dh +720 m in 0.05 yr on a 25 m cell at (-2312121, 1007756) | 10569709 |

Cache audits of the transferred states (`audit_cache`, `no_forcing_dhdt`
over the 288,090 cells thicker than 100 m): Budd extremes +8103 m/yr at the
Lambert cell above and -4233 m/yr beside it, net +21.6 Gt/yr; RC extreme
-55,214 m/yr on a 462 m grounded cell moving 43,681 m/yr at (424771,
-1800789), Victoria Land coast, with a cluster of grounded 100 to 230 m cells
at 66,000 to 101,000 m/yr around (-1546500, 831800), net +31.3 Gt/yr. A
single step under the strict contract cannot carry either; only the 1 km
Budd step-1 solve is a normal one (19 Newton iterations against the 12.9
reference; 43.6 s against 30.0 s per step, one step, no warm-up excluded).

**10-year controls (`control_transfer`, `control_native`, production
configuration, run under `MAP_CHECK_CONTROLS_AFTER_FAILED_LANE=1`).** Each
64 ranks, 128 GB, 5 h (resubmitted at that request to backfill a full
partition); cold start from the MAP through `projection.sbatch`, apparent
mass balance in balance mode and uncapped, fixed front, dt 0.05 (1 km) and
0.1 (2 km), the same tripwires as the lanes. The RC 2 km control was not
run: it would cold-start through the same gamg continuation that diverged in
the lane, and the RC snapshot is already 2.7 times the observed discharge on
that mesh.

| law | mesh | cold start | apparent MB | outcome | job |
|---|---|---|---|---|---|
| Budd | 1 km / 10 km, transferred | 8 continuation steps under gamg, 15 then 6 to 8 Newton iterations each, about 6 min | a_ref in [-8101.5, +4236.8] m/yr, net -1184.6 Gt/yr; budget SMB +2532, melt -1387, outflux -4, calving -50 Gt/yr | killed at step 16 (t = 2015.8) by the speed check, 3.84e4 m/yr at (1695321, 699618), the Lambert confluence. From step 12 the cell at (1695698, 700055) oscillated with a period of two steps and a growing amplitude, +73, -178, +282, -481, +748 m per step, flipping between floating and grounded each time; elsewhere the largest (dh/h)/dt stayed near 2 to 4 per year. Diagnostic solves of 9 to 27 s at 4 to 14 Newton iterations per 0.05 yr step, about 4 to 5 min per simulated year before that. The series has rows at 2015.0 and 2015.8 and no drift to read | 10570030 |
| Budd | 2 km / 5 km, native | 8 steps, 13 then 7 to 11 Newton iterations, about 9 min | a_ref in [-5429.4, +3082.9] m/yr, net +2157.1 Gt/yr; budget SMB +2431, melt -1242, outflux -3368 through the domain boundary at the front, calving -51 Gt/yr | killed at step 13 (t = 2016.3) by the speed check, 8.93e4 m/yr at (1687863, 702000), the Lambert confluence again. The cell at (1688440, 703000) went +247, -424, +899 m on steps 11 to 13, grounded, floating, grounded; the first ten steps were calm (largest (dh/h)/dt 0.5 to 6 per year on thin shelf cells, speed max 1.75e4 m/yr constant). Diagnostic solves of 12 to 22 s at 6 to 14 Newton iterations per 0.1 yr step, about 3 min per simulated year | 10570031 |
| RC | 1 km / 10 km, transferred | 8 steps under gamg: the first took 70 Newton iterations, 26,600 condensed iterations and 23 min, the rest 8 to 17 Newton iterations and 35 to 75 s each (the same continuation diverged at its second step on the 2 km mesh) | a_ref in [-36709.9, +55214.0] m/yr, net -1174.9 Gt/yr | killed at step 1 (t = 2015.05) by the speed check, 9.77e4 m/yr at (394334, -1806697), Victoria Land coast, the cell the cache audit named; (dh/h)/dt 17.3 per year and a 25 m cell grown to 599 m on the Antarctic Peninsula in that step | 10570032 |

The budget columns differ between the meshes by bookkeeping and not only by
state: on the buffer-0 mesh the calving front is the domain boundary, so the
front flux leaves as `outflux`; on the buffered mesh the front is inside the
domain and the fixed front removes it as `calv`.

**What the pass says.** The transfer works as a mechanism: the ring is
filled and counted, the extrapolated dofs clamped, the prior positive, the
continuation converges, the cache publishes, audits and restarts. It changes
the Budd state (discharge ratio 0.67 against 1.06) for the two reasons
above, the rebuilt friction anchor and the rebuilt grounding line, and
leaves the RC state where it is. Neither snapshot runs forward: the strict
lanes fail at step 1 on both meshes, and under the production configuration
the Budd state grows a two-step oscillation at the Lambert confluence on
both meshes, at the same cell, with the cell's grounded state flipping each
step, and is killed within 1.3 simulated years; the RC state carries
grounded cells at 66,000 to 101,000 m/yr and the speed check kills its
control at the first step, after a gamg cold start that converged on the
1 km mesh (23 min for its first step) and diverged on the 2 km mesh. The
Budd instability is therefore the snapshot's (or the Budd forward's handling
of a cell at flotation there) and not the mesh's or the transfer's, and the
resolution question gets no cost or drift row from these snapshots beyond
the setup costs (a 1 km cold start of 6 min for Budd and 30 min for RC on 64
ranks, 4 to 5 min per simulated year in the calm Budd steps; 2 km about
9 min and 3 min) and the per-step Newton counts. Both decisions wait for the
final MAPs (issue #20, issue #24). The transferred caches and the control
series stay on Quartz scratch until its purge (about 22 October 2026); the
records, scores, audits, summaries and job logs are copied under
`/N/project/ice_rheology/ISMIP7/antarctica/results/map_check_snap20260922/`.

Not measured on this pass: a clean cost per step on either mesh (every lane
stopped at step 1 or before; the controls give the production wall time per
simulated year instead), the 10-year drift, and anything on the RC 2 km
mesh beyond the score. The RC descent's diagnostic solve under `scpc_gamg`
on the buffer-0 mesh is itself a finding for the friction decision: the
production linear solver did not carry its continuation where `scpc_mumps`
did.

## Running it on Quartz

From a scratch clone of the branch (`/N/scratch/dlilien/ismip7_<branch>`,
with `sites/local.env` naming the shared `ISMIP7_DATA_ROOT` and
`ISMIP7_OBS_DATA_ROOT`), from `antarctica/`, once per law and again as each
stage settles:

```bash
make map-check MAP_CHECK_FRICTION=regularized_coulomb \
  MAP_CHECK_TARGET_MESH=/N/project/ice_rheology/ISMIP7/antarctica/mesh/antarctica_10000_1000_buffered20000.msh \
  MAP_CHECK_K_NPZ=/N/project/ice_rheology/ISMIP7/antarctica/results/issue11_melt_check/K_issue11_mesh2500.npz
```

The whole ladder for one law is about nine to ten hours of wall time and
roughly 700 core-hours, and the two laws run side by side. On a full
partition the controls' default request (the production forward's, 240 GB
and 12 h) waits for a whole free node; `MAP_CHECK_CONTROL_MEM=128G
MAP_CHECK_CONTROL_TIME=05:00:00` backfills sooner (the lanes peaked near
1 GB per rank), and a control that runs out of wall time chains itself. The mechanical
half of the question, whether the transfer works, is answerable after
`prepare_transfer`, `score_transfer` and `lane_transfer`, four to five hours
in.
