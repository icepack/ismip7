# ISMIP7 Projections-Antarctica README (draft)

Draft of 19 September 2026 against the June 2026 template (the Google
document linked from discussion #6). Submit one README per ice sheet, saved
as `README_AIS_RICE_icepack2`. Every item marked **[confirm]** needs a
decision or a name from the group before submission; everything else is
what the code on this branch does today. The forum asks that several
modelling choices be stated here rather than settled centrally; those carry
their discussion number.

Contributor names, affiliations and emails: **[confirm #38]** Andrew Hoffman
(Rice University, ah301@rice.edu) and collaborators at Indiana University and
the University of Chicago.

Date of submission: **[confirm #38]**

Ice Sheet Modeled/domain_id: AIS
Modeling group name/source_id: RICE **[confirm #38]**
Ice Sheet Model Name/ism_id: icepack2 **[confirm #38]** (icepack2 dual
shallow-shelf formulation on Firedrake 2026.4.1)

## Initialization methods

1. Data assimilation at a single epoch, no spin-up. The geometry
   (thickness, bed, surface) is BedMachine Antarctica sampled onto the
   model's cells; the basal friction and ice fluidity fields are inverted
   with the adjoint (L-BFGS-B, 200 iterations) as log-adjustments to a
   driving-stress Weertman anchor and a thermally derived fluidity prior.
   The misfit is the sigma-normalised velocity misfit against MEaSUReS
   (per-component errors, floored at 3 m/yr), ISSM's logarithmic velocity
   misfit, a thickness-tendency term against the observed 2003-2019 mean
   dH/dt (Smith et al. 2020, ISMIP7 observations kit) obtained from one
   implicit-Euler prognostic step with the model's own transport operator,
   and a term on the integrated net grounded mass balance. Tikhonov
   regularisation on both controls.
2. Targets: MEaSUReS v2 velocity where observed (the observation mask
   excludes the pole hole and unobserved cells), the Smith dH/dt over
   grounded ice, the integrated net mass balance. Outputs: log friction
   adjustment theta, log fluidity adjustment phi, and the velocity field
   consistent with them.
3. The velocity misfit (median 16-19 m/yr over observed nodes; the figure is
   mesh-specific and is to be re-measured on the chosen production mesh),
   the grounding-line discharge scored against the flux
   the observed velocity carries across the same facets
   (`antarctica/scripts/score_map.py`), a check that a forward re-solve
   reproduces the inversion's velocity, and a one-year balanced control
   whose volume above flotation and mass hold (`resid = 0`).
4. SMB climatology: RACMO2.4p1 (ERA5-forced) monthly SMB, 2000-2023 mean,
   used as the constant control climate and as the baseline the ISMIP7
   anomalies are added to. Ocean climatology: the ISMIP7 observational
   thermal forcing and salinity climatology (30_sep OI product). Calving
   fronts and ice margins are held at their BedMachine position during
   initialisation (there is no spin-up).
5. Yes. An apparent-mass-balance reference is computed once at the initial
   state with the model's own transport operator so that the initial
   thickness tendency is exactly zero; it is frozen and applied in every
   experiment. It is
   NOT reported in `acabf` (which is the SMB the model applied) but written
   alongside it as `acabf_correction` in the same units for anyone closing
   the budget.

## Projections: ice-ocean and ice-shelf fracture (AIS)

6. Ocean melt: the ISMIP7 quadratic parameterisation of Burgard et al.
   (2022), local-quadratic variant (TF_avg = TF), with one constant
   `sin(alpha)` = 5.115e-3 on every shelf. This is the ISMIP7 reference's
   mean Antarctic slope and the value the toolbox's K percentiles (K05
   4.75e-5, K50 8.5e-5, K95 1.375e-4, July 2026) were sampled with; the
   model's own draft slope is kept as an option only
   (`ISMIP7_MELT_SLOPE=local`; **[confirm #26]**, see
   `FORWARD_RUN_READINESS.md` action 5). Constants from `multimelt.constants`.
   K is dimensionless and per IMBIE basin, fitted with
   `antarctica/scripts/calibrate_melt.py` to the observed basin totals
   through the model's own melt path; the forward applies this per-basin
   field and `ISMIP7_K_SCALE` multiplies it. K* = 4.06e-5 on the 865 Gt/yr
   table (2500 m mesh) and 4.46e-5 on the July 2026 table (2 km mesh), just
   under K05; see `GEOMETRY_DISCRETIZATION.md`.
   **[confirm #42]** that every submitted run read this calibration. Thermal
   forcing (`tf`) and salinity (`so`) are read at the cell's draft from the
   ISMIP7 ocean forcing, nearest neighbour in depth and in the plane, and
   used as provided: no smoothing inside the cavities (discussion #11: the
   warm stripe through the Ross cavity is a feature of the climatology, and
   smoothing within a basin is allowed, not required). Where the forcing is
   undefined the thermal forcing is zero and the salinity 34.5. The last
   ocean year (2299 for CESM2-WACCM) is held for 2300.
   Partially floating cells: the geometry is cell-wise (DG0); a cell is
   floating when its height above flotation is negative and then receives
   the full melt, grounded cells none. No melt law acts on vertical ice
   fronts; the melt of ice that flows into an emptied marine cell is
   reported as front melt, `lifmassbf` (see the conventions below).
7. Grounding line: the flotation criterion per cell (height above
   flotation from thickness and bed); no sub-cell parameterisation.
   Basal friction is a regularised Coulomb law (`c0 = 0.5`, exact-zero on
   floating cells because the Coulomb cap is proportional to the
   effective pressure) or Budd (`N_hat = N/N_ref` with the PISM delta floor,
   gated to grounded cells by height above flotation). Effective pressure
   is the ice overburden minus the ocean pressure, so it is the migrating
   grounding line that switches friction off.
8. Calving: in the control and, by default, in the projections the calving
   front is pinned at its 2015 (BedMachine) position: ice flowing past it
   is removed and tallied as calving (`ISMIP7_FIXED_FRONT`). A level-set
   front with a von Mises calving law (Hahn, Mikula and Frolkovic 2025
   finite-volume level set; thresholds 1.0 MPa grounded, 0.15 MPa floating)
   exists but is not calibrated; **[confirm #36]** which the projections use.
   No sub-grid scheme beyond the sub-cell shed of the level-set front.
9. Ice-shelf collapse: **[confirm #10]** which of three the submitted
   projections use; the core matrix as configured today runs the first.
   `ISMIP7_FRACTURE=none`: no collapse forcing. `mask`: every FLOATING cell
   the year's ISMIP7 collapse mask flags (v2.1 for CESM2-WACCM; for
   MRI-ESM2-0, v2 under ssp585 and v1 under ssp126 and ssp370) is emptied
   and booked as calving, wherever it is. The masks flag the Ross and
   Filchner-Ronne shelves near their grounding lines first, so this opens
   holes far behind the front which the momentum balance treats as open
   ocean. `mask_front`: a flagged floating cell is
   emptied only once open water has reached it through flagged cells, so a
   shelf collapses from its front and nothing happens until the flagged
   region touches it. These are the two end-members of discussion #30
   (September 2026), where groups report 40 % to 100 % more sea level by
   2300 from the first against the second. Grounded ice is never touched,
   there is no stress condition (Lai et al. 2020), and the excess-meltwater
   and lake-property products are not used. The ssp126 and ssp585 masks of
   both ESMs end in 2299, and 2300 reads the 2299 mask. No collapse forcing
   in the historical, control or OCX runs, for which none exists
   (discussions #29 and #33); the front there is pinned at its 2015 position
   rather than following the observed fronts.
10. Tributary glaciers after a collapse: no special treatment; the front
    retreats to the new extent, the grounding line responds to the lost
    buttressing through the momentum balance, friction is unchanged.
11. The control holds the 2000-2029 climate constant, with the
    apparent-mass-balance correction. Its ocean is the ESM's own `ctrl` tree
    (`tf` and `so`, v3), read and melted as in the projections (item 6), so
    the per-basin K fitted against the 30_sep OI climatology applies
    unchanged to every ESM ocean. Its SMB is the RACMO climatology, the
    baseline the projections add the ISMIP7 anomalies to after re-referencing
    them to the 2000-2029 pool (historical 2000-2014 and ssp126 2015-2029);
    in that frame RACMO is the 2000-2029 climate. C009 and C010 differ in
    their ocean and in the historical endpoint they branch from.
    The historical runs start in 1850 from the 2015 initial state (there is
    no spin-up) and end at 1 January 2015, where the projections and the
    control branch. The ISMIP7 anomalies are relative to 1960-1989, and
    adding them to a modern baseline is what produces the jump at the start
    of a historical run that discussion #34 describes; here the anomaly's
    2000-2029 mean is removed first, no temperature forcing is applied, and
    the frozen apparent-mass-balance reference holds the initial state in
    balance, so there is no such jump. The 165 years before the initial
    state's epoch are a relaxation under that reference, not a hindcast.
    OCX (C011) runs 1979-2025. **[confirm #19]** which forcing the submitted run
    used: the ISMIP7 OCX product (RACMO2.3p2-ERA SDBN1 `acabf`, full field,
    and the expert-judgment `main` ocean, `ISMIP7_OCX_FORCING=protocol`), or
    RACMO2.4p1 actual-year SMB with the constant ocean climatology
    (`stopgap`), which is what the core ran on until September 2026. K is
    fitted to the climatology in either case, and the OCX `main` ocean
    departs from it around Mertz (discussion #48, open).

## SMB questions

19. SMB is applied as a cell-mean source in the finite-volume thickness
    transport: RACMO2.4p1 climatology plus the ISMIP7 `acabf-anomaly`
    (SDBN1 8 km, v2 for CESM2-WACCM, v1 for MRI-ESM2-0) re-referenced so the
    anomaly's mean over the control window vanishes. No surface-elevation
    feedback: neither the runoff gradient `dmrrodz` the protocol prefers nor
    `dacabfdz` is used (discussion #36), and no lapse rate (`dtsdz`), there
    being no thermal model. Precipitation is not used. Forcing is NaN
    outside the downscaled mask and is filled with zero there (discussion
    #39). CESM2-WACCM ends in 2299 in its CMIP archive; its 2300
    atmosphere files, as distributed, are the 2290-2299 mean, and are read
    as given (discussions #8 and #49). Its ocean (item 6) and collapse masks
    (item 9) hold 2299 for 2300. The organisers' reply of 23 September 2026
    on discussion #49 judges a single forcing year at the end insignificant
    and accepts a 2299 duplicate for 2300. Monthly fields are averaged to the year weighted
    by month length from each file's own time axis, and the year a file
    belongs to is taken from its name, so the differing calendars and time
    stamps of the forcing products (discussions #9 and #24) do not enter.
20. The apparent-mass-balance correction of item 5 is applied in every
    experiment, frozen at the initial state.

## GIA, bedrock and sea level

21. No bedrock adjustment.
22. BedMachine Antarctica's reference: elevations relative to the EIGEN-6C4
    geoid.
23. Not applicable.
24. Not applicable.
25. Far-field sea-level change is not included.

## Other general questions

Ice-covered area: BedMachine's mask at the initial epoch; peripheral
glaciers off the main sheet are outside the mesh; the target mask is the
initial extent, enforced at run time by the pinned front (ice past it is
removed) in the control and default projections.

PPE / ESM participation: **[confirm #38]**.

Summary paragraph: **[confirm, draft #38]** icepack2 is a finite-element
shallow-shelf model on Firedrake in its dual (velocity, membrane stress,
basal stress) formulation, with a first-order upwind finite-volume
thickness transport on the same unstructured mesh (resolution **[confirm #20]**,
pending the 1000 m inversions: 2 km at the grounding line coarsening to
180 km in the interior, or 1000 m / 10 km), an adjoint initialisation to
MEaSUReS velocities and observed thickness change, regularised Coulomb
sliding, the ISMIP7 quadratic mixed-slope ocean melt with per-basin
calibration, a pinned or level-set calving front, and the ISMIP7 collapse
masks. References: Shapero et al. 2021 (icepack); Burgard et al. 2022;
Hahn, Mikula and Frolkovic 2025; Smith et al. 2020.

## Model Characteristic Table

| Characteristic | Main suite of experiments | PPE change? |
|---|---|---|
| Mesh discretisation | Delaunay triangulation (gmsh), adaptive size field | no |
| Native grid | H: anisotropic; resolution **[confirm #20]**, pending the 1000 m inversions: the 2 km / 180 km adaptive mesh, 2 km at the grounding line and calving front to 180 km in the interior (246,677 cells), as previously run; or `antarctica_10000_1000_buffered20000`, the 1000 m / 10 km gmsh mesh (1,869,088 vertices) that has been the code default since PR #7 and on which no inversion has yet been run. V: vertically integrated (shallow shelf) | no |
| Native projection | EPSG:3031, same as BedMachine | no |
| Interpolation to diagnostic grid | conservative: exact cell-pixel overlap areas (supermesh) onto the 8 km grid; whole-pixel means for thickness, fractions and every flux, so a flux times the pixel area sums to the model's integral (`acabf` is fill outside the model domain, `libmassbffl` where no ice floats at year end, `lifmassbf` nowhere); covered-part means for elevations | no |
| Time integration | transport-first split: implicit Euler thickness transport, then the diagnostic solve at the new geometry; first order | no |
| Time step | **[confirm #20]**, pending the 1000 m inversions: 0.1 yr on the adaptive mesh, as previously run; 0.05 yr on the 1000 m / 10 km mesh, the code default since PR #7 | no |
| Advection scheme | upwind finite volume, DG0, implicit; first order | no |
| Ice flow mechanics | shallow-shelf approximation, dual finite-element formulation (CG1 velocity, DG0 membrane and basal stress) | no |
| Ice rheology | n = 3 (composite with a linear floor for thin ice) | no |
| Basal sliding | regularised Coulomb, m = 3, c0 = 0.5 (Budd available) | **[confirm #24]** |
| Basal hydrology | none | no |
| Ice-shelf fracture | **[confirm #10]** none, or the ISMIP7 collapse mask on floating cells, everywhere (`mask`) or from the front (`mask_front`); see item 9 | **[confirm]** |
| Advance and retreat | grounding line free; calving front pinned at 2015 (level-set von Mises optional) | **[confirm #36]** |
| Grounding line | flotation criterion per cell | no |
| Calving | pinned front; ice past it removed | **[confirm #36]** |
| Initial SMB | RACMO2.4p1 2000-2023 climatology | no |
| Bedrock adjustment | no | no |
| Year of initial condition | 2015 | no |
| Densities, gravity | rho_i = 917, rho_o = 1024, fresh water 1000 kg m-3, also in `params.nc` beside `CORE/`; g = 9.81 m s-2 | no |
| Variables not included | none of the mandatory set; no 3D or thermal variables (no thermal model); `hfgeoubed`, `litemp*`, `zvel*`, `thdrflf`, `deltag`, `refgeoid` absent | no |
| Days per year | 365.25: the model's year is icepack's, 31557600 s, and every model-to-SI conversion in the submitted files uses it. The forcing is converted on the way in with the tropical year, 31556926 s, a relative difference of 2e-5. The model counts time in years and has no calendar; the time axis in the files is the standard calendar (discussion #24), state at 1 January of the following year and fluxes at 1 July with bounds | no |
| Other | apparent-mass-balance correction frozen at the initial state; forcing versions cited per file in the submission | no |

## Conventions in the submitted files

Signs (discussions #16 and #22). `acabf`, `libmassbffl`, `lifmassbf` and
`licalvf` are positive for mass gained by the ice, so melt and calving are
negative.
`ligroundf` is the upwind flux across the grounding line from the velocity
the transport used, through the facets between a grounded and a floating cell
on the native mesh, divided by the area of the first FLOATING cell and
remapped conservatively. Its reference is the grounded ice sheet: positive
for grounded ice going afloat, negative where floating ice flows onto
grounded ice, at a pinning point or an ice rumple. The organisers' reply of
21 September 2026 on discussion #22 leaves the reference to each group, and
isschecker 0.5.1 bounds the field symmetrically; the group settled on this
reading on 22 September 2026. The integrated scalars carry the signs of the
fields they integrate, so `tendlicalvf`, `tendlibmassbffl` and
`tendlifmassbf` are negative and `tendligroundf` is the net grounding-line
discharge. They integrate over
true area: each native cell counts its map-plane area times af2 = (1/k)^2,
the EPSG:3031 area factor at the cell's centroid, the factor
`ismip7-scalar-processing` weights every 8 km pixel by. `topg` is not
masked to the ice, `lithk` is zero and not fill where there is no ice, and
the fill value is the finite netCDF default (discussions #10 and #19).
The fluxes are the ones the model applied: where the thickness floor held
back part of a cell's net sink, that part comes off the SMB, the melt and the
reference in proportion, so `libmassbffl` never reports melt of ice the cell
did not have. A floating cell whose base lies within 1 cm of the bed, the
checker's elevation tolerance, is written as grounded.
Sea-level estimates
(`sla20`, `slg20`, `slvaf`) are not computed by the model; **[confirm #13]** that
`ismip7-scalar-processing` was run on the gridded files.

Front melt, chosen by the group on 25 September 2026 (issue #109). The
request fills `libmassbffl` wherever no ice floats at year end, so the melt
of ice gone within the year would reach no gridded field. Most of it is
grounded ice that goes afloat into a marine cell whose shelf has gone, and
melts on arrival. The model reports that melt as `lifmassbf`, which the
request never fills: in a cell with its bed below sea level and no ice (1 m
or less) at either end of the year, the share of the year's melt that the
inflow supplied, out of what the inflow, the positive SMB and the frozen
apparent-mass-balance reference supplied together. The reference's share
stays in `libmassbffl` (issue #105), and the two fields sum to the melt the
model applied. Every yearly file names this booking, and the writer refuses
a series that mixes it with the earlier one, so a chain must not change code
versions across it. Where the melt lands depends on the time step: the
thickness floor lets a cell melt only the ice it held when a step began, so
a receiving cell ends the year holding the last step's inflow, and where
that exceeds 1 m the cell counts as ice and its melt stays in `libmassbffl`,
where the fill keeps it. In the 32 km CESM2-WACCM ssp585 at 2300 (map-plane
area, `ISMIP7_DT=0.1`), `lifmassbf` carries 603 Gt/yr, and the gridded
`libmassbffl` still leaves out 1,200 Gt/yr of the native melt: 1,094 is the
reference's share on emptied marine cells, 20 is the reference's supply
melted on emptied land cells, which the melt law counts as afloat at zero
thickness (issue #105), and 87 is shelf ice that melted away within the
year. **[confirm #109]** the production ssp585's numbers, which
the writer prints.

Compliance: isschecker 0.5.1 of 22 September 2026, which grades a range
finding by the share of values outside the bounds (discussion #46). A 32 km
CESM2-WACCM control and ssp585, 2015 to 2300, pass it with zero errors in
every test group (23 September 2026), and each submitted file set is checked
with it before upload. The bundled variable request is that release's.

Forcing versions: each run logs the product and version its readers opened,
and `core_report.py` carries those lines into the run's committed report;
**[confirm #41]** by citing them here per experiment. As audited against the
data-freeze mirror on 13 September 2026: CESM2-WACCM atmosphere SDBN1-8000m
v2, ocean v3, fracture v2.1; MRI-ESM2-0 atmosphere GEMB-SDBN1-8000m v1,
ocean v3, fracture v1; ISMIP7 ocean climatology 30_sep; observations kit
AntarcticaObsISMIP7-v1.2. Since that audit the MRI-ESM2-0 ssp585 fracture
has moved to v2 (22 September 2026, discussion #30), which the reader opens;
ssp126 and ssp370 stay at v1.
