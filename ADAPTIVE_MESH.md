# Adaptive remeshing

Adaptive remeshing after Gudmundsson et al. (2012), ported to this model in
September 2026: the desired-element-size field, `Error2EleSize`, the global
remeshing, and the field transfer between meshes.

The submission runs on the fixed gmsh pair `antarctica_10000_1000_buffered20000`
at dt 0.025 yr, chosen over the adaptive-preset mesh on 25 September 2026
(issue 20). The adaptive preset below remains available and is outside the
production path.

## Build the mesh first, do not refine mid-run

**Mid-run refinement blows up.** Refining an evolving 32 km state into 8 km
bands part-way through a run fails with all three DG0 transfer rules below in
force: the projected DG0 staircase reads as driving stress and the thickness
clamp goes from 119,000 to 256,000 Gt/yr within the first steps. This follows
from the DG0 geometry, so raising the resolution of a running model is outside
what this branch does.

**Intended use.** Build the mesh, invert on it, run forward on it.

1. Build an adaptive-preset mesh. Two sidecars are committed, both under
   `ISMIP7_LC=2000 ISMIP7_LC_COARSE=180000 ISMIP7_BUFFER_M=20000` with the adaptive
   preset, differing in where the desired sizes come from.

   Sized from the MODEL fields of the 2500 m MAP:

   ```
   adapt_mesh.py --mesh-only <2500 m MAP checkpoint> \
       --out-mesh antarctica/mesh/antarctica_ua_180000_2000.msh
   ```

   Sized from OBSERVATIONS, MEaSUReS velocity and BedMachine, where the
   scaffold supplies only the points the size field is evaluated on:

   ```
   adapt_mesh.py --mesh-only --from-obs --source-mesh <scaffold .msh> \
       --out-mesh antarctica/mesh/antarctica_ua_180000_2000_obs.msh
   ```

   `--out-mesh` fixes the name. Without it the output is
   `<reference>_adapt1.msh`, or `<reference>_<experiment>_adapt1.msh` when
   `ISMIP7_EXPERIMENT_NAME` is set, which `run_adaptive.py` does from
   `--experiment-name` so parallel experiments cannot overwrite each other in
   the shared `antarctica/mesh/`. The `.msh` files are regenerated with these
   flags; the sidecars are committed. The `_obs` mesh is what the NOTS
   inversions and the committed MAP names refer to. The new mesh's outline
   buffer is the checkpoint's `buffer_m`, else its mesh name's
   `_buffered<N>` tag, else a named `ISMIP7_BUFFER_M`, and a mesh with none
   of them is refused, since 0 and 20000 are each wrong for some legacy mesh.
   A scaffold (`--source-mesh`) takes a named `ISMIP7_BUFFER_M` first.
2. Invert on that mesh.
3. Run the forward on it, without adaptation.

**What is validated.** The identity transfer only. `adapt_mesh.py --no-remesh`
moves the state onto a fresh load of the same mesh and reproduces the run
(outflux 672 to 797 against 782 Gt/yr unadapted, VAF 57637.4 against 57637.3).
That exercises the transfer rules. `run_adaptive.py` and the segment loop are
wired and correct as far as the identity test reaches; treat a real mid-run
adaptation as unvalidated. Closed as icepack/ismip7#37, not planned for
September 2026.

Grounding-line morphing was not ported: the reference scheme carries a
mesh-deformation step with no caller, a hook commented out in the global
remeshing as "broken anyhow", and no mention in the reference documentation.
Its global remeshing is what this implements. Morphing would
be about 100 lines on top of `icepack2_tools.adapt_mesh.grounding_line_points`.

## The scheme (the reference implementation's names in brackets)

1. **Desired element size** at the nodes [`EleSizeDesired`]. Start at
   `MeshSizeMax`. For each enabled criterion
   [`ExplicitMeshRefinementCriteria`] compute a nodal error proxy `e` and map
   it [`Error2EleSize`]: `h = hMin + (e0/(e+e0))^(1/p) (hMax - hMin)` with `e0`
   the criterion's `Scale`. Take the minimum over criteria, falling back to
   `MeshSize` when none fired. Relax toward the current size with `W = 0.5` and
   clip the ratio of change to `[1/5, 5]`. Then apply the absolute bands: nodes
   within `d_i` of the grounding line get `min(h, s_i)` per row of
   `MeshAdapt.GLrange`, floored at `MeshSizeMin`, and likewise `CFrange` at the
   calving front.
2. **Global remeshing** [`explicit:global`, gmsh]: the nodal size field becomes
   gmsh's background scalar view, the outline and physical groups are rebuilt
   as `mesh_antarctica.py` builds them, and the mesh is regenerated. Then the scheme's
   element-count control: up to four rescalings of `MeshSizeMin` so the count
   lands within `[LowerLimitFactor, UpperLimitFactor] * MaxNumberOfElements`.
3. **Transfer** [`MapFbetweenMeshes`]: point-evaluation interpolation old to
   new (the reference uses FE shape functions), thickness `ThickMin` outside the old mesh,
   bed re-sampled from BedMachine with the method the MAP recorded, surface by
   flotation. Frozen anchors (`a_ref_mb`, `N_ref`, `C_w0`, `H_init`, level set)
   move too.

Criteria available, the reference list: `effective strain rates`, `effective
strain rates gradient`, `flotation` (the `DiracDelta`
`0.5 k sech^2(k (h - h_f))`), `thickness gradient`, `upper surface gradient`,
`lower surface gradient`, `|dhdt|`, `dhdt gradient`.

Two choices where the reference is nodal and this model is DG0:

- `EleSizeCurrent` is `sqrt(mean adjacent element area)`, kept despite
  being 0.66 of an edge length, because the relaxation and ratio limits are
  calibrated against it.
- Thickness transfer follows the reference's nodal surface route
  (`ISMIP7_ADAPT_GEOMETRY=bh-FROM-sBS`, the default). Under it the thickness is
  never transferred: `surface_route_thickness` moves the surface and rebuilds
  `h` against the bed re-sampled on the new mesh, so `ISMIP7_ADAPT_TRANSFER`
  cannot change the volume it lands on. The `bs-FROM-hBS` route transfers `h`
  directly and put 32 km cell-mean thicknesses onto a finely re-sampled bed,
  floating them over every trough (outflux 782 to 14,770 Gt/yr at the first
  step).

`transfer_state` consults `ISMIP7_ADAPT_TRANSFER` for three field names under
DG0 geometry: `a_ref_mb`, `thickness_dg` and `H_init`. Two are unreachable in
the shipped workflow. `a_ref_mb` is discarded whenever the checkpoint carries a
`velocity`, which every forward state checkpoint does, because the physical
divergence replaces it. `thickness_dg` is written only under CG1, where the DG0
condition fails. So the setting changes `H_init` alone, and on an initial
adaptation (`--rebuild-aref`) not even that, since `H_init` comes from the
re-sampled thickness. Everything else, the surface, the transferred physical
divergence and the `N_eff/N_ref` ratio, is interpolated unconditionally.

`project` must run on one rank: `cross_mesh_transfer` raises when the new
mesh's communicator has more than one. `run_adaptive.py` therefore takes a
separate `--adapt-launcher`, defaulting to `mpiexec -n 1`, which is safe under
either setting since the remesh is serial gmsh on rank 0 anyway. Raise it to
the forward's rank count to parallelise the transfer under the default
`interpolate`. Either way the transfer prints the volume change and the mean
front thickness before and after.

## Configuration (`ISMIP7_ADAPT_*`, defaults from the reference implementation)

| variable | reference parameter | default |
|---|---|---|
| `ISMIP7_ADAPT_MESH_SIZE` | `MeshSize` | 10 km |
| `ISMIP7_ADAPT_MESH_SIZE_MIN` / `_MAX` | `MeshSizeMin` / `MeshSizeMax` | 1 km / 10 km |
| `ISMIP7_ADAPT_GL_RANGE="5000:2000,1000:500"` | `MeshAdapt.GLrange` | none |
| `ISMIP7_ADAPT_CF_RANGE` | `MeshAdapt.CFrange` | none |
| `ISMIP7_ADAPT_CRITERIA="effective strain rates:0.01,..."` | `ExplicitMeshRefinementCriteria(I).Name/Scale[/p/EleMin/EleMax]` | none |
| `ISMIP7_ADAPT_RELAXATION_W` | `W` | 0.5 |
| `ISMIP7_ADAPT_MAX_RATIO_CHANGE` / `_MIN_` | `Max/MinRatioOfChangeInEleSizeDuringAdaptMeshing` | 5 / 0.2 |
| `ISMIP7_ADAPT_MAX_ELEMENTS` | `MaxNumberOfElements` (0 = off) | 0 |
| `ISMIP7_ADAPT_THICK_MIN` | `ThickMin` | 1 m |
| `ISMIP7_ADAPT_DIRAC_WIDTH` | `RefineDiracDeltaWidth` | 100 m |
| `ISMIP7_ADAPT_TRANSFER` | none | `interpolate` |
| `ISMIP7_ADAPT_GEOMETRY` | `MapOldToNew.Transient.Geometry` | `bh-FROM-sBS` |
| `ISMIP7_ADAPT_KEEP_CURRENT=1` | none (null test at current sizes) | off |

### The Antarctic preset: `ISMIP7_ADAPT_PRESET=ua`

The sizes the reference setup uses for Antarctica, so the adaptation is consistent with
that setup rather than with this repo's initial meshes.

| setting | value | source |
|---|---|---|
| `MeshSizeMax` (interior) | 180 km | the coupled ice-ocean pan-Antarctic mesh, GMD 18 (2025): "up to 180 km in the interior" |
| `MeshSizeMin` (grounding line) | 2 km | same paper: "adaptive refinement down to 2 km at the grounding line" |
| `MeshSize` (fallback) | 90 km | PIG-TWG example: `MeshSize = MeshSizeMax/2` |
| ice-shelf size | 10 km | inferred. The reference PIG-TWG setup uses `MeshSizeIceShelves = MeshSizeMax/5` and its pan-Antarctic runs quote 4 km, which continent-wide is about 220,000 elements on the shelves alone. 10 km keeps the 1.5 M-dof budget |
| low ground (`s < 1500 m`) | 36 km | PIG-TWG `EleSizeIndicator(s<1500) = MeshSizeMax/5`, scaled to 180 km |
| `effective strain rates` | Scale 0.001, floor 4 km | PIG-TWG Scale; pan-Antarctic "4 km in regions of high strain rate" |
| `GLrange` | 10 km: 4 km, 5 km: 2 km | inferred pan-Antarctic form of MISMIP+'s `[20000 5000; 10000 2000; 5000 500]` and the 2 km statement |
| `MaxNumberOfElements` | 250,000 | the pan-Antarctic setup: "250 000 elements" (linear, `TriNodes=3`, as here) |
| initial adaptation | up to 5 iterations | PIG-TWG `AdaptMeshMaxIterations=5` |
| remesh interval | every step in the reference | here `--adapt-every` in years, 1 yr being the practical floor |

The preset fills only what is unset, so any single `ISMIP7_ADAPT_*` variable
overrides it. Elements in these setups are linear, so a 4 km element
and a 4 km cell here resolve alike.

Timing knobs live in `run_adaptive.py`: `--adapt-every`
[`AdaptMeshTimeInterval`], `--initial-iterations` [`AdaptMeshInitial` and
`AdaptMeshMaxIterations`].

## Pieces

- `icepack_tools/adapt_mesh.py`, the shared package: the scheme itself
  (criteria, `Error2EleSize`, relaxation and ratio limits, bands, PIG-TWG
  rules, `remesh_global` with a project-supplied `build_geometry()` and the reference's
  element-count control) and the transfer helpers (`cross_mesh_transfer`,
  `preserve_front`, `physical_divergence`, `surface_route_thickness`,
  `rebuild_reference_pressure`). Nothing in it knows about Antarctica. Tests:
  `icepack_tools/test/adapt_mesh_test.py`. Import it before assembling any
  form, since it reaches `icepack2` through `icepack_tools.constants` and
  Irksome refuses to load afterwards.
- `icepack2_tools/adapt_mesh.py`: the Antarctic domain builder and
  `transfer_state`, which knows this model's checkpoint fields, frozen
  references and restart attributes.
- `antarctica/scripts/adapt_mesh.py CHK --out-checkpoint NEW [--rebuild-aref]`:
  one adaptation of a forward checkpoint, writing the `.msh`, its sidecar and a
  restartable checkpoint.
- `antarctica/scripts/run_adaptive.py`: the segment loop.
- `simulation.py`: a restart carrying `adapted_initial=1` and no `a_ref_mb`
  rebuilds the apparent-mass-balance reference on the new mesh.

## Geometry route

The reference default `MapOldToNew.Transient.Geometry = "bh-FROM-sBS"` moves the surface
and derives the thickness from it and the re-sampled bed
(`h = min(s - b, s rho_w / (rho_w - rho_i))`). The DG0 analogue of the nodal
surface is the volume-preserving CG1 lift of the cell surface, point-evaluated
at the new centroids. On the first run-step of a transient run the reference takes all
geometry from data, and `--rebuild-aref` does the same (bed and thickness from
BedMachine, `H_init = H`).

## Three things a nodal model never had to face

**The frozen mass-balance reference does not transfer.** `a_ref_mb` cancels the
old mesh's discrete flux divergence spike by spike (about 1000 m/yr at the Pine
Island grounding zone). Interpolated onto another mesh those spikes become
misplaced sources, and a same-resolution remesh blew up within a year (outflux
4,371 to 53,185 Gt/yr). The transfer carries the physical divergence instead,
`P = flux/area - a_ref`, saved as `phys_div`, and the forward rebuilds
`a_ref = flux_new/area - P` with its own upwind operator, exactly as it built
the t=0 reference.

**Budd's effective-pressure reference is frozen to the original geometry.**
`N_ref` is the t=0 effective pressure and the law reads `N_hat = N_eff/N_ref`
(cap 3, zero afloat). Re-sampling the bed changes `N_eff` cell by cell while an
interpolated `N_ref` does not follow. On a same-resolution null remesh the
grounded-cell ratio spread went from 0.91 to 1.48 (5th to 95th percentile) out
to 0.36 to 2.97, with 10% of cells losing their reference altogether, which the
law reads as triple friction; on the 8 km-band case 26% did. The transfer
carries the ratio the run had and rebuilds `N_ref = N_eff_new / ratio`.
Regularized Coulomb has no such reference, so RC forwards are immune.

**The DG0 front thickness smears under any interpolation.** The last cell at
the front is thin, a new front cell whose centroid lands in an old interior
cell inherits a thick value, and terminus traction goes as h^2 (the
`geometry.py` docstring measured a 4.7x outflux multiplier from the same effect
in the CG1 lift). The transfer prints the mean front thickness before and
after. Under the default geometry route the thickness is rebuilt from the
interpolated surface, so `preserve_front` (`ISMIP7_ADAPT_FRONT_PRESERVE`, on by
default) is what pins the boundary cells, and `ISMIP7_ADAPT_KEEP_CURRENT=1`
remeshes at current sizes to measure the transfer's own cost.

## Found on the way: Budd shelf friction was a sign test on roundoff

The identity transfer reproduced every checkpoint field to machine precision,
yet two identity checkpoints whose `surface` differed by 2e-13 m gave 672 and
1459 Gt/yr for the same year. In `build_rc_residual(fric_law="budd")` the shelf
test was `conditional(gt(N, 0), N_hat, 0)`. On a floating cell the surface is
the flotation branch, so `N = max(p_I - p_W, 0)` is a roundoff residue of
either sign; a positive residue passed the test with `N_ref` equally tiny,
where the delta floor `nhat_floor * p_I / N_ref` lifts it to the cap. That gave
triple Weertman friction on whichever shelf cells rounded positive (445 of 3515
in the 32 km control, and the 2e-13 change flipped all of them).

The gate is now HAF > 0 itself (`dual_friction.budd_nhat`), which for grounded
ice is the same statement as N > 0 (`N = rho_I g HAF` when `s = b + H`) and on
the shelf is a real negative number. An intermediate form gated on `He` alone
left 133 of the 3791 floating cells of the 32 km MAP receiving `He * nhat_cap`,
since `He` is smooth in height above flotation. Under the fixed law the two
identity checkpoints give outflux 1478 against 1475 Gt/yr and VAF 57638.204
against 57638.212 mm SLE. Census script:
`antarctica/scripts/check_budd_map.py MAP [--forward]` (old gate 418 cells at
the cap, He-only 133, production 0). The census is what calls for re-inverting
Budd MAPs. Re-solved under the fixed law, the old MAP lands at rel L2 = 0.91
from its saved velocity. A re-solve distance of this size is now known to
appear under both friction laws and for a MAP inverted under the shipped gate
(on the 2 km / 180 km adaptive mesh, 0.665 for the re-inverted Budd MAP and 0.685 for the
regularized Coulomb MAP), so it does not by itself identify the shelf gate as
the cause. Those two came from the forward's `ocean_drag` and `u_lim`, which
the 14 September inversions did not assemble; section 4 of
`antarctica/FORWARD_RUN_READINESS.md` has the measurements.
The outflux numbers quoted above (782, and the identity test's 672) were
measured with the old test in place, as was every Budd forward and Budd MAP to
date. Regularized Coulomb has a continuous `tau_cap` and never had the problem.

**Provenance of the interim form.** Between b930055 and 0cb378e the gate
carried a multiplicative `He`, which scaled grounded friction inside the
GL_WIDTH band as well as zeroing the shelf. The shipped law is HAF only, with
`He` entering the Budd branch through `tau_W`'s `exp(theta * He)` as it does
under regularized Coulomb. The re-inversions must match that form, since the
inversion and the forward assemble the same expression. The adaptive-mesh Budd MAP
inverted under the interim form is kept as
`inversion_icepack2_budd_n3_dg0_logvelnet_ua2000_He.h5`; the production
re-inversion (NOTS job 1339328) warm-starts from it under the shipped law. No
MAP inverted with the He form is to be used for a result, since a forward
records the law name alone and nothing detects the mismatch at load time.

## Level set across a remesh

`transfer_state` moves the `levelset` field as the shared package documents it:
as its P1 lift, landing back in DG0, with water (+1e6 m) outside the old mesh.
The forward's extent anchor re-solves the eikonal problem on the new mesh, the
same order the calving project uses after its remeshes. No run reads the
transferred field back yet, since a calving forward builds its level set from
the current thickness and `H_init`, so the front's sub-cell position does not
survive a remesh today. The transfer keeps the checkpoint correct for when a
forward does load it.
