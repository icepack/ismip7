# ISMIP7 Antarctica

Antarctic ISMIP7 submission on icepack2 and Firedrake. Companion to the
Greenland repository (https://github.com/dlilien/ISMIP7_Greenland_Icepack).

The pipeline is a dependency chain: data, mesh, inversion, melt calibration,
control and projections.

---

## 0. What you need, by what you want to run

Install and download the rows your column ticks. Sizes are measured.

### 0.1 Software

| | Where from | Invert | Forward run | Adapt the mesh | Submit |
|---|---|:--:|:--:|:--:|:--:|
| **Firedrake 2026.4** (brings PETSc, MUMPS, mpi4py) | firedrakeproject.org | x | x | x | x |
| **icepack2** | github.com/icepack/icepack2 | x | x | x | x |
| **icepack** (raster interpolation onto meshes) | github.com/icepack/icepack | x | x | x | x |
| **icepack_tools** (`adapt_mesh`, `levelset`, `calving`, `friction`, `grounding`) | github.com/hoffmaao/icepack_tools | | level-set front and calving laws | x | |
| **tlm_adjoint** | github.com/jrmaddison/tlm_adjoint | x | | | |
| `xarray netCDF4 scipy rasterio pyproj shapely gmsh matplotlib` (`geopandas` only to build a mesh, section 3) | pip, into the Firedrake venv | x | x | x | x |
| `earthaccess` (NSIDC), `globus-sdk` (Globus route only) | pip | x | x | x | |
| **isschecker** (`ismip7-compliance-checker`) | github.com/ismip/ISM_SimulationChecker | | | | x |

`icepack_tools` is a separate repository. Install it editable into the
same venv: `pip install -e /path/to/icepack_tools`.
`icepack2_tools/adapt_mesh.py` and `icepack2_tools/levelset.py` wrap it, and a
calving law (`ISMIP7_CALVING`) comes from its `calving` module; the rest of the
repository runs without it.

The compliance checker needs Python 3.11 or newer. Check `python -V`; a
Firedrake venv often carries an older one, in which case give the checker its
own venv and call it by absolute path.

### 0.2 Data

| | Size | How | Invert | Forward run | Adapt | Submit |
|---|---|---|:--:|:--:|:--:|:--:|
| BedMachine Antarctica v4.1, MEaSUREs velocity v2 | 8 GB | `scripts/download_data.py` (Earthdata login) | x | x | x | |
| RACMO2.4p1 SMB climatology | 2 GB | same script | | x | | |
| ISMIP7 observations MIPkit v1.2 (Smith dH/dt) | 9 GB | `scripts/download_mirror.py --product ismip7-ais-observations data/mipkit/`, landing at `ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc` (`ISMIP7_OBS_KIT` overrides). `scripts/download_forcing.py --calibration` stages the same v1.2 file in the same place over Globus | `ISMIP7_DHDT_WEIGHT` | a start before 2015 (historicals, OCX) | `--from-obs` | |
| ISMIP7 forcing per ESM and scenario: SMB anomaly 7.5 GB, SMB gradient `dacabfdz` 0.2 GB, ocean `tf` 11 GB, `so` 6.9 GB (ssp585; historical 4.3 GB) | 25 GB each | `scripts/download_mirror.py` | | x | | |
| ISMIP7 fracture (collapse mask, lake properties, excess melt) | 3 GB per scenario | same, `data/<ESM>/<scenario>/fracture/` | | `ISMIP7_FRACTURE=mask` | | |
| Ocean OI climatology and IMBIE basin numbers | 3 GB | `scripts/download_forcing.py --ocean --calibration` | | x | | |
| Whole AIS tree (all ESMs, scenarios, `ctrl`, OCX, calibration) | 313 GB | same | | | | |
| Meshes and MAP checkpoints | 15 MB, 80 MB | sections 3 and 4, or from a colleague | | x | x | |
| Melt calibration `calibration/deltaT_per_basin_1000_K6.500e-05_vertex_front.npz` (`vertex` sampling: `deltaT_per_basin_1000_K6.500e-05.npz`) | 8 kB | tracked in git (section 5) | | x | | |

Source Cooperative carries the data-freeze copy and needs no account.

```bash
# one scenario for one ESM (about 25 GB); dacabfdz is the SMB gradient the
# SMB-elevation feedback reads (section 6)
python antarctica/scripts/download_mirror.py \
    data/CESM2-WACCM/ssp585/SDBN1-8000m/acabf-anomaly/ \
    data/CESM2-WACCM/ssp585/SDBN1-8000m/dacabfdz/ \
    data/CESM2-WACCM/ssp585/ocean/tf/ data/CESM2-WACCM/ssp585/ocean/so/
# the control's ocean and SMB gradient for one ESM, which cores 9 and 10 read
# (about 18 GB)
python antarctica/scripts/download_mirror.py \
    data/CESM2-WACCM/ctrl/SDBN1-8000m/dacabfdz/ \
    data/CESM2-WACCM/ctrl/ocean/tf/ data/CESM2-WACCM/ctrl/ocean/so/
# the observations MIPkit (about 9 GB)
python antarctica/scripts/download_mirror.py --product ismip7-ais-observations data/mipkit/
# whether a local tree is current: which version a run opens
python antarctica/scripts/audit_forcing_versions.py --scenario ssp585
# whether every file is there and unchanged: read-only, exit 1 if anything is
# left to fetch or was replaced on the mirror
python antarctica/scripts/download_mirror.py --check data/CESM2-WACCM/ssp585/
```

Globus remains the archive of record and `download_forcing.py` drives it
(section 2a). Use it for anything the mirror has not synced. The MRI-ESM2-0
ssp585 collapse mask is such a case: v2 replaced v1 on Globus on 22 September
2026. `FRACTURE_MIN_VERSION` in `icepack2_tools/forcing.py` names the oldest
mask each ESM and scenario may be read at; a run under a mask mode refuses an
older one at startup, and the audit reads it `OUTDATED` and fails. Fetch it
with `download_forcing.py --scenarios --esm MRI-ESM2-0 --scenario ssp585`.

### 0.3 Short paths

**Forward run from someone else's MAP.** Firedrake, icepack2, icepack,
BedMachine, MEaSUReS, RACMO, one scenario of forcing, their `.msh` and
`inversion_*.h5`, the OI climatology and the IMBIE basin numbers. The melt
calibration is tracked (`calibration/`, section 5), so a fresh clone melts
with it and nothing is calibrated or copied.

**Inversion.** Add `tlm_adjoint`, and the MIPkit for the dH/dt term.

**Mesh adaptation.** Add `icepack_tools`. The observation-driven size field
(`adapt_mesh.py --from-obs`) needs the MIPkit and MEaSUReS.

**Submission.** Add the checker in its own Python 3.11+ venv, run with
`ISMIP7_OUTPUT=1`, then `scripts/write_ismip7_output.py`. The two halves can
sit on different machines. The writer reads the annual checkpoints, which for a
286-year experiment come to roughly 14 GB, so it belongs where the run is; it
needs the forward's own environment plus netCDF4, since it reaches icepack2
through `icepack2_tools.dual_friction` and dates its time axis with cftime.
The checker reads the written NetCDF, a few GB after zlib, so it can stay
wherever the newer Python is. Rice NOTS carries the writer's dependencies; the
checker runs here on the workstation.

### 0.4 Accounts

| For | Account | Where |
|-----|---------|-------|
| BedMachine, MEaSUREs (NSIDC) | NASA Earthdata (free) | https://urs.earthdata.nasa.gov/users/new |
| ISMIP7 forcing over Globus (the mirror needs none) | Globus + the ISMIP7 collection | https://app.globus.org |
| Submitting results | an upload folder from the ISMIP7 team | email ismip6 at gmail.com with your Globus id, `AIS`, group name and `ism_id` |

Globus transfers to this machine also need Globus Connect Personal running
locally and its endpoint UUID.

### 0.5 Running on a cluster

`antarctica/scripts/batch_runners/` is site neutral. One file per cluster in
`sites/` holds the venv path (or container image), module loads, partitions,
account, per-node limits and paths, and `submit.sh` composes the scheduler
command from it. What is yours alone on a cluster (account, a private build,
mail flags) goes in the git-ignored `sites/local.env`. Rice NOTS and IU Quartz
ship filled in, UChicago Midway is a container site waiting for its partitions
and account, `sites/template.sh` is the blank.

```bash
antarctica/scripts/batch_runners/submit.sh inversion ISMIP7_LC=2000 \
    ISMIP7_LC_COARSE=5000 \
    ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_5000_2000_buffered0.msh
ISMIP7_SITE=iu_quartz antarctica/scripts/batch_runners/submit.sh projection \
    ISMIP7_EXPERIMENT=ssp585_cesm_waccm --dry-run
```

The timing benchmark (section 7) goes through the same layer: `make timing`
and `manage_timing_campaign.py` submit `batch_runners/timing_*.script` with
`submit.sh script`, so they are the same commands at every site too.

Details, and the six steps for a first day on a cluster, in
`antarctica/scripts/batch_runners/readme.md`.

### 0.6 Repo layout

```
antarctica/
  data/            # BedMachine, velocity, RACMO                        [gitignored]
  mesh/            # *.msh, boundary_ids_antarctica_*.json, inversion_*.h5
  results/         # checkpoints, timeseries, logs                      [gitignored]
  reports/         # tracked per-core run records + MATRIX_STATUS.md
  scripts/         # entry points (sections 3 to 6)
    batch_runners/ # scheduler job scripts + sites/<cluster>.sh
ISMIP7/AIS/        # forcing tree the runtime reads                     [gitignored]
icepack2_tools/    # this repository's library
```

`FORWARD_RUN_READINESS.md` in this directory carries the data freeze, the
forcing-version audit, the OCX and control definitions, and the submission
checklist; read it before planning the core matrix.

The tracked files under `antarctica/mesh/` are the per-mesh
`boundary_ids_antarctica_*.json` sidecars (section 3). Legacy sidecars with no
mesh stem stay untracked: `boundary_ids.json` is the shared fallback that every
mesh build overwrites, so committing one would put a single build's map on the
fallback path for every clone.

---

## 1. Observational data (`download_data.py`)

```bash
cd antarctica
python scripts/download_data.py
```

| Dataset | Product | Auth | Lands in |
|---------|---------|------|----------|
| BedMachine Antarctica v4 | NSIDC-0756 | Earthdata | `<obs root>/bedmachine/` |
| MEaSUREs Ice Velocity v2 | NSIDC-0484 | Earthdata | `<obs root>/velocity/` |
| RACMO2.4p1 SMB | Zenodo `10.5281/zenodo.14217231` | none | `<obs root>/racmo/` |

The obs root is `ISMIP7_OBS_DATA_ROOT` (section 6, environment knobs), the
same directory every run reads, and it is created if missing. Existing files
are skipped, so re-running is cheap.

---

## 2. ISMIP7 forcing

Three things:

1. **The runtime tree `ISMIP7/AIS/`**, read by `icepack2_tools/forcing.py` and
   laid out per the protocol (section 2b). `ISMIP7_DATA_ROOT` overrides the
   root.
2. **`download_mirror.py`**, the Source Cooperative route: anonymous HTTPS,
   resumable, no Globus endpoint, and the only option on a machine without one.
   It takes product-relative prefixes and lands files under `--root` in the
   versioned layout `forcing.py` expects. Its module docstring documents the
   prefixes, the resume rule and the size check.
3. **`download_forcing.py`**, the Globus route, mirroring the collection
   straight into the runtime layout: climatology, bias and calibration files
   (`--ocean`, `--calibration`), and the per-(ESM, scenario) sets
   (`--scenarios`). `python scripts/preflight.py` then reports which core
   experiments the local tree can run.

### 2a. Using `download_forcing.py`

```bash
cd antarctica
# one-time: authenticate through the Globus CLI. It prints a URL to open in a
# browser and takes the resulting code back, so it works from a headless
# cluster login node; tokens are cached under ~/.globus/cli/.
python -m pip install --user globus-cli   # if `globus` is not installed
python scripts/download_forcing.py --login
python scripts/download_forcing.py --list         # legacy climatology subtree; /ISMIP7/AIS needs the Globus web app
export GLOBUS_LOCAL_ENDPOINT=<your-endpoint-uuid>
python scripts/download_forcing.py --ocean        # thetao/so/tf + climatology + bias
python scripts/download_forcing.py --calibration  # meltMIP obs melt, IMBIE2 basins, grid, topography
python scripts/download_forcing.py --scenarios    # per-(ESM, scenario) forcing (cores 1-8)
python scripts/download_forcing.py --scenarios --esm MRI-ESM2-0 --scenario historical,ssp585
python scripts/download_forcing.py --scenarios --scenario ctrl   # the control's ocean (cores 9, 10); its atmosphere comes too, unread
python scripts/download_forcing.py --scalar-processing  # ismip7-scalars grids, see "From a finished run to a submission"
python scripts/download_forcing.py --status
```

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_GLOBUS_COLLECTION` | source collection UUID | `ccc9bbd2-4091-4e35-addd-eeb639cf5332` |
| `GLOBUS_LOCAL_ENDPOINT` | your Globus Connect Personal endpoint UUID | required to transfer |

Scenario forcing lives in the collection's top-level `/ISMIP7/AIS/<ESM>/<scenario>/`
tree. `--scenarios` mirrors the minimal runtime sets (SDBN1-8000m `acabf`,
`acabf-anomaly` and the SMB gradient `dacabfdz`, ocean `tf` and `so`, fracture)
with version autodetection and checksum sync, so re-runs are completeness
checks. `--scenario ctrl` fetches the control's set, its `dacabfdz` and ocean.
Climatology, obs and calibration sets come from
`/ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/AIS`; the `OCEAN_FILES` and
`CALIBRATION_FILES` dicts in the script are the manifest. Without
`GLOBUS_LOCAL_ENDPOINT` the script prints the paths for a manual transfer in
the web app.

### 2b. The runtime tree

```
ISMIP7/AIS/
  <ESM>/<scenario>/<SDBN1|GEMB-SDBN1>-8000m/<var>/<version>/   # acabf, acabf-anomaly, ts, tas, pr
      <var>_AIS_<ESM>_<scenario>_<product>_<version>_<YEAR>.nc
  <ESM>/<scenario>/ocean/<tf|thetao|so>/<version>/
  <ESM>/<scenario>/fracture/[v*/]                     # collapse and lake masks
  meltMIP/OI_Climatology_ismip8km_60m_<tf|so|thetao>_extrap.nc
  parameterisations/ocean/imbie2/                     # basin numbers the melt offsets are stamped through
  parameterisations/{ocean,fracture}/
```

Readers pin `version=v2` for the atmosphere and `v3` for the ocean, falling
back to the highest `v<N>` present, dotted versions included, so the `v2.1`
fracture release and MRI-ESM2-0's `v1` resolve without code changes. The
atmosphere directory is whichever of `GEMB-SDBN1-8000m` and `SDBN1-8000m`
exists, so trees fetched before MRI's August 2026 rename still work. Where
both exist a run opens `GEMB-SDBN1-8000m`, the name a mirror re-sync writes,
and `audit_forcing_versions.py` audits that same directory. Fracture masks
resolve flat or versioned.

**One year bridges the end of a series.** A request for the single year after
the last one on disk reuses that year and logs it once per variable. Anything
further past the end, or a gap inside the series, raises. CESM2-WACCM's ocean
stops at 2299, so 2300 holds 2299, a duplicate the organisers accepted on 23
September 2026 (discussion #49). Its 2300 atmosphere files, withdrawn while
empty, are back as the 2290-2299 mean and are read as given; the bridge covers
them only in a tree fetched without them.

**Per-year atmosphere files hold 12 monthly slices.** The reader collapses them
to the annual mean, weighting months by length from `time_bnds`, falling back
to coordinate spacing and then to an unweighted mean with a warning. A
length-1 time axis passes through unchanged.

API: `ISMIP7Atmosphere(esm, scenario).get_smb(year, x, y, anomaly=...)`,
`ISMIP7Ocean(...).get_thermal_forcing(...)` and `.get_salinity(...)`,
`ISMIP7Fracture(...).get_collapse_mask(year, x, y)`, and
`make_forcing_callback(atm=, ocean=, fracture=, K=, K_per_basin_npz=)` which
bundles them into the callback `run_simulation` expects.

---

## 3. Build the mesh

Adaptive isotropic mesh with grounding-zone and calving-front
refinement, sized from BedMachine geometry and MEaSUREs strain rate.

```bash
cd antarctica
python scripts/mesh_antarctica.py --lc 1000 --lc-coarse 10000 --buffer-m 20000   # the production pair
ISMIP7_LC=1000 ISMIP7_LC_COARSE=10000 ISMIP7_BUFFER_M=20000 python scripts/mesh_antarctica.py
# dev mesh for inversion_icepack2.py, diagnostic_solve.py and run_eigendec.py:
python scripts/mesh_antarctica.py --lc 8000 --lc-coarse 80000 --buffer-m 20000
# mesh/antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.msh
# mesh/boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json
```

`--lc` and `--lc-coarse` are the fine (grounding line, calving front) and
coarse (interior) element sizes in metres. GL-band sizes, `shelf_size`,
`buffer_size` and the strain-rate floor scale with `lc/2500`; calving-front
decay lengths are floored at `lc` and `1.25*lc`.

`--buffer-m` pushes the outline into the ocean before meshing (default 20000,
`0` for none), letting icepack2 handle `h=0` at an interior calving front
instead of a `calving_terminus` BC. Both the `.msh` and its sidecar are named
after the exact `(COARSE, FINE, BUFFER_M)` triple, so builds never collide.

**Boundary ids.** gmsh physical groups alternate `Calving_0, Other_1, ...`,
auto-numbered from 1, so odd tags are calving and even tags are other.
`mesh_antarctica.py` writes that split to the per-mesh sidecar (named by
`scripts/mesh_naming.py`). To regenerate one without rebuilding the mesh:
`ISMIP7_BUFFER_M=<N> python scripts/make_boundary_ids.py`.

Solvers resolve the sidecar through `icepack2_tools/boundary.py`:
`ISMIP7_BNDIDS`, then the per-mesh `mesh/boundary_ids_<mesh stem>.json`, then
the shared `mesh/boundary_ids.json`. Readers hard-error on an unclassified
exterior marker or a missing id. A silent mismatch once left 95% of the ice
front without calving back-pressure. Runs print the covered front length in km
and percent.

A forward takes the sidecar name from the MAP or restart checkpoint, which
records `mesh_basename` and the `lc`, `lc_coarse` and `buffer_m` parameters, so
environment drift cannot swap a sidecar mid-trajectory.

---

## 4. Invert for basal and rheology fields (`inversion_icepack2.py`)

MAP estimate of bed friction `θ` and rheology `φ` from the diagnostic 3-field
(V x Σ x τ) system, regularized, with n=1 to 3 continuation, through
`tlm_adjoint`.

```bash
cd antarctica
ISMIP7_LC=2500 ISMIP7_LC_COARSE=64000 \
  mpiexec -n 12 python scripts/inversion_icepack2.py
# mesh/inversion_icepack2_<budd|rc>_n3_dg0_<LC>.h5
# name BOTH: the defaults are the 1000/10000 production pair, and an LC on its
# own would ask for a 2500/10000 mesh that nothing builds
```

The controls are log deviations from physical priors: `θ = log(C/C_w0)` on the
balance-friction anchor and `φ = log(A/A_prior)` on a thermomechanical fluidity
prior computed at setup and stored in the MAP. That prior reads the section 1
RACMO SMB and the section 2 `tas` climatology, so download the forcing first.
See `N3_FRAMEWORK.md`, and `ISMIP7_FLUIDITY_PRIOR=legacy` to skip it.

`ISMIP7_FRICTION` selects the law (`budd` or `regularized_coulomb`, tagged
`_budd` and `_rc` in the filename). The MAP name carries the law, the flow
exponent and the geometry space, built by `icepack2_tools/naming.py`, which the
forward, `preflight.py` and the gates all import. A MAP is valid only for its
own geometry space, since the inversion absorbs the calving-front treatment into
`θ` and `φ`. A forward that finds only a legacy untagged MAP loads it with a
loud warning. See
`../GEOMETRY_DISCRETIZATION.md`.

### Inversion solver

`ISMIP7_INVERSION_LINEAR_SOLVER` picks the linear solver of every solve
`tlm_adjoint` differentiates, and of the solve that publishes the MAP's state.
`full_mumps` factors the whole mixed Jacobian. `scpc_gamg`, the default since
3 October, is the transient's production solver: Slate eliminates `M` and `τ` cell by cell
and GAMG solves the condensed velocity system. `scpc_mumps` condenses the same
way and factors the condensed system, which makes it the exact reference for
`scpc_gamg`.

Which is faster depends on the friction law and on how an evaluation reaches
full n. On the Quartz measurements below `full_mumps` stayed the default
(issue #156), because the five-solve evaluations, then the default, were
where production inversions spent their time. The PR 155 sub-element
measurements before the closing SEP1 table are SEP2, the only scheme before
PR 158. PR 158 made SEP1 and the direct forward, one solve an evaluation, the
inversion's defaults. The older one-solve rows predate three changes that move
each side: `full_mumps` lost the SCPC zero blocks, so its seconds are upper
bounds; `scpc_gamg` keeps its solver across evaluations; both lost the prior's
per-call LU. Under SEP1 and the direct forward `scpc_gamg` is the faster at
every production resolution (the closing workstation and Quartz tables), and
it is the default since 3 October; `full_mumps` stays available as the exact
reference. Seconds a TAO iteration on 32 ranks of Quartz,
SEP2 sub-element friction, the exact front push, the mass-consistent metric:

| evaluations | mesh | `full_mumps` | `scpc_gamg` |
|---|---|---|---|
| one solve at full n (`ISMIP7_EVAL_CONTINUATION=0`; the direct forward's count) | 2 km | 227 | 175 |
| one solve at full n | 1 km | 347 | 343 |
| five solves, n from 1 (`ISMIP7_EVAL_CONTINUATION=1`, the default before PR 158) | 2 km | 982 | 2405 |
| five solves | 1 km | 1463 | 4023 on 64 ranks; on 32, evaluation 1 unfinished after 4.9 h |

Each five-solve evaluation restarts at n = 1 from the previous evaluation's
n = 3 state, and there a condensed GAMG solve needed about 625 V-cycles
against 75 at n = 3; an LU costs the same per Newton iteration whatever n is.
Under the sub-element law a condensed solve needs 103 to 125 V-cycles with one
solve an evaluation, against 31 to 64 under the cell-wise law, where
`scpc_gamg` with one solve an evaluation ran 1.8 times faster than
`full_mumps` at 2 km under today's line search, and 2.6 (2 km) and 5.8 (1 km)
times faster under the `bt` line search of the first round (the L-BFGS-B
table below). Memory decides nothing either way: the whole 1 km LU fit in
under 9 GiB a rank on 32 ranks.

`tlm_adjoint` computes the adjoint by assembling `adjoint(J)` at the recorded
state and solving it with the adjoint options (`inversion_adjoint_parameters`:
the forward's linear options, one Krylov solve, no absolute exit), so SCPC
serves the adjoint unchanged. `icepack2_tools/taped_solve.py` handles three
details of the condensed modes:

- The Newton solve holds its Jacobian at the iterate SCPC condensed
  (`frozen_linearization`), and `tlm_adjoint` refuses that callback while it
  records. The solve therefore runs with the manager paused, and the recorded
  equation confirms the converged state at SNES iteration 0. The direct
  forward (`ISMIP7_DIRECT_FORWARD`, PR 158) sends every mode this way,
  `full_mumps` with its Jacobian live, under `direct_forward_parameters`; a
  direct solve that fails leaves the state where it found it and raises, for
  the failed-trial rescue. On the slabs the direct gradient is the taped one
  to 1e-13 under the exact solvers and to 2.2e-9 under `scpc_gamg`, SEP1 as
  SEP2.
- The matrix-free adjoint solve receives no form compiler parameters, so the
  taped form carries its quadrature degree in its integrals. Without that the
  adjoint is assembled at UFL's estimated degree (up to 26) against a degree 4
  forward, and the gradient was 1e-3 off.
- One Newton solver serves the run (`StateSolverCache`). The driver builds
  its taped form once, over `(theta, phi)` on the L-BFGS-B path and over two
  Functions that TAO's fresh control copies are assigned into on the tape, so
  every evaluation reuses the solver, its SCPC context and the symbolic work
  Firedrake caches on a form. What the condensed preconditioner built at its
  first setup stays, as across the Newton iterations of one solve and the
  steps of a transient run: PETSc redoes GAMG's Galerkin products and
  Chebyshev estimates on its first interpolation, and refactors an LU on its
  first symbolic analysis. On the slabs a reused solver takes a new one's
  Newton iterations and its gradient agrees to 1.1e-12; at 32 km see the
  table closing this section.

Under `scpc_*` the taped forward keeps the transient's NLEQ-ERR line search
and solves its Newton corrections to a relative 1e-8
(`ISMIP7_INVERSION_KSP_RTOL`; the condensed solve's absolute tolerance
follows it; `ISMIP7_INVERSION_SNES_LINESEARCH` names another line search).
The transient keeps 1e-6. Each evaluation starts from the previous one's state
at new controls, at 32 km as far as ||F|| = 9e10 from the solution, and
NLEQ-ERR judges a step by the norm of the Newton correction, which a GAMG solve
to 1e-6 leaves too inexact there: at a trial point the exact solvers took in 6
Newton iterations it stalled at ||F|| = 7.9e10 for 200, at every full L-BFGS
step. Backtracking on ||F|| (`bt`) does not need exact corrections and was the
default for a day, until the synthetic slabs showed it damping nearly every
step, with an exact LU as with GAMG. Newton iterations of one taped forward:

| line search, Krylov rtol | 32 km trial point | slab, cell-wise law | slab, sub-element law | 2 km, s an evaluation |
|---|---|---|---|---|
| NLEQ-ERR, 1e-6 (the transient's) | 200, failed | 4 | 4 | |
| NLEQ-ERR, 1e-8 (the default) | 6 | 4 | 4 | 95 |
| NLEQ-ERR, 1e-10 | 6 | 4 | 4 | |
| `bt`, 1e-6 | 6 | 114 | 20 | 64 |

`tests/test_scpc_adjoint.py` holds the taped forward on both slabs to 8 Newton
iterations, which fails under `bt`.

Gradient against `full_mumps` on the synthetic slabs of
`tests/test_scpc_adjoint.py` (grounded and floating, Budd, n = 3; the
sub-element one with the exact front push on an ice-free strip and the exp
control's friction), the larger of the two controls:

| solver | relative gradient difference, cell-wise | sub-element | Taylor order |
|---|---|---|---|
| `scpc_mumps` | 2.6e-13 | 2.8e-14 | |
| `scpc_gamg` (NLEQ-ERR, 1e-8) | 4.1e-11 | 2.1e-9 | 2.00 on both |
| `scpc_gamg` at the transient's 1e-6 | 1.4e-8 | 1.6e-7 | |

The 32 km and Quartz runs below that name `bt` ran before the default moved
to NLEQ-ERR at 1e-8.

On Antarctica at 32 km (`antarctica_320000_32000_buffered0`, 6,282 vertices,
8 ranks, two runs at a time on a 16-core workstation, Budd, the `legacy` fluidity prior, a cold start,
`ISMIP7_EVAL_CONTINUATION=0`), against `full_mumps` from the same start:

| solver | optimizer | evaluations | largest relative difference in the objective | s per evaluation, forward / adjoint |
|---|---|---|---|---|
| `scpc_mumps` | L-BFGS-B | 4 | 1.1e-13 | 11 / 0.9 |
| `scpc_gamg`, `bt` | L-BFGS-B, 30 iterations | 32 | 4.3e-7 | 9.5 / 0.6 |
| `scpc_gamg`, `bt` | TAO lmvm, mass-consistent metric, bi-Laplacian prior, 20 iterations | 21 | 1.2e-8 | 19 a whole iteration |
| `scpc_gamg`, NLEQ-ERR at 1e-8 | TAO as above, with sub-element friction, the exact front push and the five-solve evaluations (`ISMIP7_SUBELEMENT_FRICTION=1`, `ISMIP7_EVAL_CONTINUATION=1`), 20 iterations | 21 | 1.6e-9 | 76 a whole iteration, against 53 |

The two L-BFGS-B runs took the same 30 iterations and 32 evaluations, from
3.0156e4 to 3.7463e3; the condensed GAMG solve averaged 51 V-cycles and the
forward 5.1 Newton iterations an evaluation. `full_mumps` took 6.5 s and 0.5 s,
which is expected at 32 km: the condensed solver pays for itself only at 1 km
on 32 ranks or more (section 7). On the TAO path both published
||F|| = 1.493 after 20 iterations; `full_mumps` failed its forward at three
trial points and took the re-ramp rescue at each, where `scpc_gamg` under
`bt` failed none. Under sub-element friction neither failed; a condensed solve
there took 10 Newton iterations and 192 V-cycles on average, against 5 and 51
under the cell-wise law. A cold start with the exp control and sub-element
friction did not climb the startup ramp at all (it diverged near n = 2.1 on
every rung, under `full_mumps` as under `scpc_gamg`), so that pair ran the log
control; the slabs cover the exp control's friction. Under SEP1 on `432c831`
the same cold start climbs the ramp on its first rung under `full_mumps`
(`runlog/test-32km-inversion-sep1-opt-exp-*`).

On Quartz (issue #156, jobs 10818443 to 10818449), from Rice's 2 km snapshot
0948, which carries controls and no state: continued on its own mesh
(`ISMIP7_MESH=checkpoint`, 925,183 vertices), and transferred onto
`antarctica_10000_1000_buffered20000` (1,869,252 vertices). Budd, the
snapshot's fluidity prior, the bi-Laplacian prior, L-BFGS-B,
`ISMIP7_EVAL_CONTINUATION=0`, the startup ramp under `scpc_mumps` in every
arm, 5 iterations at 2 km and 3 at 1 km. Medians over the evaluations after
the first; memory is sacct's AveRSS and MaxRSS a rank:

| mesh | solver | ranks | forward (s) | adjoint (s) | evaluation (s) | GiB a rank, mean / peak |
|---|---|---|---|---|---|---|
| 2 km | `full_mumps` | 32 | 102 | 52 | 168 | 2.8 / 3.6 |
| 2 km | `scpc_gamg`, `bt` | 32 | 41 | 9.3 | 64 | 3.7 / 4.0 |
| 2 km | `scpc_gamg`, `bt` | 16 | 61 | 20 | 107 | 6.1 / 6.5 |
| 2 km | `scpc_gamg`, NLEQ-ERR, Krylov rtol 1e-8 (now the default) | 32 | 70 | 12 | 95 | 3.8 / 4.0 |
| 1 km | `full_mumps` | 32 | 189 | 118 | 334 | 5.5 / 6.7 |
| 1 km | `scpc_gamg`, `bt` | 32 | 21 | 8.5 | 58 | 4.5 / 5.2 |
| 1 km | `scpc_gamg`, `bt` | 64 | 11.5 | 4.1 | 43 | 3.1 / 3.6 |

Every arm on 32 ranks ended on the objective `full_mumps` reached, to seven
digits (4.935180e4 at 2 km, 4.179544e4 at 1 km), and no forward failed. Over
every evaluation, `scpc_gamg` stayed within 3.7e-6 of `full_mumps` at 2 km
(6.5e-8 under NLEQ-ERR at 1e-8) and within 4.3e-9 at 1 km. The rest of an evaluation, 13 s at 2 km and 27 s at 1 km, is outside
the forward and the adjoint, the same under both solvers, and is most of a
64-rank evaluation; 98 % of it at 1 km was the prior's mass solve, which
`4e45164` removes (next section). The rank
count moves the first evaluation's objective (1.3e-6 at 2 km from 16 to 32
ranks, 2e-5 at 1 km from 32 to 64): a controls-only start converges the ramp
on the relative SNES test, so the state it leaves follows the partition's
rounding. Memory decides nothing at these settings: the whole 1 km LU fit in
7 GiB a rank on 32 ranks, against Rice's estimate of 410 to 480 GB for its
1 km chain (`runlog/inversion-1km.json`). Slurm samples memory every 30 s
here (`JobAcctGatherFrequency`), and a 1 km factorisation outlasts that (the
adjoint, one factorisation and its solve, took about 2 min), so its peak is in
the sample.

The production configuration (issue #156, jobs 10823630 to 10823634 and
10824069 to 10824072): the same snapshot and meshes, sub-element friction with
the exact front push, TAO with the mass-consistent metric and the bi-Laplacian
prior, `scpc_gamg` under NLEQ-ERR at 1e-8, 5 iterations at 2 km and 3 at 1 km.
Seconds a TAO iteration, the median after the first; GiB a rank, sacct's
AveRSS and MaxRSS:

| evaluations | mesh | solver | ranks | s an iteration | GiB a rank, mean / peak | V-cycles a condensed solve |
|---|---|---|---|---|---|---|
| one solve | 2 km | `full_mumps` | 32 | 227 | 2.7 / 3.4 | |
| one solve | 2 km | `scpc_gamg` | 32 | 175 | 3.6 / 3.9 | 103 |
| one solve | 1 km | `full_mumps` | 32 | 347 | 4.6 / 7.2 | |
| one solve | 1 km | `scpc_gamg` | 32 | 343 | 4.8 / 5.4 | 125 |
| five solves | 2 km | `full_mumps` | 32 | 982 | 3.6 / 4.4 | |
| five solves | 2 km | `scpc_gamg` | 32 | 2405 | 6.9 / 7.2 | 242 |
| five solves | 1 km | `full_mumps` | 32 | 1463 | 7.1 / 8.6 | |
| five solves | 1 km | `scpc_gamg` | 64 | 4023 | cancelled after 2 iterations | 403 |
| five solves | 1 km | `scpc_gamg` | 32 | none | cancelled in evaluation 1 after 4.9 h | |

Every same-rank pair that finished ended on one objective (within 1.1e-9 at
2 km and 2.1e-9 at 1 km over every iteration), the recorded solves confirmed with no
step, and the one failed trial point (the 2 km single-solve pair, iteration 3)
failed under both solvers and took the same re-ramp rescue.

The PR 155 review's two efficiency findings came after those runs. Each form
now carries the SCPC structural-zero blocks only when its own solver condenses:
the taped form under an `scpc_*` inversion solver, the startup ramp's under an
`scpc_*` lane solver. In an assembled Jacobian the blocks are 12 nonzeros a
cell, 15.8 % of the whole, and every Quartz `full_mumps` arm above ramped under
`scpc_mumps`, so it factored them, and its seconds are an upper bound. The
second finding is the one taped form and one solver for the run, described
in the list above. PR 155's head against the change, 32 km on the workstation,
8 ranks with nothing else running, Budd, the `legacy` fluidity prior, a cold
start; RSS is each rank's after the first and the last taped solve, the mean
over ranks:

| configuration | evaluations | s an evaluation, before / after | RSS a rank in MiB, before | after |
|---|---|---|---|---|
| `scpc_gamg`, L-BFGS-B, one solve an evaluation | 32 | 8.69 / 6.64 | 1125 to 1316 | 1034 to 1153 |
| `scpc_gamg`, TAO, sub-element friction, five solves an evaluation | 12 | 72.7 / 61.3 | 1027 to 1276 | 1050 to 1101 |
| `full_mumps`, L-BFGS-B, ramp under `full_mumps` | 32 | 5.56 / 4.62 | 1155 to 1011 | 998 to 1014 |
| `full_mumps`, L-BFGS-B, ramp under `scpc_mumps` | 32 | 7.36 / 4.69 | 1060 to 1073 | 962 to 1012 |

Every pair took the same iterations, with objectives within 3.5e-9 under
`scpc_gamg` and 2.2e-11 under `full_mumps`; under `scpc_gamg` the Newton
iterations matched and the V-cycles a condensed solve stayed within 2 %
(79.7 against 79.8, 207.5 against 203.4). Building an SCPC context took
0.17 s; most of what a new solver cost was the symbolic work Firedrake caches
on a form object. The blocks made a `full_mumps` evaluation 32 % slower on
PR 155's head (7.36 against 5.56 s).

RSS grew with every solver built, 6.2 MiB a rank an evaluation before and 3.8
after with one solve an evaluation, 22.6 and 4.6 with five. The rest was the
adjoint's condensed solver, which PETSc lost (issue #159). tlm_adjoint drops
its matrix-free adjoint solver without destroying it; on more than one rank
petsc4py stashes it for `PetscGarbageCleanup`, and the cleanup that destroys
it (tlm_adjoint's, in `compute_gradient`) releases the `ISMIP7SCPC` context,
whose own objects are stashed while the cleanup runs. PETSc 3.25 puts the
communicator's old garbage map back after its destroy loop, dropping the map
those stashes went into, so the condensed KSP with its GAMG hierarchy and
operator, and the weight vector, stayed at reference count 1 for good. Python
saw none of it, a later cleanup could not reach it, it scaled with the local
problem, and one rank never leaks (petsc4py destroys at once there).
`ISMIP7SCPC.destroy` now destroys the condensed KSP and operator itself and
keeps the rest of the context until the next SCPC setup, outside any cleanup.
Every evaluation in `ISMIP7_INVERSION_TIMING_JSON` carries `rss_mib`: the mean
and the maximum RSS over ranks after its adjoint (once an accepted iteration
on the TAO path), and the largest peak any rank has reached. The same 32 km
configurations rerun from 508a9be (before) and fb7c32e (after), 30 L-BFGS-B or
10 TAO iterations; MiB a rank (mean over ranks) an evaluation from that field,
over the 2nd to the 20th evaluation and over the 21st to the 32nd:

| configuration | ranks | before | after |
|---|---|---|---|
| `scpc_gamg`, L-BFGS-B, one solve an evaluation | 8 | 2.05 / 1.84 | 0.24 / -0.33 |
| `scpc_gamg`, L-BFGS-B, one solve an evaluation | 4 | 2.32 / 2.02 | 0.20 / 0.14 |
| `full_mumps`, L-BFGS-B | 8 | | 0.42 / 0.10 |
| `scpc_gamg`, TAO, sub-element friction, five solves an evaluation (an iteration) | 8 | 2.31 | 0.45 |

Each pair took the same iterations, objectives within 7.1e-11, at the same
cost (7.11 against 6.90 s an evaluation, 63.2 against 63.1 s a TAO
iteration). The L-BFGS-B rows leave out a step of 24 to 26 MiB at the 21st
evaluation, the checkpoint the driver writes at iteration 20, and every row the
first evaluation or iteration (3 to 11 MiB). On
a 2,400-cell synthetic slab on 4 ranks PETSc's own allocations grew 1.7 MiB a
rank an adjoint solve under `scpc_gamg` (0.31 under `scpc_mumps`, plus MUMPS's
factors) and 0.006 with the fix, as `full_mumps`. The lost PETSc memory grew
with the local problem, 1.5 MiB a rank an adjoint solve at 519 cells a rank and
4.6 at 2,080 (2 ranks); RSS understates it (1.7 and 4.3 MiB there, and the
32 km before arm grew alike on 8 and 4 ranks). At the slab's rate a 2 km
evaluation on 32 ranks (57,000 cells a rank) would lose about 115 MiB a rank
and a 1 km one about 230: 34 and 69 GiB a rank over 300 evaluations. Quartz did
not measure it.
`runlog/test-32km-inversion-reuse-*`, `runlog/test-32km-inversion-ramp-blocks-*`
and `runlog/test-32km-inversion-scpc-destroy-*` hold the runs, including a
first round timed beside another session's jobs.

Under `full_mumps` the 2 km inversion's RSS climbed for as long as it ran
(issue #161): job 10971250 grew 33.7 MiB a rank an evaluation over evaluations
11 to 103 and lost a rank to the OOM killer at 220G. Each evaluation builds and
drops two solvers that hold assembled matrices: the recorded solve's Newton
solver, whose mixed AIJ Jacobian is allocated even when SNES exits at
iteration 0, and the adjoint's operator, `LinearSolver` and LU. Firedrake's
`NonlinearVariationalSolver` sits in a reference cycle, so a dropped one
lives until Python's cyclic collector runs, and on more than one rank its
PETSc objects then wait for the next `PetscGarbageCleanup`. Some are never
freed: the arms below without the release ended with 25 to 37 matrices alive
in PETSc's count, against 10 (the run's own) with it, while at most four
Firedrake solvers were alive after any evaluation. That fits issue #159's
loss of objects released while a cleanup runs, reached here through the
collector; no run here isolated that path. The solvers are now destroyed when dropped
(`taped_solve.release_solver`, used by `ReleasingEquationSolver` for both
recorded equations and by `StateSolverCache` for the solver it replaces). On
Quartz, 32 ranks, 10971250's objective (rho 75 km), from its last checkpoint
or from a start that objective moves far (`runlog/test-2km-full-mumps-issue161-*`);
RSS a rank (mean) from the timing record, MiB:

| start | evaluations | without the release | with it |
|---|---|---|---|
| converged | 10 to 62 | 2,761 to 3,366 (11.6 an evaluation, still rising) | 2,657 to 2,804 (2.8; flat within 2,739 to 2,809 from evaluation 13) |
| moving (rho 750 km MAP under the rho 75 km objective) | 10 to 48 | 2,762 to 3,057 (7.8) | 2,633 to 2,781 (3.9; within 2,771 to 2,813 from evaluation 20) |
| moving, the collector left alone | 10 to 58 | 2,753 to 3,213 (9.6) | 2,635 to 2,798 (3.4; within 2,777 to 2,815 from evaluation 20) |

The peak a rank fell by about 600 MiB at 62 evaluations (5,543 to 4,915).
Objectives agree with the spread of same-code reruns, which on 32 ranks are
not bit-identical (2e-15 at evaluation 2, 1e-11 by 22, compounding through
L-BFGS-B after that), with the same Newton iterations, and the seconds an
evaluation are unchanged (136 and 128). Ruled out: lost solvers (every SNES,
KSP and PC is destroyed, here and at 32 km on the workstation); the cached forward
solver's MUMPS instance holds a constant 1,482 MB (`INFO(16)`); MUMPS's
ScaLAPACK root (`mat_mumps_icntl_13 1` changed nothing); free heap
(`malloc_trim` returned about 150 MiB a rank every evaluation and left the
trend); Python's own memory (tracemalloc flat at about 920 MiB a rank). Below
2 km the growth is small (0.5 MiB an evaluation at 32 km on 8 ranks, 4.6 at
8 km on 2, on the workstation).
`site_core.sh`'s `--kill-on-bad-exit=1` (batch_runners readme) ends a job
whose rank dies regardless.

Under SEP1 and the direct forward (PR 158's defaults; PR 155 at `432c831`), on
the workstation with other sessions sharing it: Budd, the `legacy` fluidity
prior, the log control, a cold start ramped under `full_mumps`, TAO with the
mass-consistent metric, 10 iterations. Seconds a TAO iteration, the median
after the first; the adjoint column is the iteration's time outside every
span, the adjoint and TAO's own work:

| mesh | ranks | `full_mumps` | `scpc_gamg` | adjoint, `full_mumps` / `scpc_gamg` | V-cycles a condensed solve |
|---|---|---|---|---|---|
| 32 km | 4 | 6.27 | 8.64 | 0.36 / 0.59 | 79.1 |
| 8 km | 8 | 8.22 | 9.30 | 1.35 / 0.52 | 34.5 |
| 4 km | 8 | 13.47 | 10.17 | 4.80 / 0.71 | 40.2 |

Each pair took the same Newton iterations in every direct solve and the same
objective (every iteration within 4.8e-10, gradient norms within 3.4e-9), and
no trial was lost. `scpc_gamg` comes out ahead first at 4 km, through the
adjoint: one LU factorisation of the whole mixed Jacobian under `full_mumps`,
one condensed GAMG solve under `scpc_gamg`. On all three meshes and under
both solvers 4.1 to 4.5 s of each forward lies outside its Newton solve. Peak
memory stayed between 2.4 and 2.8 GB a rank (`runlog/test-*km-inversion-sep1-*`).

On Quartz under the same defaults and the production settings (issues #156
and #157, jobs 10950098 to 10950105 on `72ac7a1`): the exp control, Rice's
2 km snapshot 0948 as the warm start with its log-velocity weight, 32 ranks,
the ramp under `scpc_mumps`, 60 iterations at 2 km and 5 at 1 km (10 for the
1 km L-BFGS-B run, job 10952185, under `scpc_gamg` only). Seconds an
evaluation, the median after the first, with the checkpoint every evaluation
writes (16 s at 2 km, 13 s at 1 km); the adjoint column holds TAO's own work
on its rows; GB a rank is sacct's MaxRSS:

| mesh | optimizer | `full_mumps` | `scpc_gamg` | forward | adjoint | GB a rank |
|---|---|---|---|---|---|---|
| 2 km | TAO | 155 | 115 | 46 / 80 | 89 / 17 | 4.8 / 2.7 |
| 2 km | L-BFGS-B, no metric | 155 | 106 | 47 / 72 | 89 / 14 | 5.1 / 2.8 |
| 2 km | L-BFGS-B, sqrt(M) coordinates | 163 | 109 | 49 / 73 | 95 / 16 | 5.0 / 2.8 |
| 1 km | TAO | 324 | 70 | 72 / 45 | 230 / 11 | 8.3 / 4.2 |
| 1 km | L-BFGS-B, no metric, 10 iterations | | 69 | 42 | 10 | 5.1 |

Each pair took the same Newton iterations in every direct solve and the same
iterates: TAO within 4.3e-7 over 61 iterations at 2 km and 1.1e-12 at 1 km,
the L-BFGS-B pairs within 3.6e-7 over their first 37 evaluations, after which
L-BFGS-B's history amplifies the difference (best objectives 5.819865e4 and
5.819830e4 without a metric). No forward failed. A condensed solve took 134 to
136 V-cycles at 2 km and 45 at 1 km. `scpc_gamg` is 1.35 to 1.50 times faster
at 2 km and 4.6 times at 1 km with half the memory, all of it in the adjoint;
its forward is the slower one at 2 km. Records
`runlog/test-*km-inversion-final-*`; the optimizer comparison on the same runs
is in `INVERSION_PRIORS.md`, issue #157 section.

### Inversion time outside the forward and the adjoint

The timing record reports `fwd_seconds`, `adj_seconds`, `total_seconds`,
`optimize_seconds` and each span as the slowest rank's duration
(`icepack2_tools/profiling.py`). Each `unspanned` duration is the slowest of
the rank-local remainders. The log prints these durations under each iteration
line. An L-BFGS-B evaluation carries `other_spans` (inside `total_seconds`) and
`before_spans` (the previous evaluation's report, term assembly, timing write
and checkpoint, and `gap`, the whole interval that holds the optimizer's own
step). A TAO iteration carries `iteration_spans` over all of its evaluations;
its `unspanned` is the adjoint and TAO itself.

Under the bi-Laplacian prior three spans held that time:

| span | path | before | since `4e45164` |
|---|---|---|---|
| `prior_solve` | both | the mass solve `M f = A θ` of each prior energy, an `EquationSolver` on the residual form: a Newton solve with a fresh MUMPS LU of `M` at every call, analysed on one rank (`ICNTL(28)=1`, the MUMPS 5.8.2 default) | `prior.BilaplacianAuxSolver`: `M` factored once (MUMPS Cholesky), then a back-substitution |
| `prior_taped` | TAO | the same solve on the tape, and one more LU of `M` in the adjoint | the same solver as a tlm_adjoint `LinearEquation`, back-substituting in the forward and the adjoint |
| `residual_norm` | TAO | the fnorm-ceiling check assembled `F` while tlm_adjoint recorded, 4 to 8 s a call at 32 km and 4 km | the check runs under `paused_manager()` |

The old solve grows with the vertex count and stays flat in the rank count,
the pattern of the 13 s and 27 s above. `scripts/probe_eval_overhead.py`
times the prior work alone on the 32 km mesh refined uniformly; seconds an
evaluation for both controls on the Mac workstation, with another 8-rank job
on its 16 cores:

| vertices | ranks | off the tape, before | parallel analysis (`ICNTL(28)=2`) | factored | TAO tape, before (forward + adjoint) | TAO tape, factored |
|---|---|---|---|---|---|---|
| 79,455 | 1 / 8 | 0.54 / 0.56 | 0.60 / 0.24 | 0.016 / 0.027 | 0.84 / 0.85 | 0.037 / 0.046 |
| 303,789 | 1 / 8 | 2.21 / 1.90 | 2.42 / 0.67 | 0.055 / 0.035 | 3.51 / 2.92 | 0.10 / 0.08 |
| 1,187,097 | 1 / 2 / 4 / 8 | 9.8 / 9.1 / 8.3 / 8.2 | 10.1 / 4.7 / 2.9 / 3.0 | 0.22 / 0.15 / 0.10 / 0.12 | 15.3 / 13.9 / 12.6 / 12.4 | 0.40 / 0.27 / 0.18 / 0.21 |

The factored solver's one factorisation took 4.1 to 4.8 s at 1.19 million
vertices. Energies agree with the old solve to 3e-16 and taped gradients to
2e-16. CG with Jacobi (22 to 26 iterations at every size) costs about what
the back-substitution does. The replicated gathers of an L-BFGS-B evaluation
(`func_to_global` four times, `global_to_func` twice) took 0.13 s on 1 rank
and 0.016 s on 8.

On Quartz the old solve is the 27 s. Jobs 10937657 (`2626c71`) and 10937658
(`4e45164`) reran job 10818449's 1 km arm (`scpc_gamg`, 64 ranks, 3
iterations) back to back on one node; seconds an evaluation, medians over
evaluations 2 to 4:

| code | evaluation | forward | adjoint | outside both | of it `prior_solve` | next largest spans |
|---|---|---|---|---|---|---|
| `2626c71` | 54.0 | 22.0 | 5.3 | 26.8 | 26.2 | `gather_gradient` 0.33, `set_controls` 0.10, `residual_norm` 0.07 |
| `4e45164` | 27.9 | 22.0 | 5.4 | 0.73 | 0.17 | the same |

That is 14 µs a vertex on Quartz against 6.9 to 8.2 on the workstation. The
factored solver's first call, which factors, took 14.8 s once per run. The
objective agreed to 7e-16 and the gradient norm to 2e-15 at every
evaluation. Between evaluations L-BFGS-B's own step, replicated on every rank
over the 3.7 million controls, took 1.0 to 1.2 s (`gap` less its spans).

In situ, base (`2626c71`, spans only) against `4e45164`, each pair back to
back: Budd, the `legacy` fluidity prior, the bi-Laplacian prior, a cold
start, `ISMIP7_EVAL_CONTINUATION=0`; `full_mumps` on 4 ranks at 32 km,
`scpc_gamg` with the ramp under `scpc_mumps` on 8 ranks at 4 km. Medians
after the first two L-BFGS-B evaluations or TAO iterations:

| mesh | optimizer | before: s an evaluation or iteration | after | spans that moved | largest relative objective difference |
|---|---|---|---|---|---|
| 32 km, 6,282 vertices | L-BFGS-B | 5.67 | 5.68 | outside the solves 0.115 to 0.080 | 1.1e-14 over 7 |
| 4 km, 117,348 vertices | L-BFGS-B | 12.5 | 13.0 | outside the solves 0.90 to 0.15 | 9.8e-14 over 5 |
| 32 km | TAO | 21.6 | 13.6 | `residual_norm` 7.86 to 0.23 | 2.2e-13 over 5 |
| 4 km | TAO | 23.0 | 14.9 | `residual_norm` 7.34 to 0.24, `prior_solve` 0.75 to 0.02, `prior_taped` 0.87 to 0.03, unspanned 1.72 to 1.30 | 9.6e-14 over 5 |

The 4 km pairs ran beside the other job (load 16 to 57), so their forward
times moved by up to 1.2 s between arms; the spans above are the comparison.
At 32 km every TAO run failed its forward at the same three trial points and
took the re-ramp rescue, now its own span (`reramp`, 2.9 s an iteration).
What an L-BFGS-B evaluation still spends outside the solves is
`residual_norm`, 0.07 s at 32 km and 0.13 s at 4 km.

### Transient (dH/dt-constrained) inversion

A velocity-only inversion fits `u` while leaving `div(h u)` unconstrained, so
the MAP can carry a flux divergence inconsistent with the observed geometry.
`ISMIP7_DHDT_WEIGHT > 0` adds one implicit-Euler prognostic step after the
diagnostic solve, using the model's own DG0 upwind operator, and scores the
resulting tendency against the observed mean dH/dt
(`icepack2_tools/obs_dhdt.py`, from the MIPkit). It needs
`ISMIP7_GEOMETRY_SPACE=dg0` and applies to grounded ice.

```bash
ISMIP7_DHDT_WEIGHT=1.0 ISMIP7_MAP_OUT=mesh/inversion_transient_2500.h5 \
  ISMIP7_LC=2500 ISMIP7_LC_COARSE=64000 \
  mpiexec -n 12 python scripts/inversion_icepack2.py
python scripts/compare_dhdt.py vel=mesh/<velocity-only>.h5 tr=mesh/<transient>.h5
```

`compare_dhdt.py` is the payoff diagnostic: both MAPs fit `u`, so the thickness
tendency is the observable that separates them. It drives the step with the
`velocity` stored in the MAP. The inversion saves that field only when the
final solve converged at the MAP's own controls and its misfit agrees with the
last accepted optimization state. MAPs written before that guard can carry a
bad one, and a score built on it means nothing. The controls are unaffected:
a forward re-solves the diagnostic from `θ` and `φ`.

`ISMIP7_DHDT_NET_SIGMA > 0` adds a term on the integrated grounded dH/dt. It is
off by default so the integrated trend stays independent validation against
IMBIE and GRACE. The per-iteration `net=` diagnostic prints either way. See
`reports/ISSUE_DRAFT_net_mass_balance_term.md`.

The MAP filename encodes friction, `LC`, geometry space and flow exponent only,
so velocity-only and transient variants collide. Give variants their own
`ISMIP7_MAP_OUT`. Every MAP records its objective as root attributes
(`misfit_norm`, `gamma_theta`, `gamma_phi`, `log_vel_weight`,
`log_vel_weight_source`, `log_vel_eps`, `dhdt_weight`, `dhdt_net_sigma`,
`mesh_basename`, `lc`, `lc_coarse`, `buffer_m`):

```bash
python -c "import h5py,sys; print(dict(h5py.File(sys.argv[1])['/'].attrs))" MAP.h5
```

A MAP without `log_vel_weight_source` predates issue 68: if a chain of several
links wrote it under `ISMIP7_LOG_VEL_WEIGHT=auto`, each link re-derived the
weight, and `log_vel_weight` is the last link's.

---

## 5. Ocean melt calibration

### The calibration every run reads

A run melts with one K and a thermal-forcing offset per IMBIE basin, the
protocol's recommendation, from the tracked file fitted under its raster
sampling (`runconfig.MELT_CALIBRATIONS`). The default sampling,
`vertex_front`, reads
`calibration/deltaT_per_basin_1000_K6.500e-05_vertex_front.npz` (issue #167,
run record `calibration-melt-1km-vertex-front`): the same K, with the offsets
refitted on the same mesh under the front-cell rule, whose melt-receiving
area is 2.4 % smaller there; they run from -0.72 K to +1.22 K. `vertex`
sampling reads `calibration/deltaT_per_basin_1000_K6.500e-05.npz`, the file
the table describes:

| | |
|---|---|
| K | 6.5e-5, the K50 of the rule-based selection below on the 1000 m / 10 km production mesh (run record `calibration-melt-toolbox-1km-rule`, IU's build), chosen by the group on 25 September 2026 (issue 26) |
| offsets | 16 IMBIE2 basins, -0.68 K to +1.20 K (Amundsen), each basin at its July 2026 total, 1067.4 Gt/yr together, fitted at that K on Rice's build of the mesh, the submission mesh (issue 20; `calibration-melt-refit-1km-rice-k50`) |
| fitted under | `ISMIP7_MELT_SLOPE=ant`, `ISMIP7_SIN_ALPHA_ANT=5.115e-3`, `ISMIP7_GEOMETRY_SPACE=dg0`, `ISMIP7_RASTER_SAMPLE=vertex` and the 30_sep OI climatology, on `antarctica_10000_1000_buffered20000`, Rice's build (1,869,252 vertices) |
| sidecar | `deltaT_per_basin_1000_K6.500e-05.source.json`: the file's sha256, the settings above, the mesh build and its vertex count, the input hashes, the jobs and the commits, with the decision fields added on promotion (below) |

Every clone carries it, so a new machine runs with it and nothing is
calibrated or copied. The offsets are stamped onto any mesh through the
IMBIE2 8 km basin grid under `ISMIP7_DATA_ROOT`. A fit records its inputs in
the npz by basename (`obs_csv`, `imbie2_nc`, `inversion`;
`calibrate_melt.input_names`), and the grid is found by that name under the
data root. The `vertex` file predates that and still names three paths on IU
Quartz; replacing it is open (issue #150).

The forward applies the melt the file was fitted to:

- it melts the cells the fit summed over, floating and holding ice on a bed
  below sea level (`forcing.melt_receiving`);
- an offsets file whose recorded slope law, slope constant or geometry
  space differs from the run's stops the run, and so does geometry sampled
  with another `raster_sample` than the file's, or a cold start that floors
  the initial thickness (`ISMIP7_H_CLAMP_INIT` above 0, the legacy friction
  law's default);
- the tracked file has to match the sha256 its sidecar records;
- `ISMIP7_K_SCALE` other than 1 is refused with an offsets file, and the
  removed `ISMIP7_K_MELT` is refused when exported;
- every run logs a `Forcing provenance:` line naming the file, its sha256
  and K, and the mesh it was fitted on beside the run's, and
  `core_report.py` carries that line into the report.

On another mesh the forward runs with these offsets and says so in that
line, so a coarse probe needs nothing more. `preflight.py` refuses a core on
a mesh other than the calibration's until `ISMIP7_DELTAT_PER_BASIN_NPZ` names
a file: offsets refitted on that mesh at the same K, or the tracked file to
run with its offsets as they are. It also refuses another build of the same
mesh name, told apart by the vertex count in the `.msh` header against the
sidecar's: the offsets were fitted on Rice's build of the production mesh
(1,869,252 vertices), and IU's has 1,869,088. A refit is a job
(`calibrate_deltaT.script`, which also reports whether the thermal forcing
rule admits the K there); `ISMIP7_INV_H5` names a checkpoint on the mesh or
the `.msh` itself:

```bash
scripts/batch_runners/submit.sh script scripts/batch_runners/calibrate_deltaT.script \
    --cd antarctica --queue debug --tasks 1 --mem 32G --time 00:40:00 \
    ISMIP7_LC=1000 ISMIP7_INV_H5=<mesh.msh> DELTAT_K=6.5e-5 DELTAT_OUT=<absolute dir>
```

The fit writes each offsets file with its sidecar, `<name>.source.json`
(`calibrate_melt.calibration_record`), and so does `select_melt_parameters.py`
for the K it selects:

| field | holds |
|---|---|
| `file`, `sha256` | the npz as written |
| `K`, `melt_slope`, `sin_alpha_ant` or `sin_alpha_cap`, `geometry_space`, `raster_sample`, `oi_version`, `rho_i`, `dt_window`, `tf_rule`, `rule_admits` | the settings of the fit and the thermal forcing rule's verdict (`null` from a parallel fit, which does not judge) |
| `mesh`, `vertices`, `cells`, `floating_cells`, `mesh_file` | the mesh's name and counts, and the sha256 and md5 of the file `ISMIP7_INV_H5` named |
| `obs_table`, `bedmachine`, `inputs_sha256` | every input, by name and sha256 |
| `melt_total_gtyr`, `melt_total_dT0_gtyr`, `unrooted` | the basins' total at their offsets as the fit summed it, the total at no offset, and the basins with no root in the window |
| `selected_as`, `refit_of`, `selection` | the K's selection: a refit at the tracked K names the tracked file and keeps its `selected_as` (K50); a selection names its toolbox commit |
| `site`, `partition`, `ranks`, `job`, `code`, `code_modified` | `ISMIP7_SITE` as `submit.sh` exports it (`null` for a fit started by hand), the Slurm partition and job, the rank count, the commit, and the tracked files modified in the checkout |

A run that names the file with `ISMIP7_DELTAT_PER_BASIN_NPZ` is then checked
against its raster sampling, and its provenance line names the mesh the
offsets were fitted on. A file written before the fits wrote sidecars is
read as given. On the 25 km rehearsal mesh a refit at the tracked K
reproduced the rehearsal's npz bit for bit, and its sidecar matched the
hand-written one in every field a fit can know (run record
`calibration-melt-sidecar-25km-check`).

Promoting a fit to the tracked calibration copies both files into
`calibration/`, adds by hand the decision fields the sidecar lacks
(`decided`, `decision`, `selected_as`, `selection`, `run_record`,
`mesh_build`) and takes a run record; a fit refuses to write into
`calibration/` itself.

`check_melt_bound.py` melts the reference geometry with the forward's own
callback and sets each basin against the total its offsets were fitted to,
and exits 1 when a basin is off by more than `--match-tol` (default 0.1
percent):

```bash
ISMIP7_INV_H5=<mesh.h5 or mesh.msh> python antarctica/scripts/check_melt_bound.py
```

On Rice's build of the production mesh (run record
`calibration-melt-refit-1km-rice-k50`, Quartz job 10649438) the forward melts
1067.389 Gt/yr against the 1067.386 its offsets were fitted to, every basin
within 0.006 Gt/yr, and no cell passes the `libmassbffl` bound (maximum 41.7
m/yr). Its earlier melt set also covered 257 836 ice-free open-ocean cells,
where it booked 156.2 Gt/yr of melt and 23.1 Gt/yr of refreezing, and 6 805
cells of bare land, where the climatology melts nothing. The two builds need
their own fits: the offsets differ by at most 0.0008 K, and the fit on IU's
build missed basins 1 and 7 on Rice's by 0.18 and 0.17 percent (job 10649416).
At 32 km (job 10644430) the offsets put the basins at 0.33 to 1.74 times their
totals and the whole at 1069.5 Gt/yr, which is why a production core runs on
the calibration's mesh.

`ISMIP7_K_PER_BASIN_NPZ` names a legacy per-basin K file (next subsection)
in place of the offsets; no run searches for one.

### The legacy per-basin K (`calibrate_melt.py`)

Solves for the Burgard quadratic-mixed-slope coefficient K, global and per
IMBIE2 basin, against integrated observed shelf melt. The target is the July
2026 table combining Paolo, Davison and Adusumilli, 1067.4 Gt/yr,
`Melt_Paolo_Davison_Adusumilli_imbie2.csv`, searched for under
`<DATA_ROOT>/meltobs/` and then `<DATA_ROOT>/parameterisations/ocean/meltobs/`.
With it in neither place the script falls back to the older Paolo and Adusumilli
table (865.0 Gt/yr) under `<DATA_ROOT>/parameterisations/ocean/meltobs/` and
prints a `[!]` line saying so; `ISMIP7_MELT_OBS_CSV` names either. Every run
prints the table it opened and its integrated target. Needs section 2 forcing
and a section 4 mesh.

The newer table comes from the Source Cooperative melt-calibration product:

```bash
python antarctica/scripts/download_mirror.py \
    --product ismip7-ais-melt-calibration data/meltobs/
```

```bash
cd antarctica
ISMIP7_LC=2500 python scripts/calibrate_melt.py
# results/calibrated_K_per_basin_<LC>.npz
```

Only the mesh is read from the MAP, so any MAP built on it serves. The default
is the section 4 name for the configured `ISMIP7_FRICTION`; `ISMIP7_INV_H5`
names a different one.

The calibration melts on the same `ISMIP7_GEOMETRY_SPACE` as the forward
(default `dg0`) and with the same `ISMIP7_MELT_SLOPE` (default `ant`, the
constant `ISMIP7_SIN_ALPHA_ANT` on every shelf): the cells the forward melts,
through the forward's own path, so the K it writes is the K the forward
applies. Its per-basin flags compare against the toolbox's July 2026 K05, K50
and K95. `ISMIP7_GEOMETRY_SPACE=cg1 ISMIP7_MELT_SLOPE=local` is the nodal
calibration the earlier K files came from, with the slope capped at 5e-3;
`ISMIP7_SIN_ALPHA_CAP` names a cap on the local slope on either geometry
(`GEOMETRY_DISCRETIZATION.md`). The file records the geometry and the slope
convention it was fitted under, and a run under another is told once at
startup.

A run reads this file only when `ISMIP7_K_PER_BASIN_NPZ` names it, and
`ISMIP7_K_SCALE` multiplies it there. `ISMIP7_K_OUT` writes it elsewhere, and
`check_melt_bound.py --npz` reads it from there:

```bash
ISMIP7_K_OUT=/scratch/check/K_2000.npz ISMIP7_LC=2000 \
    python scripts/calibrate_melt.py
```

---

### Per-basin deltaT at one toolbox K (`calibrate_deltaT.py`)

The ISMIP7 ocean-forcing recommendation calibrates one dimensionless K from
the 4-term toolbox (K05, K50, K95) and then, optionally, a thermal-forcing
offset per IMBIE basin at that K so each basin's present-day melt matches the
observed total (`optimise_deltaT` in the toolbox notebook). That offset, not
a per-basin K, is the protocol's per-basin knob.

```bash
ISMIP7_LC=2000 python antarctica/scripts/calibrate_deltaT.py          # the notebook's K05 K50 K95
ISMIP7_LC=2000 python antarctica/scripts/calibrate_deltaT.py --K 6.5e-5
```

runs the fit on the forward's own melt path (the same DG0 geometry, OI
climatology at the draft, constant mean-Antarctic slope and seawater
flotation test the forward uses, from `calibrate_melt.forward_geometry`),
solving `M_b(deltaT) = M_obs(b)` per basin by a bracketed root in plus or
minus 3 K (`melt_selection.DT_WINDOW`; the toolbox searches plus or minus
2 K and the protocol sets no window), and writes
`antarctica/results/deltaT_per_basin_<lc>_K<K>.npz` (basin ids, offsets, K,
residual, dM/dT, the slope and geometry conventions) with its
`.source.json` sidecar (section 5). A run applies it in
place of the tracked calibration with

```bash
ISMIP7_DELTAT_PER_BASIN_NPZ=antarctica/results/deltaT_per_basin_2000_K6.500e-05.npz
```

and every ocean callback (control, projections, OCX) adds the offset to TF
and melts with the file's one K. A driver refuses to start when the file is
missing or `ISMIP7_K_SCALE` is not 1, since the offsets were fitted at the
file's K.
The fit reduces every basin total across ranks, so `mpiexec -n N` writes the
same offsets as a serial run.
Measured on the 2 km MAP mesh against the July table (24 September 2026, run
record `calibration-melt-2km-1067`): K05 907, K50 1623 and K95 2625 Gt/yr
uncorrected, each brought to 1067.4 by offsets from -0.55 to +1.72 K at K05,
-0.85 to +0.87 K at K50 and -1.17 to +0.30 K at K95, every basin with a root
in the window; Amundsen (basin 9) takes the largest, +1.72 K at K05. Against
the 865 Gt/yr table the same fit needs offsets within plus or minus 1.3 K.
A file reproduces the basin totals on the mesh it was fitted on. On the
1000 m / 10 km production mesh (`calibration-melt-1km-1067`) the uncorrected
totals are 2 percent higher and the offsets differ by up to 0.12 K; the 2 km
files applied there put 1072 to 1095 Gt/yr on the fitted basins to first
order, with single basins up to 44 percent off (basin 6 at K95). The 1067.4
covers the fitted basins. Floating ice outside them (a basin the table lacks,
a coverage gap, off the 8 km grid) keeps a zero offset and still melts at the
file's K, as the toolbox applies its one K everywhere, so a run's integrated
melt exceeds 1067.4 by that amount; the first forcing step prints it
(`Melt ... Gt/yr, of which ... outside the fitted basins`).

### K05, K50 and K95 from the toolbox objective (`select_melt_parameters.py`)

Protocol section 2.3 asks for the percentiles of K to come out of the
toolbox's objective with melt computed by the model's own code on its own
grid, and prefers fitting the per-basin offsets for every K before the
objective runs. `select_melt_parameters.py` does both on the forward's DG0
cells:

```bash
ISMIP7_LC=2000 ISMIP7_INV_H5=<mesh.h5> ISMIP7_MELT_OBS_CSV=<table.csv> \
    python antarctica/scripts/select_melt_parameters.py --geometry mesh --out /abs/dir
```

For each K of the toolbox notebook's grid (120 values, 2.5e-6 to 3.0e-4, or
on its step to `--k-max`), variant `per_k` fits one offset per IMBIE2 basin
to the table with the fit above, melts the present-day climatology, the 12
ocean-model states and the
13 observed states with those offsets, and sums the toolbox's four terms on
the cells (`icepack2_tools/melt_selection.py`); variant `none` melts without
offsets, the notebook's own order. The toolbox's
`calculate_objective_function`, vendored unchanged
(`icepack2_tools/ismip7_parameter_selection_toolbox.py`, with its MIT
licence and source record), then draws the notebook's weights 100000 times
per seed; seed 0 is the headline and seeds 1 to 4 the spread. Two
departures from the notebook are deliberate: an offset is a bracketed root,
where the notebook takes the nearest point of a 0.04 K grid, so term 1 is
flat across the K whose basins all root; and each ocean state melts with its
own salinity, as notebook cells 12, 15 and 65 do (cell 63 reuses the
climatology's).

Under `per_k` the objective runs on the K that `melt_selection.TFRule`
admits. Section 2.1 asks the offsets to keep present-day thermal forcing
"not significantly below 0 degC or above 5 degC", and this repository reads
significantly as a bound on every floating cell and a bound on a share of
each basin's area. A K is admitted when, with each basin's offset applied:

| test | default | option |
|---|---|---|
| every basin reaches its observed total inside the offset window | plus or minus 3 K | `--dt-window`, `--keep-unfitted` |
| every floating cell at or above | -1.8 degC | `--tf-floor` |
| at most this share of any basin's floating area below | 25 percent below -1.0 degC | `--tf-floor-area` |
| at most this share of any basin's floating area above | 25 percent above 5.5 degC | `--tf-cap-area` |
| every floating cell at or below | no bound | `--tf-cap` |

The warm side applies the cold side's 25 percent share half a degree past
the protocol's 5 degC and bounds no single cell; `none` drops any of the
four TF tests. The toolbox searches
offsets in plus or minus 2 K and keeps a K whose fit ends at the window
edge, with its residual in term 1; the objective over every K, handled that
way inside the same 3 K window, is recorded beside the admitted one in
`selection_per_k.json`.
`--geometry notebook8km` runs the same aggregation on the notebook's 8 km
grid beside the toolbox's own `calculate_term1..4`, the check that comes
before a mesh result is read. Everything is written under `--out`:
`ensemble_<variant>.nc` (terms, offsets, residuals and TF plausibility at
every K), `selection_<variant>.json`, `tf_present_<lc>.npz` (each floating
cell's present-day TF, area and basin, enough to test another rule against
the ensemble's offsets) and, for `per_k`, `deltaT_per_basin_<lc>_K<K>.npz`
at each selected K, which `ISMIP7_DELTAT_PER_BASIN_NPZ` applies at that
file's K. On a cluster: `scripts/batch_runners/select_melt_parameters.script`.

Measured on 24 September 2026 against the July table (run records
`calibration-melt-toolbox-8km`, `-2km` and `-1km`), the evidence the
percentile was chosen from (issue 26):

| grid | offsets | K05 | K50 | K95 |
|---|---|---|---|---|
| notebook's 8 km | none | 4.75e-5 | 8.5e-5 | 1.375e-4 |
| 2 km MAP mesh | none | 4.5e-5 | 8.75e-5 | 1.40e-4 |
| 1000 m / 10 km production mesh | none | 4.5e-5 | 8.5e-5 | 1.375e-4 |
| 2 km MAP mesh | fitted for every K | 2.75e-5 | 5.75e-5 | 2.725e-4 |
| 1000 m / 10 km production mesh | fitted for every K | 2.75e-5 | 6.25e-5 | 2.70e-4 |
| 2 km MAP mesh, grid to 1e-3 | none | 4.5e-5 | 8.0e-5 | 1.375e-4 |
| 1000 m / 10 km production mesh, grid to 1e-3 | none | 4.25e-5 | 7.75e-5 | 1.375e-4 |
| 2 km MAP mesh, grid to 1e-3 | fitted for every K | 2.5e-5 | 7.5e-5 | 4.475e-4 |
| 1000 m / 10 km production mesh, grid to 1e-3 | fitted for every K | 2.5e-5 | 7.5e-5 | 4.15e-4 |
| 2 km MAP mesh, grid to 1e-3 | fitted for every K in 3 K, admitted K only | 2.5e-5 | 7.0e-5 | 3.225e-4 |
| 1000 m / 10 km production mesh, grid to 1e-3 | fitted for every K in 3 K, admitted K only | 2.5e-5 | 6.5e-5 | 2.525e-4 |

On the notebook's grid the aggregation matches the toolbox's term functions
to 1.6e-12 and reproduces the notebook's printed percentiles and totals
(877.8, 1570.7 and 2540.9 Gt/yr against 878, 1571 and 2541). With the
offsets fitted first, term 1 is flat wherever every basin roots and terms 2
to 4 place K. Below K = 4.25e-5 (4.0e-5 on the 1000 m mesh) Amundsen
(basin 9) cannot reach its total inside plus or minus 2 K, and 34 percent of
the samples (26) still land there, K05 among them; 3.8 percent (3.5) land on
the grid's top, 3.0e-4. `--k-max 1e-3` (run records ending `-wide`)
extends the grid on its step: 14 percent of the samples (13) then lie above
3.0e-4 and 0.01 percent reach 1e-3, so the tail is complete. The toolbox
divides each term by its median over the grid, so the grid's reach moves
every percentile; the table's rows to 1e-3 differ from the notebook grid's
at K50 too. At K95 the offsets put over 40 percent of the shelf area of nine
basins (eleven) below 0 degC. On Quartz the 30_sep files the
forward reads are byte-identical to the notebook's 06_nov tf v3 and so v4.

The rows marked admitted K only (25 September 2026, run records ending
`-rule`) search the offsets in 3 K and apply the thermal forcing rule above.
Every basin then fits from K = 2.5e-5, where Amundsen takes +2.96 K on the 2
km mesh (+2.88 K) and puts 12 percent of its shelf (11) above 5.5 degC, and
the admitted K run from 2.5e-5 to 3.5e-4 on the 2 km mesh and to 2.55e-4 on
the 1000 m mesh. The top is where basin 4 passes 25 percent of its area
below -1.0 degC, and that basin alone moves the top between the meshes. The
-1.8 degC floor never binds and the 5.5 degC test fails only at K of 1.0e-5
and below, so the 3 K window sets the bottom. 12 percent of the samples (11)
sit on the bottom and 3.75 percent (4.9) on the top, and the seeds agree
within one step. At the 2 km K05, K50 and K95 the present-day dM/dT is 1362,
2123 and 5305 Gt/yr per K (1386, 2087 and 4476), against 1787 to 2893 at the
notebook's percentiles; the term 3 warm-minus-cold response is 0.74, 1.56
and 5.2 times the ocean models' (0.73, 1.45 and 4.2), and at K95 46 percent
of the shelf area (40) refreezes at present day. The superseded
`-rule-cap68` runs also bounded every cell at 6.8 degC, which lifted the 2
km K05 to 2.75e-5; a 5 degC bound on every cell would leave 7.25e-5, 8.0e-5
and 3.475e-4 (7.25e-5, 7.25e-5 and 2.55e-4) with half the samples on the
bottom. The objective over every K, unfitted K kept, gives 1.75e-5, 7.75e-5
and 4.475e-4 (4.15e-4). The offsets at issue 30's K match its files within
3.7e-5 K, the root tolerance in the wider bracket.

## 6. Control and projections

Forward runs go through `scripts/simulation.py` (`setup_model` and
`run_simulation`). Drivers live in `scripts/control/`, `scripts/historical/`
and `scripts/projections/`.

```bash
cd antarctica
mpiexec -n 12 python scripts/control/run.py
# results/ctrl2015_<esm>_<lc>_{final.h5, t<year>.h5, timeseries.csv}
mpiexec -n 12 python scripts/projections/ssp585_cesm_waccm.py
# results/ssp585_cesm2_waccm_<lc>_{final.h5, timeseries.csv}
```

The other scenario drivers (ssp126 and ssp370 for both ESMs, plus `ocx.py`) are
shims over `scripts/experiment.py`. Run a historical driver first to produce
`results/hist_<esm>_<lc>_final.h5`: projections and the control both branch
from it, so they share a t=0 state and the same frozen apparent-MB correction,
and their relaxation drift cancels in projection minus control (the ISMIP6
ctrl_proj convention). Without it a projection cold-starts from BedMachine and
the control warns that it starts from a different geometry. An endpoint whose
`t_yr` is short of the branch year (a chain that stopped early) or missing is
refused (`simulation.historical_endpoint`).

Run management on the drivers: `--restart <ckpt>` or `ISMIP7_RESTART` resumes;
`ISMIP7_AUTO_RESUME=1` picks up the newest checkpoint for the experiment, which
is what lets a chained batch job continue itself, and takes precedence over the
historical endpoint so only the first link starts there; `--tag` or
`ISMIP7_RUN_TAG` suffixes the experiment name so a method line keeps and
resumes its own files; `--checkpoint-interval` sets the step-count fallback.
Checkpoints carry the mesh, geometry, inversion fields and the full `(u, M, τ)`
state, so restarts work at any rank count. A resume refuses to start when
`ISMIP7_FRICTION`, `ISMIP7_APPARENT_MB` or `ISMIP7_SMB_ELEVATION_FEEDBACK`
disagree with the checkpoint.

**SMB-elevation feedback.** Every driver adds `dacabfdz(t) (s - s_ref)` to the
SMB each step (`forcing.SMBElevationFeedback`), the form and the gradient the
SMB focus group's Atmospheric forcing README (September 2026) recommends for
Antarctica. Both surfaces are the flotation surface `max(b + h, (1 - rho_I/rho_W) h)`:
of the thickness at the start of the step, and of `H_init`, the thickness of
the chain's initial state, which every checkpoint carries. A projection or
control branched from a historical therefore measures its surface change from
the historical's initial state, and the change is exactly zero at a cold
start, where the `balance` apparent-MB reference folds the first step's SMB
into `a_ref`. The feedback enters `accum`, so the transport, the `smb_gtyr`
budget column and the submitted `acabf` all carry it. The gradient each core
reads:

| cores | gradient |
|---|---|
| 1 to 8 | the core's own ESM and scenario, same product and version as its `acabf-anomaly` |
| 9, 10 | the ESM's `ctrl`, the same in every year |
| 11 | the OCX product's, at v2 or newer in both `ISMIP7_OCX_FORCING` modes (the v1 was spatially shifted, discussion #45) |

`ISMIP7_SMB_ELEVATION_FEEDBACK=0` turns it off. With it on, a run refuses to
start on a gradient that is absent, short or below its version floor, and the
preflight reports the same cores BLOCKED. Every forward checkpoint records the
mode in its `smb_elevation_feedback` attribute (absent reads as off), and a
restart refuses a checkpoint from the other mode unless it is an adapted t=0
state: a chain carries the feedback from its cold start or not at all, so a
projection with the feedback on needs a historical run with it on. Its size at
32 km, CESM2-WACCM core 1 then core 7, as the chain with it minus the same
chain without it:

| historical start | SMB at 2015 | at 2150 | at 2300 | VAF at 2300 | record |
|---|---|---|---|---|---|
| 2003 | -0.1 Gt/yr | -22 Gt/yr | -866 Gt/yr | -43.9 mm SLE | `runlog/core07-32km-ssp585-cesm2waccm-i116y2003on.json` |
| 1850 | -3.6 Gt/yr | +39 Gt/yr | -613 Gt/yr | +17.2 mm SLE | `runlog/core07-32km-ssp585-cesm2waccm-i116on.json` |

The 1850 chain entered 2015 with a surface change of up to +790 m, 165 years
of drift, and the 2003 chain with up to +414 m. The controls from the 2003
start carry almost none of the difference: by 2301 the feedback changes VAF by
-1.3 mm SLE in core 9 and -0.9 in core 10
(`runlog/core09-32km-ctrl2015-cesm2waccm-i116y2003on.json`,
`runlog/core10-32km-ctrl2015-mriesm20-i116y2003on.json`), so core 7 minus core
9 keeps -43.7 of the -44.9 mm SLE.

**Is the run on track?**

```bash
python scripts/check_ismip6_track.py results/<exp>_timeseries.csv
# exit code 0 when no FAIL rows, so gates can chain on it
```

It audits a timeseries against IMBIE dM/dt, Rignot melt and calving, RACMO SMB
and ISMIP6-class control drift, with a runaway detector.
`scripts/compare_ismip6.py <proj.csv> <ctrl.csv>` overlays projection minus
control sea-level contribution on the ISMIP6 ensemble.

Read-only diagnostics:

| | |
|---|---|
| `compare_runs.py LABEL=results/<a> ...` | overlays budget timeseries to show where two runs part ways |
| `plot_movie.py results/<exp>` | thickness change, speed and thickness frames plus an mp4 under `figs/movie_<exp>/` |
| `region_budget.py <ckpt>.h5 [<later>.h5] [--csv <run>_timeseries.csv]` | splits the budget into grounded and floating ice, so a control that gains volume above flotation is read against the observed 2000 to 2200 Gt/yr of discharge |
| `score_map.py MAP.h5 [...] [--json OUT] [--restart STATE.h5] [--save-state DIR]` | scores an inversion by `Q(u_model)/Q(u_obs)` across its own grounding line, overall and per speed band; re-solves through `setup_model`, so periodic MAPs without a velocity work. With `ISMIP7_MESH` set the MAP is transferred first, so the same command with and without it compares a transfer; `--restart` scores a prepared state, `--json` writes the numbers with the transfer fill counts (`MAP_CHECK.md`), `--save-state` writes each solved state as `DIR/<MAP stem>_t0.h5` |
| `check_cliff_facets.py STATE.h5 [--compare OTHER.h5] [--out PREFIX]` | every facet between ice and an ice-free cell, classed by the neighbour's bed against the ice base and surface, with the push of `ISMIP7_EXACT_FRONT` versions 1 and 2, the change per IMBIE basin, and on a solved state each grounded edge cell's change against its basal drag; `--compare` takes the same MAP solved under the other version (`score_map.py --save-state` on a copy of the MAP whose `exact_front` attribute is rewritten, since a forward refuses an `ISMIP7_EXACT_FRONT` that differs from its MAP) and adds the speed change by distance to the changed facets and the grounded discharge per basin (issue #166) |
| `plot_map.py MAP.h5 [--diff B.h5]` | model and observed speed and their difference, `θ`, `C = C_w0 exp(θ)` on grounded ice, and `φ`, into `figs/maps/` |

`region_budget.py` and `score_map.py` take the run's environment, which must
match the inversion's.

### Calving front on a buffered mesh (`ISMIP7_CALVING`)

A buffered mesh has no calving sink: ice reaching the 2015 outline flows into
empty buffer cells and only shelf melt removes mass, so a control gains mass
(+1000 to +1500 Gt/yr in the September 2026 32 km case with
`ISMIP7_APPARENT_MB` unset). `ISMIP7_FIXED_FRONT` removes what crosses the
outline while leaving it fixed.

`ISMIP7_CALVING` replaces that with the shared level set
(`icepack_tools.levelset`, wrapped by `icepack2_tools/levelset.py`, the object
CalvingMIP also runs). Each step `phi` solves the eikonal problem
`|grad phi| = 1` with `phi = 0` on the facets between ice and ice-free cells of
the transport's own thickness, negative in ice and positive in water. The
calving rate `c` then retreats the front by `phi_t - c|grad phi| = 0` (Hahn,
Mikula and Frolkovic 2025, arXiv:2504.05845), linearised with the previous unit
gradient and solved by cell-centred finite volumes.

Advance needs no extra mechanism: the upwind DG0 transport fills any cell the
ice flows into and the next extent includes it. Removal conserves calved mass
in three parts, all booked to the `calv` column:

1. cells the front passed entirely (`phi > 0`) are emptied;
2. every front cell sheds `min(1, c dt L/A)` of its thickness, the mass
   `c h L dt` a front retreating at `c` loses, which carries sub-cell retreat
   between steps;
3. the sliver left in a cell that held ice when the step began.

A cell below `ISMIP7_FRONT_HMIN` (1 m) that already held ice is a retreating
front cell and is emptied. A cell that was ice-free keeps whatever the
transport put there, which is how the front advances. Outside the extent the
thickness can therefore be small and nonzero.

Momentum needs no front term: under DG0 geometry the facet term
`rho g avg(h) jump(s)` at an ice/water face is the terminus water-pressure
force for floating ice, and `ISMIP7_EXACT_FRONT` gives every ice edge inside the
mesh the free-cliff push, a grounded marine cliff and a land margin included
(version 2: the push of the face above the neighbour's bed, none against rock
above the ice surface).
The one momentum change is the drag gate: the mask is 1 only where `phi`
exceeds one cell diameter, so floor-cell ocean drag acts only in water further
than a cell from the front. Thin cells inside the t=0 extent are damped by
`h_visc_floor`, the friction law and the `ISMIP7_ALPHA_GL` collar. In the strip
a free law advances into, `C_w0` comes from the t=0 geometry where `H = 0`, so
`tau_b = 0` under both laws and the damping is `h_visc_floor` and the collar.

**The laws have one home, `icepack_tools.calving`,** beside the level set that
moves the front with their rate. The CalvingMIP project runs from the same
registry and tunes against Antarctic fronts with it (`calving/tune_greene.py`,
per-Mouginot-basin flux against the Greene et al. 2022 fronts), so a tuned
parameter means the same thing there and here. `ISMIP7_CALVING` names the law,
`ISMIP7_CALVING_PARAMS` gives its parameters as `key=value,key=value`, checked
against the law's own when the forward starts (before the MAP load), and
`ISMIP7_CALVING_MODULE` registers a law from a file first:

| `ISMIP7_CALVING` | rate `c` [m/yr] | parameters (defaults) |
|---|---|---|
| `none` | no level set; `ISMIP7_FIXED_FRONT` decides the removal | |
| `fixed` | front frozen at the t=0 extent | |
| `velocity` | `u . n`, the rate that holds a front still | `iv` (`normal`) |
| `thickness` | `max(0, 1 + (Hc - H)/Hc) |u|` where the bed is below sea level (CalvingMIP experiment 5) | `hc` (375 m), `marine_only` (true) |
| `vonmises` | `|u| sigma~ / sigma_max`, `sigma~` from the tensile principal deviatoric stresses of the solved `M` | `sigma_max_gr` (1 MPa), `sigma_max_fl` (0.15) |
| `vonmises_strain` | the same from the strain rate, `sqrt(3) B eps~^(1/n)`, `B = A^(-1/n)` from the run's fluidity (Morlighem et al. 2016) | same |
| `hfb` | `|u| min(R_xx / R_crit, ratio_max)^p`, the horizontal-force-balance threshold on the front-normal resistive stress | `sigma_max` (0), `rho_c` (`seawater`), `mode`, `stress`, `exponent`, `ratio_max` |

The two von Mises forms differ by up to `2^(1/3)` (uniaxial extension), so a
threshold tuned for one is not one for the other; `vonmises_strain` is what
`ISMIP7_CALVING=vonmises` ran before 24 September 2026. A law reads its fields
off the forward's live state (`simulation.calving_front_state`, an
`icepack_tools.calving.FrontState`): `u`, `M`, `tau` from the mixed solution,
the DG0 `h`, the bed, the height above flotation and grounded indicator under
the forward's own densities, the front normal of the level set it advances,
and the run's `A` and `n`. Its rate goes to the level set as the `prescribed`
rate. The `Calving front owner:` line and every checkpoint's `calving_law`
attribute record the law with every parameter at full precision, and
`core_report.py` lifts the line into the report.

The level set is checkpointed as `levelset` for diagnostics; a restart
rebuilds the front from the thickness. The exception is `fixed`, which anchors
on `H_init` so a resumed run does not re-freeze the front where it restarted.
The laws and the shared level set are tested in `icepack_tools`
(`test/calving_test.py`, `test/levelset_test.py`); the ISMIP7-side rules
(retreat-sliver mask, apparent-MB extent masking, the `fixed` law's t=0
anchor) are covered by `tests/`. `tests/test_levelset_laws.py` checks the
configured law driving the front against closed forms on a unit mesh, the
shed fraction and its step-size behaviour under the transport's masks, the
drag gate, and the refusal of an unknown, misconfigured or underspecified law;
`tests/test_calving_front_state.py` checks that the state a law reads is the
forward's own.

**Control and projection configurations differ.** The protocol's control is an
unforced constant-climate run with fracture, collapse, calving and GIA held at
end-of-2014 conditions (the April 2026 protocol cheat sheet in `../protocol/`),
so the control here runs with `ISMIP7_FIXED_FRONT=1` and `ISMIP7_CALVING=none`.
It also runs with `ISMIP7_APPARENT_MB`, a choice of this repository: the
protocol leaves initial conditions to each group and keeps the control to
assess drift, and whether the production runs keep the reference is a group
decision (issue #104). That is what `run_core_matrix.sh` runs and what every
control result used. `ISMIP7_CALVING=fixed` also pins the front and is a
different run: it builds a level set, so ocean drag is gated off near the front
and the retreat-sliver rule applies inside the t=0 extent, giving a slightly
different `calv` column and settled front. Both close the budget.

A law other than `fixed` is for projections. A configured law owns the front
outright, so the legacy `ISMIP7_FIXED_FRONT` mask removes nothing when a law is
set and a free law is never silently pinned. The apparent-MB reference `a_ref` is
defined only on the t=0 ice extent under every law. Under a free law it is also
cleared each step wherever the level set reports ice-free, irreversibly, so a
calved cell is not regrown and an advanced-into cell is not re-emptied. The run
log prints one `Calving front owner:` line naming the mechanism in force.

### The whole matrix in one command (`run_core_matrix.sh`)

Runs cores 1 to 11 in dependency order (historicals, controls, projections,
OCX) and finishes with the audit and ensemble comparison for each core that
reached its target year.

```bash
ISMIP7_RUN_TAG=n3 ISMIP7_LC=32000 antarctica/scripts/run_core_matrix.sh
CORES=1,2,9 antarctica/scripts/run_core_matrix.sh
```

Cores run sequentially, load-gated by `MAX_LOAD` (default cores minus 8). A
core already at its target year is skipped. Output predating the annual-mean
forcing fix is archived under `results/archive_stale_<stamp>/` and re-run.
The runner's own knobs (`CORES`, `MAX_LOAD`, `MAX_ATTEMPTS`, `FRESH`, `REUSE`,
`NRANKS`, `PROV_REF`) and the reuse rules are documented in its header.
The runner pins `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps` rather than
inheriting the development default; override it only with a solver mode that
has passed `make qualify` at the target configuration. Cluster forwards
(`batch_runners/projection.sbatch`) name their own solver, `scpc_gamg`: see
"Production configuration" in section 7.

Record each completed core with `python scripts/core_report.py --core <N>
--name <exp> --csv <timeseries.csv> --log <run.log>` (add `--ctrl-csv` for a
projection), run in that run's own shell so it captures the environment.
`--superseded "<reason>"` stamps a record when a later run replaces it.

### The relaxed initial state (`scripts/relaxation/run.py`)

An initial-state option beside the production MAP. The relaxation year rewinds
a MAP's geometry to 2014 with one year of the Smith dH/dt (the issue #117
backdating), runs it to 2015.0 at half the production step on OCX's 2014
forcing, with the apparent mass balance off and the front pinned, and an
inversion of 250 iterations re-fits the controls on the geometry it ends in.
The forwards then point `ISMIP7_INVERSION` at the relaxed MAP and start in
2003 as from any MAP. The rules are in `icepack2_tools/relaxation.py`, the
reasoning in `INVERSION_PRIORS.md` ("The relaxed re-inversion").

```bash
R=antarctica/scripts/batch_runners/submit.sh
MAP=$PWD/antarctica/results/reinvert_2km/final/<production MAP>.h5
SIZE="ISMIP7_FRICTION=<the MAP's law> ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_MESH=checkpoint"
# the year, on the MAP's own mesh: results/relax_<MAP stem>_2000_final.h5
relax=$($R projection --time <the whole year> ISMIP7_EXPERIMENT=relax ISMIP7_INVERSION=$MAP $SIZE | tail -n 1)
# the re-inversion from it, once it has finished
$R inversion --dependency afterok:${relax%%;*} --time <one link> $SIZE \
    ISMIP7_RASTER_SAMPLE=<the MAP's raster_sample> \
    ISMIP7_WARM_START=$PWD/antarctica/results/relax_<MAP stem>_2000_final.h5 \
    ISMIP7_MAXITER=250 ISMIP7_CHAIN_MAX=0 \
    ISMIP7_MAP_OUT=$PWD/antarctica/results/reinvert_2km/final/<MAP stem>_relax2014.h5
```

The re-inversion's objective settings come from the end state, which carries
the MAP's, and the strict handoff check holds them; give the run the
inversion knobs the MAP was made under, as for any chain link. That includes
`ISMIP7_RASTER_SAMPLE`, which `inversion.sbatch` sets to `vertex_front` unless
told otherwise: the sampling is an objective key, and the end state's geometry
was built under the MAP's recorded `raster_sample` (`vertex` for every MAP
inverted before issue #167), which the end state carries. Size `--time`
for one link from the MAP chain's own seconds per evaluation (about 260
evaluations for 250 iterations): `ISMIP7_MAXITER` counts per process, so a
second link would start a second 250. Give the relaxation's own submission a
`--time` that covers the whole year in one link: `afterok` waits on its first
link alone, and a re-inversion started on an unfinished end state is refused.
At 2 km the year's first solve from IU's final MAPs diverged under
`projection.sbatch`'s `scpc_gamg` (runlog
`inversion-2km-{budd,rc}-b20k-relax2014-year`); the relaxations queued on
8 October add `ISMIP7_MAP_CLIP=0` and set `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER`
to `scpc_mumps` for Budd (80 steps in 48 min) and `full_mumps` for RC, whose
`scpc_mumps` solves took 36 to 84 Newton iterations at up to 34 min a step.
A year longer than one link chains at the wall, so the re-inversion waits on
the year reaching 2015 (a trigger job) rather than on `afterok` of its first
link.

What each forward does with a relaxed MAP, recorded as `init_state` in every
checkpoint it writes:

| forward | geometry | controls | `init_state` |
|---|---|---|---|
| on the MAP's mesh | the relaxed geometry | the re-inverted ones | `relaxed` |
| on another mesh (1 km from a 2 km MAP) | that mesh's BedMachine sample | the re-inverted ones, transferred | `relaxed-controls` |

Either way the 2003 start backdates the geometry it starts from by 12 years.
The relaxed thickness stays on the MAP's mesh: carried across meshes it would
arrive as a DG0 staircase, which `../ADAPTIVE_MESH.md` measured driving the
thickness clamp from 119,000 to 256,000 Gt/yr within a few steps.

### Environment knobs (inversion)

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_MAP_OUT` | output path for the MAP, overriding the generated name. Use it for smoke tests and variants so a short run cannot replace a production MAP. A bare filename resolves under `mesh/`. A re-inversion from a relaxation's end state writes `<MAP stem>_relax<year>.h5` and stops under any other name | generated |
| `ISMIP7_MISFIT_NORM` | `sigma` divides each residual by its datum's squared error, giving a dimensionless chi^2; `none` is the legacy dimensional misfit. Selects the `ISMIP7_GAMMA_*` defaults | `sigma` |
| `ISMIP7_LOG_VEL_WEIGHT` | weight on the ISSM logarithmic velocity misfit (cost function 103). The chi^2 alone over-weights slow interior ice and leaves discharge-carrying tributaries 40 to 50% too slow; the log term is scale free. `auto` equalises it with the chi^2 term at the state the inversion starts from, and under `auto` a warm start that records a positive weight under the same `ISMIP7_MISFIT_NORM` and `ISMIP7_LOG_VEL_EPS` supplies that weight, so every link of a chain minimises one objective (issue 68). Stamped into the MAP with `log_vel_weight_source`: `requested`, `derived` or `warm_start` | `0` |
| `ISMIP7_LOG_VEL_EPS` | regularisation speed (m/yr) inside the log | `1.0` |
| `ISMIP7_WARM_START` | path to a MAP or timing-cache checkpoint used to seed `theta`/`phi` (and, when present, geometry, `fluidity_prior`, and the mixed diagnostic state). Fields are interpolated onto the live mesh, so a 1-core cache can warm-start a multi-rank invert. A relaxation's end state (`results/relax_*_final.h5`) is a warm start too: the run takes its geometry and keeps its `θ`, holds the MAP objective it carries, logs how far the relaxed geometry moved the friction anchor, and stamps the relaxed geometry's record into every checkpoint; it stops unless the state finished its year on this mesh and `ISMIP7_MAP_OUT` ends in `<MAP stem>_relax<year>.h5` | unset |
| `ISMIP7_WARM_START_THETA` | how `θ` comes over from the warm start. `θ` is a log-deviation from the friction anchor, so under a different anchor (`ISMIP7_ANCHOR_LENGTH`) the same `θ` is a different friction. `1` takes it as it is and warns when the anchors differ; `physical` rebases it on grounded ice so the friction `C_w0 exp(θ)` is the warm start's and the first solve reproduces the warm start's (from an exp-control MAP the previous anchor is its constant `friction_c_ref`); `0` starts at the new prior mean | `1` |
| `ISMIP7_WARM_START_PHI` | how `φ` comes over from the warm start. `1` takes it as it is, a deviation from this run's own fluidity prior; `physical` rebases it onto this run's prior (`ISMIP7_FLUIDITY_PRIOR`) so `A = A_prior exp(φ)` is the warm start's, e.g. a thermal-prior MAP warm-starting a Pattyn-prior run; `0` starts at the prior mean | `1` |
| `ISMIP7_WARM_START_GEOMETRY` | `1` takes thickness, bed, surface, `velocity_obs` and the mixed state from the warm start; `0` keeps this mesh's own BedMachine sample. Defaults to `1` on the same mesh, except that a MAP recording a different `lake_ice_base` or `raster_sample` does not supply its geometry, so an old MAP cannot bring the lake bowl back. `1` on a warm start recording another `raster_sample` stops the run: set `ISMIP7_RASTER_SAMPLE` to the warm start's | `1` on the same mesh, else `0` |
| `ISMIP7_WARM_START_STATE` | `1` loads the warm start's mixed state on its own mesh as the first guess when its geometry is not taken (a warm start sampled another way, issue #167), and the first evaluation's forward solves it at the full exponents instead of a cold ramp from n = 1. For RC's refit under `vertex_front` it did not work: from the ef2 state the residual norm started at 1.7e14 against the recorded 9.3e-3 and the first forward solve failed (jobs 11883148 and 11883149). `fluidity` loads instead the state of the MAP `ISMIP7_WARM_START_FLUIDITY` names, solved under the fluidity the run starts from, on that MAP's own mesh. RC's refit starts so, from Budd's refit: under Budd's fluidity RC's ef2 state started at ||F|| 2.5e21 and failed (job 11884749), Budd's state at 6.3e9, and its first forward solve converged in 75 Newton iterations (job 11884774). Refused on another mesh (`runconfig.warm_start_state`). A chain link resuming its own checkpoint drops it (`inversion.sbatch`) | `0` |
| `ISMIP7_WARM_START_FLUIDITY` | Path to another MAP: an inversion warm-started from `ISMIP7_WARM_START` takes `log_fluidity` and the `fluidity_prior` it deviates from out of this one, so the starting fluidity `A_prior exp(phi)` is that MAP's; `log_friction` and the state stay the warm start's. Fluidity does not depend on the friction law: RC's refit under the front-cell rule starts from the fluidity of Budd's refit on the same geometry (the two MAPs' fluidity priors agree to 6e-9; issue #167). The MAP's `fluidity_control` must match the run's. A chain link resuming its own checkpoint drops it (`inversion.sbatch`) | unset |
| `ISMIP7_WARM_START_FRONT_EXTEND` | `1`: a refit under a front sampling from a warm start sampled without it continues `log_friction` and `log_fluidity` harmonically over the nodes of the cells the front rule rebuilt or emptied (`geometry.front_changed_nodes`, `transfer.harmonic_extension`), and the refit sets them afresh. Every other node keeps the warm start's value, so over the band each becomes a harmonic blend of the ice upstream and the nodes seaward of the band; on a synthetic shelf the band came out as the linear interpolation between the two sides. The warm start fitted them against its own front: RC's final MAP left a band stiffer and with more friction than the rest of the ice (log fluidity -1.33 against -0.17, log friction +1.06 against -0.02), and the rule makes those cells about four times thicker (issue #167). On RC's refit it left the start worse: from the ef2 state as the first guess the residual began at ||F|| 4.1e14 with the band continued over 33,872 nodes, against 1.7e14 without it, and the first forward solve failed (job 11884485). That measurement includes the blend with the seaward nodes. RC's refit starts from the fluidity and state of Budd's refit (`ISMIP7_WARM_START_FLUIDITY`, `ISMIP7_WARM_START_STATE=fluidity`); Budd's refit started from its own controls. A `log_fluidity` taken from a MAP fitted under a front sampling (`ISMIP7_WARM_START_FLUIDITY`) is kept, and only `log_friction` is continued. The band is found by comparing this run's thickness with a vertex sample, so the warm start must be sampled with `vertex`; one sampled any other way (`cell_mean`, say) would flag interior cells and stops the run (`runconfig.front_band_extends`). The knob stays for a refit whose band is what blocks it | `0` |
| `ISMIP7_RAMP_SLIDE_FIXED` | `1` holds the sliding exponent at its target while the inversion's startup ramp climbs the flow exponent from 1 (`continuation.ramp_exponents`, `m_start`). IU tried it for RC's refit under `vertex_front` on the 2 km mesh (issue #167) and stopped it as too costly: after 2 h 29 min (job 11884073) the ramp stood at n = 1.13, its first step having taken three rungs (37, 34 and 30 min), and under `scpc_gamg` (job 11884074) its first linear solve failed as at n = m = 1. The refit starts from the fluidity and state of Budd's refit instead (`ISMIP7_WARM_START_FLUIDITY`, `ISMIP7_WARM_START_STATE=fluidity`). The ramp from n = m = 1 failed there: under full MUMPS (job 11869835) Newton ran 200 iterations at n = m = 1 in each of three rungs and diverged (residual norm 3.0e11 to 5.3e11), and under `scpc_gamg` and `scpc_mumps` (jobs 11883518 and 11883519) the first linear solve at n = m = 1 failed (1,000 Krylov iterations from a residual norm of 8.9e9) | `0` |
| `ISMIP7_ANCHOR_LENGTH` | reach (m) of the driving stress in the friction anchor `C_w0 = tau / max(|u_obs|, 1)^(1/m)`, the prior mean of the friction. `0` is the local balance, which vanishes with the surface slope and leaves ice divides with no friction in the prior. A positive length averages the driving-stress magnitude of the grounded ice over about that distance (a screened-Poisson filter with the second moment of a Gaussian of that standard deviation), so floating and ice-free cells neither add to nor dilute it. Stamped into the MAP as `friction_anchor_length`; a forward rebuilds the anchor from the MAP's value and aborts if this variable says otherwise | `0` |
| `ISMIP7_TRANSFER_FILL` | what the controls and the fluidity prior take on the dofs of the compute mesh beyond the mesh they were read from: a forward loading a MAP from another mesh, and an inversion's warm start. `extend` continues θ, φ and α harmonically from the source outline and the fluidity prior through its logarithm (`icepack2_tools/transfer.py`, `harmonic_extension`), so the bi-Laplacian prior pays almost nothing at the outline; `constant` fills θ = φ = 0 and the baseline prior `A0 * a4_factor`. Under `constant`, Rice's 2 km Budd MAP after stage 1 of issue #153 started on the 20 km buffered mesh at a smoothness cost of 1.46e5, against the 2.9e3 its last iterate recorded. A MAP records the mode of its warm start as `warm_start_fill` | `extend` |
| `ISMIP7_DRAG_GATE` | which water cells beside the ice the floor-cell ocean drag skips (`front.ocean_drag_cells`), in the inversion and the forward. `vertex` skips every cell touching the ice, at an edge or at one vertex; `facet` skips only those sharing an edge, so a cell touching the ice at one vertex drags the ice's own front node. On the 20 km buffered 2 km mesh RC's stage-1 controls gave floating ice 27 % below the observed speed under `facet`, 11 % under `vertex`, 2 % with `vertex` and a 1 m membrane floor, and 1 % with the drag off (issue #153). The MAP records `drag_gate` (`none` when no cell was dragged), and a forward runs the recorded gate and refuses a different one here; a MAP recording none, or older than the record, runs this knob, and a forward state older than the record restarts under `facet` unless it is set | `vertex` |
| `ISMIP7_EXACT_FRONT` | give every facet between ice and an ice-free cell inside the mesh the exact depth-integrated push `g (rho_I H^2 - rho_W d^2) / 2`. The DG0 facet driving stress matches it for floating ice and on flat land, falls short at a grounded marine cliff by `g D (rho_I H - rho_W D) / 2` (15 % for 1600 m of ice in 300 m of water), and on land is set mostly by the step in bed height: a median 20 times the push outward where the rock lies below the ice surface, 13 times back into the ice where it rises above (2 km buffered mesh, issue #153). `1` adds the difference (`dual_friction.front_cliff_correction` for cell-wise friction, icepack_tools' under `ISMIP7_SUBELEMENT_FRICTION=1`), which against rock above the ice surface drives the ice into the rock. `2` pushes with the face above the ice-free neighbour's bed `B` only, `g (rho_I h_e^2 - rho_W d_e^2) / 2` with `h_e = min(H, max(s - B, 0))` and `d_e = min(d, max(-B, 0))`: the free-cliff push where `B` is at or below the ice base, none where it is at or above the ice surface (issue #166); cell-wise friction only. An objective key; the MAP records the version as `exact_front` and the forward follows it, refusing an `ISMIP7_EXACT_FRONT` that differs. DG0 geometry only | `1` under the sub-element scheme, else `0` |
| `ISMIP7_LAKE_ICE_BASE` | under BedMachine's subglacial-lake mask (4, Lake Vostok) raise the bed to the ice base `s - H`. BedMachine's bed there is the lake floor, so `b + H` sits below its surface by the water column: a bowl a median 266 m and up to 916 m deep over 15,200 km2, with driving stresses near 1 MPa on its walls. Stamped into the MAP as `lake_ice_base`; a forward follows the MAP, so a MAP inverted without it keeps its own geometry | `1` |
| `ISMIP7_SKIP_CONTINUATION` | `1` skips the cold `n,m: 1→n` ramp on the initial solve, and makes each objective evaluation one taped solve at the full exponents in place of the direct forward (timing lanes). A warm start that supplies a full mixed state skips the initial ramp by itself and keeps the direct forward | `0` |
| `ISMIP7_DIRECT_FORWARD` | each objective evaluation is one untaped Newton solve at the full exponents from the last converged state, then a taped solve that starts converged, under every `ISMIP7_INVERSION_LINEAR_SOLVER` (`icepack2_tools/taped_solve.py`). A direct solve that fails goes to the failed-trial rescue, or to the line search's backtrack. `0` restores the taped 5-stage `n: 1→3` ladder in every evaluation, or one taped solve when the warm start supplied its mixed state or `ISMIP7_EVAL_CONTINUATION=0`. The MAP and the timing record carry the mode as `eval_mode` | `1` |
| `ISMIP7_DIRECT_FORWARD_MAXIT` / `_DTOL` | the direct solve's Newton iteration cap, and the residual growth (`snes_divergence_tolerance`) at which it is a lost trial | `30` / `1e6` |
| `ISMIP7_TRIAL_RESCUE_RUNGS` | rungs of the continuation ladder a failed line-search trial may climb before the trial counts as failed; `0` sends it straight to backtracking | `1` |
| `ISMIP7_EVAL_CONTINUATION` | with `ISMIP7_DIRECT_FORWARD=0` only: `0` keeps the initial `n,m: 1→n` ramp and then solves each annotated forward eval once at the full exponents, from the previous eval's state; a trial point where that solve fails takes the TAO path's re-ramp rescue. The objective is the same either way, and the MAP records the mode as `eval_continuation` | `1` |
| `ISMIP7_INVERSION_LINEAR_SOLVER` | linear solver of the annotated forwards, of the adjoint solves against them and of the publishing solve: `full_mumps` (the full mixed-Jacobian MUMPS LU), `scpc_mumps` or `scpc_gamg` (the transient's condensed modes; see "Inversion solver" in section 4). Recorded in the published (full-state) MAP as `state_solver_mode`, with its options as `state_solver_parameters`; the periodic checkpoints carry no solver record. It is outside the objective, so a chain may change it between links. The startup ramp follows `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER` | `scpc_gamg` (`full_mumps` until 3 October) |
| `ISMIP7_INVERSION_KSP_RTOL` / `ISMIP7_INVERSION_SNES_LINESEARCH` | `scpc_*` taped forwards only: outer Krylov relative tolerance and Newton line search. The direct forward inherits both through `direct_forward_parameters`; `full_mumps` keeps the shared `ISMIP7_SNES_*` settings | `1e-8` / `ISMIP7_SNES_LINESEARCH` (`nleqerr`) |
| `ISMIP7_GAMMA_THETA` / `ISMIP7_GAMMA_PHI` | Whittle-Matern prior strength on `θ` and `φ`, coupled to `ISMIP7_MISFIT_NORM` since normalising divides the misfit by about sigma^2. `scripts/lsurface.py` sweeps both on a grid and picks the L-curve corner of each (usage in its docstring) | `1e5` under `sigma`, `1e4` under `none` |
| `ISMIP7_L_REG` | prior correlation length (m) | `7.5e3` |
| `ISMIP7_MAXITER` | L-BFGS-B iteration cap | `500` |
| `ISMIP7_GRAD_PRECOND` | `none` is the raw-dof l2 metric, which is mesh dependent, so fine grounding-line cells converge slowest. `mass` optimises in `u = sqrt(M) x` under scipy, making the rate mesh independent. `mass_consistent` and `prior` run under TAO instead (scipy takes no preconditioner) with the initial inverse Hessian set to `M^-1` or to the prior covariance; `mass_consistent` is the consistent mass Riesz map; `prior`, the prior-preconditioned variant, was not used and is experimental here. `none` measured fastest: at 2 km it reached TAO's 60-iteration objective in 14 evaluations against TAO's 79 (`INVERSION_PRIORS.md`, issue #157 section), and the chains take it since 3 October | `none` |
| `ISMIP7_PRIOR_FORM` | `laplacian` uses `A = delta*M + gamma*K` as the prior precision; `bilaplacian` uses `A M^-1 A`, the squared-operator prior of Villa et al. (2021), the operator that a 2-D Whittle-Matern field needs to be function-valued. Different priors, not two spellings of one: their gammas are not convertible and their MAPs are not comparable, so the MAP stamps `prior_form` | `laplacian` |
| `ISMIP7_PRIOR_SIGMA_THETA` / `_PHI`, `ISMIP7_PRIOR_RHO` | `bilaplacian` only: the log-deviation scale and correlation length (m), converted to `(delta, gamma)` by the closed forms of Villa et al. (2021), `sigma^2 = 1/(4 pi gamma delta)`, `rho = sqrt(8 gamma/delta)`. The un-squared form has no such closed form, which is why its gamma can only be tuned | `0.3` / `0.3` / `ISMIP7_L_REG` |
| `ISMIP7_PRECOND_STEP0` | TAO metrics only: the largest change the FIRST step may make to a control, in that control's units, applied per control block. L-BFGS's first step is `-H_0 g` at unit length with no curvature pair to rescale it, and one evaluation outside the region where the forward has a solution returns NaN that every later trial point inherits | `0.15` |
| `ISMIP7_CHECKPOINT_EVERY_IT` | inversion checkpoint cadence: `ISMIP7_MAP_OUT` is rewritten every N accepted TAO iterates (every N evaluations on the scipy path). A chain link resumes from the last one, so a long interval repeats work after a wall-clock kill | `1` |
| `ISMIP7_GRAD_CHECK` | `1` runs a Taylor test of the taped objective's gradient at the start controls, prints the spread of the gradient and of the metric-scaled first step over grounded and floating nodes, and exits without writing a MAP or a `.done` marker. TAO metrics only (`ISMIP7_GRAD_PRECOND=mass_consistent` or `prior`); on a scipy metric the driver refuses it at startup | `0` |
| `ISMIP7_GRAD_CHECK_SEED` | the Taylor test's largest perturbation, relative to the controls' l-inf norm | `1e-3` |
| `ISMIP7_GTOL` | TAO metrics only: `tao_gatol` on the prior-metric gradient norm `sqrt(g' A^-1 g)`, which is mesh independent unlike the raw l2 norm the scipy path prints. `0` leaves stopping to `ISMIP7_FTOL` and `ISMIP7_MAXITER` | `0` |
| `ISMIP7_FTOL` | the relative-decrease stopping rule, on both optimizer paths: stop once `(J_old - J_new) / max(|J_old|, |J_new|, 1) <= ftol` (scipy's L-BFGS-B `ftol` is the same expression). `1e-10` converges a production inversion; the L-surface sweeps pass `1e-4`. `0` disables it | `1e-10` |
| `ISMIP7_MIN_ITER` | TAO metrics only: iterations before `ISMIP7_FTOL` may stop the run | `3` |
| `ISMIP7_SIGMA_U_FLOOR` | floor on the per-component MEaSUREs error (m/yr), so near-zero errors cannot let a few nodes dominate | `1.0` |
| `ISMIP7_SIGMA_U_UNOBS` | sigma (m/yr) where MEaSUREs reports no error. Those nodes carry a zero-filled `u_obs`, so they need a large sigma when `ISMIP7_OBS_MASK=0` | `1e4` |
| `ISMIP7_OBS_MASK` | `0` drops the velocity-observation mask | `1` |
| `ISMIP7_DHDT_WEIGHT` | weight on the dH/dt chi^2; `0` disables the transient constraint. Needs `ISMIP7_GEOMETRY_SPACE=dg0` | `0` |
| `ISMIP7_DHDT_SIGMA` | assumed dH/dt uncertainty (m/yr), hand set since the MIPkit ships no uncertainty field | `0.1` |
| `ISMIP7_DHDT_DT` | timestep of the prognostic step (yr) | `1.0` |
| `ISMIP7_DHDT_VAR` | `dhdt_smith` (firn corrected, 2003 to 2019) or `dhdt_cpom` (uncorrected, so not interchangeable) | `dhdt_smith` |
| `ISMIP7_DHDT_MELT` | `0` drops ocean melt from the prognostic source, which otherwise melts with the forward's calibration (section 5); melt is zero on grounded ice | `1` |
| `ISMIP7_DHDT_CLIM_START` / `_END` | RACMO climatology window for that source | `2003` / `2019` |
| `ISMIP7_DHDT_REACH` | pixel-to-cell reach as a multiple of `sqrt(area)`, rejecting pixels outside the mesh that nearest-centroid assignment would snap onto boundary cells | `0.75` |
| `ISMIP7_DHDT_NET_SIGMA` | sigma (Gt/yr) on the integrated grounded dH/dt; `0` disables the net term. Active only with `ISMIP7_DHDT_WEIGHT > 0` | `0` |
| `ISMIP7_OBS_KIT` | path to `AntarcticaObsISMIP7-v*.nc`. The kit is needed only to build the dH/dt cache rasters in `<ISMIP7_OBS_DATA_ROOT>/dhdt_cache/`; with those staged it may be absent. A path that does not exist is a hard error | newest under `<DATA_ROOT>/obs/mipkit` |

### Environment knobs (forward runs)

`icepack2_tools/runconfig.py` owns the run-shaping knobs (`ISMIP7_LC`,
`ISMIP7_LC_COARSE`, `ISMIP7_FRICTION`, `ISMIP7_GEOMETRY_SPACE`,
`ISMIP7_N_FLOW`), so an unset knob cannot mean one resolution to the inversion
and another to the preflight.

PETSc and recovery defaults have the same single-owner rule:
`icepack2_tools/solverconfig.py` supplies the runtime dictionaries, the timing
metadata, and the core-report provenance block. `simulation.py` must not
redeclare those literals.

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_LC` / `ISMIP7_LC_COARSE` | fine and coarse mesh resolution tags, selecting mesh and MAP | `1000` / `10000`, the production pair (`2500` / `64000` until 2026-09-19) |
| `ISMIP7_BUFFER_M` | outline buffer (m) the mesh is built with and named by (`runconfig.BUFFER_M_DEFAULT`, which the outline extraction and the mesh and sidecar names all read) | `20000` |
| `ISMIP7_MESH_FRONT` | the ice edge whose marine front the production mesh's nodes follow. `bm` names the `_frontbm` build (`mesh_antarctica.py --front bm`), whose nodes and edges lie on BedMachine's marine front with lc-sized cells along it; built and measured for issue #167 and kept as an option, since under `vertex_front` it carries the current meshes' front flux within 3 % at 6 % (2 km) and 8 % (1 km) more cells (run records `test-i167-frontbm-meshes`, `test-i167-front-probes`). Adaptation refuses a front mesh | `none` |
| `ISMIP7_MESH` | mesh path for the inversion and tools. A forward takes its mesh from the checkpoint unless this names another mesh, in which case the MAP is transferred onto it; a file with the checkpoint mesh's name and another triangulation is refused (`ISMIP7_MESH_BUILD_CHECK`). `checkpoint` means the mesh embedded in the MAP or restart file: `site_env.sh` always exports a derived path, so this is how a job submitted through `submit.sh projection` runs MAP-native. For the inversion, `checkpoint` is the mesh inside `ISMIP7_WARM_START`, under the `mesh_basename` that file records, so a MAP released without its .msh can be continued on its own mesh. Its recorded `lc`, `lc_coarse` and `buffer_m`, or values derived from a standard basename, supply the output MAP and timing provenance and the derived output names. A contradiction or incomplete mesh identity is refused | `mesh/antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.msh` |
| `ISMIP7_RASTER_SAMPLE` | how BedMachine lands on a DG0 cell. `vertex_front` (issue #167) is `vertex` with the marine front rebuilt from BedMachine's mask (`geometry.front_cells`): a cell holds ice when half its samples are BedMachine ice, an ice-free cell holding ocean is emptied, and the ice cells beside it take BedMachine's mean thickness and bed over their ice. Under `vertex` alone the cells the front crosses hold a fraction of the front's thickness on a buffered mesh (2 km front band 39.5 m against BedMachine's 163 m, 179 against 998 Gt/yr out of the front under `velocity_obs`). `vertex` projects the CG1 vertex interpolant; `cell_mean` takes the raster's true cell mean. `cell_mean` measured rougher: neighbouring cells share two of three vertex samples, so `vertex` damps jumps by construction. Cell means raised interior surface jumps 6% and bed and thickness jumps 35%, and at 2 km the momentum solve did not converge within 60 minutes. It does classify flotation better (32 km misclassification 9.1% to 3.2%), so the knob stays. Stamped into the MAP and read back by a forward on the MAP's mesh; a transfer rebuilds under the run's own. Each sampling melts with the calibration fitted under it (`runconfig.MELT_CALIBRATIONS`). Reproduce with `probe_raster_sampling.py` and `probe_front_cells.py` | `vertex_front` |
| `ISMIP7_INVERSION` | explicit MAP path for a forward or preflight. The forward checks the MAP's recorded `friction`, `n_flow` and `geometry_space` against the run and aborts on a mismatch, warning only when the MAP predates those attributes; `preflight.py` checks that the file exists. Use it to A/B MAPs on one mesh, or, with `ISMIP7_MESH` also set (the timing matrix, `make map-check`), to run a MAP on a different mesh: its continuous fields are then interpolated onto `ISMIP7_MESH` by strict point location (`icepack2_tools/transfer.py`), and a target dof outside the MAP's outline takes a stated fill (0 for the log controls, the constant baseline for the fluidity prior, the raster sample for `velocity_obs`), counted and printed as `Transfer fill:` lines (`MAP_CHECK.md`) | derived |
| `ISMIP7_MAP_CLIP` | bound on the absolute value of a MAP's log controls: a forward clips theta and phi to it when it loads a MAP and prints the count as `MAP clip:`, and an inversion clips its warm start's theta to it. `0` disables it. Forwards from IU's final 2 km MAPs (the exact_front version 2 Budd and RC MAPs on the 20 km buffered mesh, issues #153 and #166) run `0`, with `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER` `scpc_mumps` for Budd and `full_mumps` for RC, passed at submission, so the forward runs the controls the inversion fitted: RC's theta reaches 13.0, and 856 floating nodes of Budd's phi lie beyond 10. From Budd's version 1 final, `scpc_gamg` diverged under the clip from a backdated geometry, and without the backdate its step took about 25 times as long as one under `scpc_mumps`; from RC's version 2 final, `scpc_mumps` took 36 to 84 Newton iterations a solve and 13 to 34 min a step (run record `test-2km-budd-b20k-final-forward-diagnostics`). IU expects the regularization to set the solver a MAP needs, so Rice's more strongly regularized MAPs may run `scpc_mumps` or `scpc_gamg`; this is untested, and the membrane floor also moves the start-up cost (a 1 m floor cold start from RC's stage 2 MAP took 4.8 h under `scpc_gamg`, `test-2km-rc-b20k-vgate-floor1`) | `10`, or `6` at `ISMIP7_N_FLOW=4` |
| `ISMIP7_CALVING` | `none`, or a law in `icepack_tools.calving`: `fixed`, `velocity`, `thickness`, `vonmises`, `vonmises_strain`, `hfb` (see above) | `none` |
| `ISMIP7_CALVING_PARAMS` | the law's parameters, `key=value,key=value` (e.g. `sigma_max_fl=0.2,sigma_max_gr=1`), checked against the law at startup. Replaces `ISMIP7_CALVING_SIGMA_MAX_GROUNDED` / `_FLOATING`, which are refused when a law is set | the law's defaults |
| `ISMIP7_CALVING_MODULE` | a Python file whose `@icepack_tools.calving.register` laws join the registry before `ISMIP7_CALVING` is looked up | unset |
| `ISMIP7_FRACTURE` | `mask` applies the ISMIP7 collapse forcing to every floating cell it flags, booked as calving; `mask_front` only to the flagged cells open water has reached, so no hole opens behind a standing front (the two end-members of discussion #30). Masks exist for the SSPs only, so the control, historicals and OCX abort on either. Needs DG0. Every run prints its mode once (`Ice-shelf collapse forcing: ISMIP7_FRACTURE=...`, `none` included). Under a mask mode the timeseries columns `collapse_flagged_cells`, `collapse_removed_cells` and `collapse_held_cells` count the flagged floating cells, the ones the mode has emptied and the ones it leaves standing (always 0 under `mask`), all ranks summed; the budget lines, a closing log line and the core report repeat them | `none` |
| `ISMIP7_OCX_FORCING` | what core 11 runs on. `protocol` is the ISMIP7 OCX product (RACMO2.3p2-ERA SDBN1 `acabf`, expert-judgment ocean `tf`/`so`), and the run refuses to start without it. `stopgap` is RACMO2.4p1 actual-year SMB with the constant OI ocean climatology, what the core ran on before the product was readable here. K is fitted to the climatology, so read `check_melt_bound.py --ocx` first (discussion #48) | `protocol` |
| `ISMIP7_OCX_OCEAN` | which expert-judgment OCX ocean scenario to read: `main` (the core one), `cold`, `warm` or `vary`. A member other than `main` writes to `ocx_<member>` | `main` |
| `ISMIP7_SMB_ELEVATION_FEEDBACK` | the SMB-elevation feedback, `dacabfdz` times the surface change since the chain's initial state, added to the SMB every step (section 6). `1` or unset turns it on; `0` or an empty value turns it off, and the value set is closed. With it on a run refuses to start without the gradient it reads, and a resume refuses a checkpoint written in the other mode | `1` |
| `ISMIP7_OUTPUT` | `1` records the ISMIP7 yearly fields and scalars (`<exp>_<lc>_ismip7_annual_<year>.h5`, `<exp>_<lc>_ismip7_scalars.csv`), regridded afterwards by `write_ismip7_output.py`. The value set is closed, so a typo is rejected at startup. `projection.sbatch` and `run_core_matrix.sh` default it to `1`, and `=0` or an empty value turns it off there. A chained projection must export it on every link; a link that cold-starts mid-year logs the gap and begins at the next 1 January. Resuming continues a series, and a cold start into a populated series is refused | unset (`1` under the core-experiment runners) |
| `ISMIP7_BNDIDS` | boundary-id JSON override | per-mesh sidecar, else `mesh/boundary_ids.json` |
| `ISMIP7_GEOMETRY_SPACE` | `dg0` (one thickness for terminus force and mass flux) or `cg1` (legacy, A/B only). Selects the MAP. See `../GEOMETRY_DISCRETIZATION.md` | `dg0` |
| `ISMIP7_DATA_ROOT` | forcing tree root | `<repo>/ISMIP7/AIS` |
| `ISMIP7_OBS_DATA_ROOT` | BedMachine, MEaSUREs velocity, RACMO and the dH/dt cache observational-data root; name it in a site file when these files do not live beside the code. Also a write target: with the MIPkit present `obs_dhdt` builds `<root>/dhdt_cache/` here, so staged cache tifs belong under this root, wherever it points | `<repo>/antarctica/data` |
| `ISMIP7_T_END` / `ISMIP7_DT` | end time and timestep (yr). `t=Y.0` is 1 January of year Y, so a run covering 2015 to 2300 ends at `2301` and a historical covering 2003 to 2014 ends at `2015`. Each driver owns its end (historical `2015`, ssp370 `2101`, other projections and control `2301`, OCX `2026`). The step defaults to `runconfig.DT_DEFAULT`, the production step, which `projection.sbatch` also exports | driver's own / `0.025` |
| `ISMIP7_GEOMETRY_BACKDATE` | years of the Smith et al. (2020) mean dH/dt a cold start undoes on grounded ice before it runs (issue #117). Unset: `2015 - ISMIP7_T_START` for a start from 2003 up to 2015 (the historicals and OCX start in 2003), none from 2015 on, and a start before 2003 is refused. The friction anchors stay on the 2015 geometry; floating ice keeps its 2015 thickness. `0` turns it off | driver's start |
| `ISMIP7_RELAX_START` / `ISMIP7_RELAX_DT` / `ISMIP7_RELAX_FORCING` | the relaxation year of the relaxed initial state (`scripts/relaxation/run.py`): its start (from 2003, before 2015), its step, which has to divide the year, and its forcing, `ocx` (OCX's own for the year, as `ISMIP7_OCX_FORCING` selects) or `none` (no SMB or melt, a smoke test needing no forcing data). It refuses `ISMIP7_APPARENT_MB`, `ISMIP7_OUTPUT=1`, a calving law and an unpinned front, and `projection.sbatch` defaults the first two off for it | `2014` / half of `ISMIP7_DT` / `ocx` |
| `ISMIP7_FRICTION` | `budd`, `regularized_coulomb` or `budd_legacy`; selects the MAP. The set is closed, so a misspelling is rejected at startup | `budd` |
| `ISMIP7_OUTPUT_INTERVAL` | budget log line every N steps; the timeseries gets a row every step | `10` |
| `ISMIP7_CHECKPOINT_EVERY_YR` / `ISMIP7_KEEP_CHECKPOINTS` | checkpoint cadence in model years, and how many to keep besides `_final.h5` | `5` / `3` |
| `ISMIP7_RESTART` | restart checkpoint | `hist_<esm>[_<tag>]_<lc>_final.h5` if present; refused when it is short of the branch year |
| `ISMIP7_MESH_BUILD_CHECK` | refuse a MAP or restart whose mesh has the `ISMIP7_MESH` file's name and another triangulation (vertex and cell counts and two centroid sums, `transfer.meshes_match`), the way two sites' builds of the production mesh differ; `0` transfers across them on purpose | `1` |
| `ISMIP7_AUTO_RESUME` | resume from this experiment's newest checkpoint when no `ISMIP7_RESTART` is given. An integer flag, `=0` disables it, since the runners export it unconditionally and `--export=ALL` cannot unset. `projection.sbatch` refuses to chain when it is off | unset |
| `ISMIP7_RUN_TAG` | experiment-name suffix for a parallel method line | unset |
| `ISMIP7_WALL_STOP_MIN` | wall-clock budget in minutes from process start, checked before each step against the longest step so far, so the run writes its final checkpoint and exits with `t_yr` short of `t_end` for a chained job to resume. `projection.sbatch` derives it from the job's own TimeLimit, holding back `ISMIP7_WALL_MARGIN_MIN` minutes (25, sized for a 1 km final write; a one-hour 2 km link can pass less). `0` disables | `0` |
| `ISMIP7_EXPERIMENT_NAME` | the run's identity, used by `adapt_mesh.py` to name adapted meshes and sidecars so parallel experiments cannot overwrite each other. Set by `run_adaptive.py --experiment-name`. See `../ADAPTIVE_MESH.md` | unset |
| `ISMIP7_APPARENT_MB` | `1` or `balance` zeroes the t=0 thickness tendency; `div` cancels only the flux divergence; `0`, `off`, `none` and empty disable it | unset |
| `ISMIP7_FIXED_FRONT` | hold the calving front at the t=0 extent, tallying inflow beyond it as calving. `=0` disables. Ignored whenever an `ISMIP7_CALVING` law is configured | unset |
| `ISMIP7_FRONT_ADVANCE` | `free` (default): a level-set law's front advances wherever the transport carries ice into an empty cell the law does not remove; `none`: retreat-only, nothing advances past the t=0 extent, ice reaching it is removed each step and booked as calving (the ISMIP6 retreat-only scheme; the ISMIP7 submission's front). Under `none` a cell a front rule has emptied (beyond the law's front, shed whole, or a retreat sliver in a front cell) stays ice-free for the rest of the run, its inflow booked as calving; the mask is saved in checkpoints as `calved_cells` (`front.held_calved`, `tests/test_front_advance.py`). Needs `ISMIP7_CALVING`. |
| `ISMIP7_TRIPWIRE_U_MAX` / `ISMIP7_TRIPWIRE_H_MAX` / `ISMIP7_TRIPWIRE_DH_RATE` / `ISMIP7_TRIPWIRE_HMIN` | runaway tripwire: fail the step when max speed exceeds `U_MAX` [m/yr], max thickness exceeds `H_MAX` [m], or a cell that entered the step at least `HMIN` thick thickens at a relative rate `(dh/h)/dt` above `DH_RATE` [1/yr] (a rate so every dt scores the same physics alike; thinner cells are reported, never tripped: buffer cells fill by more than their own thickness); every step prints a `tripwire step-k:` line with the worst cells; unset = off (timing lanes export 2e4 / 5000 / 20 / 100) | _(unset)_ |
| `ISMIP7_LEGACY_TRANSPORT` | restore the pre-July-2026 CG-projection transport (needs `cg1`) | unset |
| `ISMIP7_SNES_TYPE` / `ISMIP7_SNES_MAXIT` | diagnostic Newton type and iteration cap | `newtonls` / `200` |
| `ISMIP7_SNES_LINESEARCH` | line search used by `newtonls` | `nleqerr` |
| `ISMIP7_SNES_RTOL` / `ISMIP7_SNES_ATOL` / `ISMIP7_SNES_STOL` | initial nonlinear relative, absolute and step tolerances. After the initial/restart solve, the absolute tolerance follows the self-scaled policy below | `1e-8` / `1e-50` / `0` |
| `ISMIP7_SNES_DIVERGENCE_TOL` | residual-growth divergence threshold; PETSc's `-3` (`PETSC_UNLIMITED`) disables this test (`-1` means `PETSC_DETERMINE`, restoring the default `1e4`) | `-3` |
| `ISMIP7_SNES_ATOL_SCALE` / `ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE` | persistent absolute tolerance after a converged setup solve (`scale * achieved norm`) / after accepting a loaded hard-era state (`scale * loaded-state norm`) | `100` / `1e-6` |
| `ISMIP7_SNES_KSP_EW` | enable PETSc Eisenstat-Walker variable inner tolerance for an A/B test | `0` |
| `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER` | `schur_gamg` or `schur_mumps`: legacy PETSc `selfp` approximation; `scpc_gamg` or `scpc_mumps`: exact cell-local Slate elimination and an assembled velocity solve; `full_mumps`: complete mixed-Jacobian reference. Legacy `iterative`/`mumps` aliases mean `schur_gamg`/`schur_mumps` | `scpc_gamg` for cluster forwards (`batch_runners/projection.sbatch`); `full_mumps` for a forward driver run by hand, the core runner, the workstation launchers and the inversion's startup ramp (its taped solves follow `ISMIP7_INVERSION_LINEAR_SOLVER`); the timing Makefile's `TIMING_SOLVER` (`scpc_mumps`) |
| `ISMIP7_MUMPS_ANALYSIS` | `parallel`: MUMPS factorizations of the `scpc_*` and `schur_mumps` velocity block use the distributed PT-Scotch analysis (ICNTL 28=2, 29=1) where PETSc was built with PT-Scotch, and MUMPS's sequential analysis on a build without it, where the distributed one makes the condensed preconditioner return NaN; `sequential` always uses MUMPS's own analysis and ordering. Any other value is an error. The cache fingerprint (`timing_campaign.solver_configuration_fingerprint`) normalises these keys, so a cache prepared under one analysis validates under the other | `parallel`; the inversion sets `sequential` in its own process unless the environment names one (its taped solves and direct forward under `ISMIP7_INVERSION_LINEAR_SOLVER=scpc_mumps`, its startup ramp under a `scpc_mumps` or `schur_mumps` lane solver) |
| `ISMIP7_FREEZE_LINEARIZATION` | `scpc_*` only (their Jacobian is matrix-free): build it on a copy of the state that is refreshed only when SNES re-forms the Jacobian, so the NLEQ-ERR line search's simplified-Newton solve sees the Jacobian SCPC condensed instead of one that has followed the state to the trial point. `0` restores the live state every lane before 2026-09-19 ran with; records carry `solver_configuration.linearization_state` (`frozen`/`live`/`assembled`) | `1` |
| `ISMIP7_KSP_RTOL` / `ISMIP7_KSP_MAXIT` | outer FGMRES relative tolerance / iteration limit for the iterative diagnostic mode | `1e-6` / `1000` |
| `ISMIP7_CONDENSED_KSP_TYPE` / `ISMIP7_CONDENSED_KSP_ATOL_FACTOR` / `ISMIP7_CONDENSED_KSP_RTOL` / `ISMIP7_CONDENSED_KSP_RESTART` | `scpc_gamg` only: the Krylov method that iterates on SCPC's assembled condensed velocity system around GAMG; its **absolute** tolerance as a fraction of `ISMIP7_KSP_RTOL` (FGMRES hands the preconditioner unit vectors and the elimination is exact, so the outer relative residual after one iteration is the inner absolute residual: this is the loosest inner solve that leaves the outer FGMRES one iteration); its relative tolerance, parked out of reach; and the GMRES restart. `preonly` restores one V-cycle per outer iteration | `fgmres` / `0.5` / `1e-12` / `100` |
| `ISMIP7_CONDENSED_NEAR_NULLSPACE` | `scpc_gamg` only: `none` leaves PETSc's default (the two translations); `rigid_body` adds the in-plane rotation of the condensed velocity space. Measured on 2500/25000 × 16: 2 % fewer V-cycles, each 14 % dearer | `none` |
| `ISMIP7_CONDENSED_PETSC_OPTIONS` | `scpc_gamg` only: further unprefixed `name=value` options (or bare flags) of the condensed solve for a tuning rung, e.g. `"pc_gamg_threshold=0.02 mg_levels_ksp_max_it=4"`; applied last, recorded in the lane's `diagnostic_petsc_options` and printed in its log | empty |
| `ISMIP7_SNES_MONITOR` / `ISMIP7_SNES_LOG` | enable diagnostic SNES/KSP monitors and optionally route them to a file; KSP output includes short and true residual lines per iteration, with a header before each mixed solve | `0` / stdout |
| `ISMIP7_SOLVER_VIEW` | emit `snes_view`, outer `ksp_view`, and the SCPC condensed `ksp_view`; enabled by `make debug` and `make reference` to expose the actual block sizes and hierarchy | `0` |
| `ISMIP7_TRANSPORT_KSP_RTOL` / `ISMIP7_TRANSPORT_KSP_MAXIT` | GMRES relative tolerance / iteration limit for the persistent DG0 transport solver (`ismip7_transport_` PETSc prefix) | `1e-10` / `500` |
| `ISMIP7_MASS_RESIDUAL_TOL_GT` | fail-loud absolute tolerance for both the discrete transport identity and the complete step mass budget | `5e-5` Gt |
| `ISMIP7_RESCUE_ENABLED` | permit a failed direct transient diagnostic solve to enter the continuation/trust-region/subcycle rescue ladder; set to `0` for strict timestep qualification | `1` |
| `ISMIP7_DELTAT_PER_BASIN_NPZ` | the melt calibration, one K and a TF offset per basin (`select_melt_parameters.py`, `calibrate_deltaT.py`); every ocean callback melts with its K. Refused with `ISMIP7_K_SCALE` other than 1, and when fitted under another slope or geometry than the run's (section 5) | the tracked file of the run's raster sampling (`runconfig.MELT_CALIBRATIONS`), `calibration/deltaT_per_basin_1000_K6.500e-05_vertex_front.npz` under the default |
| `ISMIP7_K_PER_BASIN_NPZ` | a legacy per-basin K file from `calibrate_melt.py`, read in place of the offsets; refused together with `ISMIP7_DELTAT_PER_BASIN_NPZ` | unset |
| `ISMIP7_K_SCALE` | multiplies a legacy per-basin K; refused with an offsets file | `1` |
| `ISMIP7_K_MELT` | removed with the tracked calibration, and refused when exported | |
| `ISMIP7_MELT_OBS_CSV` | per-basin melt observation table read by `scripts/calibrate_melt.py`; columns are located by header name, so either published table serves | `Melt_Paolo_Davison_Adusumilli_imbie2.csv` under `<DATA_ROOT>/meltobs/`, else under `<DATA_ROOT>/parameterisations/ocean/meltobs/`, else the older Paolo and Adusumilli table with a `[!]` line |
| `ISMIP7_MELT_SLOPE` | the draft slope the quadratic melt law sees, in the forward and in `scripts/calibrate_melt.py` and `scripts/calibrate_deltaT.py`: `ant` is one constant `sin(alpha)` on every shelf, the protocol's reference ("mean Antarctic slope, no slope dependency"); `local` is this mesh's draft slope. A K or deltaT file records the convention it was fitted under and a run under the other is told once | `ant` |
| `ISMIP7_SIN_ALPHA_ANT` | the constant under `ant`. The toolbox's K percentiles were sampled with 5.1117e-3, which its own gamma_T conversion gives and the notebook's slope recipe on the 8 km BedMap3 v3 topography reproduces; the default rounds it up by 0.065 percent | `5.115e-3` |
| `ISMIP7_SIN_ALPHA_CAP` | `local` slope only: cap on `sin(alpha)` in `scripts/calibrate_melt.py`; the forward applies none | none under `dg0`, `5e-3` under `cg1` |
| `ISMIP7_K_OUT` | output path for `scripts/calibrate_melt.py`, overriding the generated name. A bare filename resolves under `results/` | `results/calibrated_K_per_basin_<lc>.npz` |
| `ISMIP7_ESM` | ESM for the control | `CESM2-WACCM` |
| `ISMIP7_CLIM_SCENARIO` / `_START` / `_END` | reference-climate pool: the scenario pooled with `historical`, and the window, shared by the control's SMB climatology and the projections' aSMB re-reference through `icepack2_tools/climatology.py`. A partial pool warns | `ssp126` / `2000` / `2029` |
| `ISMIP7_H_CLAMP` | thickness floor (m) | `0` |
| `ISMIP7_NO_CALVING_TERMINUS` | drop the calving-terminus BC | unset |
| `ISMIP7_SUBCYCLES` / `ISMIP7_RESCUE_MAXIT` | dt-subcycle rescue ladder, and the Newton cap on its rungs | `1,4,16` / `600` |
| `ISMIP7_SUBSTEP_ADAPT` | `1` splits each macro step into `m` substeps chosen from a backward-Euler error estimate on the DG0 thickness, in place of the `ISMIP7_SUBCYCLES` retry list; forcing, output, checkpoints and budgets keep the macro step. The rules are in `icepack2_tools/substep.py`; regression test `tests/test_substep.py` | `0` |
| `ISMIP7_SUBSTEP_TOL` / `_INIT` / `_MAX` / `_QUIET` / `_HMIN` | the adaptive controller's tolerance (m), starting and largest substep count, quiet macro steps before `m` halves, and the thickness (m) below which a cell is left out of the estimate | `1.0` / `1` / `64` / `20` / `10` |
| `ISMIP7_FSSA_THETA` | weight of the free-surface stabilization of the lagged thickness-velocity coupling (`icepack2_tools/fssa.py`): the momentum balance carries the gravity load at the surface the step is about to produce, as a bulk viscosity on `div(u - u_ref)` with `u_ref` set by `ISMIP7_FSSA_REFERENCE`. `1`, the default since 6 October 2026, makes the lagged step stable at any size; `0` leaves the term out of the residual. `u_ref` and the step travel in checkpoints, and a restart from a checkpoint stepped without the stabilization keeps it off unless this knob is set (`solverconfig.forward_fssa_theta`); a prepared timing or map-check cache starts like a cold start. The forward logs the weight it steps with, and `core_report.py` lifts that line; regression test `tests/test_fssa.py` | `1` |
| `ISMIP7_FSSA_REFERENCE` | the velocity the stabilization measures the surface change from. `start`: the velocity the apparent mass balance was built from (the starting velocity when there is none), so t=0 reproduces the MAP velocity. `step`: the velocity of the last advance, with that advance's thickness tendency as a load (saved in checkpoints as `fssa_tendency`), so the load vanishes at every steady state of the transport. `auto`: `start` under `ISMIP7_APPARENT_MB`, `step` without it. A restart stops when the reference resolves to another one than its checkpoint was stepped with; a checkpoint that records only `fssa_tau` counts as `start`, and one whose `fssa_tau` is 0 was never stepped and loads under either | `auto` |
| `ISMIP7_SUBSTEP_TRACE` | `x,y;x,y` in km: after every transport advance, print thickness, source, flux divergence, front flags and calving rate for each cell within 2.3 km of each point | unset |
| `ISMIP7_STAGE_READS` / `ISMIP7_STAGE_WRITES` | a single-node job reads its MAP or restart from a node-local copy, and writes checkpoints and yearly output on node-local disk before renaming them into place (`icepack2_tools/staging.py`); over NFS a 2 km MAP took 45 min to read and a checkpoint 9 min to write on NOTS. `0` uses the networked path directly | `1` under Slurm, `0` off it |
| `ISMIP7_H_OCEAN` / `ISMIP7_K_LIM` | front backstops read by `scripts/simulation.py`: the thickness (m) at which the ice-free ocean drag ramps to zero, and the speed-limiter coefficient the rescue ladder raises for a rescue solve | `10.0` / `1e-3` |
| `ISMIP7_ALPHA_GL` | grounding-line coercivity, Budd only, read by `scripts/simulation.py` and `scripts/inversion_icepack2.py` | `0.5` (`0` for RC) |
| `ISMIP7_RC_HVISC_FLOOR` | membrane-only thickness floor (m) under both laws (`h_visc_floor`): a cell thinner than this carries this much ice in the membrane term alone. The first water row beside the ice carries it too, coupling the front to the ocean drag one cell out: with the vertex gate on the 20 km buffered 2 km mesh, floating ice started 11 % slow at 10 m, 4 % at 2.5 m and 2 % at 1 m, and 2.5 m costs 14 % more forward time at 1 km (GEOMETRY_DISCRETIZATION.md, issue #153). The inversion records it in the MAP as `h_visc_floor`; a forward runs the recorded floor and refuses a different one here, and a MAP older than the record runs this knob | `2.5` in the inversion; in a forward the MAP's, and `10.0` for a MAP older than the record |
| `ISMIP7_RC_CW0_FLOOR` | RC `C_w0` floor, read by `scripts/simulation.py` and `scripts/inversion_icepack2.py` | `0.0` |
| `ISMIP7_M_SLIDE` | sliding exponent, read by `scripts/inversion_icepack2.py`, `scripts/simulation.py`, `scripts/thermo_prior.py`, `scripts/plot_map.py`, `scripts/run_eigendec.py` | `3.0` |

> **dt guidance.** Production runs on the 1000 m mesh use `ISMIP7_DT=0.025`,
> the step the group chose with the mesh on 25 September 2026 (issue 20):
> `runconfig.DT_DEFAULT`, which `projection.sbatch` also exports. The timing
> matrix ran 0.05 there, and at 0.05 a 1 km control from a transferred 2 km
> Budd MAP diverged at the Lambert confluence (`MAP_CHECK.md`); at 0.025
> Rice's 1 km historicals ran 78 and 66 model years with no rescue step.
> At 2500 m and coarser use `0.1`; `0.25` is
> acceptable there when 10 steps per year is too costly, under the production
> closure only: under the matrix's strict contract `0.125` passed and `0.25`
> ran away at both 2500 m and 5000 m, so the stable step does not grow with
> the mesh. `dt=1.0` mis-melts per step and resurrects clamped cells.

---

## 7. Timing benchmark (`make timing`)

Resolution vs. core-count wall-clock benchmark for the transient solver. The
campaign is a staged state machine: it prepares an exact-mesh state, runs one
strict scout per mesh, and submits scaling lanes only after the corresponding
scout passes. Every job is submitted through `batch_runners/submit.sh`, so the
commands below are the same on any cluster with a site file (§0.5): the site
supplies the account and partitions, `SLURM_QUEUE=short|long|debug` picks a
partition by class, and `SLURM_PARTITION=` / `SLURM_CONSTRAINT=` name one
outright. Lanes are single-node everywhere; one that this site's nodes cannot
hold is recorded `not_runnable` (`exceeds_site_cores`, `exceeds_site_mem`) and
shown as such in the matrix, so matrices from different sites stay comparable.
Each record's `host` block names the site and node that measured it. Run from
`antarctica/` (needs §1 data, Firedrake, and Slurm):

```bash
cd antarctica
make timing
# after the currently active stage settles, run the same command again
make timing
# synchronize from the cluster with a trailing-slash rsync source
# (REMOTE_RESULTS=host:path/to/antarctica/results/; defaults to IU Quartz)
make sync-results
make matrix
# → TIMING_MATRIX.md, leading with wall time per simulated year and per
#   285-year projection, then 20 configured lanes plus NOT PLANNED cells
#   (MATRIX_OUTPUT=TIMING_MATRIX_QUARTZ_SCPC_MUMPS.md keeps one committed
#   file per site and solver)
# → results/timing/timing_<tag>_<LC>_<LC_coarse>_<ncores>.json
```

`make timing` first reuses a valid five-step qualification of the campaign's
solver (`TIMING_SOLVER`, default `scpc_mumps`; or runs it once), creates the
rank-independent improved inversion if needed, and
then advances only stages whose prerequisites already exist. It never
resubmits an active job or an accepted result. Failed lanes remain failed for
inspection; use `FORCE_TIMING=1` only for an intentional retry.

**The solver is a campaign parameter.** `TIMING_SOLVER=scpc_mumps|scpc_gamg`
names the diagnostic solver the lanes time: the same exact cell-local
condensation, with the condensed velocity system factored by MUMPS or
preconditioned by GAMG. It leads the campaign tag
(`scpc_gamg_10step_dt0p125at2500_…`), so the two campaigns keep separate
records, status stamps and matrices and neither report accepts the other's
lanes. The prepared caches are **not** per solver: they are always
`scpc_mumps` states (`timing_campaign.CACHE_SOLVER_MODE`, spelled into the
cache filenames), a lane checks its cache against the fingerprint of the
solver that *prepared* it rather than its own, and a GAMG lane therefore
starts from exactly the state the MUMPS lane of the same mesh started from
— the linear solver is the only difference between the two matrices. A GAMG
campaign on a site that already holds the caches needs no prepare or invert:

```bash
make qualify TIMING_SOLVER=scpc_gamg        # 2-step then 5-step gate, 2.5 km, 16 ranks
make timing-scout TIMING_SOLVER=scpc_gamg   # one scout per mesh, from the existing caches
make timing-scale TIMING_SOLVER=scpc_gamg   # once scouts have passed
make matrix TIMING_SOLVER=scpc_gamg MATRIX_OUTPUT=TIMING_MATRIX_QUARTZ_SCPC_GAMG.md
```

(`make timing TIMING_SOLVER=scpc_gamg` is the same thing as one re-runnable
command.) `make matrix` renders the campaign `make timing` last launched
unless the command line names one — the tag, or any parameter of it such as
`TIMING_SOLVER`. The matrix's *Solver work* table gives Newton iterations per
step × iterations on the condensed system per Newton iteration for every
accepted lane, from SCPC's own count: SNES's `linear_iterations` leaves out the
solve the NLEQ-ERR line search makes for its simplified Newton step, which is a
back-substitution under MUMPS and most of the Krylov work under GAMG.

**The line search must see the Jacobian that was condensed.** A matrix-free
Jacobian is the form's action at whatever the state Function holds, and
Firedrake writes every point the residual is evaluated at into it. NLEQ-ERR
evaluates the residual at its trial point and then solves, with the same KSP,
for the simplified Newton step `J(x_k)⁻¹F(x_trial)`: the operator had become
`J(x_trial)` while SCPC's condensed system was still `x_k`'s. Under exact
condensed MUMPS the Newton-step solves took one outer iteration and the
line-search solves 7–36 (1134 of a 2500/25000 × 16 lane's 1246 outer
iterations, each a mixed-Jacobian action and three Slate sweeps), and the step
the line search judged was not the one the method defines. Since 2026-09-19 the
`scpc_*` Jacobian is built on a copy of the state refreshed only when SNES
re-forms it (`preconditioners.frozen_linearization`; `make solver-smoke` holds
an exact condensed solve to one outer iteration per solve and checks the root
against the live one). `full_mumps` assembles its Jacobian and was always
frozen. Lanes accepted before that date ran live; `ISMIP7_FREEZE_LINEARIZATION=0`
reproduces them.

`scpc_gamg` iterates where iterations are cheap. The first configuration ran
one V-cycle per outer FGMRES iteration on the matrix-free mixed system, so each
iteration paid a mixed-Jacobian action and three Slate sweeps (0.19 s against
0.03 s on the assembled system). The mode now solves the assembled condensed
system with FGMRES + GAMG to an absolute tolerance just under the outer
relative one, which is exactly what leaves the outer FGMRES the single
iteration it has under MUMPS. Quartz, 2500/25000 × 16, frozen linearization,
identical 107-iteration Newton path and final state in every lane: single
V-cycle 79.0 s/step, inner solve to a relative 1e-7 27.6, to the absolute
tolerance 19.3 (19 V-cycles a solve), `scpc_mumps` 13.2. An exact coarse solve
changes nothing (19.1); point-block Jacobi smoothing and the rigid-body
near-nullspace lose. The knobs (`ISMIP7_CONDENSED_*`, table above) touch only
`scpc_gamg`'s options, never the `scpc_mumps` fingerprint the caches are held
to. Tune on one mesh with probe lanes, which leave the campaign's records
alone; each job's log names its condensed options and per-solve
`condensed_its`:

```bash
ISMIP7_CONDENSED_PETSC_OPTIONS="pc_gamg_threshold=0.02" \
  make timing-probe TIMING_SOLVER=scpc_gamg TIMING_ONLY_MESH=2500/25000 \
  SLURM_QUEUE=debug SLURM_TIME=01:00:00 FORCE_TIMING=1
```

**Production configuration.** The two Quartz matrices of 2026-09-19,
[`TIMING_MATRIX_QUARTZ_SCPC_MUMPS.md`](TIMING_MATRIX_QUARTZ_SCPC_MUMPS.md) and
[`TIMING_MATRIX_QUARTZ_SCPC_GAMG.md`](TIMING_MATRIX_QUARTZ_SCPC_GAMG.md), are
the same campaign from the same caches with every accepted lane under the
frozen linearization, so they differ in the linear solver alone. They set the
production forward configuration: the **1000 m / 10 km mesh, `scpc_gamg`, up to
64 ranks**, and the group chose the mesh for the submission on 25 September
2026 with a **step of `0.025` yr** (issue 20), half the matrix's 0.05 (see the
dt guidance above). Minutes of transient loop per simulated year on that mesh
at the matrix's 0.05 (and the 285-year extrapolation):

| 1000/10000 | 16 ranks | 32 ranks | 64 ranks |
|---|---|---|---|
| `scpc_mumps` | 34.0 (6.7 d) | 23.4 (4.6 d) | 23.2 (4.6 d) |
| `scpc_gamg` | 32.3 (6.4 d) | 15.7 (3.1 d) | 10.0 (47.5 h) |

At 0.025 a simulated year takes twice the steps. Rice measured 22.5 min per
simulated year on 32 Cascade Lake ranks under `scpc_gamg`, 62 simulated years
per 24 h job. The 64-rank Quartz figure at 0.025 is not yet measured; at the
matrix's 30 s a step it is about 20 min, about 4 days per 285 years.

The factorization stops scaling past 32 ranks and the V-cycles do not. The two
64-rank lanes take the same nonlinear path (129 Newton iterations and 286
condensed solves over the ten steps) and agree on the total outflux to 9e-13
relative, and the GAMG lane peaks at 1.0 GiB a rank against 1.5. At 2000 m
and coarser `scpc_mumps` is still the faster on 16 ranks, by 15 to 33 %, and
on 32 it is never faster: level at 2000/20000 and 2500/25000, 6 % slower at
2000/40000 and 16 % at 2500/50000. A coarse run on 16 ranks is the one case
for naming it.

Neither matrix has a 500 m timing, and that is not for want of trying: under
`scpc_mumps` both 32-rank lanes ran at dt 0.025 and tripped the runaway
tripwire at step 1, while the 64-rank pair was either not run or blocked by its
scout; under `scpc_gamg` none were run. 500 m is an open stability question,
not merely an untried one. Closed as icepack/ismip7#22, not planned for
September 2026.

`batch_runners/site_env.sh` names the mesh, `projection.sbatch` the solver and
the step, and each `sites/<name>.sh` the rank count (64 on Quartz).
`ISMIP7_DIAGNOSTIC_LINEAR_SOLVER` owns these transient lanes and the
inversion's startup ramp. `ISMIP7_INVERSION_LINEAR_SOLVER` separately selects
the inversion's taped forwards, adjoints and publishing solve, with
`scpc_gamg` as its default. Section 4 records the inversion measurements.

**Starting state.** The production MAP is inverted on the production mesh: a
1 km / 10 km Budd inversion, warm-started from the 2 km Budd snapshot 0241,
was running at Rice on 25 September (issue #24), on Rice's build of the
mesh (1,869,252 vertices), which is the submission mesh. A site with its own
build (IU's has 1,869,088 vertices) uses Rice's `.msh` (release
`maps-2km-snap-2026-09-24`, md5 `5d318c0a`) with the MAP: a forward
interpolates a MAP onto whatever `ISMIP7_MESH` names, and it refuses a file
with the MAP mesh's name and another triangulation (`ISMIP7_MESH_BUILD_CHECK=0`
allows that transfer on purpose). On Quartz the default name has held Rice's
build since 26 September, and IU's is kept as
`mesh/antarctica_10000_1000_buffered20000_iubuild.msh`.

Until that MAP exists a forward starts from a coarse MAP by transfer, the way
the matrix's own lanes do. Name the MAP and let the forward interpolate it onto
the mesh `site_env.sh` exports: `simulation.py` keeps the MAP's own mesh as the
interpolation source and uses `ISMIP7_MESH` only for the target spaces, which is
the same path the 500 m lanes take from the 2.5 km MAP. A target dof outside
the MAP's mesh takes a stated fill and is counted (`MAP_CHECK.md`).

```bash
MAP=$ISMIP7_REPO/antarctica/results/timing/inversion
MAP=$MAP/inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5
submit.sh projection ISMIP7_EXPERIMENT=control \
  ISMIP7_FRICTION=budd ISMIP7_INVERSION=$MAP
```

`make map-check MAP_CHECK_FRICTION=regularized_coulomb|budd` runs one released
2 km MAP through that transfer and its checks, on the MAP's own mesh and on
this one, and prints where each stage stands (`MAP_CHECK.md`).

`ISMIP7_FRICTION` has to come with it: the campaign source is a Budd MAP, this
section's default is regularized Coulomb, and the forward aborts on a MAP whose
recorded law disagrees with the run. Without `ISMIP7_INVERSION` the runner falls
back to `ISMIP7_MAP_DEFAULT`, which names an RC MAP at the run's own resolution
that has never been inverted, and warns at submission that the file is absent.
The MAP on the production mesh, and closing the Budd/RC gap, are both open
under the inversions; icepack/ismip7#21 was closed as a duplicate on
22 September. (issue #24)

The stages and contracts are:

1. **Inversion provenance** — the only timing input is the campaign source
   MAP (`TIMING_INVERSION`). Since 2026‑09‑17 that is
   `results/timing/inversion/inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5`,
   the 250‑iteration re‑inversion of the imported `dg0_logvelnet` MAP on the
   2500/25000 campaign mesh under the HAF‑gated Budd law: the invert stage's
   own output for that mesh, so the manager never re‑inverts or re‑prepares
   the source mesh (its cache is the MAP itself, as the invert job published
   it). The imported old‑gate MAP
   (`mesh/inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5`) remains
   the solver‑qualification input only: under the fixed law its transferred
   controls run away at step 3 on every mesh. Cache manifests record the
   source checksum and lane records its basename (`initial_state_source`),
   so a cache or a current‑campaign record from any other source is never
   reused. The redistribution launcher has distinct input/output arguments;
   it never rewrites the source in place.
2. **Prepared caches (`make timing-prepare`)** — one job for each of the ten
   `(LC, LC_coarse)` meshes constructs `bed` and `thickness` as cell averages
   of BedMachine on that exact target mesh, recomputes the hydrostatic surface,
   performs the adaptive cold continuation, and saves the complete mixed state
   at 2015.0 without a transport step. Only continuous inversion products
   (controls, fluidity prior, and observed velocity) are transferred from the
   imported MAP; direct DG0-to-DG0 interpolation is forbidden because it
   aliases source cells at target centroids into artificial surface jumps. It
   includes geometry, velocity, membrane and basal stress, controls, the exact
   velocity observation field, physical priors,
   frozen reference fields, mesh identity, solver configuration, and source
   inversion checksum. Cache schema v3 requires an explicit inventory of these
   fields, including `velocity_obs`, plus the target-geometry source and
   construction method. The job then repacks the cache on one rank and
   publishes the HDF5 file and JSON manifest atomically. Mesh, inversion,
   physics, solver, or cache-schema changes invalidate it. The cache tag and
   the campaign tag moved from `v3` to `v4` on 2026‑09‑16 when the Budd shelf
   gate was fixed (height above flotation instead of the sign of roundoff
   `N`, ported from hoffmaao/antarctica e602705): the manifest records
   `friction_gate: haf` and validation requires it, so a v3 cache is never
   reused, and v3 records are the old-law archive.
3. **Per-mesh short inversion (`make timing-inversion`, optional)** — off by
   default since 2026‑09‑15: `make timing` runs lanes from the transferred
   prepare state (`TIMING_INITIAL_STATE=prepare`); `TIMING_INVERT=1` or
   `TIMING_INITIAL_STATE=invert` re-inverts each mesh first and records those
   lanes under the `…_reinverted` tag. Evidence for the default: every
   strict-lane runaway sits in a floating margin cell, which the
   grounded-only dH/dt misfit never sees, so 5- and 250-iteration inverts
   moved the runaway (Amundsen → Pine Island) without curing it. When it
   runs: after a valid prepared cache exists, one L-BFGS job per mesh
   re-inverts with the
   log-velocity + dH/dt + net-balance objective
   (`ISMIP7_LOG_VEL_WEIGHT=auto`, `ISMIP7_DHDT_WEIGHT=1`,
   `ISMIP7_DHDT_NET_SIGMA=10`). Default length is
   `TIMING_INVERSION_MAXITER=250` (override with e.g. `=5` for a debug pass);
   walltime defaults to `TIMING_INVERSION_TIME=24:00:00`. The job warm-starts
   from the pristine copy of that mesh's prepared cache
   (`…_buffered20000.prepare.h5`: controls, fluidity prior, geometry, and
   mixed diagnostic state) and skips the cold `1→n` continuation. Ranks are 32
   for LC &lt; 2500 m and 16 otherwise; memory follows `INVERSION_MEMORY_BY_LC`,
   sized for the full mixed-Jacobian MUMPS factorisation plus the adjoint
   tape, which leaves room under `scpc_gamg`, the default since 3 October. **The 500 m meshes are not
   re-inverted** (their invert needs ≈430 GB on one node); their lanes start
   from the prepared cache, i.e. the transferred 2.5 km MAP, and
   `TIMING_MATRIX.md` says so per mesh under "Initial states". After L-BFGS
   the job publishes the mixed state at the returned controls with a
   self-scaled absolute tolerance (`snes_atol = ISMIP7_SNES_ATOL_SCALE ×` the
   residual the last accepted forward reached, bounded at
   `ISMIP7_FINAL_SNES_MAXIT=50` iterations, reason always printed): when the
   returned controls are the last evaluated point that solve confirms the
   state at iteration 0 instead of grinding a converged residual against a
   relative test. The parallel MAP is rewritten on one rank to
   `results/timing/inversion/…_{lc}_{lc_coarse}_{N}iter.h5` with the full
   mixed diagnostic state, a profiling JSON (including the final-solve
   outcome) is written under `results/timing/`, then that checkpoint is
   published as the timing cache (no second cold prepare) so scout/scale
   provenance points at the short invert. The cache manifest keeps two
   solver facts apart: `diagnostic_solver_mode` (`scpc_mumps`, the mode every
   campaign cache is held to, whatever solver the lanes time) and
   `state_solver` (the selected inversion mode that actually produced the
   state). The inversion applies the forward's floor-cell stabilizers
   (`ISMIP7_OCEAN_DRAG`, `ISMIP7_H_OCEAN`, `ISMIP7_U_LIM`, owned by
   `runconfig.residual_stabilizers`), so its mixed state is a solution of the
   residual the forward assembles at restart — before 2026‑09‑14 it was not
   (`||F||` 1e1 in the inversion vs 1e10 in the forward on the same state).
   Every state checkpoint, prepared cache or inversion MAP, records
   `full_state_residual`; a restart accepts a mixed state without a solve
   only when its residual is within `ISMIP7_SNES_ATOL_SCALE ×` that record,
   and otherwise re-solves from the loaded guess (bounded, step-size exit
   live) before the first step. `FOLLOW_PREPARE=1` queues each invert behind
   its active prepare job; a running invert is never resubmitted, even with
   `FORCE_TIMING=1`.
4. **Cache audit / contract probe** — optional diagnostics on a prepared cache:
   ```console
   make timing-cache-audit TIMING_ONLY_MESH=2500/25000 \
     SLURM_QUEUE=debug SLURM_TIME=00:15:00
   ```

   This submits a read-only assembly of the exact initial DG0 upwind
   `div(h*u)` in one cache. It writes a JSON record under `results/timing/`
   with global extrema, the fixed-front-masked equivalent, and a ranked
   **hotspot table** (`--top`, default 20) scored by `|a_ref|·dt/h` — the
   fraction of a cell's thickness the no-forcing tendency would move in one
   step — with grounded/floating/buffer flags; it performs no diagnostic or
   transport solve. On 2500/25000 the top rows are the cells where every
   strict lane has run away.

   To run one mesh under a different physics contract without touching the
   campaign records (the experiment ladder):

   ```console
   make timing-probe TIMING_ONLY_MESH=2500/25000 TIMING_CONTRACT=divfront \
     TIMING_SCOUT_MONITOR=1 SLURM_QUEUE=debug SLURM_TIME=01:00:00
   ```

   Contracts: `strict` (no apparent MB, no calving sink — the campaign
   default), `div` (the exact uncapped `ISMIP7_APPARENT_MB=div` correction
   built from the pristine cached state; `make timing-amb-probe` is its
   alias), `front` (`ISMIP7_FIXED_FRONT=1`: ice advected beyond the 2015
   extent is removed each step and tallied as calving), `divfront` (both —
   the production closure). Probe records carry the `…_probe` suffix. Under
   `div` and no forcing the corrected state is stationary away from
   positivity-limited initially ice-free cells, so a passing probe diagnoses
   the runaway but is not a representative timing measurement.

5. **Scouts (`make timing-scout`)** — one lowest-retained-core lane per mesh:
   32 ranks at 500 m and 16 ranks elsewhere. A scout passes only after five
   direct steps, the exact final year, no negative solve reason or rescue
   label, matching cache provenance (the per-mesh invert MAP, or the
   transferred source MAP at 500 m), and transport and complete-step mass
   residuals no larger than `5e-5 Gt`.
6. **Scaling (`make timing-scale`)** — higher-core lanes are submitted only
   for meshes whose scout passes. A failed scout stamps its scaling lanes
   `BLOCKED BY SCOUT`; the jobs are not submitted. `FOLLOW_INVERT=1` lets
   scouts and scaling lanes queue behind an active invert job
   (`--dependency=afterok`). `TIMING_INITIAL_STATE=prepare` runs **control
   lanes** from the transferred prepare cache with no invert gate (the cache
   must have been published by `make timing-prepare`, not by an invert);
   their records and stamps carry the `…_transferred` tag, so they never
   count as campaign lanes, and `make matrix TIMING_TAG=<campaign tag>_transferred`
   renders them separately. Inversion records written before the publish
   gate (2026‑09‑14) are rejected; re-run the invert with `FORCE_TIMING=1`.
7. **Strict transient timing** — all matrix lanes use the campaign's
   `TIMING_SOLVER` (a lane whose tag and solver disagree refuses to start),
   disable rescue, and restrict subcycles to `1`. The interval is
   `MATRIX_STEPS` steps of `MATRIX_DT_2500 × LC / 2500` years (defaults 10 and
   0.125 since the 2026-09-16 dt ladder: 5 × 0.25 runs away at step 2 while
   10 × 0.125 and 20 × 0.0625 complete the 1.25 yr window with flat speed and
   thickness; both are written into the campaign tag after the solver, e.g.
   `scpc_mumps_10step_dt0p125at2500_…`, so a dt-ladder rung never mixes with
   another) and the physics contract is `TIMING_CONTRACT` (default `strict`;
   a non-strict contract suffixes the tag). The scaled step is capped at an
   absolute `MATRIX_DT_MAX` = 0.125 yr (2026-09-17: both 5000 m scouts ran
   away at dt 0.25 just as 2500/25000 does, a near-flotation cell flipping
   between grounded and floating every step, so the stable step does not grow
   with the mesh); the 5000 m lanes therefore cover 1.25 yr, not 2.5. The cap
   is a constant, not part of the tag, and applies from campaign v4 on. Every lane runs the runaway
   **tripwire**: the first step whose maximum speed exceeds
   `ISMIP7_TRIPWIRE_U_MAX` (2e4 m/yr), whose maximum thickness exceeds
   `ISMIP7_TRIPWIRE_H_MAX` (5000 m), or in which a cell at least
   `ISMIP7_TRIPWIRE_HMIN` (100 m) thick thickens at a relative rate
   `(dh/h)/dt` above `ISMIP7_TRIPWIRE_DH_RATE` (20/yr; a rate rather than a
   per-step fraction, so each rung of a dt ladder is judged alike) fails at
   once with the cell's coordinates, entry/exit thickness and
   grounded/floating/buffer flags (`RUNAWAY TRIPWIRE step-k: …`, category
   `runaway_tripwire`), instead of three steps later when the transport
   budget finally breaks.
   Thinner cells never trip (a 0.3 m buffer cell filling at 70 m/yr is not
   a runaway) but every step's `tripwire step-k:` line names the worst
   relative and absolute thickness change so a seed is visible before it
   trips. The first diagnostic, transport, tripwire or
   mass-budget failure ends the lane. The primary timer starts immediately
   before the step loop, after cache loading, solver construction, and
   transport setup; `setup_seconds` is reported separately. JSON records are
   atomically written on catchable success or failure and include solver
   histories, global mesh counts, extrema, residuals, phase, and completed
   steps. Slurm-only termination (including exit 137) remains distinct from a
   recorded numerical failure. The manager pins
   `ISMIP7_SNES_DIVERGENCE_TOL=-3`, PETSc's `PETSC_UNLIMITED` sentinel, so the
   solver reaches its real convergence or iteration-limit result; the legacy
   value `-1` means `PETSC_DETERMINE` and accidentally restores the default
   `1e4` growth cutoff.
8. **Reporting and synchronization** — `make sync-results` uses
   `quartz:/N/project/ice_rheology/ISMIP7/antarctica/results/` and the local
   canonical `results/` with trailing slashes so the directory contents merge
   at the correct level. Override `QUARTZ_RESULTS` if the remote checkout
   moves. `make matrix` refuses to run if `antarctica/results/` is nested
   beneath this directory. It labels numerical failures, OOM/external exits,
   incomplete output, scout-blocked lanes, and `NOT PLANNED` separately.

The selected matrix retains both coarse-mesh ratios and contains 20 lanes:

| Fine LC | Retained ranks per mesh |
|---:|---|
| 500 m | 32, 64 |
| 1000 m | 16, 32, 64 |
| 2000 m | 16, 32 |
| 2500 m | 16, 32 |
| 5000 m | 16 |

This is intentionally aggressive. Preliminary runs showed that 64 ranks gave
only about 7–20% improvement over 32 on 2000–5000 m meshes; the
5000/100000 case showed none, while 16→32 remained useful at 2000–2500 m.
The slow/OOM-prone 500 m 16-rank lanes and both 5000 m higher-rank lanes are
therefore omitted. Both coarse ratios remain because their observed scaling
differs. Timesteps remain resolution-scaled from the measured safe 2.5 km
value: 0.05, 0.10, 0.20, 0.25, and 0.50 yr for 500–5000 m, respectively.
Revisit them only if strict scouts using the corrected inversion and caches
fail.

Individual stages can be run separately: `make meshes`, `make redistribute`,
`make qualify`, `make timestep-probe`, `make timing-prepare`,
`make timing-inversion`, `make timing-cache-audit`, `make timing-amb-probe`,
`make timing-scout`,
`make timing-scale`, `make sync-results`, and `make matrix`.
`make timing-dry-run` prints the staged commands without submitting jobs.
`make map-check` and `make map-check-dry-run` are the same shape for one
released MAP (`MAP_CHECK.md`).
`make transient` routes matrix work through the scout gate;
qualification/debug calls still use its direct compatibility path.
Use `TIMING_ONLY_MESH=2500/25000` to prepare, invert, audit, probe, or scout one exact
mesh, and add `TIMING_SCOUT_MONITOR=1` to write the scout/probe SNES/KSP
monitor under `results/logs/`. Before treating a live-looking status as active,
the campaign
manager checks both `squeue` and `sacct`; a cancelled, timed-out, OOM, or
otherwise terminated allocation is recorded as an external failure instead of
being silently resubmitted.
`make matrix-legacy` regenerates the archived 30-lane report at
`TIMING_MATRIX_scpc_mumps_5step_dt0p25at2500.md`, using the copied Slurm logs
to retain the distinction between numerical failures and incomplete output.
`make solver-smoke` is the cheap setup check: on a 2×2 mesh and two MPI ranks
it exercises the real mixed field shapes, all new condensation/MUMPS modes,
and two coefficient updates through the persistent transport solver. It does
not replace the Antarctic qualification probes.

For a one-hour probe on the Slurm debug queue, use `make debug`. It submits one
transient job on the same `antarctica_64000_2500.msh` / DG0-logvelnet MAP pair
with `CORES=16`, using `--partition=debug --time=01:00:00 --mem=64G`. It uses
exact Slate local condensation plus MUMPS/PT-Scotch on the velocity system,
enables diagnostic SNES/KSP monitoring, and writes
`results/logs/timing_scpc_mumps_debug_lc2500_lcc64000_n<cores>_<stamp>.log`.
`make reference` uses that same input pair with the memory-heavy full
mixed-Jacobian MUMPS baseline. Both are diagnostic probes, not substitutes for
`make qualify`. Set `CORES=8` for an 8-core local probe.

Slurm scripts live in `scripts/batch_runners/` (`timing_prepare.script`,
`timing_redistribute.script`, and `timing_transient.script`), following the
Quartz module-load pattern used by the other batch scripts.

---

## Outputs

Per experiment in `results/`:

- `<exp>_final.h5`, a self-contained restart checkpoint: mesh, geometry,
  inversion fields, the full `(u, M, τ)` state, the frozen apparent-MB
  reference, and `levelset` when a calving law is configured. Under `dg0` the
  saved `thickness` is the prognostic state; a `cg1` run also saves
  `thickness_dg`. The `geometry_space` and `mesh_basename` attributes let a
  restart resolve the same sidecar, and `dt_yr` records the step. A resume
  that continues its own series at another step, as a run taken past a crash
  may need, prints a `WARNING: step change on resume` line and continues.
- `<exp>_t<year>.h5`, periodic checkpoints, keeping the
  `ISMIP7_KEEP_CHECKPOINTS` most recently written.
- `<exp>_timeseries.csv`, one row per step, the year at six decimals:
  `year, vaf_mm_sle, mass_gt, smb_gtyr, melt_gtyr, outflux_gtyr, calv_gt,
  clamp_gt, resid_gt, amb_gtyr`. The residual must close to 0.00. SMB and
  melt are what the advances applied: no forcing acts on open ocean or on
  cells a front rule holds ice-free (`front.unforced_cells`), so `clamp` is
  only the positivity limit on thin ice and `calv` only ice that crossed the
  front. How the ISMIP7 fields book this forcing, the outflux and the
  apparent-MB reference is in `ISMIP7_README_AIS_RICE_icepack2.md`.

VAF is in mm of sea-level equivalent, mass in Gt, both over map-plane area.
The ISMIP7 scalars of a run with `ISMIP7_OUTPUT=1`
(`<exp>_<lc>_ismip7_scalars.csv`) integrate over true area, map-plane area
times af2 = (1/k)^2 of EPSG:3031, as the organisers' tool does, so their mass
sits about 2.6 % above `mass_gt`.

### From a finished run to a submission

The chain has been run end to end on a workstation, with a two-year 32 km
control standing in for a production run (22 September 2026):

```bash
# 1. the run itself banks the yearly fields and the native scalars
ISMIP7_OUTPUT=1 mpiexec -n 8 python antarctica/scripts/control/run.py

# 2. onto the 8 km grid, in the request's files and names
python antarctica/scripts/write_ismip7_output.py \
    antarctica/results/<exp>_<lc>_ismip7_annual.h5 --out-dir submission \
    --esm CESM2-WACCM --scenario ctrl --exp C009 \
    --source-id RICE --ism-id icepack2 \
    --scalars antarctica/results/<exp>_<lc>_ismip7_scalars.csv

# 3. the compliance checker over what came out
python -m isschecker --variable-list ismip7 \
    --source-path submission/AIS/RICE/icepack2/CORE/C009

# 4. the model's densities into the upload, then the sea-level scalars
#    nothing here computes, from the tools' own venv
ismip7-scalars-set-params --region AIS --group RICE --model icepack2 \
    --rhoi 917 --rhow 1024 --rhof 1000 --modelpath submission/AIS
python antarctica/scripts/download_forcing.py --scalar-processing
ismip7-scalars --region AIS --group RICE --model icepack2 \
    --experiment ctrl --modelid m001 --esm CESM2-WACCM --forcingid f001 \
    --configid C009 --exp-group CORE --hist ctrl --refyear 2016 \
    --datapath ISMIP7/Output-Processing/Data/AIS \
    --modelpath submission/AIS --outpath scalars

# 5. the tool's scalars against the model's own
python antarctica/scripts/compare_scalars.py \
    --submission submission/AIS/RICE/icepack2/CORE/C009 \
    --tool scalars/nc/AIS/RICE/icepack2/CORE/C009 \
    --datapath ISMIP7/Output-Processing/Data/AIS \
    --params submission/AIS/RICE/icepack2/params.nc --refyear 2016 --native-af2 \
    --native-csv antarctica/results/<exp>_<lc>_ismip7_scalars.csv \
    --out-csv scalars/comparison.csv --out-md scalars/comparison.md
```

On a cluster, steps 4 and 5 are one serial job,
`scripts/batch_runners/scalar_processing.script`, whose header carries the
`submit.sh` line and the knobs.

The tools' venv, from any Python 3.11 to 3.14 (on Quartz, the
`python/3.14.5` module). Neither package is on PyPI, so both come from their
tags:

```bash
python3 -m venv <venv>
<venv>/bin/pip install \
    "isschecker @ git+https://github.com/ismip/ISM_SimulationChecker@0.5.1" \
    "ismip7-scalars @ git+https://github.com/ismip/ismip7-scalar-processing@3f36eb3"
```

Things that are easy to get wrong:

- **Both tools need Python 3.11 to 3.14**, in an interpreter of their own:
  the workstation's Firedrake environment is 3.10, and on a cluster the
  Firedrake environment's `PYTHONPATH` would shadow the venv's packages (the
  job script scrubs it). Where conda is unavailable, `nix` provides one, and
  the nix interpreter then needs the shared libraries it cannot see: gcc's C++
  runtime, zlib, expat and udunits, on `LD_LIBRARY_PATH`, with
  `UDUNITS2_XML_PATH` set. On Quartz every dependency of the checker installs
  as a wheel, and conda-forge carries its release too. Put the venv and pip's
  cache on scratch, since the home file quota is small.
- **The scalar tool needs four auxiliary grids** per region (the area factor,
  the extended Rignot basins, the glacier and ice-cap area factor and the
  maximum-extent mask), which live on Globus under
  `/ISMIP7/Output-Processing/Data` rather than with the forcing.
  `--scalar-processing` fetches them; `--scalar-resolution` picks the grid.
  Copying `/ISMIP7/Output-Processing/Data` with the Globus web app into
  `ISMIP7/Output-Processing/Data/` of the checkout, the same path, needs no
  CLI login (IU Quartz, 23 September 2026).
- **`params.nc` carries the model's densities, in the upload**:
  `AIS/RICE/icepack2/params.nc`, beside `CORE/`, where the organisers' run of
  the tool reads it and where the tool looks by default, so it runs without
  `--params-path`. Give it ours: 917 ice and 1024 seawater, the pair
  `simulation.py` builds the surface with, and 1000 fresh water. The tool
  defaults to 1027 seawater, and the comparison refuses any other densities.
  `scalar_processing.script` reads the tree's file and stops without one.
- **The tool stops without a historical run.** A run on its own names itself
  as `--hist` with the stamped year of its first state as `--refyear`. A
  projection paired with its historical names `--hist historical
  --hist-configid C001` (C002 for MRI-ESM2-0) and no year: the tool looks a
  year up in the historical and then in the projection, so a stray one quietly
  becomes the reference.
- **`--refyear` takes the stamped year.** A state variable is
  stamped 1 January of the following year, so a run starting in 2015 has 2016
  as its first state, and a reference of 2015 is not found.
- **The tool's files carry the submission's own names** (`lim_AIS_RICE_...`),
  so `--outpath` stays outside the upload tree. The upload carries the model's
  scalars; the tool's copies of them fail the checker. `params.nc` sits in the
  tree above `CORE/`, where the checker, which reads one set-counter directory,
  never sees it.
- **`--native-af2` says the model's scalars integrate over true area**, as a
  current forward writes them. The writer checks every year of the scalars
  CSV against the annual files, refuses a series that mixes true-area and
  map-plane years, and stamps the scalar files with the one it found
  (`scalar_area`). The comparison exits 2 when the switch disagrees with the
  stamp; leave it off for a run from before the change.
- **Core 11 goes through `scripts/isschecker_ocx.py`.** `--scenario ocx`
  names its forcing field `ERA5`, as the organisers' conventions document
  does (issue #18).
  Stock isschecker 0.5.1 has no `ocx` experiment, so over a C011 directory it
  reports one naming error and checks no file, and its field 5 takes CMIP
  models only. The script runs 0.5.1 with both added in-process, and takes
  the checker's arguments:

  ```bash
  <venv>/bin/python antarctica/scripts/isschecker_ocx.py \
      --source-path submission/AIS/RICE/icepack2/CORE/C011 --variable-list ismip7
  ```

  The `ocx` row is the Protocol Overview's C011 window, a start from 1990 to
  2015 and an end in 2025, so a start before 1990 or an end before 2025
  draws time errors. The log's version line says the checker was patched.
  The script refuses any release but 0.5.1, and stops once a release has an
  `ocx` row of its own.
- **A tree written before the whole-pixel flux means needs `--overlap`.**
  Its flux files carry no `flux_pixel_mean`, and `acabf` there is a mean over
  the covered part of a pixel, which the comparison undoes with the writer's
  cached `<annual>.overlap.npz`. A current tree needs nothing undone.

Measured on that rehearsal: `isschecker` 0.5.1 over 31 files reports zero
errors in variable presence, naming, numerical, spatial, consistency and
attribute tests, and 93 time errors, which are the three-per-file
experiment-length checks a two-year run cannot satisfy. `ismip7-scalars` 0.1.0
then wrote `sla20`, `slg20` and `slvaf`, each with its glacier and ice-cap
variant, in NetCDF and CSV. On three full-length 32 km ssp585 runs (IU Quartz,
23 September 2026) every identity `compare_scalars.py` checks holds, and
`reports/scalar_comparison_32km.md` sets out what differs and why. The tool's
slvaf runs above the model's as shelves thin and the grounding line retreats,
by about as much on a 1 km mesh as on a 32 km one (issue #99).
`scripts/grid_vaf_attribution.py` measures it for any state or year of annual
output, and `reports/grid_vaf_resolution.md` has the September 2026 numbers.

At full length, 2015 to 2300, a 32 km control and ssp585 pass 0.5.1 with zero
errors in every test group, the length checks included (23 September 2026, run
records `core09-32km-ctrl2015-cesm2waccm-p4` and
`core07-32km-ssp585-cesm2waccm-p4`). The first full-length pass failed on two
things a short run does not reach, both since fixed in the output: melt booked
before the positivity limiter, and cells within 1 cm of flotation written as
floating (`antarctica/FORWARD_RUN_READINESS.md`, action 10).

### The run log (`build_runlog.py`)

None of the files above are in git: a checkpoint is hundreds of MB and a
timeseries is regenerated. What reaches the repository is a record. Every
simulation, of any kind, leaves one JSON file in `antarctica/runlog/`, and

```bash
make -C antarctica runlog                                # render reports/SIMULATIONS.md
make -C antarctica runlog RUNLOG_CSV=results/runlog.csv  # and the records flat
make -C antarctica runlog-check                          # exits 1 on a stale render
```

renders them into `antarctica/reports/SIMULATIONS.md`, grouped by task kind,
with a summary table per group and every recorded field below it. A record is
reviewed in a pull request beside the code its run used, which is what makes it
a trace rather than a note; `runlog-check` and `tests/test_runlog.py` refuse a
stale render or a malformed record.

The group's shared progress sheet imports `RUNLOG_CSV`, one row per record and
one column per field in `build_runlog.py` order, so a run's configuration and
outcome are typed only in its record.

`id`, `task`, `title`, `status` and `institution` are required and everything else is
optional, so a planned run is five lines and a finished one carries its whole
provenance: mesh, MAP and iteration count, what it branched from, forcing
versions, the melt and front configuration, site, ranks, job ids, code SHA,
cost per model year, the audit verdict, whether the ISMIP7 output was written,
the checker version and result, the scalars, and the upload. `build_runlog.py`
owns the field list, so a field it does not know is an error rather than a
typo that renders as nothing.

A core experiment still gets its full report from `core_report.py`, with the
budget rows and the ensemble overlay. The run log is the index across all of
them, and across the runs that are not core experiments: the inversions, the
calibrations and the forward tests. A superseded run keeps its record with
`status` set to `superseded` and the reason in `notes`. `institution` names the
institution that owns the run, or `unassigned` for a planned run; `results` is
relative to the checkout on the site the record names, and never an account or
home-directory path.

---

## Known issues

**Diagnostic-Newton wall on hard projection geometries.** The forward blow-ups
are fixed (balanced apparent-MB init plus persistent DG0 thickness state), and
the Newton can still stall on evolved geometries. `ISMIP7_SNES_TYPE` and
`ISMIP7_SNES_MAXIT` are the knobs. The aSMB-forced walls seen so far predate
the annual-mean forcing fix and may be forcing induced; see
`reports/MATRIX_STATUS.md`.

**Wall retry.** When the rescue ladder is exhausted the run saves and stops
short of its target year. Relaunching from that state has cleared the wall in 3
of 3 observed cases (ssp585-CESM at 2096.7, CTRL-CESM at 2268, CTRL-MRI at
2250, both CTRLs then reaching 2300), because a fresh process re-runs the n=1
to n continuation at the loaded geometry. `run_core_matrix.sh` does this
automatically on a workstation. On a cluster the relaunch is yours: resubmit,
and `ISMIP7_AUTO_RESUME` picks the run up. A checkpoint written with
`stalled=1` makes `projection.sbatch` report and exit 1 without a successor, so
an unattended chain cannot spend days re-attempting the same years.

**Mesh and sidecar naming is per `(COARSE, FINE, BUFFER_M)`.** Older meshes
built before this convention need renaming or rebuilding.

**`icepack2_tools/coupled.py`** sketches ice-plume coupling against a
`PlumeModel` that does not exist in this tree. Nothing runs it.

**Timing runs are deliberately strict.** A matrix lane does not enter the
diagnostic rescue ladder or retry with transport subcycles. The first failure
is written to its timing JSON and status stamp, and a failed scout blocks that
mesh's scaling lanes. Inspect the recorded category and the Slurm log before
using `FORCE_TIMING=1` for an intentional retry.

---

## References

- Burgard et al. 2022, *The Cryosphere*, basal-melt parameterisation assessment.
- multimelt reference implementation: https://github.com/ClimateClara/multimelt
- ISMIP7 ocean forcing pipeline: https://github.com/ismip/ismip7-antarctic-ocean-forcing
- Greenland companion: https://github.com/dlilien/ISMIP7_Greenland_Icepack
