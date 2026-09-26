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
| **icepack_tools** (`adapt_mesh`, `levelset`, `friction`, `grounding`) | github.com/hoffmaao/icepack_tools | | level-set front | x | |
| **tlm_adjoint** | github.com/jrmaddison/tlm_adjoint | x | | | |
| `xarray netCDF4 scipy rasterio pyproj shapely gmsh matplotlib` (`geopandas` only to build a mesh, section 3) | pip, into the Firedrake venv | x | x | x | x |
| `earthaccess` (NSIDC), `globus-sdk` (Globus route only) | pip | x | x | x | |
| **isschecker** (`ismip7-compliance-checker`) | github.com/ismip/ISM_SimulationChecker | | | | x |

`icepack_tools` is a separate repository. Install it editable into the
same venv: `pip install -e /path/to/icepack_tools`.
`icepack2_tools/adapt_mesh.py` and `icepack2_tools/levelset.py` wrap it; the
rest of the repository runs without it.

The compliance checker needs Python 3.11 or newer. Check `python -V`; a
Firedrake venv often carries an older one, in which case give the checker its
own venv and call it by absolute path.

### 0.2 Data

| | Size | How | Invert | Forward run | Adapt | Submit |
|---|---|---|:--:|:--:|:--:|:--:|
| BedMachine Antarctica v4.1, MEaSUREs velocity v2 | 8 GB | `scripts/download_data.py` (Earthdata login) | x | x | x | |
| RACMO2.4p1 SMB climatology | 2 GB | same script | | x | | |
| ISMIP7 observations MIPkit v1.2 (Smith dH/dt) | 9 GB | `scripts/download_mirror.py --product ismip7-ais-observations data/mipkit/`, landing at `ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc` (`ISMIP7_OBS_KIT` overrides). `scripts/download_forcing.py --calibration` stages the same v1.2 file in the same place over Globus | `ISMIP7_DHDT_WEIGHT` | | `--from-obs` | |
| ISMIP7 forcing per ESM and scenario: SMB anomaly 7.5 GB, ocean `tf` 11 GB, `so` 6.9 GB (ssp585; historical 4.3 GB) | 25 GB each | `scripts/download_mirror.py` | | x | | |
| ISMIP7 fracture (collapse mask, lake properties, excess melt) | 3 GB per scenario | same, `data/<ESM>/<scenario>/fracture/` | | `ISMIP7_FRACTURE=mask` | | |
| Ocean OI climatology and IMBIE basin numbers | 3 GB | `scripts/download_forcing.py --ocean --calibration` | | x | | |
| Whole AIS tree (all ESMs, scenarios, `ctrl`, OCX, calibration) | 313 GB | same | | | | |
| Meshes and MAP checkpoints | 15 MB, 80 MB | sections 3 and 4, or from a colleague | | x | x | |
| Melt calibration `calibration/deltaT_per_basin_1000_K6.500e-05.npz` | 8 kB | tracked in git (section 5) | | x | | |

Source Cooperative carries the data-freeze copy and needs no account.

```bash
# one scenario for one ESM (about 25 GB)
python antarctica/scripts/download_mirror.py \
    data/CESM2-WACCM/ssp585/SDBN1-8000m/acabf-anomaly/ \
    data/CESM2-WACCM/ssp585/ocean/tf/ data/CESM2-WACCM/ssp585/ocean/so/
# the control's ocean for one ESM, which cores 9 and 10 read (about 18 GB)
python antarctica/scripts/download_mirror.py \
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
tree. `--scenarios` mirrors the minimal runtime sets (SDBN1-8000m `acabf` and
`acabf-anomaly`, ocean `tf` and `so`, fracture) with version autodetection and
checksum sync, so re-runs are completeness checks. Climatology, obs and
calibration sets come from `/ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/AIS`; the
`OCEAN_FILES` and `CALIBRATION_FILES` dicts in the script are the manifest.
Without `GLOBUS_LOCAL_ENDPOINT` the script prints the paths for a manual
transfer in the web app.

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
atmosphere directory is whichever of `SDBN1-8000m` and `GEMB-SDBN1-8000m`
exists, `SDBN1` first, so trees fetched before MRI's August 2026 rename still
work. Fracture masks resolve flat or versioned.

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
(`misfit_norm`, `gamma_theta`, `gamma_phi`, `log_vel_weight`, `log_vel_eps`,
`dhdt_weight`, `dhdt_net_sigma`, `mesh_basename`, `lc`, `lc_coarse`,
`buffer_m`):

```bash
python -c "import h5py,sys; print(dict(h5py.File(sys.argv[1])['/'].attrs))" MAP.h5
```

---

## 5. Ocean melt calibration

### The calibration every run reads

A run melts with one K and a thermal-forcing offset per IMBIE basin, the
protocol's recommendation, from a tracked file,
`calibration/deltaT_per_basin_1000_K6.500e-05.npz`:

| | |
|---|---|
| K | 6.5e-5, the K50 of the rule-based selection below on the 1000 m / 10 km production mesh (run record `calibration-melt-toolbox-1km-rule`), chosen by the group on 25 September 2026 (issue 26) |
| offsets | 16 IMBIE2 basins, -0.68 K to +1.20 K (Amundsen), each basin at its July 2026 total, 1067.4 Gt/yr together |
| fitted under | `ISMIP7_MELT_SLOPE=ant`, `ISMIP7_SIN_ALPHA_ANT=5.115e-3`, `ISMIP7_GEOMETRY_SPACE=dg0`, `ISMIP7_RASTER_SAMPLE=vertex` and the 30_sep OI climatology, on `antarctica_10000_1000_buffered20000` |
| sidecar | `deltaT_per_basin_1000_K6.500e-05.source.json`: the file's sha256, the settings above, the input hashes, the job and the commits |

Every clone carries it, so a new machine runs with it and nothing is
calibrated or copied. The offsets are stamped onto any mesh through the
IMBIE2 8 km basin grid under `ISMIP7_DATA_ROOT`.

The forward applies the melt the file was fitted to:

- it melts the cells the fit summed over, floating and holding ice
  (`forcing.melt_receiving`);
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
run with its offsets as they are.

```bash
ISMIP7_LC=<lc> ISMIP7_INV_H5=<mesh.h5> python antarctica/scripts/calibrate_deltaT.py --K 6.5e-5
```

`check_melt_bound.py` melts the reference geometry with the forward's own
callback and sets each basin against the total its offsets were fitted to,
and exits 1 when a basin is off by more than `--match-tol` (default 0.1
percent):

```bash
ISMIP7_INV_H5=<mesh.h5> python antarctica/scripts/check_melt_bound.py
```

On the production mesh (run record `calibration-melt-forward-1km-k50`, Quartz
job 10644378) the forward melts 1067.390 Gt/yr against the 1067.386 its offsets
were fitted to, every basin within 0.006 Gt/yr, and no cell passes the
`libmassbffl` bound (maximum 41.7 m/yr). Its earlier melt set also covered
257 687 ice-free open-ocean cells, where it booked 156.3 Gt/yr of melt and
23.1 Gt/yr of refreezing, and 6 666 cells of bare land, where the climatology
melts nothing. At 32 km (job 10644430) the same offsets put the basins
at 0.33 to 1.74 times their totals and the whole at 1069.5 Gt/yr, which is why
a production core runs on the calibration's mesh.

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
residual, dM/dT, the slope and geometry conventions). A run applies it in
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
the control warns that it starts from a different geometry.

Run management on the drivers: `--restart <ckpt>` or `ISMIP7_RESTART` resumes;
`ISMIP7_AUTO_RESUME=1` picks up the newest checkpoint for the experiment, which
is what lets a chained batch job continue itself, and takes precedence over the
historical endpoint so only the first link starts there; `--tag` or
`ISMIP7_RUN_TAG` suffixes the experiment name so a method line keeps and
resumes its own files; `--checkpoint-interval` sets the step-count fallback.
Checkpoints carry the mesh, geometry, inversion fields and the full `(u, M, τ)`
state, so restarts work at any rank count. A resume refuses to start when
`ISMIP7_FRICTION` or `ISMIP7_APPARENT_MB` disagree with the checkpoint.

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
| `score_map.py MAP.h5 [...] [--json OUT] [--restart STATE.h5]` | scores an inversion by `Q(u_model)/Q(u_obs)` across its own grounding line, overall and per speed band; re-solves through `setup_model`, so periodic MAPs without a velocity work. With `ISMIP7_MESH` set the MAP is transferred first, so the same command with and without it compares a transfer; `--restart` scores a prepared state, `--json` writes the numbers with the transfer fill counts (`MAP_CHECK.md`) |
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
force. The one momentum change is the drag gate: the mask is 1 only where `phi`
exceeds one cell diameter, so floor-cell ocean drag acts only in water further
than a cell from the front. Thin cells inside the t=0 extent are damped by
`h_visc_floor`, the friction law and the `ISMIP7_ALPHA_GL` collar. In the strip
a free law advances into, `C_w0` comes from the t=0 geometry where `H = 0`, so
`tau_b = 0` under both laws and the damping is `h_visc_floor` and the collar.

`vonmises` is Morlighem et al. 2016 verbatim:
`c = |u| sqrt(3) B eps~^(1/n) / sigma_max`, with `eps~` from the tensile
principal strain rates, `B = A^(-1/n)`, and separate grounded and floating
thresholds. Those thresholds are the tuning targets: a 2015 control should hold
the observed front (the obs kit's 24 yearly Greene masks, 1997 to 2021) and
discharge about 1300 Gt/yr. The level set is checkpointed as `levelset` for
diagnostics; a restart rebuilds the front from the thickness. The exception is
`fixed`, which anchors on `H_init` so a resumed run does not re-freeze the
front where it restarted. The shared implementation's tests are
`icepack_tools/test/levelset_test.py`; the ISMIP7-side rules (retreat-sliver
mask, apparent-MB extent masking, the `fixed` law's t=0 anchor) are covered by
`tests/`. `tests/test_levelset_laws.py` checks each law against closed forms on
a unit mesh: the `vonmises` rate under both thresholds, the shed fraction and
its step-size behaviour under the transport's masks, the drag gate, and the
refusal of an unknown or underspecified law.

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

`vonmises` is for projections. A configured law owns the front outright, so the
legacy `ISMIP7_FIXED_FRONT` mask removes nothing when a law is set and
`vonmises` is never silently pinned. The apparent-MB reference `a_ref` is
defined only on the t=0 ice extent under every law. Under a free law it is also
cleared each step wherever the level set reports ice-free, irreversibly, so a
calved cell is not regrown and an advanced-into cell is not re-emptied. The run
log prints one `Calving front owner:` line naming the mechanism in force.

### A calving law from hoffmaao/calving (`forward_calving.py`)

The laws developed for CalvingMIP (github.com/hoffmaao/calving: `fixed`,
`velocity`, `position`, `thickness`, `vonmises`, `vonmises_strain`, `hfb`
and any `--law-module` plugin) run in the forward unchanged. A law reads
seven fields from its model, and `simulation.LiveCalvingState` supplies
them from the live dual state: `u`, `M`, `tau` from the mixed solution,
the DG0 `h`, `haf` and the grounded indicator as UFL on the cells, and the
front normal from the level set the forward advances. The law's rate goes
to the shared level set as its `prescribed` law, so the front is retreated
and the shed mass is tallied exactly as under `ISMIP7_CALVING=vonmises`.

```bash
CALVING_DIR=/path/to/calving mpiexec -n 16 python antarctica/scripts/forward_calving.py \
    --law hfb --law-param sigma_max=0.15 --experiment control --tag hfb_test
```

wraps the experiment driver (`control`, `ocx`, `ssp126|ssp370|ssp585` with
`--esm`, `hist`), which keeps its own forcing, restart and auto-resume;
leave `ISMIP7_CALVING` unset. `--tag` (or `ISMIP7_RUN_TAG`) is required, so a
law-driven run never resumes from or overwrites a stock run's checkpoints. The
law's `describe()` is written to every checkpoint's `calving_law` attribute and
to the `Calving front owner:` line, which `core_report.py` lifts into the
report. Any object with `rate(model, t)` returning a UFL rate on the cells and
`describe()` can be placed in `ctx["calving_law"]` before `run_simulation` the
same way. A threshold is fitted on Antarctica
with `calving/tune_greene.py --experiment vonmises|hfb --state <forward
checkpoint>` (per-Mouginot-basin flux against the Greene et al. 2022 fronts,
on the same level-set normal), which is what makes a tuned parameter mean
the same thing in both places.

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

### Environment knobs (inversion)

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_MAP_OUT` | output path for the MAP, overriding the generated name. Use it for smoke tests and variants so a short run cannot replace a production MAP. A bare filename resolves under `mesh/` | generated |
| `ISMIP7_MISFIT_NORM` | `sigma` divides each residual by its datum's squared error, giving a dimensionless chi^2; `none` is the legacy dimensional misfit. Selects the `ISMIP7_GAMMA_*` defaults | `sigma` |
| `ISMIP7_LOG_VEL_WEIGHT` | weight on the ISSM logarithmic velocity misfit (cost function 103). The chi^2 alone over-weights slow interior ice and leaves discharge-carrying tributaries 40 to 50% too slow; the log term is scale free. `auto` equalises it with the chi^2 term at the warm-start state. Stamped into the MAP | `0` |
| `ISMIP7_LOG_VEL_EPS` | regularisation speed (m/yr) inside the log | `1.0` |
| `ISMIP7_WARM_START` | path to a MAP or timing-cache checkpoint used to seed `theta`/`phi` (and, when present, geometry, `fluidity_prior`, and the mixed diagnostic state). Fields are interpolated onto the live mesh, so a 1-core cache can warm-start a multi-rank invert | unset |
| `ISMIP7_SKIP_CONTINUATION` | `1` skips the cold `n,m: 1→n` ramp on the initial solve and inside each annotated forward eval (single solve at full exponents). Auto-enabled when the warm start supplies a full mixed state | `0` |
| `ISMIP7_GAMMA_THETA` / `ISMIP7_GAMMA_PHI` | Whittle-Matern prior strength on `θ` and `φ`, coupled to `ISMIP7_MISFIT_NORM` since normalising divides the misfit by about sigma^2. `scripts/lsurface.py` sweeps both on a grid and picks the L-curve corner of each (usage in its docstring) | `1e5` under `sigma`, `1e4` under `none` |
| `ISMIP7_L_REG` | prior correlation length (m) | `7.5e3` |
| `ISMIP7_MAXITER` | L-BFGS-B iteration cap | `500` |
| `ISMIP7_GRAD_PRECOND` | `none` is the raw-dof l2 metric, which is mesh dependent, so fine grounding-line cells converge slowest. `mass` optimises in `u = sqrt(M) x` under scipy, making the rate mesh independent. `mass_consistent` and `prior` run under TAO instead (scipy takes no preconditioner) with the initial inverse Hessian set to `M^-1` or to the prior covariance; `mass_consistent` is the consistent mass Riesz map; `prior`, the prior-preconditioned variant, was not used and is experimental here. Defaults to `none` to keep runs comparable with everything measured so far | `none` |
| `ISMIP7_PRIOR_FORM` | `laplacian` uses `A = delta*M + gamma*K` as the prior precision; `bilaplacian` uses `A M^-1 A`, the squared-operator prior of Villa et al. (2021), the operator that a 2-D Whittle-Matern field needs to be function-valued. Different priors, not two spellings of one: their gammas are not convertible and their MAPs are not comparable, so the MAP stamps `prior_form` | `laplacian` |
| `ISMIP7_PRIOR_SIGMA_THETA` / `_PHI`, `ISMIP7_PRIOR_RHO` | `bilaplacian` only: the log-deviation scale and correlation length (m), converted to `(delta, gamma)` by the closed forms of Villa et al. (2021), `sigma^2 = 1/(4 pi gamma delta)`, `rho = sqrt(8 gamma/delta)`. The un-squared form has no such closed form, which is why its gamma can only be tuned | `0.3` / `0.3` / `ISMIP7_L_REG` |
| `ISMIP7_PRECOND_STEP0` | TAO metrics only: the largest change the FIRST step may make to a control, in that control's units, applied per control block. L-BFGS's first step is `-H_0 g` at unit length with no curvature pair to rescale it, and one evaluation outside the region where the forward has a solution returns NaN that every later trial point inherits | `0.15` |
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
| `ISMIP7_BUFFER_M` | outline buffer (m) in the default mesh and sidecar names | `20000` |
| `ISMIP7_MESH` | mesh path for the inversion and tools. A forward takes its mesh from the checkpoint unless this names another mesh, in which case the MAP is transferred onto it. `checkpoint` means the mesh embedded in the MAP or restart file: `site_env.sh` always exports a derived path, so this is how a job submitted through `submit.sh projection` runs MAP-native | `mesh/antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.msh` |
| `ISMIP7_RASTER_SAMPLE` | how BedMachine lands on a DG0 cell. `vertex` projects the CG1 vertex interpolant; `cell_mean` takes the raster's true cell mean. `cell_mean` measured rougher: neighbouring cells share two of three vertex samples, so `vertex` damps jumps by construction. Cell means raised interior surface jumps 6% and bed and thickness jumps 35%, and at 2 km the momentum solve did not converge within 60 minutes. It does classify flotation better (32 km misclassification 9.1% to 3.2%), so the knob stays. Stamped into the MAP and read back by the forward. Reproduce with `probe_raster_sampling.py` | `vertex` |
| `ISMIP7_INVERSION` | explicit MAP path for a forward or preflight. The forward checks the MAP's recorded `friction`, `n_flow` and `geometry_space` against the run and aborts on a mismatch, warning only when the MAP predates those attributes; `preflight.py` checks that the file exists. Use it to A/B MAPs on one mesh, or, with `ISMIP7_MESH` also set (the timing matrix, `make map-check`), to run a MAP on a different mesh: its continuous fields are then interpolated onto `ISMIP7_MESH` by strict point location (`icepack2_tools/transfer.py`), and a target dof outside the MAP's outline takes a stated fill (0 for the log controls, the constant baseline for the fluidity prior, the raster sample for `velocity_obs`), counted and printed as `Transfer fill:` lines (`MAP_CHECK.md`) | derived |
| `ISMIP7_CALVING` | `none`, `fixed` or `vonmises` (see above) | `none` |
| `ISMIP7_CALVING_SIGMA_MAX_GROUNDED` / `_FLOATING` | von Mises thresholds (MPa) | `1.0` / `0.15` |
| `ISMIP7_FRACTURE` | `mask` applies the ISMIP7 collapse forcing to every floating cell it flags, booked as calving; `mask_front` only to the flagged cells open water has reached, so no hole opens behind a standing front (the two end-members of discussion #30). Masks exist for the SSPs only, so the control, historicals and OCX abort on either. Needs DG0. Every run prints its mode once (`Ice-shelf collapse forcing: ISMIP7_FRACTURE=...`, `none` included). Under a mask mode the timeseries columns `collapse_flagged_cells`, `collapse_removed_cells` and `collapse_held_cells` count the flagged floating cells, the ones the mode has emptied and the ones it leaves standing (always 0 under `mask`), all ranks summed; the budget lines, a closing log line and the core report repeat them | `none` |
| `ISMIP7_OCX_FORCING` | what core 11 runs on. `protocol` is the ISMIP7 OCX product (RACMO2.3p2-ERA SDBN1 `acabf`, expert-judgment ocean `tf`/`so`), and the run refuses to start without it. `stopgap` is RACMO2.4p1 actual-year SMB with the constant OI ocean climatology, what the core ran on before the product was readable here. K is fitted to the climatology, so read `check_melt_bound.py --ocx` first (discussion #48) | `protocol` |
| `ISMIP7_OCX_OCEAN` | which expert-judgment OCX ocean scenario to read: `main` (the core one), `cold`, `warm` or `vary`. A member other than `main` writes to `ocx_<member>` | `main` |
| `ISMIP7_OUTPUT` | `1` records the ISMIP7 yearly fields and scalars (`<exp>_<lc>_ismip7_annual_<year>.h5`, `<exp>_<lc>_ismip7_scalars.csv`), regridded afterwards by `write_ismip7_output.py`. The value set is closed, so a typo is rejected at startup. `projection.sbatch` and `run_core_matrix.sh` default it to `1`, and `=0` or an empty value turns it off there. A chained projection must export it on every link; a link that cold-starts mid-year logs the gap and begins at the next 1 January. Resuming continues a series, and a cold start into a populated series is refused | unset (`1` under the core-experiment runners) |
| `ISMIP7_BNDIDS` | boundary-id JSON override | per-mesh sidecar, else `mesh/boundary_ids.json` |
| `ISMIP7_GEOMETRY_SPACE` | `dg0` (one thickness for terminus force and mass flux) or `cg1` (legacy, A/B only). Selects the MAP. See `../GEOMETRY_DISCRETIZATION.md` | `dg0` |
| `ISMIP7_DATA_ROOT` | forcing tree root | `<repo>/ISMIP7/AIS` |
| `ISMIP7_OBS_DATA_ROOT` | BedMachine, MEaSUREs velocity, RACMO and the dH/dt cache observational-data root; name it in a site file when these files do not live beside the code. Also a write target: with the MIPkit present `obs_dhdt` builds `<root>/dhdt_cache/` here, so staged cache tifs belong under this root, wherever it points | `<repo>/antarctica/data` |
| `ISMIP7_T_END` / `ISMIP7_DT` | end time and timestep (yr). `t=Y.0` is 1 January of year Y, so a run covering 2015 to 2300 ends at `2301` and a historical covering 2003 to 2014 ends at `2015`. Each driver owns its end (historical `2015`, ssp370 `2101`, other projections and control `2301`, OCX `2026`) | driver's own / `1.0` |
| `ISMIP7_GEOMETRY_BACKDATE` | years of the Smith et al. (2020) mean dH/dt a cold start undoes on grounded ice before it runs (issue #117). Unset: `2015 - ISMIP7_T_START` for a start from 2003 up to 2015 (the historicals and OCX start in 2003), none from 2015 on, and a start before 2003 is refused. The friction anchors stay on the 2015 geometry; floating ice keeps its 2015 thickness. `0` turns it off | driver's start |
| `ISMIP7_FRICTION` | `budd`, `regularized_coulomb` or `budd_legacy`; selects the MAP. The set is closed, so a misspelling is rejected at startup | `budd` |
| `ISMIP7_OUTPUT_INTERVAL` | timeseries row every N steps | `10` |
| `ISMIP7_CHECKPOINT_EVERY_YR` / `ISMIP7_KEEP_CHECKPOINTS` | checkpoint cadence in model years, and how many to keep besides `_final.h5` | `5` / `3` |
| `ISMIP7_RESTART` | restart checkpoint | `hist_<esm>[_<tag>]_<lc>_final.h5` if present |
| `ISMIP7_AUTO_RESUME` | resume from this experiment's newest checkpoint when no `ISMIP7_RESTART` is given. An integer flag, `=0` disables it, since the runners export it unconditionally and `--export=ALL` cannot unset. `projection.sbatch` refuses to chain when it is off | unset |
| `ISMIP7_RUN_TAG` | experiment-name suffix for a parallel method line | unset |
| `ISMIP7_WALL_STOP_MIN` | wall-clock budget in minutes from process start, checked before each step against the longest step so far, so the run writes its final checkpoint and exits with `t_yr` short of `t_end` for a chained job to resume. `projection.sbatch` derives it from the job's own TimeLimit, holding back 25 minutes. `0` disables | `0` |
| `ISMIP7_EXPERIMENT_NAME` | the run's identity, used by `adapt_mesh.py` to name adapted meshes and sidecars so parallel experiments cannot overwrite each other. Set by `run_adaptive.py --experiment-name`. See `../ADAPTIVE_MESH.md` | unset |
| `ISMIP7_APPARENT_MB` | `1` or `balance` zeroes the t=0 thickness tendency; `div` cancels only the flux divergence; `0`, `off`, `none` and empty disable it | unset |
| `ISMIP7_FIXED_FRONT` | hold the calving front at the t=0 extent, tallying inflow beyond it as calving. `=0` disables. Ignored whenever an `ISMIP7_CALVING` law is configured | unset |
| `ISMIP7_TRIPWIRE_U_MAX` / `ISMIP7_TRIPWIRE_H_MAX` / `ISMIP7_TRIPWIRE_DH_RATE` / `ISMIP7_TRIPWIRE_HMIN` | runaway tripwire: fail the step when max speed exceeds `U_MAX` [m/yr], max thickness exceeds `H_MAX` [m], or a cell that entered the step at least `HMIN` thick thickens at a relative rate `(dh/h)/dt` above `DH_RATE` [1/yr] (a rate so every dt scores the same physics alike; thinner cells are reported, never tripped: buffer cells fill by more than their own thickness); every step prints a `tripwire step-k:` line with the worst cells; unset = off (timing lanes export 2e4 / 5000 / 20 / 100) | _(unset)_ |
| `ISMIP7_LEGACY_TRANSPORT` | restore the pre-July-2026 CG-projection transport (needs `cg1`) | unset |
| `ISMIP7_SNES_TYPE` / `ISMIP7_SNES_MAXIT` | diagnostic Newton type and iteration cap | `newtonls` / `200` |
| `ISMIP7_SNES_LINESEARCH` | line search used by `newtonls` | `nleqerr` |
| `ISMIP7_SNES_RTOL` / `ISMIP7_SNES_ATOL` / `ISMIP7_SNES_STOL` | initial nonlinear relative, absolute and step tolerances. After the initial/restart solve, the absolute tolerance follows the self-scaled policy below | `1e-8` / `1e-50` / `0` |
| `ISMIP7_SNES_DIVERGENCE_TOL` | residual-growth divergence threshold; PETSc's `-3` (`PETSC_UNLIMITED`) disables this test (`-1` means `PETSC_DETERMINE`, restoring the default `1e4`) | `-3` |
| `ISMIP7_SNES_ATOL_SCALE` / `ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE` | persistent absolute tolerance after a converged setup solve (`scale * achieved norm`) / after accepting a loaded hard-era state (`scale * loaded-state norm`) | `100` / `1e-6` |
| `ISMIP7_SNES_KSP_EW` | enable PETSc Eisenstat-Walker variable inner tolerance for an A/B test | `0` |
| `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER` | `schur_gamg` or `schur_mumps`: legacy PETSc `selfp` approximation; `scpc_gamg` or `scpc_mumps`: exact cell-local Slate elimination and an assembled velocity solve; `full_mumps`: complete mixed-Jacobian reference. Legacy `iterative`/`mumps` aliases mean `schur_gamg`/`schur_mumps` | `scpc_gamg` for cluster forwards (`batch_runners/projection.sbatch`); `full_mumps` for a forward driver run by hand, the core runner and the workstation launchers, and the inversion's own linear solve whatever is set; the timing Makefile's `TIMING_SOLVER` (`scpc_mumps`) |
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
| `ISMIP7_DELTAT_PER_BASIN_NPZ` | the melt calibration, one K and a TF offset per basin (`select_melt_parameters.py`, `calibrate_deltaT.py`); every ocean callback melts with its K. Refused with `ISMIP7_K_SCALE` other than 1, and when fitted under another slope or geometry than the run's (section 5) | `calibration/deltaT_per_basin_1000_K6.500e-05.npz` |
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
| `ISMIP7_H_OCEAN` / `ISMIP7_K_LIM` | front backstops read by `scripts/simulation.py`: the thickness (m) at which the ice-free ocean drag ramps to zero, and the speed-limiter coefficient the rescue ladder raises for a rescue solve | `10.0` / `1e-3` |
| `ISMIP7_ALPHA_GL` | grounding-line coercivity, Budd only, read by `scripts/simulation.py` and `scripts/inversion_icepack2.py` | `0.5` (`0` for RC) |
| `ISMIP7_RC_HVISC_FLOOR` / `ISMIP7_RC_CW0_FLOOR` | RC viscous-thickness and `C_w0` floors, read by `scripts/simulation.py` and `scripts/inversion_icepack2.py` | `10.0` / `0.0` |
| `ISMIP7_M_SLIDE` | sliding exponent, read by `scripts/inversion_icepack2.py`, `scripts/simulation.py`, `scripts/thermo_prior.py`, `scripts/plot_map.py`, `scripts/run_eigendec.py` | `3.0` |

> **dt guidance.** Production projections on the 1000 m mesh run
> `ISMIP7_DT=0.05`, the step the timing matrix ran there and the one
> `projection.sbatch` defaults to. At 2500 m and coarser use `0.1`; `0.25` is
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
64 ranks, `dt = 0.05` yr**. Minutes of transient loop per simulated year on
that mesh (and the 285-year extrapolation):

| 1000/10000 | 16 ranks | 32 ranks | 64 ranks |
|---|---|---|---|
| `scpc_mumps` | 34.0 (6.7 d) | 23.4 (4.6 d) | 23.2 (4.6 d) |
| `scpc_gamg` | 32.3 (6.4 d) | 15.7 (3.1 d) | 10.0 (47.5 h) |

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
**Inversions are not part of this:** `inversion_icepack2.py` factors the
complete mixed Jacobian with MUMPS, because `tlm_adjoint` differentiates
through that solve, and no setting changes it. That is also why the switch is
made in the forward runner and not in `solverconfig`'s default, which the
inversion reads to stamp its MAP.

**Starting state, unsettled.** No MAP has been inverted on the 1000 m mesh, and
the plan for now is not to invert one: transfer the coarse MAP instead, the way
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
Re-inverting on the production mesh, and closing the Budd/RC gap, are both
open under the 2 km inversions; icepack/ismip7#21 was closed as their
duplicate on 22 September. (issue #24)

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
   tape rather than the condensed transient. **The 500 m meshes are not
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
   campaign cache is held to, whatever solver the lanes time) and `state_solver` (`full_mumps`, what actually produced the
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
  restart resolve the same sidecar.
- `<exp>_t<year>.h5`, periodic checkpoints, keeping the
  `ISMIP7_KEEP_CHECKPOINTS` most recently written.
- `<exp>_timeseries.csv`, one row per `OUTPUT_INTERVAL` steps:
  `year, vaf_mm_sle, mass_gt, smb_gtyr, melt_gtyr, outflux_gtyr, calv_gt,
  clamp_gt, resid_gt, amb_gtyr`. The residual must close to 0.00. SMB and
  melt are what the advances applied: no forcing acts on open ocean or on
  cells a front rule holds ice-free (`front.unforced_cells`), so `clamp` is
  only the positivity limit on thin ice and `calv` only ice that crossed the
  front. The ISMIP7 `acabf` and `libmassbffl` fields book the same forcing.

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
`reports/scalar_comparison_32km.md` sets out what differs and why.

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
