# ISMIP7 Antarctica — data, setup, and how to run

Antarctic ISMIP7 submission built on **icepack2 / Firedrake**. This document is
the end-to-end guide: what to install, how to download every input (including the
Globus forcing), and the exact command sequence to reproduce a run. Companion to
David Lilien's Greenland repo (https://github.com/dlilien/ISMIP7_Greenland_Icepack).

> New here? Read top to bottom once. The pipeline is a strict dependency chain —
> each step consumes the previous step's output (data → mesh → inversion →
> calibration → control/projection).

---

## 0. Prerequisites

### Software

- **Firedrake** with **icepack2** (the mixed/3-field formulation) and
  `tlm_adjoint` for the inversions. Firedrake brings its own PETSc/petsc4py/MPI.
  Activate that environment before anything else, e.g.:
  ```bash
  source ~/venv-firedrake/bin/activate
  ```
- Python packages used by the scripts (install into the Firedrake venv):
  ```bash
  pip install globus-sdk earthaccess xarray netCDF4 scipy shapely gmsh colorcet
  ```
- **gmsh** (the `gmsh` Python module above is sufficient for meshing).

### Accounts / credentials

| For | Account | Where |
|-----|---------|-------|
| BedMachine, MEaSUREs velocity (NSIDC) | NASA Earthdata (free) | https://urs.earthdata.nasa.gov/users/new |
| ISMIP7 ocean/atmosphere forcing | Globus + access to the ISMIP6/7 collection | https://app.globus.org |

To pull data **to this machine** over Globus you also need **Globus Connect
Personal** running locally and its endpoint UUID
(https://www.globus.org/globus-connect-personal).

### Repo layout (what lives where)

```
antarctica/
  data/            # observational inputs (BedMachine, velocity, RACMO)  [gitignored]
  mesh/            # *.msh + boundary_ids_antarctica_*.json + inversion_*.h5  [gitignored except boundary_ids*.json]
  results/         # checkpoints (*.h5), timeseries (*.csv), logs        [gitignored]
  scripts/         # all entry points (see §3–§6)
ISMIP7/AIS/        # ISMIP7 forcing tree the *runtime* reads             [gitignored]
icepack2_tools/    # reusable library (mesh, forcing, eikonal, grounding, regrid)
```

Everything large is gitignored. The only tracked files under `antarctica/mesh/`
are the tiny `boundary_ids_antarctica_*.json` (the gmsh physical-line →
calving/other map; see §3).

---

## 1. Observational data (`download_data.py`)

BedMachine geometry, MEaSUREs velocity, and RACMO SMB. RACMO is public (Zenodo);
the NSIDC products need an Earthdata login (handled interactively by
`earthaccess`).

```bash
cd antarctica
python scripts/download_data.py
```

| Dataset | Product | Auth | → lands in |
|---------|---------|------|-----------|
| BedMachine Antarctica v4 | NSIDC-0756 | Earthdata | `data/bedmachine/` |
| MEaSUREs Ice Velocity v2 | NSIDC-0484 | Earthdata | `data/velocity/` |
| RACMO2.4p1 SMB | Zenodo `10.5281/zenodo.14217231` | none | `data/racmo/` |

The scripts skip files that already exist, so re-running is cheap.

---

## 2. ISMIP7 forcing via Globus

There are **two distinct things** here:

1. **The runtime tree `ISMIP7/AIS/`** — this is what the simulation actually
   reads at run time (via `icepack2_tools/forcing.py`). Its layout follows the
   official ISMIP7 protocol (see §2b). Point the code at it with
   `ISMIP7_DATA_ROOT` (defaults to `<repo>/ISMIP7/AIS`).
2. **`download_forcing.py`** — a helper that mirrors the Globus collection
   straight into `ISMIP7/AIS/`, matching the runtime layout - no rename or
   reorganization step: climatology/bias/calibration files (`--ocean` /
   `--calibration`) plus the per-(ESM, scenario) runtime forcing sets from
   the top-level `/ISMIP7/AIS` tree (`--scenarios`). After downloading,
   `python scripts/preflight.py` reports which core experiments the local
   tree can actually run.

### 2a. Using `download_forcing.py`

```bash
cd antarctica

# one-time: authenticate through the Globus CLI. On a remote Quartz session it
# prints a URL; open it in your browser and paste the resulting code back.
# The CLI caches credentials under ~/.globus/cli/.
python -m pip install --user globus-cli   # if `globus` is not installed
python scripts/download_forcing.py --login

# browse the legacy climatology subtree to sanity-check paths
# (the /ISMIP7/AIS scenario tree is not listed here; use the Globus web app)
python scripts/download_forcing.py --list

# tell the script where to *put* the files (your Globus Connect Personal endpoint)
export GLOBUS_LOCAL_ENDPOINT=<your-endpoint-uuid>

# download (no mode flag runs --ocean + --calibration; --dry-run to preview)
python scripts/download_forcing.py --ocean        # CESM2-WACCM thetao/so/tf + climatology + bias
python scripts/download_forcing.py --calibration  # meltMIP obs melt, IMBIE2 basins, grid, topography
python scripts/download_forcing.py --scenarios    # per-(ESM, scenario) runtime forcing (cores 1-8)
python scripts/download_forcing.py --scenarios --esm MRI-ESM2-0 --scenario historical,ssp585
python scripts/download_forcing.py --status        # what's present locally
```

Knobs:

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_GLOBUS_COLLECTION` | source collection UUID | `ccc9bbd2-4091-4e35-addd-eeb639cf5332` |
| `GLOBUS_LOCAL_ENDPOINT` | **your** local Globus Connect Personal endpoint UUID | _(required to transfer)_ |

Remote base paths on the collection: the per-year scenario forcing lives in
the top-level `/ISMIP7/AIS/<ESM>/<scenario>/` tree (both CESM2-WACCM and
MRI-ESM2-0, scenarios historical/ssp126/ssp370/ssp585 and more - everything
cores 1-8 need). `--scenarios` mirrors the minimal runtime sets from it
(SDBN1-8000m acabf and acabf-anomaly, ocean tf and so, fracture) as
recursive-directory transfers with version autodetection (highest `v<N>` on
the share) and checksum-level sync, so re-runs are idempotent completeness
checks. The climatology/bias/obs/calibration sets still come from
`/ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/AIS`; their exact file lists are the
`OCEAN_FILES` and `CALIBRATION_FILES` dicts at the top of
`scripts/download_forcing.py` — that file is the authoritative manifest, and
its status comment records the collection layout.

No `GLOBUS_LOCAL_ENDPOINT`? The script prints the remote/local paths so you can do
the transfer by hand in the Globus web app instead.

### 2b. The runtime tree `ISMIP7/AIS/` (what `forcing.py` reads)

`icepack2_tools/forcing.py` resolves data as follows (override the root with
`ISMIP7_DATA_ROOT`):

```
ISMIP7/AIS/
  <ESM>/<scenario>/SDBN1-8000m/<var>/<version>/      # atmosphere (acabf, acabf-anomaly, ts, tas, pr, ...)
      <var>_AIS_<ESM>_<scenario>_SDBN1-8000m_<version>_<YEAR>.nc
  <ESM>/<scenario>/ocean/<tf|thetao|so>/<version>/   # ocean thermal forcing, salinity, temp
  <ESM>/<scenario>/fracture/[v*/]                     # ice-shelf collapse / lake masks (flat or versioned)
  meltMIP/OI_Climatology_ismip8km_60m_<tf|so|thetao>_extrap.nc   # CTRL climatology
  parameterisations/ocean/imbie2/                     # IMBIE2 basin numbers (per-basin K)
  parameterisations/ocean/{bfrns,meltobs,shelfmask,floatingmasks,...}
  parameterisations/fracture/
```

Defaults baked into the readers: atmosphere pins `version=v2` at `SDBN1-8000m`,
ocean pins `version=v3`; when the pinned version directory is absent the readers
fall back to the highest `v<N>` subdir present, so MRI-ESM2-0 `v1` and future
re-releases resolve without code changes. Fracture masks are found both flat in
`fracture/` and inside `fracture/v*/` (highest version wins). The forcing API:

- `ISMIP7Atmosphere(esm, scenario).get_smb(year, x, y, anomaly=…)`
- `ISMIP7Ocean(esm, scenario).get_thermal_forcing(...)` / `.get_salinity(...)`
- `ISMIP7Fracture(esm, scenario).get_collapse_mask(year, x, y)`
- `make_forcing_callback(atm=, ocean=, fracture=, K=…, K_per_basin_npz=…)` — bundles
  all three into the per-step callback `run_simulation` expects.

---

## 3. Build the mesh

Adaptive isotropic mesh with Ua-style grounding-zone + calving-front refinement,
sized from BedMachine geometry and MEaSUREs strain rate (needs §1 data).

```bash
cd antarctica
python scripts/mesh_antarctica.py --lc 2500 --lc-coarse 64000 --buffer-m 20000
# or equivalently via env vars (defaults match the rest of the pipeline):
ISMIP7_LC=2500 ISMIP7_LC_COARSE=64000 ISMIP7_BUFFER_M=20000 python scripts/mesh_antarctica.py
# dev mesh used by inversion_icepack2.py / diagnostic_solve.py / run_eigendec.py:
python scripts/mesh_antarctica.py --lc 8000 --lc-coarse 80000 --buffer-m 20000
# → mesh/antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.msh
#    (e.g. antarctica_64000_2500_buffered20000.msh)
# → mesh/boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json
#    (e.g. boundary_ids_antarctica_64000_2500_buffered20000.json)
```

`--lc` / `--lc-coarse` select the fine (grounding-line/calving-front) and
coarse (interior) element sizes in meters. The GL-band element sizes,
`shelf_size`, `buffer_size`, and strain-rate floor all scale with `lc/2500`;
calving-front decay lengthscales are floored at `lc` and `1.25*lc` so they
never fall below the mesh's own fine resolution.

The mesh outline is pushed `--buffer-m` / `ISMIP7_BUFFER_M` meters into the ocean before
meshing (default `20000`; pass `0` for no buffer), which lets icepack2 handle
`h=0` at the (now-interior) calving front instead of needing a
`calving_terminus` BC. Because this changes the boundary topology, both the
`.msh` and its sidecar are named after the exact `(COARSE, FINE, BUFFER_M)`
combination used to build them — building a different resolution or buffer
never collides with or silently invalidates a previous build.

**Boundary IDs.** The gmsh physical groups come out alternating
`Calving_0, Other_1, Calving_2, …`, auto-numbered `1,2,3,…`, so **odd tag =
calving, even tag = other**. `mesh_antarctica.py` automatically writes that
split to `mesh/boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json`
(read by every solver via `ISMIP7_BNDIDS`, which itself defaults using the
same `(COARSE, FINE, BUFFER_M)`-based naming — see `scripts/mesh_naming.py`).
These small JSONs are the *only*
tracked files in `mesh/`. To regenerate one for an existing mesh without
rebuilding it, run `ISMIP7_BUFFER_M=<N> python scripts/make_boundary_ids.py`
(or pass explicit `ISMIP7_MESH`/`ISMIP7_BNDIDS` paths) — it parses the mesh's
`$PhysicalNames` block and writes `{"calving":[odd…], "other":[even…]}`.

---

## 4. Invert for basal/rheology fields (`inversion_icepack2.py`)

MAP estimate of the bed friction `θ` and rheology `φ` from the diagnostic
3-field (V × Σ × τ) system, regularized, n=1→3 continuation, via `tlm_adjoint`.

```bash
cd antarctica
ISMIP7_LC=2500 mpiexec -n 12 python scripts/inversion_icepack2.py
# → mesh/inversion_icepack2_budd_n3_<LC>.h5   (the MAP checkpoint every forward run loads)
```

The controls are log-deviations from physical prior means: `θ = log(C/C_w0)`
on the balance-friction anchor and `φ = log(A/A_prior)` on a thermomechanical
fluidity prior the inversion computes at setup (and stores in the MAP). That
prior solve reads the §1 RACMO SMB **and the §2 ISMIP7 `tas` climatology**, so
download the §2 forcing before inverting; see `N3_FRAMEWORK.md` for the
method (and `ISMIP7_FLUIDITY_PRIOR=legacy` to skip it).

The friction law is selected with `ISMIP7_FRICTION` (`budd`, the default, or
`regularized_coulomb`); the MAP checkpoint name carries a matching `_budd` /
`_rc` tag, and the forward runs load the checkpoint for whichever law they are
started with. This `antarctica-n3` branch also appends an `_n3` flow-exponent
tag (from `map_n_tag()`) so n=3 and n=4 MAPs coexist on disk; see
`N3_FRAMEWORK.md`. See the script header for the full set of regularization /
iteration options.

---

## 5. Calibrate ocean melt (`calibrate_melt.py`)

Solves for the Burgard quadratic-mixed-slope coefficient **K** (global and
per-IMBIE2-basin) by matching integrated observed shelf melt
(Paolo/Adusumilli ≈ 865 Gt/yr). Needs §2 forcing + §4 inversion mesh.

```bash
cd antarctica
ISMIP7_LC=2500 python scripts/calibrate_melt.py
# → results/calibrated_K_per_basin_<LC>.npz
```

The control run (§6) requires this `.npz`. Projections can either use it
(`K_per_basin_npz=`) or a scalar `ISMIP7_K_MELT`.

---

## 6. Run control & projections

All forward runs go through `scripts/simulation.py` (`setup_model` +
`run_simulation`). Drivers live in `scripts/control/` and `scripts/projections/`.

```bash
cd antarctica

# Control (CTRL2015): fixed 2000–2029 SMB climatology + OI ocean climatology,
# per-basin calibrated K, melt recomputed each step from evolving geometry.
mpiexec -n 12 python scripts/control/run.py
# → results/ctrl2015_<esm>_<lc>_{final.h5, t<year>.h5, timeseries.csv}

# Core Experiment 7: SSP5-8.5 / CESM2-WACCM, 2015–2300
mpiexec -n 12 python scripts/projections/ssp585_cesm_waccm.py
# → results/ssp585_cesm2_waccm_<lc>_{final.h5, timeseries.csv}
```

Other scenario drivers in `scripts/projections/` (ssp126/ssp370 × CESM2-WACCM /
MRI-ESM2-0, plus `ocx.py`) are thin shims over `scripts/experiment.py` and
follow the same pattern. Historical spin-up drivers are in `scripts/historical/`;
run one first to produce `results/hist_<esm>_<lc>_final.h5`, which the
projections pick up automatically via `ISMIP7_RESTART` (otherwise they
cold-start from BedMachine).

Restart / run-management flags on the control driver: `--restart <ckpt>`
(or `ISMIP7_RESTART`) resumes from a checkpoint; `ISMIP7_AUTO_RESUME=1`
picks up the newest checkpoint for the experiment unattended; `--tag`
(or `ISMIP7_RUN_TAG`) suffixes the experiment name so a tagged run keeps -
and resumes - its own output files; `--checkpoint-interval` sets the
step-count fallback cadence. Checkpoints are self-contained (mesh, geometry,
inversion fields, and the full `(u, M, τ)` solver state), so restarts are
seamless at any MPI rank count. Each checkpoint also records the friction law
and whether an apparent-MB correction was active; a resume refuses to start
(with a message naming the fix) if `ISMIP7_FRICTION` or `ISMIP7_APPARENT_MB`
doesn't match the checkpoint, since a silent mismatch would run cleanly but
produce wrong physics.

**Is the run on track?** Audit any timeseries CSV against observed Antarctic
budget envelopes (IMBIE dM/dt, Rignot melt/calving, RACMO SMB, ISMIP6-class
control drift) plus a runaway detector:

```bash
python scripts/check_ismip6_track.py results/<exp>_timeseries.csv
# exit code 0 iff no FAIL rows, so launch gates can chain on it
```

### Environment knobs (all forward runs)

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_LC` | fine mesh resolution tag (selects mesh + inversion h5) | `2500` |
| `ISMIP7_LC_COARSE` | coarse mesh tag | `64000` |
| `ISMIP7_BUFFER_M` | outline buffer (m) used to resolve the default mesh/boundary-id filenames (see §3) | `20000` |
| `ISMIP7_MESH` | override mesh `.msh` path | `mesh/antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.msh` |
| `ISMIP7_BNDIDS` | override boundary-id JSON | `mesh/boundary_ids_antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.json` |
| `ISMIP7_INVERSION` | override inversion checkpoint (θ/φ MAP); fields are interpolated onto `ISMIP7_MESH` when set | `mesh/inversion_icepack2_budd_n3_<LC>.h5` |
| `ISMIP7_DATA_ROOT` | ISMIP7 forcing tree root | `<repo>/ISMIP7/AIS` |
| `ISMIP7_T_END` / `ISMIP7_DT` | end year / timestep (yr) | `2300` / `1.0` |
| `ISMIP7_FRICTION` | friction law (`budd`, `regularized_coulomb`) - selects the MAP h5 | `budd` |
| `ISMIP7_OUTPUT_INTERVAL` | write a timeseries/log row every N steps | `10` |
| `ISMIP7_CHECKPOINT_EVERY_YR` | checkpoint cadence in model years (`0` = use step count) | `5` |
| `ISMIP7_KEEP_CHECKPOINTS` | periodic checkpoints kept on disk (plus `_final.h5`) | `3` |
| `ISMIP7_RESTART` | restart checkpoint | `hist_<esm>_<lc>_final.h5` if present |
| `ISMIP7_AUTO_RESUME` | set to resume from the newest own checkpoint unattended | _(unset)_ |
| `ISMIP7_APPARENT_MB` | apparent-mass-balance init: `1`/`balance` zeroes the t=0 thickness tendency (ISMIP6 ctrl_proj-style), `div` cancels only the flux divergence | _(unset)_ |
| `ISMIP7_FIXED_FRONT` | set to hold the calving front at the t=0 extent (inflow beyond it tallied as calving) | _(unset)_ |
| `ISMIP7_LEGACY_TRANSPORT` | set to restore the pre-Jul-2026 CG-projection transport scheme | _(unset)_ |
| `ISMIP7_SNES_TYPE` / `ISMIP7_SNES_MAXIT` | diagnostic Newton type / max iterations | `newtonls` / `200` |
| `ISMIP7_K_MELT` | scalar Burgard K (projections) | `1.15e-4` (Burgard K50) |
| `ISMIP7_K_PER_BASIN_NPZ` | per-basin K file (control) | `results/calibrated_K_per_basin_<lc>.npz` |
| `ISMIP7_ESM` | ESM for control (`CESM2-WACCM`, `MRI-ESM2-0`) | `CESM2-WACCM` |
| `ISMIP7_H_CLAMP` | thickness floor (m) | `0` |
| `ISMIP7_NO_CALVING_TERMINUS` | set to drop the calving-terminus BC | _(unset)_ |

> **dt guidance** (from a dt-convergence sweep): use `ISMIP7_DT=0.1` for
> production projections; `0.25` is acceptable if 10 steps/yr is too costly for a
> 285-yr run. `dt=1.0` over/under-melts per step and resurrects clamped cells.

---

## 7. Timing benchmark (`make timing`)

Resolution vs. core-count wall-clock benchmark for the transient solver. Run
from `antarctica/` (needs §1 data, Firedrake, and Slurm on the cluster):

```bash
cd antarctica
make timing
# → TIMING_MATRIX.md
# → results/timing/timing_<LC>_<LC_coarse>_<ncores>.json  (one per cell)
```

The pipeline has five idempotent stages (each skips work already done):

1. **Meshes** — build 10 meshes at `buffer=20000 m`: LC ∈ {500, 1000, 2000,
   2500, 5000} m with LC_coarse = 10×LC and 20×LC.
2. **Inversion** — one inversion at LC=2500 / LC_coarse=25000 (via
   `scripts/batch_runners/timing_inversion.script` on Slurm, or `mpiexec`
   locally).
3. **Redistribute** — rewrite the complete tagged MAP checkpoint on a single
   core (`scripts/redistribute_checkpoint.py`) so any rank count can load it;
   fields, root metadata, and mesh provenance are retained.
4. **Transient** — 30 short runs (10 mesh combos × 16/32/64 cores): 5 years
   (`2015`–`2020`, `dt=1.0`), zero SMB/melt forcing. Every resolution
   loads its own mesh via `ISMIP7_MESH` and warm-starts θ/φ and the physical
   prior from the single inversion via cross-mesh interpolation
   (`ISMIP7_INVERSION=mesh/inversion_icepack2_budd_n3_2500.h5`).
5. **Matrix** — aggregate JSON timing records into [`TIMING_MATRIX.md`](TIMING_MATRIX.md).

Individual stages can be run separately: `make meshes`, `make inversion`,
`make redistribute`, `make transient`, `make matrix`.

Slurm scripts live in `scripts/batch_runners/` (`timing_inversion.script`,
`timing_meshes.script`, `timing_redistribute.script`, and
`timing_transient.script`), following the Quartz module-load pattern used by
the other batch scripts. On Quartz, run `make timing` from the `antarctica/`
directory; mesh generation and checkpoint redistribution are submitted as
single-rank jobs, while the transient jobs use 16, 32, and 64 ranks.

---

## Outputs

Per experiment in `results/`:
- `<exp>_final.h5` — final state checkpoint (Firedrake `CheckpointFile`),
  self-contained for restart: mesh, geometry, inversion fields, DG0 thickness
  (`thickness_dg`, the prognostic state), the full `(u, M, τ)` solver state,
  and the frozen apparent-MB reference when one is active.
- `<exp>_t<year>.h5` — periodic checkpoints (every `ISMIP7_CHECKPOINT_EVERY_YR`
  model years, default 5; only the newest `ISMIP7_KEEP_CHECKPOINTS` are kept).
- `<exp>_timeseries.csv` — one row per `OUTPUT_INTERVAL` steps with columns
  `year, vaf_mm_sle, mass_gt, smb_gtyr, melt_gtyr, outflux_gtyr, calv_gt,
  clamp_gt, resid_gt, amb_gtyr`: the mass-budget audit (SMB, shelf melt,
  boundary outflux, fixed-front calving, clamp/limiter corrections, the
  apparent-MB source, and the budget residual, which must close to 0.00).

VAF is reported in mm of sea-level equivalent; mass in Gt.

---

## Known issues (read before trusting a long run)

- **Diagnostic-Newton wall on hard projection geometries.** The earlier
  forward blow-ups are fixed (balanced apparent-MB init + persistent DG0
  thickness state; the 32 km CTRL and ssp585 both audit ON TRACK via
  `check_ismip6_track.py`), but the diagnostic Newton can still stall on
  evolved projection geometries (first seen at the 2024.5 ssp585 state).
  `ISMIP7_SNES_TYPE` / `ISMIP7_SNES_MAXIT` are the knobs for experimenting;
  a robust fix is the next work item.
- **Upstream forcing moved (resolved 2026-07-19).** The per-year scenario
  forcing was not withdrawn - it moved to the top-level `/ISMIP7/AIS` tree
  during the collection reorganization. Mirror it with
  `scripts/download_forcing.py --scenarios` (§2a) and run
  `scripts/preflight.py` to see what can run locally.
- **Mesh/boundary-id naming is now per-`(COARSE, FINE, BUFFER_M)`.** Mesh and
  sidecar filenames are tagged with the exact resolution and outline buffer
  used to build them (`antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.msh` /
  `boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json`, see §3),
  so a sidecar can no longer silently mismatch a mesh built with a different
  resolution or buffer. If you have older meshes/sidecars built before this
  naming convention, rename them to match or rebuild via `mesh_antarctica.py`.
- **`icepack2_tools/coupled.py`** is a WIP sketch of ice↔plume coupling and
  references a `PlumeModel` that does not yet exist in this tree — not wired into
  any run.
- **Timing runs inherit any diagnostic-Newton convergence failure.** If a
  transient job stops in its rescue ladder, rerun that cell after inspecting
  its `timing_*.txt` / `timing_*.err` files; the matrix builder leaves missing
  cells as `—` rather than inventing a timing value.

---

## References

- Burgard et al. 2022, *The Cryosphere* — basal-melt parameterisation assessment.
- multimelt (the reference implementation): https://github.com/ClimateClara/multimelt
- ISMIP7 ocean forcing pipeline: https://github.com/ismip/ismip7-antarctic-ocean-forcing
- Greenland companion: https://github.com/dlilien/ISMIP7_Greenland_Icepack
