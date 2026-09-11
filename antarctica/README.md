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
  reports/         # tracked per-core run records + MATRIX_STATUS.md
  scripts/         # all entry points (see §3–§6)
ISMIP7/AIS/        # ISMIP7 forcing tree the *runtime* reads             [gitignored]
icepack2_tools/    # reusable library (mesh, forcing, eikonal, grounding, regrid)
```

Everything large is gitignored. The only tracked files under `antarctica/mesh/`
are the tiny per-mesh `boundary_ids_antarctica_*.json` (the gmsh physical-line →
calving/other map; see §3). The legacy sidecars that carry no mesh stem
(`boundary_ids.json`, `_2500`, `_aniso`, `_buffered`) are deliberately **not**
tracked: they predate the per-mesh convention, and `boundary_ids.json` is the
shared fallback that every mesh build overwrites - committing it would put one
build's map on the fallback path for every clone.

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
`fracture/` and inside `fracture/v*/` (highest version wins).

**Per-year atmosphere files are monthly.** Each `<var>_..._<YEAR>.nc` holds 12
slices (`time` = days since `<YEAR>-01-15`), so the reader collapses the time
axis to that year's **annual mean**, weighting months by length from
`time_bnds` (falling back to the time-coordinate spacing, then to an unweighted
mean with a warning). A length-1 time axis passes through unchanged, so annual
files are unaffected. The forcing API:

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
(the name is built by `scripts/mesh_naming.py`, the single owner of the mesh /
sidecar filename convention). These small JSONs are the *only* tracked files
in `mesh/`. To regenerate one for an existing mesh without rebuilding it, run
`ISMIP7_BUFFER_M=<N> python scripts/make_boundary_ids.py`
(or pass explicit `ISMIP7_MESH`/`ISMIP7_BNDIDS` paths) — it parses the mesh's
`$PhysicalNames` block and writes `{"calving":[odd…], "other":[even…]}`.

Every solver resolves the sidecar through `icepack2_tools/boundary.py`, in
order: `ISMIP7_BNDIDS` if set, else the **per-mesh**
`mesh/boundary_ids_<mesh stem>.json`, else the shared `mesh/boundary_ids.json`.
The per-mesh name is what makes a stale sidecar impossible to pick up by
accident: the shared file is overwritten by every mesh build, and a sidecar
built for another mesh is not valid here (the physical-group count depends on
resolution *and* on the outline buffer). Every reader **hard-errors** when an
exterior marker on the mesh is unclassified or a sidecar id is absent from it -
a mismatch used to be legal and silent, and left ~95% of the ice front with no
calving back-pressure. Runs print the covered front length in km and percent.

A forward run does not take this name from its environment: the MAP / restart
checkpoint records the `.msh` it was written on (`mesh_basename`, plus the
`lc` / `lc_coarse` / `buffer_m` parameters to rebuild the name), and that
record wins, so `ISMIP7_LC` / `ISMIP7_LC_COARSE` / `ISMIP7_BUFFER_M` drifting
between the inversion and the forward cannot swap the sidecar underneath a
trajectory.

---

## 4. Invert for basal/rheology fields (`inversion_icepack2.py`)

MAP estimate of the bed friction `θ` and rheology `φ` from the diagnostic
3-field (V × Σ × τ) system, regularized, n=1→3 continuation, via `tlm_adjoint`.

```bash
cd antarctica
ISMIP7_LC=2500 mpiexec -n 12 python scripts/inversion_icepack2.py
# → mesh/inversion_icepack2_budd_n3_dg0_<LC>.h5   (the MAP checkpoint every forward run loads)
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
`N3_FRAMEWORK.md`. Since Aug 2026 the name additionally carries the **geometry
space** it was inverted under (`_dg0` for the `ISMIP7_GEOMETRY_SPACE=dg0`
default, untagged for `cg1`): the inversion absorbs the calving-front treatment
into `θ`/`φ`, so a MAP is only valid for its own geometry space. The whole name
is built by one helper (`icepack2_tools/naming.py`), which the forward,
`preflight.py` and the launch gates all import. A forward that finds only the
legacy untagged MAP loads it and warns loudly - runnable as a smoke test, not a
result. See `../GEOMETRY_DISCRETIZATION.md`. The regularization,
misfit-normalization and iteration knobs are tabulated under "Environment knobs
(inversion)" below.

### Transient (dH/dt-constrained) inversion

A velocity-only inversion fits `u` but never constrains `div(h u)`, so the MAP
can carry a flux divergence inconsistent with the observed geometry and the
forward then drifts (which `ISMIP7_APPARENT_MB` has been masking).
`ISMIP7_DHDT_WEIGHT > 0` adds **one implicit-Euler prognostic step** after the
diagnostic solve - the model's own DG0 upwind FV operator, so the inversion is
penalised for the divergence its *own* transport scheme produces - and scores
the resulting tendency against an observed mean dH/dt map
(`icepack2_tools/obs_dhdt.py`, from the ISMIP7 observations MIPkit). It
requires `ISMIP7_GEOMETRY_SPACE=dg0` and is restricted to grounded ice; see the
`obs_dhdt.py` module docstring for the target field, its 2003-2019 protocol
asterisk, and why the pixels are binned rather than point-sampled.

```bash
ISMIP7_DHDT_WEIGHT=1.0 ISMIP7_MAP_OUT=mesh/inversion_transient_2500.h5 \
  ISMIP7_LC=2500 mpiexec -n 12 python scripts/inversion_icepack2.py
python scripts/compare_dhdt.py vel=mesh/<velocity-only>.h5 tr=mesh/<transient>.h5
```

`compare_dhdt.py` is the payoff diagnostic: the t=0 *velocity* misfit cannot
tell a velocity-only MAP from a transient one (both fit `u`), so the thickness
tendency is the observable that can. It drives that step with the `velocity`
stored in the MAP, so its score is only as good as that field. The inversion
refuses to save a final-solve velocity whose misfit grossly disagrees with the
last accepted optimization state - it warns and leaves the field out - so a MAP
from a failed final solve has no `velocity` to score. MAPs written before that
guard can still carry one, and a dH/dt score built on it is meaningless. The
controls are unaffected either way: a forward re-solves the diagnostic from
`θ`/`φ` and never reads the stored velocity.

The per-cell dH/dt term above still leaves the *integrated* mass trend free.
`ISMIP7_DHDT_NET_SIGMA > 0` adds a second term penalising the net
grounded+observed dH/dt integral directly. It is **off by default, on purpose**:
the integrated trend is the one number every downstream assessment
(IMBIE/GRACE consistency) checks, so assimilating it forfeits it as
*independent* validation. The per-iteration `net=` diagnostic prints either
way, so the bias stays visible without being penalised. See
`reports/ISSUE_DRAFT_net_mass_balance_term.md`.

Because the MAP filename encodes only friction, `LC`, geometry space and flow
exponent, a velocity-only MAP and a transient one land on the **same path**.
Give variants their own `ISMIP7_MAP_OUT`. Every MAP checkpoint also records the
objective that produced it as root attributes - `misfit_norm`, `gamma_theta`,
`gamma_phi`, `dhdt_weight`, `dhdt_net_sigma` (the *resolved* value: 0 whenever
the term was not actually built), alongside `mesh_basename` and the
`lc`/`lc_coarse`/`buffer_m` mesh parameters - so a MAP already on disk can be
identified:

```bash
python -c "import h5py,sys; print(dict(h5py.File(sys.argv[1])['/'].attrs))" MAP.h5
```

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
run one first to produce `results/hist_<esm>_<lc>_final.h5`: the projections
AND the control both branch from it automatically, so they share the same t=0
state and their shared relaxation drift (and the identical frozen apparent-MB
correction) cancels in projection-minus-control (the ISMIP6 ctrl_proj
convention). Without it a projection cold-starts from BedMachine, and the
control cold-starts with a loud warning that it starts from a DIFFERENT
geometry than the hist-branched projections, so projection-minus-CTRL will
not cleanly isolate the forced response.

Restart / run-management flags on the control driver: `--restart <ckpt>`
(or `ISMIP7_RESTART`) resumes from a checkpoint; `ISMIP7_AUTO_RESUME=1`
picks up the newest checkpoint for the experiment unattended; `--tag`
(or `ISMIP7_RUN_TAG`, honored by every forward driver, not just the control)
suffixes the experiment name so a tagged method line (e.g. the n=3 matrix)
keeps - and resumes - its own output files, with the historical → projection
/ CTRL restart chain staying within that line; `--checkpoint-interval` sets
the step-count fallback cadence. Checkpoints are self-contained (mesh, geometry,
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

For the forced response, `scripts/compare_ismip6.py <proj.csv> <ctrl.csv>` overlays our projection-minus-CTRL sea-level contribution on the ISMIP6 ensemble to check broad consistency (auto-selects the scenario pool; `--exps` to override).

### The whole matrix in one command (`run_core_matrix.sh`)

`scripts/run_core_matrix.sh` runs core experiments 1-11 end to end in protocol
dependency order (both historicals first, since the CTRLs and the projections
branch from the historical endpoint, then the CTRLs, the projections, OCX) and
finishes with the observational audit and the ensemble comparison above for
each core it brought to its target year.

```bash
ISMIP7_RUN_TAG=n3 ISMIP7_LC=32000 antarctica/scripts/run_core_matrix.sh
CORES=1,2,9 antarctica/scripts/run_core_matrix.sh    # a subset
```

Cores run **sequentially** and are load-gated (`MAX_LOAD`, default
cores - 8), so the machine stays usable for other work. A core whose
timeseries already reaches its target year is skipped; output that predates
the annual-mean atmosphere-forcing fix is archived under
`results/archive_stale_<stamp>/` and re-run rather than reused. Resolution,
flow exponent, run tag, dt, friction and rank count come from the `ISMIP7_*`
knobs below; the runner's own knobs (`CORES`, `MAX_LOAD`, `MAX_ATTEMPTS`,
`FRESH`, `REUSE`, `NRANKS`, `PROV_REF`) and the exact reuse/dependency rules
are documented in its header, which is their authoritative reference. Its one
non-obvious behavior, wall retry, is described under "Known issues" below.
The runner explicitly pins `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps`; it
does not inherit the development default. Override it only with a solver mode
that has passed `make qualify` at the target configuration.

Each completed core is recorded with
`python scripts/core_report.py --core <N> --name <exp> --csv <timeseries.csv>
--log <run.log>` (add `--ctrl-csv` for a projection), run **in the run's own
shell** so it captures that run's `ISMIP7_*` environment. It writes the
tracked markdown record under `reports/`; `reports/MATRIX_STATUS.md` carries
the matrix-wide status. `--superseded "<reason>"` stamps an existing record
with a validity banner when a later run replaces it, so a superseded result
cannot be read as current.
The report resolves unset solver defaults and embeds the complete diagnostic
and transport PETSc dictionaries, continuation/rescue settings, tolerance
policy, and requested-versus-canonical solver mode.

### Environment knobs (inversion)

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_MAP_OUT` | full output path for the MAP h5, overriding the generated name. Use it for smoke tests and variant inversions so a short run cannot replace a converged production MAP. A bare filename resolves under `mesh/`; the directory is created and probed for writability at startup | _(generated name)_ |
| `ISMIP7_MISFIT_NORM` | `sigma`: divide each residual by its own datum's squared error, making the misfit a dimensionless chi^2 so terms of different units can be traded off. `none`: legacy dimensional misfit. **Selects the `ISMIP7_GAMMA_*` defaults** (see below) | `sigma` |
| `ISMIP7_GAMMA_THETA` / `ISMIP7_GAMMA_PHI` | Whittle-Matern prior strength on `θ` / `φ`. Default is coupled to `ISMIP7_MISFIT_NORM`, because normalizing divides the misfit by ~sigma^2 and would otherwise weaken the prior by the same factor | `1e5` under `sigma`, `1e4` under `none` |
| `ISMIP7_L_REG` | prior correlation length (m) | `7.5e3` |
| `ISMIP7_MAXITER` | L-BFGS-B iteration cap | `500` |
| `ISMIP7_GRAD_PRECOND` | optimization metric. `none` is the raw-dof Euclidean l2 metric, which is **mesh-dependent**: a gradient entry scales with its dof's cell area, so the fine grounding-line cells converge slowest. `mass` optimizes in `u = sqrt(M) x` (M = lumped mass), i.e. steepest descent in L2, which makes the convergence rate mesh-independent. Defaults to `none` **deliberately**, so runs stay comparable with everything measured so far; flip after the current A/B. Any other value aborts at startup | `none` |
| `ISMIP7_SIGMA_U_FLOOR` | floor on the per-component MEaSUREs velocity error (m/yr); without it the near-zero errors let a few nodes dominate the functional | `1.0` |
| `ISMIP7_SIGMA_U_UNOBS` | sigma (m/yr) assigned where MEaSUREs reports no error at all. Those nodes also have a zero-filled `u_obs`, so they must be given a *large* sigma, not the floor, or they would carry maximal weight on a fabricated zero velocity when `ISMIP7_OBS_MASK=0` | `1e4` |
| `ISMIP7_OBS_MASK` | `0` drops the velocity-observation mask (unobserved nodes re-enter the misfit) | `1` |
| `ISMIP7_DHDT_WEIGHT` | weight on the dH/dt chi^2 term; `0` disables the transient constraint entirely (requires `ISMIP7_GEOMETRY_SPACE=dg0` when on) | `0` |
| `ISMIP7_DHDT_SIGMA` | assumed dH/dt uncertainty (m/yr). A hand-set scalar: the MIPkit ships **no** uncertainty field for either dH/dt product | `0.1` |
| `ISMIP7_DHDT_DT` | timestep of the single prognostic step (yr) | `1.0` |
| `ISMIP7_DHDT_VAR` | observed field: `dhdt_smith` (firn-corrected, 2003-2019 mean) or `dhdt_cpom` (**not** firn corrected, so not interchangeable) | `dhdt_smith` |
| `ISMIP7_DHDT_MELT` | `0` drops ocean melt from the prognostic step's source (SMB only). On by default so the step matches the forward's forcing. The per-basin K comes from `ISMIP7_K_PER_BASIN_NPZ`, else `results/calibrated_K_per_basin_<lc>.npz`, else the 2500 m file; if none exists this warns and falls back to SMB-only rather than aborting, since melt is zero on grounded ice and the misfit is grounded-only | `1` |
| `ISMIP7_DHDT_CLIM_START` / `_END` | RACMO SMB climatology window for that source | `2003` / `2019` |
| `ISMIP7_DHDT_REACH` | pixel-to-cell reach as a multiple of the cell scale `sqrt(area)`; rejects raster pixels lying outside the mesh that nearest-centroid assignment would otherwise snap onto boundary cells | `0.75` |
| `ISMIP7_DHDT_NET_SIGMA` | sigma (Gt/yr) on the *integrated* grounded+observed dH/dt; `0` disables the net mass-balance term. Off by default on purpose - see §4 - and only ever active when `ISMIP7_DHDT_WEIGHT > 0` | `0` |
| `ISMIP7_OBS_KIT` | path to `AntarcticaObsISMIP7-v*.nc` | newest under `<DATA_ROOT>/obs/mipkit` |

---

### Environment knobs (all forward runs)

The run-shaping knobs - `ISMIP7_LC`, `ISMIP7_LC_COARSE`, `ISMIP7_FRICTION`
and `ISMIP7_GEOMETRY_SPACE` below, plus `ISMIP7_N_FLOW` (see
`../COMPOSITE_RHEOLOGY.md`) - have exactly one owner in code,
`icepack2_tools/runconfig.py`. Every driver, probe and gate reads them through
it, so an unset knob cannot mean one resolution to the inversion and another
to the preflight. A run that wants something else exports it, which is also
how it reaches the core report.

PETSc and recovery defaults have the same single-owner rule:
`icepack2_tools/solverconfig.py` supplies the runtime dictionaries, the timing
metadata, and the core-report provenance block. `simulation.py` must not
redeclare those literals.

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_LC` | fine mesh resolution tag (selects mesh + inversion h5) | `2500` |
| `ISMIP7_LC_COARSE` | coarse mesh tag | `64000` |
| `ISMIP7_BUFFER_M` | outline buffer (m) used to resolve the default mesh/boundary-id filenames (see §3) | `20000` |
| `ISMIP7_MESH` | mesh `.msh` path (inversion and tools). A forward takes its mesh from the MAP/restart checkpoint, which records its own mesh basename and parameters, so here it only names the boundary sidecar for a legacy checkpoint that carries no such record | `mesh/antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.msh` |
| `ISMIP7_BNDIDS` | override boundary-id JSON | `mesh/boundary_ids_antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.json` if present, else `mesh/boundary_ids.json` |
| `ISMIP7_INVERSION` | override inversion checkpoint (θ/φ MAP); fields are interpolated onto `ISMIP7_MESH` when set | `mesh/inversion_icepack2_budd_n3_<LC>.h5` |
| `ISMIP7_OBS_DATA_ROOT` | BedMachine, MEaSUREs velocity, and RACMO observational-data root; useful when these files live on an external volume | `<repo>/antarctica/data` |
| `ISMIP7_GEOMETRY_SPACE` | space for `h`/`s`/`b` (`dg0`: one thickness for the terminus force and the mass flux; `cg1`: legacy, for A/B only) - also selects the MAP h5 (see `../GEOMETRY_DISCRETIZATION.md`) | `dg0` |
| `ISMIP7_DATA_ROOT` | ISMIP7 forcing tree root | `<repo>/ISMIP7/AIS` |
| `ISMIP7_T_END` / `ISMIP7_DT` | end year / timestep (yr) | `2300` / `1.0` |
| `ISMIP7_FRICTION` | friction law (`budd`, `regularized_coulomb`) - selects the MAP h5 | `budd` |
| `ISMIP7_OUTPUT_INTERVAL` | write a timeseries/log row every N steps | `10` |
| `ISMIP7_CHECKPOINT_EVERY_YR` | checkpoint cadence in model years (`0` = use step count) | `5` |
| `ISMIP7_KEEP_CHECKPOINTS` | periodic checkpoints kept on disk (plus `_final.h5`) | `3` |
| `ISMIP7_RESTART` | restart checkpoint | `hist_<esm>[_<tag>]_<lc>_final.h5` if present |
| `ISMIP7_AUTO_RESUME` | set to resume from the newest own checkpoint unattended | _(unset)_ |
| `ISMIP7_RUN_TAG` | experiment-name suffix for a parallel method line (see run-management flags above) | _(unset)_ |
| `ISMIP7_APPARENT_MB` | apparent-mass-balance init: `1`/`balance` zeroes the t=0 thickness tendency (ISMIP6 ctrl_proj-style), `div` cancels only the flux divergence | _(unset)_ |
| `ISMIP7_FIXED_FRONT` | set to hold the calving front at the t=0 extent (inflow beyond it tallied as calving) | _(unset)_ |
| `ISMIP7_LEGACY_TRANSPORT` | set to restore the pre-Jul-2026 CG-projection transport scheme (requires `ISMIP7_GEOMETRY_SPACE=cg1`) | _(unset)_ |
| `ISMIP7_SNES_TYPE` / `ISMIP7_SNES_MAXIT` | diagnostic Newton type / max iterations | `newtonls` / `200` |
| `ISMIP7_SNES_LINESEARCH` | line search used by `newtonls` | `nleqerr` |
| `ISMIP7_SNES_RTOL` / `ISMIP7_SNES_ATOL` / `ISMIP7_SNES_STOL` | initial nonlinear relative, absolute and step tolerances. After the initial/restart solve, the absolute tolerance follows the self-scaled policy below | `1e-8` / `1e-50` / `0` |
| `ISMIP7_SNES_DIVERGENCE_TOL` | residual-growth divergence threshold; negative disables this PETSc test | `-1` |
| `ISMIP7_SNES_ATOL_SCALE` / `ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE` | persistent absolute tolerance after a converged setup solve (`scale * achieved norm`) / after accepting a loaded hard-era state (`scale * loaded-state norm`) | `100` / `1e-6` |
| `ISMIP7_SNES_KSP_EW` | enable PETSc Eisenstat-Walker variable inner tolerance for an A/B test | `0` |
| `ISMIP7_DIAGNOSTIC_LINEAR_SOLVER` | `schur_gamg` or `schur_mumps`: legacy PETSc `selfp` approximation; `scpc_gamg` or `scpc_mumps`: exact cell-local Slate elimination and an assembled velocity solve; `full_mumps`: complete mixed-Jacobian reference. Legacy `iterative`/`mumps` aliases mean `schur_gamg`/`schur_mumps` | `full_mumps` for forward drivers and the core runner; `scpc_mumps` in the timing Makefile |
| `ISMIP7_KSP_RTOL` / `ISMIP7_KSP_MAXIT` | outer FGMRES relative tolerance / iteration limit for the iterative diagnostic mode | `1e-6` / `1000` |
| `ISMIP7_SNES_MONITOR` / `ISMIP7_SNES_LOG` | enable diagnostic SNES/KSP monitors and optionally route them to a file; KSP output includes short and true residual lines per iteration, with a header before each mixed solve | `0` / stdout |
| `ISMIP7_SOLVER_VIEW` | emit `snes_view`, outer `ksp_view`, and the SCPC condensed `ksp_view`; enabled by `make debug` and `make reference` to expose the actual block sizes and hierarchy | `0` |
| `ISMIP7_TRANSPORT_KSP_RTOL` / `ISMIP7_TRANSPORT_KSP_MAXIT` | GMRES relative tolerance / iteration limit for the persistent DG0 transport solver (`ismip7_transport_` PETSc prefix) | `1e-10` / `500` |
| `ISMIP7_MASS_RESIDUAL_TOL_GT` | fail-loud absolute tolerance for both the discrete transport identity and the complete step mass budget | `5e-5` Gt |
| `ISMIP7_K_MELT` | scalar Burgard K (projections) | `1.15e-4` (Burgard K50) |
| `ISMIP7_K_PER_BASIN_NPZ` | per-basin K file (control) | `results/calibrated_K_per_basin_<lc>.npz` |
| `ISMIP7_ESM` | ESM for control (`CESM2-WACCM`, `MRI-ESM2-0`) | `CESM2-WACCM` |
| `ISMIP7_CLIM_SCENARIO` / `ISMIP7_CLIM_START` / `_END` | reference-climate pool: the scenario pooled with `historical`, and the window, used both for the control's SMB climatology and for the projections' aSMB re-reference. `ssp126` is the protocol pool (cheat sheet, April 2026); the two uses share one owner (`icepack2_tools/climatology.py`) because a disagreement makes projection-minus-control difference two unrelated baselines. A partial pool warns rather than refusing, and the coverage line reaches the core report | `ssp126` / `2000` / `2029` |
| `ISMIP7_H_CLAMP` | thickness floor (m) | `0` |
| `ISMIP7_NO_CALVING_TERMINUS` | set to drop the calving-terminus BC | _(unset)_ |
| `ISMIP7_RESCUE_ENABLED` | permit a failed direct transient diagnostic solve to enter the continuation/trust-region/subcycle rescue ladder; set to `0` for strict timestep qualification | `1` |
| `ISMIP7_SUBCYCLES` | dt-subcycle rescue ladder: a step that fails the rescue solves rewinds its own advance and retries at `dt/m` for each `m` in this list | `1,4,16` |
| `ISMIP7_RESCUE_MAXIT` | Newton iteration cap on the rescue rungs (hard-era steps converge linearly and need the extra patience) | `600` |

> **dt guidance** (from a dt-convergence sweep): use `ISMIP7_DT=0.1` for
> production projections; `0.25` is acceptable if 10 steps/yr is too costly for a
> 285-yr run. `dt=1.0` over/under-melts per step and resurrects clamped cells.

---

## 7. Timing benchmark (`make timing`)

Resolution vs. core-count wall-clock benchmark for the transient solver. The
campaign is a staged state machine: it prepares an exact-mesh state, runs one
strict scout per mesh, and submits scaling lanes only after the corresponding
scout passes. Run from `antarctica/` (needs §1 data, Firedrake, and Slurm):

```bash
cd antarctica
make timing
# after the currently active stage settles, run the same command again
make timing
# synchronize from Quartz with trailing-slash source and destination paths
make sync-results
make matrix
# → TIMING_MATRIX.md, showing 20 configured lanes plus NOT PLANNED cells
# → results/timing/timing_<tag>_<LC>_<LC_coarse>_<ncores>.json
```

`make timing` first reuses a valid five-step `scpc_mumps` qualification (or
runs it once), creates the rank-independent improved inversion if needed, and
then advances only stages whose prerequisites already exist. It never
resubmits an active job or an accepted result. Failed lanes remain failed for
inspection; use `FORCE_TIMING=1` only for an intentional retry.

The stages and contracts are:

1. **Inversion provenance** — the only timing input is
   `inversion_icepack2_budd_n3_dg0_logvelnet_2500_1core.h5`, a serial repack of
   the imported improved `dg0_logvelnet` MAP. The redistribution launcher has
   distinct input/output arguments; it never rewrites the source in place.
   The campaign tag contains `dg0_logvelnet_cached_strict`, so an old-MAP
   record cannot satisfy this campaign.
2. **Prepared caches (`make timing-prepare`)** — one job for each of the ten
   `(LC, LC_coarse)` meshes performs the adaptive cold continuation and saves
   the complete mixed state at 2015.0 without a transport step. It includes
   geometry, velocity, membrane and basal stress, controls, physical priors,
   frozen reference fields, mesh identity, solver configuration, and source
   inversion checksum. The job then repacks the cache on one rank and publishes
   the HDF5 file and JSON manifest atomically. Mesh, inversion, physics, solver,
   or cache-schema changes invalidate it.
3. **Scouts (`make timing-scout`)** — one lowest-retained-core lane per mesh:
   32 ranks at 500 m and 16 ranks elsewhere. A scout passes only after five
   direct steps, the exact final year, no negative solve reason or rescue
   label, matching cache provenance, and transport and complete-step mass
   residuals no larger than `5e-5 Gt`.
4. **Scaling (`make timing-scale`)** — higher-core lanes are submitted only
   for meshes whose scout passes. A failed scout stamps its scaling lanes
   `BLOCKED BY SCOUT`; the jobs are not submitted.
5. **Strict transient timing** — all matrix lanes use `scpc_mumps`, disable
   rescue, and restrict subcycles to `1`. The first diagnostic, transport, or
   mass-budget failure ends the lane. The primary timer starts immediately
   before the five-step loop, after cache loading, solver construction, and
   transport setup; `setup_seconds` is reported separately. JSON records are
   atomically written on catchable success or failure and include solver
   histories, global mesh counts, extrema, residuals, phase, and completed
   steps. Slurm-only termination (including exit 137) remains distinct from a
   recorded numerical failure.
6. **Reporting and synchronization** — `make sync-results` uses
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
`make timing-scout`, `make timing-scale`, `make sync-results`, and
`make matrix`. `make timing-dry-run` prints the staged commands without
submitting jobs. `make transient` routes matrix work through the scout gate;
qualification/debug calls still use its direct compatibility path.
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
- `<exp>_final.h5` — final state checkpoint (Firedrake `CheckpointFile`),
  self-contained for restart: mesh, geometry, inversion fields, the full
  `(u, M, τ)` solver state, and the frozen apparent-MB reference when one is
  active. Under the `dg0` geometry default the saved `thickness` **is** the
  prognostic transport state; a `cg1` run additionally saves the separate DG0
  carrier as `thickness_dg`, since there the CG1 `thickness` is only its lift.
  The `geometry_space` and `mesh_basename` attributes record the discretization
  and the `.msh` the trajectory started on, so a restart resolves the same
  boundary sidecar; restarting into a different geometry space projects and
  warns loudly (see `../GEOMETRY_DISCRETIZATION.md`).
- `<exp>_t<year>.h5` — periodic checkpoints (every `ISMIP7_CHECKPOINT_EVERY_YR`
  model years, default 5; only the `ISMIP7_KEEP_CHECKPOINTS` most recently
  *written* are kept - by write time, not by highest year, so a re-run that
  rewinds to the historical endpoint keeps its own states).
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
  thickness state), but the diagnostic Newton can still stall on
  evolved projection geometries. `ISMIP7_SNES_TYPE` / `ISMIP7_SNES_MAXIT` are
  the knobs for experimenting; a robust fix is the next work item. The
  aSMB-forced walls seen so far are suspect: they predate the annual-mean
  atmosphere-forcing fix and may be forcing-induced rather than a solver
  limit - see `reports/MATRIX_STATUS.md` for which runs still stand.
  **Wall retry.** When the in-run rescue ladder is exhausted the run saves and
  stops short of its target year, and simply relaunching from that saved state
  clears the wall: a fresh process re-runs the n=1→n continuation at the
  loaded geometry, which the in-run ladder cannot do (3 of 3 observed walls
  resumed - ssp585-CESM at 2096.7, CTRL-CESM at 2268, CTRL-MRI at 2250 - and
  both CTRLs then reached 2300). `run_core_matrix.sh` does this automatically,
  relaunching from the newest checkpoint at or before the timeseries' last
  year while each attempt keeps advancing, and giving up on a stall.
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
- **Timing runs are deliberately strict.** A matrix lane does not enter the
  diagnostic rescue ladder or retry with transport subcycles. The first
  failure is written to its timing JSON/status stamp, and a failed scout blocks
  that mesh's scaling lanes. Inspect the recorded category and Slurm log before
  using `FORCE_TIMING=1` for an intentional retry.

---

## References

- Burgard et al. 2022, *The Cryosphere* — basal-melt parameterisation assessment.
- multimelt (the reference implementation): https://github.com/ClimateClara/multimelt
- ISMIP7 ocean forcing pipeline: https://github.com/ismip/ismip7-antarctic-ocean-forcing
- Greenland companion: https://github.com/dlilien/ISMIP7_Greenland_Icepack
