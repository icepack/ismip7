r"""compare_scalars.py: the organisers' scalar processing set against the
model's own scalars (issue #13).

A 6 x 5 grid over three years stands in for a submission. The model's scalars
are built the way the writer's conservative remap makes them hold: a grid sum
of a flux is the mesh sum, since every flux is a whole-pixel mean. A tree
written before the whole-pixel means (``whole_pixel=False``) holds acabf as a
mean over the covered part of a pixel and libmassbffl over the floating part.
On top of that the model books a known amount of melt where the grid holds no
floating ice and keeps a known mass in cells thinner than lithk shows, and
books front melt as lifmassbf in a pixel with no floating ice (issue #109),
where no fill touches it. The
tool's output comes from a transcription of its own expressions
(ismip7_scalars 0.1.0, scalars.py and slc/), including A2020's step by step
accumulation, and from the tool itself where it is installed.
"""
import csv
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "antarctica" / "scripts"))
import compare_scalars as cs  # noqa: E402

NX, NY, DX = 6, 5, 8000.0
X = -20000.0 + DX * np.arange(NX)
Y = -16000.0 + DX * np.arange(NY)
YEARS = (2015, 2016, 2017)
UNITS = "days since 1850-01-01"
TAG = "AIS_RICE_icepack2_m001_CESM2-WACCM_f001_ssp585_C007_2015-2017"
FILL = netCDF4.default_fillvals["f4"]
EXTRA_MELT = -2.5e5      # kg/s the model books off ice floating at year end
THIN = 3.0e9             # m3 of ice in cells of 1 m or less
FRONT_MELT = -4.0e-6     # kg m-2 s-1 of lifmassbf in a grounded pixel's empty cells
FL_VARS = ("acabf", "libmassbfgr", "libmassbffl", "licalvf", "lifmassbf", "ligroundf", "dlithkdt")
# the stand-in for af2_AIS_08000m_v1.nc
AF2 = (1.02 + 0.01 * np.arange(NX)[None, :] + 0.002 * np.arange(NY)[:, None]).astype("f4")


def stamp(year, month, day):
    import cftime
    return float(cftime.date2num(cftime.datetime(year, month, day, calendar="standard"),
                                 UNITS, calendar="standard"))


def st_times():
    return np.array([stamp(y + 1, 1, 1) for y in YEARS], dtype="f4")


def fl_times():
    return np.array([stamp(y, 7, 1) for y in YEARS], dtype="f4")


def write_gridded(folder, var, cube, flux, attrs=None):
    with netCDF4.Dataset(folder / f"{var}_{TAG}.nc", "w") as ds:
        for k, v in (attrs or {}).items():
            ds.setncattr(k, v)
        ds.createDimension("time", None); ds.createDimension("y", NY); ds.createDimension("x", NX)
        ds.createVariable("x", "f8", ("x",))[:] = X
        ds.createVariable("y", "f8", ("y",))[:] = Y
        t = ds.createVariable("time", "f4", ("time",))
        t.units, t.calendar, t.long_name = UNITS, "standard", "time"
        t[:] = fl_times() if flux else st_times()
        if flux:
            ds.createDimension("nv", 2)
            t.bounds = "time_bnds"
            ds.createVariable("time_bnds", "f4", ("time", "nv"))[:] = [
                [stamp(y, 1, 1), stamp(y + 1, 1, 1)] for y in YEARS]
        v = ds.createVariable(var, "f4", ("time", "y", "x"), fill_value=np.float32(FILL))
        arr = np.asarray(cube, dtype="f4").copy()
        arr[~np.isfinite(arr)] = FILL
        v[:] = arr


def write_series(folder, var, values, flux, time_dtype="f4", tag=TAG, attrs=None):
    with netCDF4.Dataset(folder / f"{var}_{tag}.nc", "w") as ds:
        for k, v in (attrs or {}).items():
            ds.setncattr(k, v)
        ds.createDimension("time", None)
        t = ds.createVariable("time", time_dtype, ("time",))
        t.units, t.calendar, t.long_name = UNITS, "standard", "time"
        t[:] = fl_times() if flux else st_times()
        ds.createVariable(var, time_dtype if time_dtype == "f8" else "f4", ("time",))[:] = values


def write_grid(folder, name, arr, dtype, flip=False):
    with netCDF4.Dataset(folder / cs.GRID_FILES.get(name, "iaf2_GIC_AIS_{res}000m_v0.nc")
                         .format(res="08"), "w", format="NETCDF3_CLASSIC") as ds:
        ds.createDimension("x", NX); ds.createDimension("y", NY)
        ds.createVariable("x", "f8", ("x",))[:] = X
        ds.createVariable("y", "f8", ("y",))[:] = Y[::-1] if flip else Y
        ds.createVariable(name, dtype, ("y", "x"))[:] = arr


def write_params(path, rhow=1024.0):
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as ds:
        for k, v in (("rhoi", 917.0), ("rhow", rhow), ("rhof", 1000.0), ("oarea", 3.625e14)):
            ds.createVariable(k, "f8")[()] = v


def model_fields(moving_bed=False, front_melt=True):
    r"""Three years of a small ice sheet: a grounded interior, a shelf with
    one half-floating pixel, one half-covered pixel at the domain edge, and a
    column outside the domain (topg filled there). ``front_melt`` books
    lifmassbf in one pixel with no floating ice, as the forward has since
    issue #109; without it lifmassbf is zero, the booking before."""
    rng = np.random.default_rng(7)
    topg = np.full((NY, NX), 200.0)
    topg[:, 3:] = -600.0
    topg[:, 5] = np.nan                               # outside the domain
    cov = np.ones((NY, NX)); cov[:, 5] = 0.0; cov[2, 4] = 0.5
    f = {v: [] for v in ("lithk", "topg", "sftgrf", "sftflf") + FL_VARS}
    for k, yr in enumerate(YEARS):
        h = np.zeros((NY, NX))
        h[:, :3] = 2000.0 - 50.0 * k + 10.0 * rng.random((NY, 3))
        h[:, 3] = 800.0 - 20.0 * k                    # grounded on the marine bed
        h[:, 4] = 300.0 - 5.0 * k                     # floating
        h[2, 4] *= cov[2, 4]                          # whole-pixel mean of a half-covered pixel
        b = topg.copy()
        if moving_bed and k == 2:
            b[0, 0] += 5.0
        hf = np.maximum(-np.nan_to_num(b), 0.0) * 1024.0 / 917.0
        grounded = (h > 1.0) & (h > hf)
        floating = (h > 1.0) & ~grounded
        sftflf = floating.astype(float)
        sftflf[1, 4] = 0.5                            # half the pixel floats at year end
        f["lithk"].append(h); f["topg"].append(b)
        f["sftgrf"].append(grounded.astype(float)); f["sftflf"].append(sftflf)
        acabf = np.where(cov > 0, 1e-5 * (1 + rng.random((NY, NX))), np.nan)
        f["acabf"].append(acabf)
        f["libmassbffl"].append(np.where(sftflf > 0, -2e-5 * (1 + rng.random((NY, NX))), np.nan))
        f["libmassbfgr"].append(np.where(grounded, 0.0, np.nan))
        f["licalvf"].append(np.where(sftflf > 0, -1e-6 * rng.random((NY, NX)), 0.0))
        lif = np.zeros((NY, NX))
        if front_melt:
            lif[3, 3] = FRONT_MELT * (1 + k)             # grounded pixel, sftflf 0
        f["lifmassbf"].append(lif)
        lg = np.zeros((NY, NX)); lg[:, 4] = 3e-6; lg[0, 3] = -1e-6
        f["ligroundf"].append(lg)
        prev = f["lithk"][k - 1] if k else h
        f["dlithkdt"].append((h - prev) / cs.SECONDS_PER_YEAR)
    return {v: np.array(c, dtype="f4").astype(np.float64) for v, c in f.items()}, cov


def native_scalars(f, cov, weight=1.0):
    r"""What the model's mesh sums give, from the float32 grid as written.
    ``weight`` is the area factor a forward with true-area scalars carries."""
    A = DX * DX
    w = np.asarray(weight, dtype=np.float64)
    out = {s: [] for s in cs.ST_SCALARS + tuple(s for s, _ in cs.FL_SCALARS)}
    for k in range(len(YEARS)):
        h, b = f["lithk"][k], np.nan_to_num(f["topg"][k])
        hf = np.maximum(-b, 0.0) * 1024.0 / 917.0
        out["lim"].append(917.0 * (np.sum(h * w) * A + THIN))
        out["limnsw"].append(917.0 * np.sum(np.maximum(h - hf, 0.0) * w) * A)
        out["iareagr"].append(np.sum(f["sftgrf"][k] * w) * A)
        out["iareafl"].append(np.sum(f["sftflf"][k] * w) * A)
        z = {v: np.nan_to_num(f[v][k]) for v in FL_VARS}
        out["tendacabf"].append(np.sum(z["acabf"] * cov * w) * A)
        out["tendlibmassbfgr"].append(0.0)
        out["tendlibmassbffl"].append(np.sum(z["libmassbffl"] * f["sftflf"][k] * w) * A + EXTRA_MELT)
        out["tendlicalvf"].append(np.sum(z["licalvf"] * w) * A)
        out["tendlifmassbf"].append(np.sum(z["lifmassbf"] * w) * A)
        out["tendligroundf"].append(np.sum(z["ligroundf"] * w) * A)
    return out


def tool_scalars(sub, datapath, params, refyear=2016):
    r"""The tool's numbers, transcribed from ismip7_scalars 0.1.0 (scalars.py
    run, compute_st_series, compute_slc_series, _cumulative_a2020,
    run_fl_scalars; slc_vaf, slc_G2020, slc_A2020) for a run that is its own
    reference, reading the files the way it does: whole cubes, masked."""
    def cube(var):
        with netCDF4.Dataset(next(sub.glob(f"{var}_*.nc"))) as ds:
            return ds.variables[var][:, :, :]
    with netCDF4.Dataset(datapath / "af2_AIS_08000m_v1.nc") as ds:
        af2 = ds.variables["af2"][:, :]
    with netCDF4.Dataset(datapath / "maxmask1_AIS_08000m_v0.nc") as ds:
        maxmask1 = ds.variables["maxmask1"][:, :]
    with netCDF4.Dataset(datapath / "iaf2_GIC_AIS_08000m_v0.nc") as ds:
        gic = np.ones_like(ds.variables["iaf2"][:, :])
    with netCDF4.Dataset(params) as ds:
        RHOI, RHOSW, RHOFW = (float(ds.variables[k][()]) for k in ("rhoi", "rhow", "rhof"))
    AO, area_m2 = 3.625e14, (float("08") * 1000.0) ** 2
    lithk, topg, sftgrf, sftflf = cube("lithk"), cube("topg"), cube("sftgrf"), cube("sftflf")
    region = maxmask1 * 0 + 1
    A = region * af2 * area_m2
    out = {s: [] for s in cs.SCALARS}
    for n in range(len(YEARS)):
        H = lithk[n, :, :] * maxmask1
        B = topg[n, :, :]
        hf = np.maximum(-B, 0) * RHOSW / RHOI
        out["lim"].append(np.sum(H * A) * RHOI)
        out["limnsw"].append(np.sum(np.maximum(H - hf, 0) * A) * RHOI)
        out["iareagr"].append(np.sum(sftgrf[n, :, :] * A))
        out["iareafl"].append(np.sum(sftflf[n, :, :] * A))
    W = region * (af2 * area_m2)
    for s, v in cs.FL_SCALARS:
        out[s] = list(np.einsum("nyx,yx->n", np.ma.filled(cube(v), 0.0), W))

    ref = [y + 1 for y in YEARS].index(refyear)
    H0 = lithk[ref, :, :] * maxmask1 * gic
    B0 = topg[ref, :, :]
    S0 = topg[ref, :, :] * 0.0

    def vaf(H, B, S):
        hf = np.maximum(S - B, 0.0) * RHOSW / RHOI
        return np.sum(np.maximum(H - hf, 0.0) * A)

    def g20(H0, H, B0, B):
        def v(H, B):
            return np.sum(np.maximum(H + np.minimum(B, 0.0) * RHOSW / RHOI, 0.0) * A)

        def pov(B):
            return np.sum((np.maximum(-B, 0.0) * A))

        def den(H):
            return np.sum((H * (RHOI / RHOFW - RHOI / RHOSW) * A))
        af = -(v(H, B) / AO * RHOI / RHOSW - v(H0, B0) / AO * RHOI / RHOSW)
        return af - (pov(B) / AO - pov(B0) / AO) - (den(H) / AO - den(H0) / AO)

    def a20(H0, H, B0, B, S0, S):
        def masks(H, B, S):
            F = H - RHOSW / RHOI * (S - B)
            O = np.zeros_like(H); O[F < 0] = 1
            L = 1 - O
            I = np.zeros_like(H); I[H > 0] = 1
            return I, L, I * L
        I0, L0, G0 = masks(H0, B0, S0)
        I, L, G = masks(H, B, S)
        Hn = RHOSW / RHOI * np.maximum(S - B, 0)
        Hn0 = RHOSW / RHOI * np.maximum(S0 - B0, 0)
        HF, HF0 = G * (H - Hn), G0 * (H0 - Hn0)
        HM = (H - H0) * L0 * L + (HF - HF0) * (1 - L0 * L)
        HV = (1 - RHOFW / RHOSW) * ((H - H0) - (HF - HF0)) * (1 - L0 * L)
        return -RHOI / RHOFW * np.sum((HM + HV) * A) / AO

    acc, cumul = 0.0, [0.0]
    Hp, Bp = lithk[0, :, :] * maxmask1 * gic, topg[0, :, :]
    for n in range(1, len(YEARS)):
        Hh, Bh = lithk[n, :, :] * maxmask1 * gic, topg[n, :, :]
        acc += a20(Hp, Hh, Bp, Bh, S0, S0)
        cumul.append(acc)
        Hp, Bp = Hh.copy(), Bh.copy()
    for n in range(len(YEARS)):
        H, B = lithk[n, :, :] * maxmask1 * gic, topg[n, :, :]
        out["slvaf"].append(-(vaf(H, B, S0) / AO * RHOI / RHOFW - vaf(H0, B0, S0) / AO * RHOI / RHOFW))
        out["slg20"].append(g20(H0, H, B0, B))
        out["sla20"].append(cumul[n] - cumul[ref])
    return {s: np.asarray(v, dtype=np.float64) for s, v in out.items()}


def build(tmp_path, *, rhow=1024.0, flip=False, moving_bed=False, shift_tool=False,
          perturb_sftgrf=False, csv_off=False, ground_near=False, whole_pixel=True,
          native_af2=False, stamp=True, native_scale=1.0, front_melt=True,
          perturb_front=False):
    r"""A submission, its grids, params.nc, the tool's output and the
    writer's overlap cache; returns the argument list for compare_scalars.
    ``whole_pixel=False`` writes the tree the way the writer did before its
    flux means were whole-pixel means. ``native_af2`` gives the model's
    scalars the area factor, ``stamp=False`` leaves the writer's stamp of
    their area off, and ``native_scale`` puts them that factor off.
    ``front_melt=False`` writes lifmassbf as zero, and ``perturb_front`` puts
    the model's tendlifmassbf one percent off the grid's in one year."""
    sub = tmp_path / "tree" / "AIS" / "RICE" / "icepack2" / "CORE" / "C007"
    tool = tmp_path / "out" / "tool" / "nc" / "AIS" / "RICE" / "icepack2" / "CORE" / "C007"
    data = tmp_path / "grids"
    for d in (sub, tool, data):
        d.mkdir(parents=True)
    f, cov = model_fields(moving_bed=moving_bed, front_melt=front_melt)
    N = native_scalars(f, cov, weight=AF2 if native_af2 else 1.0)
    N = {s: [v * native_scale for v in vals] for s, vals in N.items()}
    if perturb_front:
        N["tendlifmassbf"][1] *= 1.01
    grid = {v: f[v].copy() for v in ("lithk", "topg", "sftgrf", "sftflf") + FL_VARS}
    if whole_pixel:
        # the same fluxes as whole-pixel means: fill stays fill (NaN * 0)
        grid["acabf"] = grid["acabf"] * cov
        grid["libmassbffl"] = grid["libmassbffl"] * f["sftflf"]
    if perturb_sftgrf:
        grid["sftgrf"][1, 0, 0] = 0.0
    if ground_near:
        # the writer's rule for a cell a few millimetres afloat: grounded in
        # the masks, so no floating melt there, while the model's own areas
        # still count it floating
        grid["sftflf"][1, 0, 4], grid["sftgrf"][1, 0, 4] = 0.0, 1.0
        grid["libmassbffl"][1, 0, 4] = np.nan
    for v in ("lithk", "topg", "sftgrf", "sftflf"):
        write_gridded(sub, v, grid[v], flux=False)
    for v in FL_VARS:
        write_gridded(sub, v, grid[v], flux=True,
                      attrs={cs.FLUX_MEAN_ATTR: cs.WHOLE_PIXEL} if whole_pixel else None)
    # the writer takes each scalar from the CSV's seven digits into float32,
    # and stamps the area it found them to integrate over
    rows = [{"year": yr, **{s: f"{N[s][k]:.6e}" for s in N}} for k, yr in enumerate(YEARS)]
    area = {cs.SCALAR_AREA_ATTR: cs.TRUE_AREA if native_af2 else cs.MAP_PLANE} if stamp else None
    for s in N:
        write_series(sub, s, np.array([float(r[s]) for r in rows], dtype="f4"),
                     flux=s not in cs.ST_SCALARS, attrs=area)
    if csv_off:
        rows[1]["lim"] = f"{float(rows[1]['lim']) * 1.001:.6e}"
    native_csv = tmp_path / "run_ismip7_scalars.csv"
    with open(native_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    mm = np.ones((NY, NX), dtype="i4"); mm[4, 1] = 0          # one ice pixel outside the mask
    write_grid(data, "af2", AF2, "f4", flip=False)
    write_grid(data, "maxmask1", mm, "i4", flip=flip)
    write_grid(data, "iaf2", np.ones((NY, NX), dtype="f4"), "f4")
    params = tmp_path / "tree" / "AIS" / "RICE" / "icepack2" / "params.nc"   # in the upload
    write_params(params, rhow=rhow)

    T = tool_scalars(sub, data, params)
    for s in cs.SCALARS:
        vals = np.roll(T[s], 1) if shift_tool and s == "lim" else T[s]
        write_series(tool, s, vals, flux=s not in cs.ST_SCALARS + cs.SEA_LEVEL, time_dtype="f8")

    rows_npz = np.flatnonzero(cov.ravel() > 0)
    counts = np.zeros(NX * NY, dtype=int); counts[rows_npz] = 1
    np.savez(tmp_path / "run.overlap.npz", data=cov.ravel()[rows_npz] * DX * DX,
             indices=np.arange(len(rows_npz)), indptr=np.concatenate([[0], np.cumsum(counts)]),
             shape=np.array([NX * NY, len(rows_npz)]))
    return ["--submission", str(sub), "--tool", str(tool), "--datapath", str(data),
            "--params", str(params), "--refyear", "2016", "--native-csv", str(native_csv),
            "--overlap", str(tmp_path / "run.overlap.npz"),
            "--out-csv", str(tmp_path / "cmp.csv"), "--out-md", str(tmp_path / "cmp.md")], f, cov


def rows_of(tmp_path):
    with open(tmp_path / "cmp.csv") as fh:
        return {(r["scalar"], int(r["year"])): {k: (float(v) if k not in ("scalar", "year") else v)
                                                 for k, v in r.items()} for r in csv.DictReader(fh)}


def gate_line(tmp_path, name):
    return next(ln for ln in (tmp_path / "cmp.md").read_text().splitlines()
                if ln.startswith(f"| {name} |"))


def test_a_consistent_submission_passes_and_every_difference_is_named(tmp_path):
    args, f, cov = build(tmp_path)
    assert cs.main(args) == 0, (tmp_path / "cmp.md").read_text()
    rows = rows_of(tmp_path)
    A = DX * DX
    for (s, yr), r in rows.items():
        parts = r["T_minus_R"] + r["d_area"] + r["d_mm"] + r["d_fill"] + r["d_resid"]
        assert parts == pytest.approx(r["T_minus_N"], rel=1e-9, abs=1e-12 * max(1.0, abs(r["N"])))
    k = 1
    # the area factor on lim, by hand: rho_i sum(h * maxmask1 * (af2 - 1)) dx^2
    mm = np.ones((NY, NX)); mm[4, 1] = 0
    want = 917.0 * np.sum(f["lithk"][k] * mm * (AF2.astype(np.float64) - 1.0)) * A
    assert rows[("lim", 2016)]["d_area"] == pytest.approx(want, rel=1e-6)
    assert rows[("lim", 2016)]["d_mm"] == pytest.approx(-917.0 * f["lithk"][k][4, 1] * A, rel=1e-6)
    # the thin-cell mass lithk leaves out, and the melt booked off the floating ice
    lim = rows[("lim", 2016)]
    assert lim["d_resid"] == pytest.approx(-917.0 * THIN, abs=cs.CSV_DIGITS * abs(lim["N"]))
    assert rows[("tendlibmassbffl", 2016)]["d_resid"] == pytest.approx(-EXTRA_MELT, rel=1e-3)
    # whole-pixel means are the tool's own convention: nothing to undo
    for s in ("tendacabf", "tendlibmassbffl"):
        assert all(rows[(s, yr)]["d_fill"] == 0.0 for yr in YEARS)
    # sea level: zero at the reference, and the native counterpart from lim and limnsw
    assert rows[("slvaf", 2015)]["T"] == 0.0 and rows[("slvaf", 2015)]["N"] == 0.0
    assert rows[("sla20", 2017)]["T"] == pytest.approx(rows[("slg20", 2017)]["T"], abs=1e-9)
    md = (tmp_path / "cmp.md").read_text()
    assert "FAIL" not in md and "Exit status 0." in md
    assert "- flux means: whole pixel" in md and "| tendacabf against N | pass |" in md


def test_a_tree_written_before_the_whole_pixel_means_still_compares(tmp_path):
    r"""No ``flux_pixel_mean`` in the flux files: acabf is a mean over the
    covered part of a pixel and libmassbffl over the floating part, and the
    comparison undoes both."""
    args, f, cov = build(tmp_path, whole_pixel=False)
    assert cs.main(args) == 0, (tmp_path / "cmp.md").read_text()
    rows = rows_of(tmp_path)
    A, k = DX * DX, 1
    # the fill conventions: the half-covered acabf pixel and the half-floating melt pixel
    acabf = np.nan_to_num(f["acabf"][k])
    assert rows[("tendacabf", 2016)]["d_fill"] == pytest.approx(np.sum(acabf * (1 - cov)) * A, rel=1e-6)
    bmb = np.nan_to_num(f["libmassbffl"][k])
    assert rows[("tendlibmassbffl", 2016)]["d_fill"] == pytest.approx(
        np.sum(bmb * (1 - f["sftflf"][k])) * A, rel=1e-6)
    assert rows[("tendlibmassbffl", 2016)]["d_resid"] == pytest.approx(-EXTRA_MELT, rel=1e-3)
    md = (tmp_path / "cmp.md").read_text()
    assert "covered part of a pixel" in md and "| tendacabf against N | pass |" in md


def test_flux_files_that_disagree_on_their_pixel_means_are_refused(tmp_path, capsys):
    args, _, _ = build(tmp_path)
    melt = next((tmp_path / "tree").rglob("libmassbffl_*.nc"))
    with netCDF4.Dataset(melt, "a") as ds:
        ds.delncattr(cs.FLUX_MEAN_ATTR)
    assert cs.main(args) == 2
    assert "disagree" in capsys.readouterr().err


def test_a_tool_series_a_year_out_fails_the_replay(tmp_path):
    args, _, _ = build(tmp_path, shift_tool=True)
    assert cs.main(args) == 1
    assert "| FAIL |" in gate_line(tmp_path, "T against the replay R")


def test_the_tool_run_with_other_densities_is_caught(tmp_path):
    args, _, _ = build(tmp_path, rhow=1027.0)
    assert cs.main(args) == 1
    assert "| FAIL |" in gate_line(tmp_path, "densities")
    assert "| pass |" in gate_line(tmp_path, "T against the replay R")


def test_one_grounded_pixel_off_breaks_the_area_identity(tmp_path):
    args, _, _ = build(tmp_path, perturb_sftgrf=True)
    assert cs.main(args) == 1
    assert "| FAIL |" in gate_line(tmp_path, "forbidden-policy sums against N")


def test_floating_area_written_as_grounded_trades_between_the_areas(tmp_path):
    args, _, _ = build(tmp_path, ground_near=True)
    assert cs.main(args) == 0, (tmp_path / "cmp.md").read_text()
    rows = rows_of(tmp_path)
    moved = rows[("iareagr", 2016)]["d_resid"]
    assert moved == pytest.approx(DX * DX, rel=1e-6)
    assert rows[("iareafl", 2016)]["d_resid"] == pytest.approx(-moved, rel=1e-6)
    assert "near flotation: 1 years" in (tmp_path / "cmp.md").read_text()


def test_a_moving_bed_is_caught(tmp_path):
    args, _, _ = build(tmp_path, moving_bed=True)
    assert cs.main(args) == 1
    assert "| FAIL |" in gate_line(tmp_path, "topg constant in time")


def test_a_csv_that_disagrees_with_the_scalar_files_fails(tmp_path):
    args, _, _ = build(tmp_path, csv_off=True)
    assert cs.main(args) == 1
    assert "| FAIL |" in gate_line(tmp_path, "scalar files against the CSV")


def test_a_grid_on_other_axes_is_refused(tmp_path, capsys):
    args, _, _ = build(tmp_path, flip=True)
    assert cs.main(args) == 2
    assert "array position" in capsys.readouterr().err


def test_the_submission_folder_is_refused_as_tool_output(tmp_path, capsys):
    args, _, _ = build(tmp_path)
    args[args.index("--tool") + 1] = args[args.index("--submission") + 1]
    assert cs.main(args) == 2
    assert "own output tree" in capsys.readouterr().err


def test_a_reference_year_that_is_not_there_exits_2(tmp_path):
    args, _, _ = build(tmp_path)
    args[args.index("--refyear") + 1] = "2030"
    assert cs.main(args) == 2


def test_without_the_overlap_cache_an_old_tree_s_acabf_is_reported_and_not_gated(tmp_path):
    args, _, _ = build(tmp_path, whole_pixel=False)
    i = args.index("--overlap")
    del args[i:i + 2]
    assert cs.main(args) == 0
    assert "| not checked |" in gate_line(tmp_path, "tendacabf against N")
    assert "no --overlap" in (tmp_path / "cmp.md").read_text()


def test_whole_pixel_means_need_no_overlap_cache(tmp_path):
    args, _, _ = build(tmp_path)
    i = args.index("--overlap")
    del args[i:i + 2]
    assert cs.main(args) == 0, (tmp_path / "cmp.md").read_text()
    assert "| pass |" in gate_line(tmp_path, "tendacabf against N")
    assert "no --overlap" not in (tmp_path / "cmp.md").read_text()


def test_native_scalars_over_true_area_leave_no_area_term(tmp_path):
    r"""Issue #97: with the area factor on both sides the area term is zero,
    and every gate still holds."""
    args, _, _ = build(tmp_path, native_af2=True)
    assert cs.main(args + ["--native-af2"]) == 0, (tmp_path / "cmp.md").read_text()
    assert all(r["d_area"] == 0.0 for r in rows_of(tmp_path).values())
    md = (tmp_path / "cmp.md").read_text()
    assert "- native scalars: over true area" in md and "the area term is zero" in md


def test_a_switch_that_contradicts_the_writer_s_stamp_is_refused(tmp_path, capsys):
    args, _, _ = build(tmp_path / "true", native_af2=True)
    assert cs.main(args) == 2
    assert "stamped the native scalars true_area" in capsys.readouterr().err
    args, _, _ = build(tmp_path / "plane")
    assert cs.main(args + ["--native-af2"]) == 2
    assert "stamped the native scalars map_plane" in capsys.readouterr().err


def test_unstamped_true_area_scalars_without_the_switch_fail_the_sums(tmp_path):
    args, _, _ = build(tmp_path, native_af2=True, stamp=False)
    assert cs.main(args) == 1
    assert "| FAIL |" in gate_line(tmp_path, "forbidden-policy sums against N")


def test_a_native_factor_five_percent_off_fails_under_the_switch(tmp_path):
    r"""The switch allows af2's largest change between neighbouring pixels
    (1 % on this small grid, 5.9e-4 on the 8 km one); a missing or doubled
    factor exceeds it."""
    args, _, _ = build(tmp_path, native_af2=True, native_scale=1.05)
    assert cs.main(args + ["--native-af2"]) == 1
    assert "| FAIL |" in gate_line(tmp_path, "forbidden-policy sums against N")


def test_strict_fails_on_the_named_differences(tmp_path):
    args, _, _ = build(tmp_path)
    assert cs.main(args + ["--strict"]) == 1


def test_the_constants_are_the_model_s():
    pytest.importorskip("firedrake")
    from icepack2_tools import ismip7_output
    import write_ismip7_output as wio
    assert cs.RHO_I == ismip7_output.RHO_I
    assert cs.SECONDS_PER_YEAR == ismip7_output.SECONDS_PER_YEAR
    assert (cs.FLUX_MEAN_ATTR, cs.WHOLE_PIXEL) == (wio.FLUX_MEAN_ATTR, wio.FLUX_MEAN)
    assert ((cs.SCALAR_AREA_ATTR, cs.TRUE_AREA, cs.MAP_PLANE)
            == (wio.SCALAR_AREA_ATTR, wio.TRUE_AREA, wio.MAP_PLANE))


def test_it_imports_without_firedrake():
    code = ("import sys; sys.modules['firedrake'] = None; "
            f"sys.path.insert(0, {str(REPO / 'antarctica' / 'scripts')!r}); import compare_scalars")
    subprocess.run([sys.executable, "-c", code], check=True)


def test_the_tool_itself_agrees_with_the_replay(tmp_path):
    r"""The tool finds params.nc in the upload tree, with no --params-path,
    as the organisers' run of it will."""
    scalars = pytest.importorskip("ismip7_scalars.scalars")
    args, _, _ = build(tmp_path)
    out = tmp_path / "real"
    rc = scalars.main(["--region", "AIS", "--group", "RICE", "--model", "icepack2",
                       "--modelid", "m001", "--esm", "CESM2-WACCM", "--forcingid", "f001",
                       "--experiment", "ssp585", "--configid", "C007", "--hist", "ssp585",
                       "--refyear", "2016", "--datapath", args[args.index("--datapath") + 1],
                       "--modelpath", str(tmp_path / "tree" / "AIS"), "--outpath", str(out)])
    assert rc == 0
    args[args.index("--tool") + 1] = str(out / "nc" / "AIS" / "RICE" / "icepack2" / "CORE" / "C007")
    assert cs.main(args) == 0, (tmp_path / "cmp.md").read_text()


def test_front_melt_in_a_pixel_with_no_floating_ice_sums_to_the_model(tmp_path):
    r"""lifmassbf is a whole-pixel mean under ``forbidden`` (issue #109), so
    its grid sum is the model's tendlifmassbf wherever the melt sits, and the
    forbidden-policy gate holds it there; a tree from before, with lifmassbf
    zero, passes the same gate."""
    args, _, _ = build(tmp_path)
    assert cs.main(args) == 0, (tmp_path / "cmp.md").read_text()
    rows = rows_of(tmp_path)
    for yr in YEARS:
        r = rows[("tendlifmassbf", yr)]
        assert r["N"] < 0.0
        assert abs(r["d_resid"]) <= 1e-6 * abs(r["N"])
    assert "tendlifmassbf: the front melt booked" in (tmp_path / "cmp.md").read_text()

    old = tmp_path / "old"
    old.mkdir()
    args, _, _ = build(old, front_melt=False)
    assert cs.main(args) == 0, (old / "cmp.md").read_text()
    assert "tendlifmassbf: zero in every year" in (old / "cmp.md").read_text()


def test_a_model_front_melt_off_the_grid_fails(tmp_path):
    args, _, _ = build(tmp_path, perturb_front=True)
    assert cs.main(args) == 1
    assert "| FAIL |" in gate_line(tmp_path, "forbidden-policy sums against N")
