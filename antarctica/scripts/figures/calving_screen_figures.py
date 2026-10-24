#!/usr/bin/env python3
r"""Figures of a calving screen, in the layout of the ISMIP6 Antarctica
projection paper (Seroussi et al. 2020) and the house style (figstyle.py).

A screen is a set of historical forwards that differ only in the calving
law. The pinned-front member (``fixed``) plays the part the control run
plays in the projections: every "relative" quantity is member minus
pinned front, and sea-level contributions are positive for mass loss.

Inputs, under ``--data``:
  hist_<esm>_<tag>_<member>_<lc>_timeseries.csv   run_simulation's budget rows
  score_<tag>_<member>.json                       score_front.py output (areas)
  fields_<tag>_<member>.npz                       extract_screen_fields.py (maps, regions)
  mesh_<tag>.npz                                  the triangulation (extract_screen_fields.py)

Figures (Seroussi et al. 2020 analogue in brackets):
  fig_budget        SMB, basal melt, calving flux and VAF of every member [Fig. 1]
  fig_extent        the front's rms distance to the observed coast per basin,
                    and the ice area against the Greene et al. (2022) record [Figs. 2, 3]
  fig_contribution  sea-level contribution relative to the pinned front, thin
                    lines per member, thick per law, the spread at the end [Fig. 7]
  fig_regional      the regional split of that contribution [Figs. 5, 8]
  fig_response      mean and spread across the free members of the thickness
                    and speed change relative to the pinned front [Fig. 6]

    python calving_screen_figures.py --data screen_data --tag calvdg1h --out figures/calvdg1h
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.tri import Triangulation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs  # noqa: E402

GT_PER_MM = 361.8          # Gt of ice per mm of sea-level equivalent
RHO_ICE, RHO_WATER = 917.0, 1028.0
# Rignot et al. (2013), Science 341: Antarctic calving flux 1265 +/- 139 Gt/yr
OBS_CALVING = (1265.0, 139.0)

# ── member grammar: one colour + marker per law, line style per parameter ──
LAWS = {  # key prefix -> (label, colour, marker)
    "fixed":  ("pinned front",   fs.BLACK,  "D"),
    "nocalv": ("no calving",     fs.GREY,   "x"),
    "hfb":    ("resistive stress", fs.RED,  "o"),
    "vm":     ("von Mises",       fs.BLUE,  "s"),
    "eig":    ("eigencalving",    fs.GREEN, "^"),
    "thk":    ("thickness",       fs.PURPLE, "v"),
}
LINESTYLES = ["-", "--", ":", "-."]
# Mouginot et al. (2017) basins of the observations kit, grouped as ISMIP6 does
REGIONS = {"Antarctic Peninsula": "AP", "West Antarctica": "WAIS", "East Antarctica": "EAIS"}


def law_of(member):
    for key in LAWS:
        if member.startswith(key):
            return key
    raise KeyError(f"no law grammar for member {member!r}")


def law_params(described):
    r"""The parameters of a law as its ``describe()`` names them, e.g.
    ``hfb(sigma_max=0.086, ...)`` -> {"sigma_max": 0.086, ...}."""
    m = re.match(r"\w+\((.*?)\)", described)
    if not m:
        raise ValueError(f"not a law description: {described!r}")
    params = {}
    for item in filter(None, (p.strip() for p in m.group(1).split(","))):
        k, v = item.split("=", 1)
        try:
            params[k] = float(v)
        except ValueError:
            params[k] = v
    return params


def member_label(member, law, described=None):
    r"""``hfb086`` run as ``hfb(sigma_max=0.086, ...)`` -> "resistive stress,
    0.086 MPa". The parameter comes from the law the checkpoint recorded
    (``described``), because the member names use more than one scale; a
    member without that record is labelled by its name."""
    base = LAWS[law][0]
    if law in ("fixed", "nocalv"):
        return base
    if described is None:
        return f"{base} ({member})"
    p = law_params(described)
    if law == "hfb":
        return f"{base}, {p['sigma_max']:g} MPa"
    if law == "vm":
        grounded = p["sigma_max_gr"] == p["sigma_max_fl"]
        return f"{base}, {p['sigma_max_fl']:g} MPa" + (" (grounded too)" if grounded else "")
    if law == "eig":
        return f"{base}, K = " + f"{p['K']:.3g}".replace("e+0", "e").replace("e+", "e")
    if law == "thk":
        return f"{base}, {p['hc']:g} m"
    return f"{base} ({member})"


class Member:
    def __init__(self, name, csv_path, score_path=None, fields_path=None):
        self.name, self.law = name, law_of(name)
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        self.ts = {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if rows[0][k] not in ("",)}
        y = self.ts["year"]
        self.t0 = round(y[0] - (y[1] - y[0]), 3) if len(y) > 1 else y[0]
        self.dt = np.diff(np.concatenate([[self.t0], y]))
        self.score = json.load(open(score_path)) if score_path and os.path.exists(score_path) else None
        self.fields = np.load(fields_path) if fields_path and os.path.exists(fields_path) else None

    def described(self):
        r"""The law the member's checkpoint recorded, or ``None`` without one."""
        if self.fields is None or "calving_law" not in self.fields:
            return None
        return str(self.fields["calving_law"]) or None

    # time series helpers -------------------------------------------------
    def years(self):
        return self.ts["year"]

    def calving_rate(self):
        r"""Calving flux [Gt/yr] per step (the CSV books Gt per step)."""
        return self.ts["calv_gt"] / self.dt

    def yearly(self, v):
        r"""Block means over model years (year t covers (t, t+1]); returns (year, mean)."""
        yr = np.floor(self.ts["year"] - 1e-6).astype(int)
        out = [(y, np.average(v[yr == y], weights=self.dt[yr == y])) for y in np.unique(yr)]
        return np.array([o[0] + 0.5 for o in out]), np.array([o[1] for o in out])

    def dvaf_mm(self):
        r"""VAF change since the start [mm SLE], positive for gain."""
        return self.ts["vaf_mm_sle"] - self.ts["vaf_mm_sle"][0] + self.dt[0] * 0.0

    def style(self, k=0):
        label, colour, marker = LAWS[self.law]
        return dict(color=colour, marker=marker, ls=LINESTYLES[k % len(LINESTYLES)],
                    markevery=max(len(self.ts["year"]) // 12, 1), lw=fs.LW, ms=fs.MS,
                    markerfacecolor="white" if k else colour, markeredgewidth=0.7, label=member_label(self.name, self.law, self.described()))


def load_members(data, tag):
    pat = re.compile(rf"hist_(\w+?)_{re.escape(tag)}_(\w+?)_(\d+)_timeseries\.csv$")
    members = {}
    for p in sorted(glob.glob(os.path.join(data, f"hist_*_{tag}_*_timeseries.csv"))):
        m = pat.search(os.path.basename(p))
        if not m:
            continue
        name = m.group(2)
        members[name] = Member(name, p, os.path.join(data, f"score_{tag}_{name}.json"),
                               os.path.join(data, f"fields_{tag}_{name}.npz"))
    if "fixed" not in members:
        sys.exit(f"no pinned-front member ({tag}_fixed) under {data}")
    # order: pinned, no calving, then the laws in grammar order
    order = list(LAWS)
    return dict(sorted(members.items(), key=lambda kv: (order.index(kv[1].law), kv[0])))


def styled(members):
    r"""Assign a line style per member within its law, and return [(member, style)]."""
    count, out = {}, []
    for m in members.values():
        k = count.get(m.law, 0)
        count[m.law] = k + 1
        out.append((m, m.style(k)))
    return out


def relative(member, fixed, key):
    r"""Member minus pinned front on the member's time axis."""
    return np.interp(member.years(), fixed.years(), fixed.ts[key])


def contribution_mm(member, fixed):
    r"""Sea-level contribution relative to the pinned front [mm SLE], positive = loss."""
    dv = member.dvaf_mm() - (relative(member, fixed, "vaf_mm_sle") - fixed.ts["vaf_mm_sle"][0])
    return -dv


def headroom(ax, frac=0.18):
    r"""Room above the data for the panel letter (the letters always sit top left)."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + frac * (hi - lo))


def legend_across_top(fig, handles, y=1.0):
    if not handles:
        return None
    return fs.top_legend(fig, handles, ncol=min(len(handles), 4), y=y)


# ── Fig. 1 analogue: the budget of every member ─────────────────────────────
def fig_budget(members, out):
    fig, axs = plt.subplots(2, 2, figsize=(fs.FULL, 0.72 * fs.FULL), sharex=True)
    axs = axs.ravel()
    handles = []
    for m, st in styled(members):
        y = m.years()
        axs[0].plot(y, m.ts["smb_gtyr"], **st)
        axs[1].plot(y, -np.abs(m.ts["melt_gtyr"]), **st)
        yc, c = m.yearly(m.calving_rate())
        axs[2].plot(yc, -c, **{**st, "markevery": 1})
        axs[3].plot(y, m.dvaf_mm(), **st)
        handles.append(Line2D([], [], **{k: v for k, v in st.items() if k != "markevery"}))
    lo, hi = -(OBS_CALVING[0] + OBS_CALVING[1]), -(OBS_CALVING[0] - OBS_CALVING[1])
    axs[2].axhspan(lo, hi, color="0.85", zorder=0, lw=0)
    axs[2].text(0.98, 0.04, "Rignot et al. (2013)", transform=axs[2].transAxes, ha="right", va="bottom",
                fontsize=fs.FS_TICK, color="0.35")
    for ax, letter, lab in zip(axs, "abcd", ["surface mass balance (Gt yr$^{-1}$)", "basal melt (Gt yr$^{-1}$)",
                                            "calving flux (Gt yr$^{-1}$)", "VAF change (mm SLE)"]):
        fs.grid(ax)
        fs.panel_letter(ax, letter)
        ax.set_ylabel(lab)
        ax.margins(x=0.02)
        headroom(ax)
    for ax in axs[2:]:
        ax.set_xlabel("year")
    fig.align_ylabels(axs[[0, 2]])
    fig.align_ylabels(axs[[1, 3]])
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.09, top=0.88, wspace=0.28, hspace=0.1)
    legend_across_top(fig, handles, y=0.9)
    fs.save(fig, out / "fig_budget")


# ── Fig. 7 analogue: contribution relative to the pinned front ──────────────
def fig_contribution(members, out):
    fixed = members["fixed"]
    fig = plt.figure(figsize=(fs.FULL, 0.5 * fs.FULL))
    gs = fig.add_gridspec(1, 2, width_ratios=[5, 1], wspace=0.05, left=0.09, right=0.99, bottom=0.14, top=0.86)
    ax, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    handles, finals = [], {}
    by_law = {}
    for m, st in styled(members):
        if m.name == "fixed":
            continue
        y, c = m.years(), contribution_mm(m, fixed)
        thin = {**st, "lw": 0.6 * fs.LW, "ms": 0.8 * fs.MS, "alpha": 0.9}
        ax.plot(y, c, **thin)
        handles.append(Line2D([], [], **{k: v for k, v in thin.items() if k != "markevery"}))
        if y[-1] >= fixed.years()[-1] - 0.5:   # the spread at the end is for members that reach it
            finals[m.name] = (c[-1], st["color"], st["marker"], st["markerfacecolor"])
        by_law.setdefault(m.law, []).append((y, c))
    # thick per-law means where a law has several members
    for law, runs in by_law.items():
        if len(runs) < 2:
            continue
        end = min(y[-1] for y, _ in runs)   # the mean stops where the first member of the law stops
        common = np.arange(fixed.years()[0], end + 1e-9, 0.25)
        series = [np.interp(common, y, c) for y, c in runs]
        ax.plot(common, np.mean(series, axis=0), color=LAWS[law][1], lw=2.2 * fs.LW, zorder=4, solid_capstyle="round")
    ax.axhline(0, color=fs.BLACK, lw=0.6)
    fs.grid(ax)
    fs.panel_letter(ax, "a")
    ax.set_xlabel("year")
    ax.set_ylabel("sea-level contribution\nrelative to the pinned front (mm SLE)")
    ax.margins(x=0.02)
    # the spread at the end, one marker per member on a shared axis
    for name, (v, colour, marker, face) in finals.items():
        axb.plot([0], [v], marker=marker, color=colour, markerfacecolor=face, markeredgewidth=0.7, ms=fs.MS + 1,
                 ls="none", zorder=3)
    headroom(ax)
    axb.set_ylim(ax.get_ylim())
    axb.set_xticks([0])
    axb.set_xticklabels([f"{int(fixed.years()[-1])}"])
    axb.yaxis.tick_right()
    axb.tick_params(labelsize=fs.FS_TICK)
    axb.axhline(0, color=fs.BLACK, lw=0.6)
    fs.grid(axb)
    fs.panel_letter(axb, "b", x=0.3)
    legend_across_top(fig, handles, y=0.88)
    fs.save(fig, out / "fig_contribution")


# ── Figs. 2/3 analogue: the front against the record ──────────────────────
BASINS = ["A-A'", "A'-B", "B-C", "C-C'", "C'-D", "D-D'", "D'-E", "E-E'", "E'-F'", "F'-G", "G-H", "H-H'",
          "H'-I", "I-I''", "I''-J", "J-J''", "J''-K", "K-A"]   # the kit's Mouginot basins 1..18


def scored_states(m):
    return [s for s in m.score["states"] if s.get("epoch") is not None] if m.score else []


def fig_extent(members, mesh, out):
    r"""(a) the rms distance from each member's front to the observed coast per basin at the
    year of the field files (or the last scored state), (b) the ice area against the record.
    At 2 km the members' extents differ by a cell at the fronts, so a map of where they
    disagree (the count raster in the mesh file) reads as blank at continental scale; the
    per-basin misfit is the figure that carries the information."""
    scored = [(m, st) for m, st in styled(members) if scored_states(m)]
    if not scored:
        print("fig_extent: no front scores, skipped")
        return
    year = float(mesh["t_yr"]) if mesh is not None and "t_yr" in mesh else None
    fig = plt.figure(figsize=(fs.FULL, 0.5 * fs.FULL))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.22, left=0.07, right=0.99, bottom=0.17, top=0.82)
    ax, axa = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    handles, rec, years_used = [], None, set()
    for m, st in scored:
        st_ = {k: v for k, v in st.items() if k != "markevery"}
        states = scored_states(m)
        pick = min(states, key=lambda s: abs(s["t_yr"] - year)) if year is not None else states[-1]
        years_used.add(round(pick["t_yr"]))
        rms = np.array([pick["basins"].get(str(b), {}).get("front", {}).get("rms_m", np.nan) for b in range(1, 19)]) / 1e3
        ax.plot(np.arange(1, 19), rms, **{**st_, "ls": "none", "ms": fs.MS + 1})
        yrs = np.array([s["t_yr"] for s in states])
        area = np.array([s["ice_area_km2"] for s in states])
        if rec is None:
            rec = (yrs, np.array([s["obs_ice_area_km2"] for s in states]))
        axa.plot(yrs, (area - area[0]) / 1e3, **st_)
        handles.append(Line2D([], [], **st_))
    # the whole-front rms of each member as a thin line at the right edge
    ax.set_xticks(np.arange(1, 19))
    ax.set_xticklabels(BASINS, rotation=60, ha="right", fontsize=fs.FS_TICK - 0.5)
    ax.set_ylabel("front rms distance to the observed coast (km)")
    ax.set_ylim(bottom=0)
    fs.grid(ax)
    fs.panel_letter(ax, "a")
    headroom(ax)
    fs.title(ax, ", ".join(f"{y}" for y in sorted(years_used)) + " against the Greene et al. (2022) coast")
    axa.plot(rec[0], (rec[1] - rec[1][0]) / 1e3, color=fs.BLACK, lw=1.8 * fs.LW, marker="*", ms=fs.MS + 2, zorder=5)
    handles.append(Line2D([], [], color=fs.BLACK, lw=1.8 * fs.LW, marker="*", ms=fs.MS + 2, label="Greene et al. (2022)"))
    fs.grid(axa)
    fs.panel_letter(axa, "b")
    axa.set_xlabel("year")
    axa.set_ylabel("ice area change (10$^3$ km$^2$)")
    axa.margins(x=0.03)
    headroom(axa)
    legend_across_top(fig, handles, y=0.9)
    fs.save(fig, out / "fig_extent")


# ── Figs. 5/8 analogue: regional split ───────────────────────────────────────
def fig_regional(members, out):
    fixed = members["fixed"]
    if fixed.fields is None or "region_dvaf_gt" not in fixed.fields:
        print("fig_regional: no regional totals in the field files, skipped")
        return
    regions = [str(r) for r in fixed.fields["region_names"]]
    rows = [(m, st) for m, st in styled(members) if m.name != "fixed" and m.fields is not None]
    nb, nr = len(rows), len(regions)
    fig, ax = plt.subplots(figsize=(fs.FULL, 0.42 * fs.FULL))
    width = 0.8 / nb
    handles = []
    for i, (m, st) in enumerate(rows):
        dv = -(m.fields["region_dvaf_gt"] - fixed.fields["region_dvaf_gt"]) / GT_PER_MM
        xs = np.arange(nr) - 0.4 + (i + 0.5) * width
        ax.bar(xs, dv, width=width * 0.92, color=st["color"], edgecolor=fs.BLACK, lw=0.4,
               hatch=None if st["ls"] == "-" else "///", zorder=3)
        if "region_smb_grounded_gt" in m.fields:
            smb = -(m.fields["region_smb_grounded_gt"] - fixed.fields["region_smb_grounded_gt"]) / GT_PER_MM
            ax.plot(xs, smb, ls="none", marker="D", color=fs.BLACK, markerfacecolor="white", ms=fs.MS, zorder=4)
        handles.append(Line2D([], [], color=st["color"], marker="s", ls="none", ms=fs.MS + 2, label=st["label"]))
    if "region_smb_grounded_gt" in fixed.fields:
        handles.append(Line2D([], [], color=fs.BLACK, marker="D", markerfacecolor="white", ls="none", ms=fs.MS,
                              label="grounded SMB anomaly"))
    ax.axhline(0, color=fs.BLACK, lw=0.6)
    ax.set_xticks(np.arange(nr))
    ax.set_xticklabels(regions)
    ax.set_ylabel(f"contribution {fixed.t0:.0f}-{float(fixed.fields['t_yr']):.0f} (mm SLE)")
    fs.grid(ax)
    ax.xaxis.grid(False)
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.82)
    legend_across_top(fig, handles, y=0.84)
    fs.save(fig, out / "fig_regional")


# ── Fig. 6 analogue: mean and spread of the response ────────────────────────
def cell_to_vertex(mesh, v):
    r"""Area-weighted average of a cell field onto the vertices (for filled contours)."""
    tri, area = mesh["cells"], mesh["area"]
    num = np.bincount(tri.ravel(), weights=np.repeat(v * area, 3), minlength=len(mesh["x"]))
    den = np.bincount(tri.ravel(), weights=np.repeat(area, 3), minlength=len(mesh["x"]))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def fig_response(members, mesh, out):
    fixed = members["fixed"]
    free = [m for m in members.values() if m.fields is not None and m.law not in ("fixed", "nocalv")]
    if fixed.fields is None or not free or mesh is None:
        print("fig_response: field files missing, skipped")
        return
    tri = Triangulation(mesh["x"], mesh["y"], mesh["cells"])
    ice = (fixed.fields["h"] > 1.0) | (fixed.fields["h0"] > 1.0)
    # the window: the ice at the start, padded by 2% of its width
    cx, cy = mesh["cx"][ice], mesh["cy"][ice]
    pad = 0.02 * (cx.max() - cx.min())
    box = (cx.min() - pad, cx.max() + pad, cy.min() - pad, cy.max() + pad)
    dh = np.array([(m.fields["h"] - m.fields["h0"]) - (fixed.fields["h"] - fixed.fields["h0"]) for m in free])
    du = np.array([m.fields["speed"] - fixed.fields["speed"] for m in free])
    panels = [("thickness change (m)", dh.mean(0), dh.std(0)), ("speed change (m yr$^{-1}$)", du.mean(0), du.std(0))]
    aspect = (box[3] - box[2]) / (box[1] - box[0])
    fig, axs, caxs = fs.map_grid(2, [1, 1], aspect, left=0.35, right=0.1)
    letters = iter("abcd")
    for r, (label, mean, sd) in enumerate(panels):
        for c, (v, kind, name) in enumerate([(mean, "difference", "mean"), (sd, "spread", "standard deviation")]):
            vv = cell_to_vertex(mesh, np.where(ice, v, np.nan))
            ax = axs[r][c]
            tri.set_mask(~ice)
            m = fs.tri_field(ax, tri, np.nan_to_num(vv), kind, Normalize(*fs.limits(v[ice], kind, percentile=2.0)))
            fs.frame(ax, box, scale=(r == 1 and c == 0))
            fs.colorbar(fig, m, caxs[r][c], f"{name}, {label}")
            fs.panel_letter(ax, next(letters))
            if r == 0:
                fs.title(ax, name)
        fs.row_label(axs[r][0], label.split(" (")[0])
    n = len(free)
    fig.suptitle(f"{n} free member{'s' if n != 1 else ''} minus the pinned front, {fixed.t0:.0f}-{float(fixed.fields['t_yr']):.0f}",
                 y=0.995, fontsize=fs.FS)
    fs.save(fig, out / "fig_response")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default=None, help="comma-separated figure names")
    args = ap.parse_args()
    fs.apply_rc()
    out = Path(args.out)
    members = load_members(args.data, args.tag)
    print("members:", ", ".join(members))
    mesh_path = os.path.join(args.data, f"mesh_{args.tag}.npz")
    mesh = np.load(mesh_path) if os.path.exists(mesh_path) else None
    todo = args.only.split(",") if args.only else ["budget", "contribution", "extent", "regional", "response"]
    for name in todo:
        {"budget": lambda: fig_budget(members, out),
         "contribution": lambda: fig_contribution(members, out),
         "extent": lambda: fig_extent(members, mesh, out),
         "regional": lambda: fig_regional(members, out),
         "response": lambda: fig_response(members, mesh, out)}[name]()


if __name__ == "__main__":
    main()
