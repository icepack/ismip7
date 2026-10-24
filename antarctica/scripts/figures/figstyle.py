r"""House figure style for this repository's figures.

The group's figures emulate Morrow, Ghelichkhan & Kramer (2026), whose
figure code is github.com/sghelichkhani/morrow2026 (@ 832e49f). This module
restates that grammar so every figure built here follows it; change a value
here, not in a figure. Sources, by file in that repository:

``parallel_scaling/plot_paper_figures.py``, ``reported.py``
    Times serif text with STIX math; series in a fixed colour + marker (red
    #d62728 circles, blue #3182bd squares, black diamonds, green #2ca02c
    triangles, purple #9467bd), a marker on every data point, line width
    1.9 / marker size 8.5 at their drawing size; ``grid(which="both",
    alpha=0.3)``; panel letters in light-grey circles at the top left of each
    axes; short headings in light-yellow rounded boxes; one framed legend
    across the top of the figure, shared by every panel.
``verification/tracy/plot_solution.py``, ``verification/vauclin/plot_solution.py``
    field panels attached with no gaps and no ticks; 10 discrete colour bands
    in panels and colour bars alike; a rotated, yellow-boxed row label on the
    left; one colour bar per row on the right; time stamps in light-blue boxes.
``parallel_scaling/figstyle.py``
    ``save``: figures written without creation timestamps, so regenerating an
    unchanged figure gives an empty git diff.

Morrow et al. draw at 9-12 inches with 15-17 pt text and scale the figures
down on the page. Here they are drawn at print size (Copernicus full width
17.8 cm), with every length scaled by the same ~0.58, so the PDF needs no
rescaling. Map panels, which Morrow et al. leave bare, carry a scale bar
(polar stereographic, so no north arrow). Panel letters are lower case, as
in their mass-conservation figure and the Copernicus caption convention.
Findings go in captions, not titles.

    import figstyle as fs
    fs.apply_rc()
    fig, axs, caxs = fs.map_grid(2, [2, 1], aspect)
    m = fs.tri_field(axs[0][0], tri, v, "difference", norm)
    fs.colorbar(fig, m, caxs[0][0], label)
    fs.panel_letter(axs[0][0], "a"); fs.row_label(axs[0][0], "thickness change")
    fs.save(fig, path)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import LogFormatterMathtext, MaxNLocator

# ── series grammar (reported.py) ───────────────────────────────────────────
RED, BLUE, BLACK, GREEN, PURPLE = "#d62728", "#3182bd", "#000000", "#2ca02c", "#9467bd"
GREY = "gray"
LW, MS = 1.1, 4.5                 # their 1.9 / 8.5, scaled to print size

# ── boxes (plot_paper_figures.py, tracy/plot_solution.py) ──────────────────
LETTER_BOX = dict(boxstyle="circle,pad=0.25", facecolor="lightgrey", edgecolor="black", lw=0.6)
TITLE_BOX = dict(boxstyle="round,pad=0.3", facecolor="lightyellow", edgecolor="black", lw=0.6)
TIME_BOX = dict(boxstyle="round,pad=0.3", facecolor="lightblue", edgecolor="black", lw=0.6)

# ── page and type ──────────────────────────────────────────────────────────
FULL = 17.8 / 2.54                 # inches
DPI = 300
FS, FS_LABEL, FS_TICK = 8.5, 9.0, 7.5
N_LEVELS = 10                      # discrete colour bands

# ── field styles: colormap, symmetric about zero ───────────────────────────
STYLES = {
    "theta":      ("RdBu_r", True),
    "difference": ("RdBu_r", True),
    "residual":   ("RdBu_r", True),
    "gradient":   ("RdBu_r", True),
    "speed":      ("magma",  False),
    "ratio":      ("Blues_r", False),
    "elevation":  ("cividis", False),
    "thickness":  ("Blues",   False),
    "spread":     ("Purples", False),
    "count":      ("viridis", False),
}
MINUS = "\N{MINUS SIGN}"
HALO_LIGHT = [pe.Stroke(linewidth=1.9, foreground="white"), pe.Normal()]
HALO_DARK = [pe.Stroke(linewidth=1.9, foreground=BLACK), pe.Normal()]


def apply_rc():
    """Morrow et al. rcParams (Times, STIX math), at print size; call once before plotting."""
    plt.rcParams.update({
        "font.family": "serif",
        # STIX General (Times-based, TrueType, bundled with matplotlib) stands in where
        # Times New Roman is not installed; the OpenType Nimbus Roman embeds as a CFF
        # font that PDF readers flag as mismatched.
        "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS_LABEL,
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK, "legend.fontsize": FS,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "lines.linewidth": LW, "lines.markersize": MS, "axes.unicode_minus": True,
        "legend.edgecolor": "black", "legend.fancybox": False, "legend.framealpha": 1.0,
        "patch.linewidth": 0.6,
        "figure.dpi": 100, "savefig.dpi": DPI, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def grid(ax):
    """plot_paper_figures.py: major and minor grid, alpha 0.3, under the data."""
    ax.grid(True, which="both", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)


def panel_letter(ax, letter, x=0.075, y=0.91):
    """Grey-circle panel letter at the top left (plot_paper_figures.panel_letter)."""
    ax.text(x, y, letter, transform=ax.transAxes, ha="center", va="center", fontsize=FS,
            zorder=6, bbox=LETTER_BOX)


def title(ax, text, pad=6):
    """A short heading in a yellow box above the axes."""
    ax.set_title(text, bbox=TITLE_BOX, pad=pad, fontsize=FS)


def row_label(ax, text, x=-0.06):
    """Rotated, yellow-boxed row label outside the left edge (tracy/plot_solution.py)."""
    ax.text(x, 0.5, text, transform=ax.transAxes, rotation=90, ha="right", va="center",
            fontsize=FS, bbox=TITLE_BOX)


def top_legend(fig, handles, ncol=None, y=1.0):
    """One framed legend across the top of the figure, shared by every panel."""
    return fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, y),
                      ncol=ncol or len(handles), frameon=True, handlelength=1.8,
                      columnspacing=1.4, borderpad=0.4)


def num(x, spec):
    """Format a number for figure text with a true minus sign, as on the tick labels."""
    return format(x, spec).replace("-", MINUS)


def sci(x, digits=0):
    r"""Mathtext scientific notation, e.g. 5e6 -> $5 \times 10^{6}$."""
    e = int(np.floor(np.log10(abs(x))))
    return rf"${x / 10 ** e:.{digits}f} \times 10^{{{e}}}$"


# ── colour limits ──────────────────────────────────────────────────────────
def nice(x):
    """Round up to 1, 2, 2.5 or 5 times a power of ten, so colour-bar ends are round."""
    if x <= 0:
        return 1.0
    p = 10 ** np.floor(np.log10(x))
    return float(min(m * p for m in (1, 2, 2.5, 5, 10) if m * p >= x * (1 - 1e-9)))


def limits(a, kind, percentile=1.0):
    """Robust colour limits: percentile-clipped, symmetric where the sign matters, rounded."""
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if STYLES[kind][1]:
        m = nice(np.percentile(np.abs(a), 100 - percentile))
        return -m, m
    return 0.0, nice(np.percentile(a, 100 - percentile))


# ── attached map panels (tracy/plot_solution.py layout) ────────────────────
def map_grid(nrows, groups, aspect, width=FULL, left=0.30, right=0.55, top=0.30, bottom=0.06,
             cbar_w=0.10, cbar_pad=0.05, group_gap=0.62, hgap=0.04):
    r"""Rows of attached map panels, ``groups`` panels per colour bar (e.g. [2, 1]: two
    panels sharing one colour bar, then a third with its own), sized so each panel shows
    the map at ``aspect`` = height / width with no distortion. Lengths are in inches.

    Returns the figure, ``axs[r][c]`` and ``caxs[r][g]`` (one colour-bar axes per group)."""
    n = sum(groups)
    fixed = left + right + len(groups) * (cbar_pad + cbar_w) + (len(groups) - 1) * group_gap
    w = (width - fixed) / n
    h = w * aspect
    height = top + nrows * h + (nrows - 1) * hgap + bottom
    fig = plt.figure(figsize=(width, height))
    axs, caxs = [], []
    for r in range(nrows):
        y0 = height - top - (r + 1) * h - r * hgap
        x, row, crow = left, [], []
        for g in groups:
            for _ in range(g):
                ax = fig.add_axes([x / width, y0 / height, w / width, h / height])
                ax.set_xticks([])
                ax.set_yticks([])
                row.append(ax)
                x += w
            x += cbar_pad
            crow.append(fig.add_axes([x / width, (y0 + 0.06 * h) / height, cbar_w / width,
                                      0.88 * h / height]))
            x += cbar_w + group_gap
        axs.append(row)
        caxs.append(crow)
    return fig, axs, caxs


def frame(ax, box, scale=False):
    """Fix a map panel to ``box`` = (x0, x1, y0, y1) in metres, bare as in Morrow et al."""
    ax.set(xlim=box[:2], ylim=box[2:])
    ax.set_xticks([])
    ax.set_yticks([])
    if scale:
        scalebar(ax)


def scalebar(ax, loc=(0.07, 0.06), color=BLACK):
    """Kilometre scale bar sized to a round fifth of the panel width."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    raw = abs(x1 - x0) / 5e3
    mag = 10 ** np.floor(np.log10(raw))
    L = float(min([1, 2, 5, 10], key=lambda m: abs(m * mag - raw)) * mag)
    bx, by = x0 + loc[0] * (x1 - x0), y0 + loc[1] * (y1 - y0)
    ax.plot([bx, bx + L * 1e3], [by, by], color=color, lw=2.0, solid_capstyle="butt",
            zorder=6, path_effects=HALO_LIGHT)
    ax.annotate(f"{L:g} km", (bx + L * 5e2, by), textcoords="offset points", xytext=(0, 3),
                ha="center", fontsize=FS_TICK, color=color, zorder=6,
                path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])


def tri_field(ax, tri, v, kind, norm=None):
    r"""Field in ``N_LEVELS`` discrete bands on ``tri``, rasterized under vector overlays.
    ``norm`` defaults to :func:`limits`; a LogNorm gives log-spaced bands."""
    cmap = STYLES[kind][0]
    if norm is None:
        norm = Normalize(*limits(v, kind))
    log = isinstance(norm, LogNorm)
    if log:
        lv = np.geomspace(norm.vmin, norm.vmax, N_LEVELS + 1)
        v = np.maximum(v, 0.5 * norm.vmin)
    else:
        lv = np.linspace(norm.vmin, norm.vmax, N_LEVELS + 1)
    lo, hi = np.nanmin(v) < lv[0], np.nanmax(v) > lv[-1]
    extend = {(1, 1): "both", (1, 0): "min", (0, 1): "max"}.get((int(lo), int(hi)), "neither")
    m = ax.tricontourf(tri, v, lv, cmap=cmap, norm=norm if log else None, extend=extend, zorder=0)
    ax.set_rasterization_zorder(0.5)
    return m


def colorbar(fig, mappable, cax, label):
    """Discrete-band colour bar with its label rotated alongside (tracy/plot_solution.py)."""
    cb = fig.colorbar(mappable, cax=cax)
    cb.set_label(label, rotation=90, labelpad=4)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(width=0.5, length=2.5, labelsize=FS_TICK)
    lv = np.asarray(mappable.levels)
    if np.all(lv > 0) and lv[-1] / lv[0] > 50:
        cb.set_ticks([10.0 ** k for k in range(int(np.ceil(np.log10(lv[0]))),
                                               int(np.floor(np.log10(lv[-1]))) + 1)])
        cb.ax.yaxis.set_major_formatter(LogFormatterMathtext())
    else:
        t = MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]).tick_values(lv[0], lv[-1])
        eps = 1e-9 * (lv[-1] - lv[0])
        cb.set_ticks(t[(t >= lv[0] - eps) & (t <= lv[-1] + eps)])
        cb.ax.yaxis.set_major_formatter(lambda v, _: num(float(f"{v:.3g}"), "g"))
    return cb


def grounding_line(ax, tri, haf, color=BLACK, lw=0.8, ls="-"):
    """Grounding line = zero contour of height above flotation, haloed so it reads on any fill:
    a black line gets a white halo, a coloured one a black halo."""
    cs = ax.tricontour(tri, haf, levels=[0.0], colors=color, linewidths=lw, linestyles=ls)
    cs.set(path_effects=HALO_LIGHT if color == BLACK else HALO_DARK)
    return cs


def verification(ax, items):
    """A checklist panel: a check, a claim, and the evidence under it in grey."""
    ax.axis("off")
    yy, step = 0.80, 0.80 / max(len(items), 1)
    for head, sub in items:
        kw = dict(transform=ax.transAxes, va="center", in_layout=False)
        ax.text(0.02, yy, "\N{CHECK MARK}", color=GREEN, fontsize=FS + 1, fontfamily="DejaVu Sans",
                **kw)   # Times has no check mark
        ax.text(0.10, yy, head, fontsize=FS, color=BLACK, **kw)
        ax.text(0.10, yy - 0.42 * step, sub, fontsize=FS_TICK, color="0.35", **kw)
        yy -= step


# ── output (parallel_scaling/figstyle.py) ──────────────────────────────────
_PDF_METADATA = {"CreationDate": None, "Producer": None, "Creator": None}
_PNG_METADATA = {"Software": ""}


def save(fig, path, pad_inches=0.04):
    """PNG at 300 dpi and a vector PDF (fills rasterized), with no embedded timestamp."""
    path.parent.mkdir(parents=True, exist_ok=True)
    for ext, meta in (("png", _PNG_METADATA), ("pdf", _PDF_METADATA)):
        fig.savefig(str(path.with_suffix(f".{ext}")), bbox_inches="tight",
                    pad_inches=pad_inches, metadata=dict(meta))
    plt.close(fig)
    print(f"Saved {path.with_suffix('')}.{{png,pdf}}", flush=True)
