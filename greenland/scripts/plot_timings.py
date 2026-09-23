#! /usr/bin/env python3
"""Plot diagnostic solve timings from test_logs/timings_*.txt files."""

import argparse
import csv
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median


DEFAULT_LOG_DIR = Path("test_logs")
DEFAULT_FIG_DIR = Path(__file__).resolve().parent.parent / "figs" / "timings"
MESH_COLORS = {
    "detailed": "tab:blue",
    "promice": "tab:orange",
    "simple": "tab:green",
    "buffered": "tab:red",
}
TIMING_RE = re.compile(r"^timings_(?P<mesh_kind>.+)_n(?P<nproc>\d+)\.txt$")


@dataclass(frozen=True)
class TimingRow:
    """One timing measurement from a diagnostic solve log."""

    mesh_kind: str
    nproc: int
    lc: int
    seconds: float
    damping: float
    delta0: float

def configure_matplotlib_cache() -> None:
    """Point Matplotlib and fontconfig at writable cache directories."""
    cache_root = Path(os.environ.get("TMPDIR", "/tmp")) / "greenland_plot_cache"
    matplotlib_cache = cache_root / "matplotlib"
    font_cache = cache_root / "fontconfig"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    font_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(font_cache))
    os.environ.setdefault("MPLBACKEND", "Agg")



def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot diagnostic solve timings from test_logs/timings_*.txt."
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Directory containing timings_<mesh_kind>_n<nproc>.txt files.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=DEFAULT_FIG_DIR,
        help="Directory where timing figures will be saved.",
    )
    parser.add_argument(
        "--output-name",
        default="diagnostic_timings.png",
        help="Output figure filename.",
    )
    return parser.parse_args()


def parse_float(value: str) -> float:
    """Parse a float, including nan values from timing logs."""
    value = value.strip()
    if value.lower() == "nan":
        return math.nan
    return float(value)


def read_timing_file(path: Path) -> list[TimingRow]:
    """Read one diagnostic timing file."""
    match = TIMING_RE.match(path.name)
    if match is None:
        return []

    mesh_kind = match.group("mesh_kind")
    nproc = int(match.group("nproc"))
    rows = []
    with path.open(newline="") as fin:
        reader = csv.DictReader(fin, skipinitialspace=True)
        for row in reader:
            if not row.get("time (s)"):
                continue
            rows.append(
                TimingRow(
                    mesh_kind=mesh_kind,
                    nproc=nproc,
                    lc=int(row["lc (m)"]),
                    seconds=float(row["time (s)"]),
                    damping=parse_float(row["LS damping"]),
                    delta0=parse_float(row["TR delta0"]),
                )
            )
    return rows


def read_timings(log_dir: Path) -> list[TimingRow]:
    """Read all diagnostic timing files from a log directory."""
    rows = []
    for path in sorted(log_dir.glob("timings_*_n*.txt")):
        rows.extend(read_timing_file(path))
    return rows


def aggregate_timings(rows: list[TimingRow]) -> dict[tuple[str, int], list[tuple[int, float]]]:
    """Aggregate repeated timing rows by median runtime."""
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.mesh_kind, row.nproc, row.lc)].append(row.seconds)

    series = defaultdict(list)
    for (mesh_kind, nproc, lc), seconds in grouped.items():
        series[(mesh_kind, nproc)].append((lc, median(seconds)))

    return {key: sorted(values) for key, values in series.items()}


def plot_timings(rows: list[TimingRow], fig_path: Path) -> None:
    """Plot diagnostic solve runtimes by resolution and mesh kind."""
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    series = aggregate_timings(rows)
    if not series:
        raise ValueError("No timing rows found to plot")

    fig, ax = plt.subplots(figsize=(9, 6))
    nprocs = sorted({nproc for _, nproc in series})
    linestyles = ["-", "--", ":", "-."]
    markers = ["o", "s", "^", "D", "v", "P"]

    for mesh_kind in MESH_COLORS:
        for nproc in nprocs:
            values = series.get((mesh_kind, nproc))
            if not values:
                continue
            x, y = zip(*values)
            ax.plot(
                x,
                y,
                color=MESH_COLORS[mesh_kind],
                linestyle=linestyles[nprocs.index(nproc) % len(linestyles)],
                marker=markers[nprocs.index(nproc) % len(markers)],
                label=f"{mesh_kind}, n={nproc}",
            )

    for (mesh_kind, nproc), values in sorted(series.items()):
        if mesh_kind in MESH_COLORS:
            continue
        x, y = zip(*values)
        ax.plot(x, y, color="0.4", marker="o", label=f"{mesh_kind}, n={nproc}")

    ax.set_xscale("log") # , base=2)
    ax.set_yscale("log")
    ax.set_xlabel("Fine target length scale, lc (m)")
    ax.set_ylabel("Diagnostic solve time (s)")
    ax.set_title("Greenland diagnostic solve timings")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()

    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)


def main() -> None:
    """Plot diagnostic timing logs from command-line arguments."""
    args = parse_args()
    fig_path = args.fig_dir.expanduser() / args.output_name
    rows = read_timings(args.log_dir.expanduser())
    plot_timings(rows, fig_path)
    print(f"Saved {fig_path}")


if __name__ == "__main__":
    main()
