#!/usr/bin/env python3
"""
Plot generated Antarctic meshes using Firedrake's triplot.

Uses the same mesh visualization approach as the icepack tutorials.

Usage:
    python scripts/plot_meshes.py
"""

import os
import re
import glob
import firedrake
import matplotlib.pyplot as plt
from icepack.plot import subplots

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")

_BASENAME_RE = re.compile(r"antarctica_(\d+)_(\d+)(?:_buffered(\d+))?$")


def parse_mesh_basename(basename):
    """Split 'antarctica_<COARSE>_<FINE>[_buffered<BUFFER_M>]' into ints.

    The `_buffered<BUFFER_M>` tag postdates most of the meshes on disk, so it
    is optional; buffer_m is None for a name that doesn't carry it.

    Returns (coarse_m, fine_m, buffer_m).
    """
    m = _BASENAME_RE.match(basename)
    if not m:
        raise ValueError(
            f"Mesh filename '{basename}' doesn't match "
            "antarctica_<COARSE>_<FINE>[_buffered<BUFFER_M>] (see mesh_naming.py)"
        )
    coarse_m, fine_m = int(m.group(1)), int(m.group(2))
    buffer_m = None if m.group(3) is None else int(m.group(3))
    return coarse_m, fine_m, buffer_m


def buffer_label(buffer_m):
    """Buffer tag for a plot title; legacy names carry no buffer."""
    if buffer_m is None:
        return "buffer unknown"
    return f"buffer {buffer_m // 1000} km"


def report_skipped(skipped):
    """Report every unparseable name once, rather than aborting on the first."""
    if not skipped:
        return
    print(f"Skipped {len(skipped)} mesh file(s) with unrecognized names:")
    for basename in skipped:
        print(f"  {basename}")


def plot_single_mesh(mesh_fn, ax, title):
    """Plot a single mesh on the given axes using firedrake.triplot."""
    mesh = firedrake.Mesh(mesh_fn)
    firedrake.triplot(
        mesh, axes=ax, interior_kw={"linewidth": 0.1}, boundary_kw={"linewidth": 0.8}
    )
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    # Find all mesh files, sorted by resolution (coarsest first)
    mesh_files = sorted(glob.glob(os.path.join(MESH_DIR, "antarctica_*.msh")))
    if not mesh_files:
        print("No mesh files found. Run mesh_antarctica.py first.")
        return

    plots, skipped = [], []
    for mesh_fn in mesh_files:
        basename = os.path.splitext(os.path.basename(mesh_fn))[0]
        try:
            coarse_m, fine_m, buffer_m = parse_mesh_basename(basename)
        except ValueError:
            skipped.append(basename)
            continue
        plots.append((mesh_fn, basename, coarse_m, fine_m, buffer_m))

    if not plots:
        print(f"Found {len(mesh_files)} mesh file(s), none with a recognized name.")
        report_skipped(skipped)
        return

    n_meshes = len(plots)
    print(f"Found {n_meshes} mesh file(s)")

    # Multi-panel figure with all resolutions
    fig, axes = subplots(1, n_meshes, figsize=(5 * n_meshes, 5))
    if n_meshes == 1:
        axes = [axes]

    for ax, (mesh_fn, basename, coarse_m, fine_m, buffer_m) in zip(axes, plots):
        title = (
            f"{fine_m // 1000}–{coarse_m // 1000} km ({buffer_label(buffer_m)})"
        )

        print(f"  Plotting {basename}...")
        plot_single_mesh(mesh_fn, ax, title)

    fig.suptitle("Antarctica Adaptive Meshes (EPSG:3031)", fontsize=14, y=1.02)
    fig.tight_layout()

    out_fn = os.path.join(FIG_DIR, "antarctica_meshes.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    print(f"Saved: {out_fn}")
    plt.close(fig)

    # Also make individual plots for each mesh (higher detail)
    for mesh_fn, basename, coarse_m, fine_m, buffer_m in plots:
        fine_km, coarse_km = fine_m // 1000, coarse_m // 1000

        fig, ax = subplots(1, 1, figsize=(10, 10))
        title = (
            f"Antarctica Mesh: {fine_km}–{coarse_km} km "
            f"({buffer_label(buffer_m)})"
        )
        plot_single_mesh(mesh_fn, ax, title)

        out_fn = os.path.join(FIG_DIR, f"{basename}.png")
        fig.savefig(out_fn, dpi=200, bbox_inches="tight")
        print(f"Saved: {out_fn}")
        plt.close(fig)

    report_skipped(skipped)


if __name__ == "__main__":
    main()
