#!/usr/bin/env python3
"""
Rewrite an inversion checkpoint on a single core for partition-independent reuse.

Firedrake checkpoint files written under one MPI partition layout may not load
cleanly on a different number of ranks. This script loads the full checkpoint
(mesh + all saved fields) serially and rewrites it atomically so downstream
multi-rank jobs can load it on any core count.

Usage:
    python scripts/redistribute_checkpoint.py --lc 2500
    python scripts/redistribute_checkpoint.py --input mesh/in.h5 --output mesh/out.h5
"""

import argparse
import os

import firedrake as fd
from firedrake import COMM_WORLD
from firedrake.petsc import PETSc

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MESH_DIR = os.path.join(_ROOT, "mesh")


def default_checkpoint(lc):
    """Return the MAP filename used by the current inversion defaults."""
    friction = os.environ.get("ISMIP7_FRICTION", "budd")
    friction_tag = {"regularized_coulomb": "_rc", "budd": "_budd"}.get(
        friction, ""
    )
    n_flow = float(os.environ.get("ISMIP7_N_FLOW", "3.0"))
    n_tag = "" if abs(n_flow - 4.0) < 1e-9 else f"_n{int(round(n_flow))}"
    return os.path.join(
        MESH_DIR, f"inversion_icepack2{friction_tag}{n_tag}_{lc}.h5"
    )

# Fields saved by inversion_icepack2.py and simulation.py. Optional fields are
# skipped so older checkpoints can still be rewritten, but every field needed
# by a current MAP or restart is preserved when present.
_CHECKPOINT_FIELDS = (
    "log_friction",
    "log_fluidity",
    "velocity_obs",
    "obs_mask",
    "thickness",
    "bed",
    "surface",
    "fluidity_prior",
    "velocity",
    "membrane_stress",
    "basal_stress",
    "H_init",
    "phi_eff",
    "C_w0",
    "N_ref",
    "a_ref_mb",
    "thickness_dg",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lc",
        type=int,
        default=int(os.environ.get("ISMIP7_LC", "2500")),
        help="LC tag used to resolve the default checkpoint path",
    )
    parser.add_argument(
        "--input",
        default=None,
        help=("Input checkpoint .h5 (default: the MAP path selected by "
              "ISMIP7_FRICTION and ISMIP7_N_FLOW)"),
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output checkpoint .h5 (default: same as --input, in-place rewrite)",
    )
    return parser.parse_args()


def main():
    if COMM_WORLD.size != 1:
        raise RuntimeError(
            "redistribute_checkpoint.py must run on a single MPI rank "
            f"(got {COMM_WORLD.size}). Run without mpiexec/srun."
        )

    args = parse_args()
    in_fn = args.input or default_checkpoint(args.lc)
    out_fn = args.output or in_fn

    if not os.path.isfile(in_fn):
        raise FileNotFoundError(f"Checkpoint not found: {in_fn}")

    PETSc.Sys.Print(f"Loading checkpoint: {in_fn}")
    with fd.CheckpointFile(in_fn, "r") as chk:
        mesh = chk.load_mesh()
        root_attrs = {
            key: value
            for key, value in chk.attributes("/").items()
            if key != "dmplex_storage_version"
        }
        loaded = {}
        for name in _CHECKPOINT_FIELDS:
            try:
                loaded[name] = chk.load_function(mesh, name=name)
            except (KeyError, RuntimeError, ValueError):
                PETSc.Sys.Print(f"  (field '{name}' not present, skipping)")

    if "log_friction" not in loaded or "log_fluidity" not in loaded:
        raise ValueError(
            f"Checkpoint {in_fn} is missing log_friction/log_fluidity "
            "(was it produced by inversion_icepack2.py?)"
        )

    tmp_fn = out_fn + ".tmp"
    PETSc.Sys.Print(f"Writing redistributed checkpoint: {out_fn}")
    with fd.CheckpointFile(tmp_fn, "w") as chk:
        chk.save_mesh(mesh)
        for name, fn in loaded.items():
            chk.save_function(fn, name=name)
        for key, value in root_attrs.items():
            chk.set_attr("/", key, value)

    os.replace(tmp_fn, out_fn)
    PETSc.Sys.Print(f"Done: {out_fn}")


if __name__ == "__main__":
    main()
