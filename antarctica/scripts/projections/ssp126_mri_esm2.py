#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 6: SSP1-2.6 with MRI-ESM2-0 (2015-2300).

Usage:
    mpiexec -n 24 python scripts/projections/ssp126_mri_esm2.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiment import run_core_experiment

if __name__ == "__main__":
    run_core_experiment(
        core=6, title="SSP1-2.6 with MRI-ESM2-0 (2015-2300)", name="ssp126_mri_esm2_0",
        esm="MRI-ESM2-0", scenario="ssp126",
        t_start_default=2015.0, t_end_default=2300.0,
        restart_from_hist=True,
    )
