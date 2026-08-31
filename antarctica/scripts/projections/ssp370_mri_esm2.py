#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 4: SSP3-7.0 with MRI-ESM2-0 (2015-2100).

Usage:
    mpiexec -n 24 python scripts/projections/ssp370_mri_esm2.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiment import run_core_experiment

if __name__ == "__main__":
    run_core_experiment(
        core=4, title="SSP3-7.0 with MRI-ESM2-0 (2015-2100)", name="ssp370_mri_esm2_0",
        esm="MRI-ESM2-0", scenario="ssp370",
        t_start_default=2015.0, t_end_default=2100.0,
        restart_from_hist=True,
    )
