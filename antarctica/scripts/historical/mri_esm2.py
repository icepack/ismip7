#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 2: Historical with MRI-ESM2-0 (1850-2014).

Usage:
    mpiexec -n 24 python scripts/historical/mri_esm2.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiment import run_core_experiment

if __name__ == "__main__":
    run_core_experiment(
        core=2, title="Historical with MRI-ESM2-0 (1850-2014)", name="hist_mri_esm2_0",
        esm="MRI-ESM2-0", scenario="historical",
        t_start_default=1850.0, t_end_default=2014.0,
        restart_from_hist=False,
    )
