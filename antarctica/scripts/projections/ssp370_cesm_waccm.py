#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 3: SSP3-7.0 with CESM2-WACCM (2015-2100).

Usage:
    mpiexec -n 24 python scripts/projections/ssp370_cesm_waccm.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiment import run_core_experiment

if __name__ == "__main__":
    run_core_experiment(
        core=3, title="SSP3-7.0 with CESM2-WACCM (2015-2100)", name="ssp370_cesm2_waccm",
        esm="CESM2-WACCM", scenario="ssp370",
        t_start_default=2015.0, t_end_default=2100.0,
        restart_from_hist=True,
    )
