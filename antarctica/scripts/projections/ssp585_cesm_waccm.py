#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 7: SSP5-8.5 with CESM2-WACCM (2015-2300).

Usage:
    mpiexec -n 24 python scripts/projections/ssp585_cesm_waccm.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiment import run_core_experiment

if __name__ == "__main__":
    run_core_experiment(
        core=7, title="SSP5-8.5 with CESM2-WACCM (2015-2300)", name="ssp585_cesm2_waccm",
        esm="CESM2-WACCM", scenario="ssp585",
        t_start_default=2015.0, t_end_default=2300.0,
        restart_from_hist=True,
    )
