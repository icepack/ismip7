r"""Reference-climate window and pool scenario for the ISMIP7 runs.

The control's SMB climatology and the projections' aSMB re-reference pool must
be built from the SAME (window, scenario) pair. If they disagree, projection
minus control differences two different baselines - the same class of error as
branching a control from a cold start instead of the historical.

The control run, the projection driver, the preflight gate and the per-core
report all need these values, so they are owned here instead of being copied
into each reader. A copy that drifts silently reintroduces the baseline
mismatch, and the report - the only committed record of a run - would then
document an effective value the run never used.

The pool scenario is ssp126, per the ISMIP7 protocol cheat sheet (April 2026):
"The climatology (2000-2029) should be created using a combination of the
historical and the SSP126 simulations" (discussion #28: last 15 yr historical
+ first 15 yr ssp126).

Pure Python: importable without Firedrake, so the preflight and the report
stay fast.
"""

import os

CLIM_START_DEFAULT = "2000"
CLIM_END_DEFAULT = "2029"
CLIM_SCENARIO_DEFAULT = "ssp126"


def clim_start():
    r"""First year of the reference-climate window."""
    return int(os.environ.get("ISMIP7_CLIM_START", CLIM_START_DEFAULT))


def clim_end():
    r"""Last year of the reference-climate window (inclusive)."""
    return int(os.environ.get("ISMIP7_CLIM_END", CLIM_END_DEFAULT))


def clim_scenario():
    r"""Scenario pooled with ``historical`` to build the climatology."""
    return os.environ.get("ISMIP7_CLIM_SCENARIO", CLIM_SCENARIO_DEFAULT)


CLIM_POOL_MARKER = "SMB climatology pool:"


def clim_pool_missing(years):
    r"""Years of the reference-climate window absent from ``years``."""
    return sorted(set(range(clim_start(), clim_end() + 1)) - set(years))


def describe_clim_pool(years, source):
    r"""One greppable line describing the pool a run ACTUALLY built.

    A run whose pool covers 15 of the window's 30 years re-references its
    anomalies to a different mean than one with full coverage, while both are
    differenced against the same control. A partial pool is a degraded run,
    not necessarily an invalid one, so the runtimes warn rather than refuse -
    which only works if the effective coverage survives into the record.
    ``core_report.py`` lifts every line carrying CLIM_POOL_MARKER out of the
    run log and into the report, so the two are afterwards distinguishable.

    The status and the numerator both come from ``clim_pool_missing``, so a
    caller that passes unfiltered years cannot get a COMPLETE line printed
    directly above the PARTIAL warning the same year list triggers.
    """
    want = clim_end() - clim_start() + 1
    missing = clim_pool_missing(years)
    in_window = want - len(missing)
    have = sorted(set(years))
    span = f"{have[0]}-{have[-1]}" if have else "no years"
    status = "COMPLETE" if not missing else "EMPTY" if not in_window else "PARTIAL"
    return (f"{CLIM_POOL_MARKER} {status} {in_window}/{want} yr, {span} "
            f"(historical+{clim_scenario()}, window "
            f"{clim_start()}-{clim_end()}, {source})")
