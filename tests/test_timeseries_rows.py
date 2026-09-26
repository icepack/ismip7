r"""The timeseries year column at the production step, 0.025 yr (issue 20).

``run_simulation`` writes one row per step, trims the series back to the
resume year on a warm restart, and the audits read the step back from the
year column. With the year written to one decimal, three rows in four repeat
a year at 0.025: the track audit inferred a step of 0 and divided by it, and
the trim dropped the rows of the step a chain link resumed at. These are the
rules in ``icepack2_tools.timeseries``, exercised on the cases that failed.

Serial, no Firedrake, no data files.
"""
import os
import re
import subprocess
import sys

import pytest

from icepack2_tools.timeseries import (
    format_year, rows_kept_on_resume, step_from_years,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRACK = os.path.join(REPO, "antarctica", "scripts", "check_ismip6_track.py")
HEADER = "year,vaf_mm_sle\n"
DT = 0.025


def years(start, stop, dt):
    n = int(round((stop - start) / dt))
    return [start + k * dt for k in range(1, n + 1)]


def series(start, stop, dt, fmt=format_year):
    return [HEADER] + [f"{fmt(t)},0.0\n" for t in years(start, stop, dt)]


def test_every_row_at_the_production_step_has_its_own_year():
    written = [format_year(t) for t in years(2015.0, 2016.0, DT)]
    assert len(written) == 40 and len(set(written)) == 40
    assert written[0] == "2015.025000" and written[-1] == "2016.000000"


@pytest.mark.parametrize("t_start, t_written", [
    (2015.675, 2015.75),   # a wall-clock stop mid-year, rows past it on disk
    (2020.0, 2020.1),      # the five-yearly checkpoint, the job ran on after it
])
def test_a_resume_keeps_every_row_to_its_own_year_and_none_after(t_start, t_written):
    kept = rows_kept_on_resume(series(2015.0, t_written, DT), HEADER, t_start, DT)
    assert kept[0] == HEADER
    got = [float(line.split(",")[0]) for line in kept[1:]]
    assert got == pytest.approx(years(2015.0, t_start, DT))
    assert got[-1] == pytest.approx(t_start)


def test_a_resume_supplies_a_missing_header_and_drops_unreadable_rows():
    lines = ["garbled\n", "2015.025000,0.0\n", "not a row\n", "2015.050000,0.0\n",
             "2015.075000,0.0\n"]
    kept = rows_kept_on_resume(lines, HEADER, 2015.05, DT)
    assert kept == [HEADER, "2015.025000,0.0\n", "2015.050000,0.0\n"]
    assert rows_kept_on_resume([], HEADER, 2015.05, DT) == [HEADER]


def test_the_step_is_read_from_a_full_precision_series():
    col = [float(format_year(t)) for t in years(1850.0, 2015.0, DT)]
    assert step_from_years(col) == pytest.approx(DT, rel=1e-9)


@pytest.mark.parametrize("dt", [0.05, 0.025])
def test_the_step_is_read_from_a_legacy_one_decimal_series(dt):
    r"""The median of neighbouring differences reads such a series as 0.1 at
    dt 0.05 and 0 at dt 0.025; the span over the step count does not."""
    col = [float(f"{t:.1f}") for t in years(2015.0, 2025.0, dt)]
    assert step_from_years(col) == pytest.approx(dt, rel=1e-2)


@pytest.mark.parametrize("col", [[], [2015.025], [2015.0, 2015.0, 2015.0]])
def test_a_series_that_carries_no_step_is_refused(col):
    with pytest.raises(ValueError):
        step_from_years(col)


def _budget_csv(path, dt, fmt):
    cols = ("year,vaf_mm_sle,mass_gt,smb_gtyr,melt_gtyr,outflux_gtyr,"
            "calv_gt,clamp_gt,resid_gt,amb_gtyr\n")
    with open(path, "w") as f:
        f.write(cols)
        for k, t in enumerate(years(2015.0, 2020.0, dt)):
            mass = 2.4e7 - 100.0 * k * dt
            f.write(f"{fmt(t)},57000.0,{mass:.2f},2500.0,1100.0,1200.0,"
                    f"{1200.0 * dt:.4f},0.0,0.0,0.0\n")


@pytest.mark.parametrize("dt, fmt, rel", [
    (0.025, format_year, 1e-9),
    # one decimal, as every series before this change: the median read 0.1
    (0.05, lambda t: f"{t:.1f}", 2e-2),
    # and 0 here, which the audit divided by
    (0.025, lambda t: f"{t:.1f}", 2e-2),
])
def test_the_track_audit_reads_the_step_without_being_told(tmp_path, dt, fmt, rel):
    r"""``core_report.py`` runs the audit with no ``--dt``."""
    csv_fn = tmp_path / "ctrl_1000_timeseries.csv"
    _budget_csv(csv_fn, dt, fmt)
    r = subprocess.run([sys.executable, TRACK, str(csv_fn)],
                       capture_output=True, text=True)
    assert "Traceback" not in r.stderr, r.stderr
    assert r.returncode in (0, 1), r.stdout + r.stderr
    (read,) = re.findall(r"dt=(\S+) yr", r.stdout)
    assert float(read) == pytest.approx(dt, rel=rel)
