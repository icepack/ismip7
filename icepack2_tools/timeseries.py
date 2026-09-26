r"""The year column of a run's budget timeseries, written and read back.

``run_simulation`` appends one row per step to ``<exp>_<lc>_timeseries.csv``
and, on a warm restart, keeps the rows up to the resume year before it
appends. The audits (``check_ismip6_track.py``, ``compare_runs.py``, and
``core_report.py`` through the first) recover the step from the year column,
because ``calv_gt``, ``clamp_gt`` and ``resid_gt`` are Gt per step and become
rates only once divided by it.

The year is written to six decimals. At one decimal a step under 0.1 yr
repeats years across rows (three rows in four at the production step of
0.025), so the median step read back is 0 or 0.1, and the resume trim, which
compares years to half a step, drops rows of the step it resumes at and keeps
rows the restart writes again. Pure Python, so these rules are tested without
Firedrake.
"""

__all__ = ["YEAR_DECIMALS", "format_year", "rows_kept_on_resume", "step_from_years"]

YEAR_DECIMALS = 6


def format_year(t):
    r"""The year field of a timeseries row."""
    return f"{t:.{YEAR_DECIMALS}f}"


def rows_kept_on_resume(lines, header, t_start, dt):
    r"""The lines a warm restart at ``t_start`` keeps of an existing series.

    ``lines`` is the file as read, header first; ``header`` stands in when
    the file has none, so the first entry of the result is always the header
    the resumed series carries. A row is kept when its year is at or before
    ``t_start``, to half a step; the rows after it are the restart's to write
    again. A row whose year does not parse is dropped.
    """
    kept = [lines[0]] if lines and lines[0].startswith("year") else [header]
    for line in lines[1:]:
        try:
            if float(line.split(",", 1)[0]) <= t_start + 0.5 * dt:
                kept.append(line)
        except (ValueError, IndexError):
            pass
    return kept


def step_from_years(years):
    r"""The time step of a series, from its year column: the span over the
    number of steps.

    Exact for a series written one row per step at full precision, and right
    to a rounding of the span for an older series written to one decimal,
    where neighbouring rows differ by 0 or 0.1 and their median is wrong.
    Raises ``ValueError`` for a series too short or too flat to carry a step.
    """
    n = len(years)
    if n < 2:
        raise ValueError(f"a step needs at least two rows, and the series has {n}")
    span = float(years[-1]) - float(years[0])
    if not span > 0.0:
        raise ValueError(f"the series spans {span:g} yr, so it carries no step")
    return span / (n - 1)
