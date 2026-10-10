r"""relaxation.py - the relaxed initial state.

A production MAP is rewound one year with the Smith et al. (2020) mean dH/dt
(the backdating of issue #117), run forward to the geometry year at half the
production step on OCX's forcing for that year, with the apparent mass
balance off and the front pinned (``antarctica/scripts/relaxation/run.py``),
and re-inverted from the state the year ends in (``inversion_icepack2.py``
with that state as ``ISMIP7_WARM_START``). The MAP that comes out carries the
relaxed geometry, which serves the re-inversion alone: a forward on any mesh,
the MAP's own included, takes the MAP's controls onto that mesh's BedMachine
geometry, as it does for any MAP, so every resolution starts from the same
2015 geometry, which the 2003 start of the historicals and OCX backdates.

The rules the driver, the inversion and the forward share live here, free of
Firedrake, so the suite tests them directly.
"""
import os

from .handoff import OBJECTIVE_KEYS, OBJECTIVE_RECORD_KEYS
from .runconfig import DHDT_WINDOW, GEOMETRY_YEAR, dt as _dt

RELAX_START_DEFAULT = "2014"
RELAX_FORCINGS = ("ocx", "none")
RELAX_FORCING_DEFAULT = "ocx"

# geometry_source_method of a MAP re-inverted on a relaxed geometry
GEOMETRY_METHOD_RELAXED = "relaxed-forward-v1"
# Set only on the relaxation driver's own final checkpoint: the state a
# re-inversion starts from. Forward checkpoints never carry it on.
END_STATE_ATTR = "relaxation_end_state"
# How the relaxed geometry was made. Written on the end state, copied into the
# MAP re-inverted from it, and carried by every forward checkpoint after that.
RELAX_KEYS = (
    "relax_source_map", "relax_source_sha256", "relax_t_start", "relax_t_end",
    "relax_dt", "relax_backdate_years", "relax_forcing",
    "relax_source_objective_total", "relax_fssa_theta", "relax_fssa_reference",
)
# what a relaxed MAP adds: the end state it was re-inverted from
RELAX_MAP_KEYS = RELAX_KEYS + ("relax_state", "relax_state_sha256")

# init_state, which a forward records in every checkpoint it writes:
#   observed          a MAP inverted on BedMachine geometry
#   relaxed           a relaxed MAP on its own mesh, starting from its geometry:
#                     forwards cold-started before 10 October 2026 record it,
#                     and their restarts keep it; no cold start takes it now
#   relaxed-controls  a relaxed MAP's controls on the forward mesh's BedMachine
#                     geometry, on any mesh
INIT_STATES = ("observed", "relaxed", "relaxed-controls")
INIT_STATE_ATTR = "init_state"


def _value(value):
    r"""An attribute as a plain Python value: h5py hands back bytes and
    numpy scalars."""
    if isinstance(value, bytes):
        return value.decode()
    item = getattr(value, "item", None)
    return item() if callable(item) else value


def relax_window():
    r"""``(t_start, t_end)`` of the relaxation year: ``ISMIP7_RELAX_START``
    (default 2014) to the geometry year. A start inside the Smith window and
    before the geometry year; any other is refused."""
    t_start = float(os.environ.get("ISMIP7_RELAX_START", RELAX_START_DEFAULT))
    if not DHDT_WINDOW[0] <= t_start < GEOMETRY_YEAR:
        raise ValueError(
            f"ISMIP7_RELAX_START={t_start:g}: the relaxation starts before the "
            f"{GEOMETRY_YEAR:g} geometry and no earlier than {DHDT_WINDOW[0]:g}, "
            f"where the Smith mean dH/dt that rewinds it begins")
    return t_start, GEOMETRY_YEAR


def relax_dt(t_start, t_end):
    r"""Step of the relaxation year [yr]: ``ISMIP7_RELAX_DT``, by default half
    of ``ISMIP7_DT``, so it follows the production step. It has to divide the
    window, or the year would not end at the geometry year a re-inversion
    requires."""
    raw = os.environ.get("ISMIP7_RELAX_DT", "").strip()
    step = float(raw) if raw else 0.5 * _dt()
    if not step > 0.0:
        raise ValueError(f"the relaxation step must be positive, got {step:g}")
    n = (t_end - t_start) / step
    if abs(n - round(n)) > 1e-6:
        raise ValueError(
            f"a step of {step:g} yr does not divide {t_start:g} to {t_end:g}; "
            f"set ISMIP7_RELAX_DT to one that does")
    return step


def relax_backdate_years(t_start):
    r"""Years of the Smith mean dH/dt the relaxation undoes: the geometry year
    minus the start. ``ISMIP7_GEOMETRY_BACKDATE`` may repeat that number and
    may not change it."""
    years = GEOMETRY_YEAR - float(t_start)
    named = os.environ.get("ISMIP7_GEOMETRY_BACKDATE", "").strip()
    if named and abs(float(named) - years) > 1e-9:
        raise ValueError(
            f"ISMIP7_GEOMETRY_BACKDATE={named} but a relaxation from "
            f"{t_start:g} undoes {years:g} years; unset it")
    return years


def relax_forcing():
    r"""``ISMIP7_RELAX_FORCING``: ``ocx`` (default), OCX's own forcing for the
    year as ``ISMIP7_OCX_FORCING`` selects it, or ``none``, no SMB and no melt,
    a smoke test that needs no forcing data."""
    value = os.environ.get("ISMIP7_RELAX_FORCING", RELAX_FORCING_DEFAULT).strip().lower()
    if value not in RELAX_FORCINGS:
        raise ValueError(
            f"ISMIP7_RELAX_FORCING must be one of {RELAX_FORCINGS}, got {value!r}")
    return value


def relax_environment_problems(*, apparent_mb, output, calving, fixed_front):
    r"""What in the resolved run settings the relaxation year cannot run
    under, as sentences; empty when it can run. ``apparent_mb`` is
    ``runconfig.apparent_mb_mode()``, ``output`` ``ismip7_output()``,
    ``calving`` ``calving_law()`` and ``fixed_front`` ``fixed_front()``."""
    problems = []
    if apparent_mb is not None:
        problems.append(
            f"ISMIP7_APPARENT_MB={apparent_mb} would hold the year at the "
            f"rewound state; set ISMIP7_APPARENT_MB=0")
    if output:
        problems.append(
            "ISMIP7_OUTPUT=1 would write ISMIP7 files for a run that is no "
            "experiment; set ISMIP7_OUTPUT=0")
    if calving != "none":
        problems.append(
            f"ISMIP7_CALVING={calving}: the front stays pinned through the "
            f"year so the re-inversion keeps the MAP's extent; set "
            f"ISMIP7_CALVING=none")
    if not fixed_front:
        problems.append(
            "the front stays pinned through the year; set ISMIP7_FIXED_FRONT=1")
    return problems


def relax_experiment_name(source_map, tag=""):
    r"""The relaxation's experiment name, from the MAP it relaxes:
    ``relax_<MAP stem>``, with ``_<tag>`` for ``ISMIP7_RUN_TAG``."""
    stem = os.path.splitext(os.path.basename(source_map))[0]
    return f"relax_{stem}" + (f"_{tag}" if tag else "")


def relaxed_map_name(source_map, t_start):
    r"""Basename of the MAP re-inverted from a relaxation of ``source_map``
    started at ``t_start``: the source's stem with ``_relax<year>``."""
    stem = os.path.splitext(os.path.basename(source_map))[0]
    return f"{stem}_relax{int(round(float(t_start)))}.h5"


def is_relaxed(attrs):
    r"""Whether ``attrs`` (a checkpoint's attributes) describe a relaxed
    geometry: a relaxation's end state, or a MAP re-inverted on one."""
    return (bool(int(_value(attrs.get(END_STATE_ATTR, 0)) or 0))
            or str(_value(attrs.get("geometry_source_method", ""))) == GEOMETRY_METHOD_RELAXED)


def relaxation_attrs(*, source_map, source_sha256, source_attrs, t_start, t_end,
                     dt, backdate_years, forcing):
    r"""The attributes the relaxation writes on its final checkpoint.

    ``source_attrs`` are the source MAP's attributes. Its objective settings
    (``handoff.OBJECTIVE_KEYS``) go across under their own names, so the
    re-inversion holds the MAP's objective, its log-velocity weight included,
    exactly as a chained link does. The value the MAP recorded at its iterate
    stays behind: it belongs to the MAP's geometry, and a warm start without
    one skips the handoff gap check. It travels as
    ``relax_source_objective_total`` for the log. A MAP that is relaxed
    already is refused.
    """
    if is_relaxed(source_attrs):
        raise ValueError(
            f"{source_map} is relaxed already; relax the MAP it came from")
    out = {key: _value(source_attrs[key]) for key in OBJECTIVE_KEYS
           if key in source_attrs}
    out.update({
        END_STATE_ATTR: 1,
        "relax_source_map": os.path.realpath(source_map),
        "relax_source_sha256": str(source_sha256),
        "relax_t_start": float(t_start),
        "relax_t_end": float(t_end),
        "relax_dt": float(dt),
        "relax_backdate_years": float(backdate_years),
        "relax_forcing": str(forcing),
    })
    if "objective_total" in source_attrs:
        out["relax_source_objective_total"] = float(_value(source_attrs["objective_total"]))
    assert not set(out) & set(OBJECTIVE_RECORD_KEYS)
    return out


def relax_fssa_attrs(theta, reference):
    r"""The free-surface stabilization the year ran under, for the end state's
    record: the weight and the resolved reference a forward's ``setup_model``
    returns (``fssa_theta``, ``fssa_reference``). A weight of 0 is off, and
    its reference is recorded as ``none``."""
    theta = float(theta)
    return {"relax_fssa_theta": theta,
            "relax_fssa_reference": str(reference) if theta > 0 else "none"}


def end_state_problems(attrs, *, same_mesh, geometry_taken, map_out):
    r"""Why a relaxation end state cannot seed a re-inversion, as sentences;
    empty when it can. The year has to have reached its end, unstalled, and
    the re-inversion has to run on its mesh and take its geometry: the
    controls it fits belong to that geometry. ``map_out`` is the MAP the
    re-inversion writes, which has to carry ``relaxed_map_name``, so it never
    replaces the production MAP a forward loads by default."""
    problems = []
    t_end = float(_value(attrs.get("relax_t_end", GEOMETRY_YEAR)))
    t_yr = attrs.get("t_yr")
    if t_yr is None or abs(float(_value(t_yr)) - t_end) > 1e-6:
        problems.append(
            f"it stopped at t={_value(t_yr)}, before the end of its year "
            f"{t_end:g}; finish the relaxation first")
    if int(_value(attrs.get("stalled", 0)) or 0):
        problems.append("its solver stalled")
    if not same_mesh:
        problems.append(
            "it is on another mesh; run the re-inversion on the relaxation's "
            "mesh (ISMIP7_MESH=checkpoint)")
    if not geometry_taken:
        problems.append(
            "its geometry is not taken (ISMIP7_WARM_START_GEOMETRY=0, or a "
            "lake_ice_base or a raster_sample other than the run's)")
    want = relaxed_map_name(str(_value(attrs["relax_source_map"])),
                            _value(attrs["relax_t_start"]))
    if os.path.basename(map_out) != want:
        problems.append(
            f"it would write {os.path.basename(map_out)}; a MAP re-inverted "
            f"from it is named {want}, so set ISMIP7_MAP_OUT to a path ending "
            f"in {want}")
    return problems


def inherited_geometry(warm_attrs, *, geometry_taken, warm_basename, warm_sha256=None):
    r"""The geometry record an inversion writes when its warm start's geometry
    is relaxed: ``geometry_source``, ``geometry_source_method`` and the
    ``RELAX_MAP_KEYS``. None when the warm start's geometry is not relaxed,
    and the inversion keeps its BedMachine record.

    From an end state the record names that state (``warm_basename``,
    ``warm_sha256``) and copies its ``RELAX_KEYS``; from a relaxed MAP or a
    re-inversion's own checkpoint (the next link of a chain) it is copied
    whole. A relaxed warm start whose geometry the run does not take is
    refused: its controls were fitted to the relaxed geometry alone.
    """
    if not is_relaxed(warm_attrs):
        return None
    if not geometry_taken:
        raise RuntimeError(
            f"the warm start {warm_basename} carries a relaxed geometry and "
            f"its controls were fitted to it, but this run does not take that "
            f"geometry; run on its mesh with ISMIP7_MESH=checkpoint and leave "
            f"ISMIP7_WARM_START_GEOMETRY unset")
    if int(_value(warm_attrs.get(END_STATE_ATTR, 0)) or 0):
        out = {key: _value(warm_attrs[key]) for key in RELAX_KEYS if key in warm_attrs}
        out.update({
            "geometry_source": str(warm_basename),
            "geometry_source_method": GEOMETRY_METHOD_RELAXED,
            "relax_state": str(warm_basename),
            "relax_state_sha256": str(warm_sha256),
        })
        return out
    out = {key: _value(warm_attrs[key]) for key in RELAX_MAP_KEYS if key in warm_attrs}
    out["geometry_source"] = str(_value(warm_attrs.get("geometry_source", warm_basename)))
    out["geometry_source_method"] = GEOMETRY_METHOD_RELAXED
    return out


def init_state(map_attrs):
    r"""The ``INIT_STATES`` entry of a forward cold-started from a MAP with
    attributes ``map_attrs``. A relaxed MAP gives its controls alone, on its
    own mesh as on any other, so its forwards start from the geometry a
    forward from the MAP it was relaxed from starts from."""
    if str(_value(map_attrs.get("geometry_source_method", ""))) != GEOMETRY_METHOD_RELAXED:
        return "observed"
    return "relaxed-controls"


def describe_relaxation(attrs):
    r"""One line naming how a relaxed geometry was made, from its record."""
    def _get(key, default="?"):
        return _value(attrs.get(key, default))
    if "relax_fssa_theta" not in attrs:
        fssa = "FSSA unrecorded"
    elif float(_get("relax_fssa_theta")) > 0:
        fssa = (f"FSSA theta {float(_get('relax_fssa_theta')):g}, "
                f"{_get('relax_fssa_reference')} reference")
    else:
        fssa = "FSSA off"
    return (f"relaxed {float(_get('relax_t_start', 'nan')):g} to "
            f"{float(_get('relax_t_end', 'nan')):g} at dt "
            f"{float(_get('relax_dt', 'nan')):g} yr ({_get('relax_forcing')} forcing, "
            f"{fssa}) "
            f"from {os.path.basename(str(_get('relax_source_map')))}")


def anchor_ratio_counts(c_new, c_old, grounded, floor=1e-8):
    r"""This rank's share of the anchor change a relaxed geometry makes:
    ``ln R = ln(C_w0 relaxed / C_w0 at the MAP's geometry)`` over the grounded
    entries, as ``(n_grounded, n_over_0.1, n_over_0.3, sum |ln R|, min ln R,
    max ln R)``; the inversion reduces them over ranks. Under a kept theta the
    friction at the first evaluation is the MAP's times R on these cells."""
    import numpy as np
    g = np.asarray(grounded, dtype=bool)
    lnr = np.log(np.maximum(np.asarray(c_new, dtype=float)[g], floor)
                 / np.maximum(np.asarray(c_old, dtype=float)[g], floor))
    a = np.abs(lnr)
    return (int(g.sum()), int((a > 0.1).sum()), int((a > 0.3).sum()),
            float(a.sum()),
            float(lnr.min()) if lnr.size else np.inf,
            float(lnr.max()) if lnr.size else -np.inf)
