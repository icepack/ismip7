r"""The production setting a run gets with no knob set (issue 20).

The group chose the submission mesh and step on 25 September 2026: the
1000 m / 10 km gmsh pair with a 20 km buffer,
``antarctica_10000_1000_buffered20000``, at dt 0.025 yr. The shell side of the
mesh is pinned by ``test_site_core.py`` and the runner's step by
``test_projection_chain.py``; this pins what a script or a driver run by hand
resolves through ``icepack2_tools.runconfig`` and ``mesh_naming``.

Serial, no Firedrake, no data files.
"""
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

from icepack2_tools.runconfig import dt, lc, lc_coarse  # noqa: E402
from mesh_naming import (  # noqa: E402
    DEFAULT_BUFFER_M, bndids_filename, get_buffer_m, mesh_basename,
)

PRODUCTION_MESH = "antarctica_10000_1000_buffered20000"


@pytest.fixture
def unset(monkeypatch):
    for name in ("ISMIP7_DT", "ISMIP7_LC", "ISMIP7_LC_COARSE", "ISMIP7_BUFFER_M"):
        monkeypatch.delenv(name, raising=False)
    return monkeypatch


def test_the_default_mesh_is_the_production_mesh(unset):
    assert get_buffer_m() == DEFAULT_BUFFER_M
    assert mesh_basename(lc_coarse(), lc(), get_buffer_m()) == PRODUCTION_MESH


def test_the_production_mesh_carries_a_tracked_sidecar(unset):
    sidecar = bndids_filename(lc_coarse(), lc(), get_buffer_m())
    assert os.path.basename(sidecar) == f"boundary_ids_{PRODUCTION_MESH}.json"
    subprocess.run(["git", "-C", REPO, "ls-files", "--error-unmatch", sidecar],
                   check=True, capture_output=True)


def test_the_default_step_is_the_production_step(unset):
    assert dt() == 0.025


def test_a_named_step_wins(unset):
    unset.setenv("ISMIP7_DT", "0.1")
    assert dt() == 0.1


# --- one outline buffer -------------------------------------------------------
# The outline extraction defaulted to 0 while every mesh name defaulted to
# 20000, so a bare call built an unbuffered outline for a buffered name.

from icepack2_tools.runconfig import (  # noqa: E402
    BUFFER_M_DEFAULT, buffer_m, mesh_build_check,
)
from mesh_naming import (  # noqa: E402
    buffer_from_name, mesh_stem, resolve_outline_buffer,
)


def test_the_names_and_the_outline_read_one_buffer_default(unset):
    assert buffer_m() == DEFAULT_BUFFER_M == float(BUFFER_M_DEFAULT) == 20000.0
    unset.setenv("ISMIP7_BUFFER_M", "0")
    assert buffer_m() == get_buffer_m() == 0.0


@pytest.mark.parametrize("name, buffer", [
    (PRODUCTION_MESH, 20000.0),
    (f"antarctica/mesh/{PRODUCTION_MESH}.msh", 20000.0),
    ("antarctica_5000_2000_buffered0", 0.0),
    (f"{PRODUCTION_MESH}_adapt1", 20000.0),
    ("antarctica_64000_2500_buffered", None),   # legacy, no size
    ("antarctica_ua_180000_2000_obs", None),
    ("antarctica_320000_32000", None),
])
def test_a_mesh_name_records_its_buffer_or_says_nothing(name, buffer):
    assert buffer_from_name(name) == buffer


def test_a_new_mesh_takes_the_first_buffer_that_is_known():
    r"""adapt_mesh.py: a checkpoint's record, else its name's tag, else a
    named ISMIP7_BUFFER_M, else nothing, which it refuses."""
    assert resolve_outline_buffer(0.0, PRODUCTION_MESH, "5000") == 0.0
    assert resolve_outline_buffer(None, PRODUCTION_MESH, "5000") == 20000.0
    assert resolve_outline_buffer(None, "antarctica_ua_180000_2000", "5000") == 5000.0
    assert resolve_outline_buffer(None, "antarctica_ua_180000_2000") is None


def test_the_outline_is_buffered_by_the_shared_default(unset):
    r"""A synthetic 400 km ice square at 4 km pixels: with no buffer named
    the outline grows by 20 km on each side, and an explicit 0 keeps it."""
    pytest.importorskip("gmsh")
    import numpy as np
    from icepack2_tools.mesh import extract_ice_outline
    x = y = np.arange(0.0, 1.0e6, 4.0e3)
    mask = np.zeros((y.size, x.size), dtype=np.int8)
    mask[50:150, 50:150] = 2                     # grounded ice
    bare = extract_ice_outline(mask, x, y, buffer_m=0.0)
    grown = extract_ice_outline(mask, x, y)
    side = np.sqrt(bare.area)
    assert side == pytest.approx(4.0e5, rel=0.05)
    assert np.sqrt(grown.area) == pytest.approx(side + 2 * 2.0e4, rel=0.02)
    unset.setenv("ISMIP7_BUFFER_M", "0")
    assert extract_ice_outline(mask, x, y).area == pytest.approx(bare.area)


def test_mesh_names_compare_without_directory_or_extension():
    assert mesh_stem(f"/a/b/{PRODUCTION_MESH}.msh") == mesh_stem(PRODUCTION_MESH)
    assert mesh_stem(f"{PRODUCTION_MESH}.msh") == PRODUCTION_MESH


# --- the mesh build check's flag --------------------------------------------

@pytest.mark.parametrize("value, check", [
    (None, True), ("", True), ("0", False), ("1", True),
])
def test_the_build_check_flag_reads_like_auto_resume(unset, value, check):
    if value is None:
        unset.delenv("ISMIP7_MESH_BUILD_CHECK", raising=False)
    else:
        unset.setenv("ISMIP7_MESH_BUILD_CHECK", value)
    assert mesh_build_check() is check


def test_a_flag_that_is_not_an_integer_is_an_error(unset):
    unset.setenv("ISMIP7_MESH_BUILD_CHECK", "no")
    with pytest.raises(ValueError, match="ISMIP7_MESH_BUILD_CHECK"):
        mesh_build_check()


def test_mesh_naming_imports_from_outside_the_repo(tmp_path):
    r"""``diagnostic_solve.py`` and any script run from ``antarctica/`` import
    ``mesh_naming`` with only the scripts directory on ``sys.path``."""
    scripts = os.path.join(REPO, "antarctica", "scripts")
    code = (
        "import sys\n"
        "sys.path[:] = [p for p in sys.path[1:]\n"
        "               if 'site-packages' not in p and not p.startswith(%r)]\n"
        "sys.path.insert(0, %r)\n"
        "import mesh_naming\n"
        "print(mesh_naming.DEFAULT_BUFFER_M)\n" % (REPO, scripts)
    )
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    r = subprocess.run([sys.executable, "-c", code], cwd=tmp_path, env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert float(r.stdout) == DEFAULT_BUFFER_M
