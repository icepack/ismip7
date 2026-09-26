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
