r"""Two builds of one mesh name are told apart before a checkpoint is moved
across them (issue 20).

Rice's and IU's builds of ``antarctica_10000_1000_buffered20000`` have
1,869,252 and 1,869,088 vertices. A forward names its mesh with
``ISMIP7_MESH`` and interpolates the MAP onto it, so a site's own build ran
the MAP through a transfer without a word; ``setup_model`` now refuses the
same name with another triangulation (``ISMIP7_MESH_BUILD_CHECK=0`` allows
it). ``meshes_match`` is the test: counts, then two centroid sums, since
counts alone and ``|c|^2`` alone both take a flipped diagonal for the same
mesh.
"""
import pytest

firedrake = pytest.importorskip("firedrake")
from firedrake import RectangleMesh, UnitSquareMesh  # noqa: E402

from icepack2_tools.transfer import meshes_match  # noqa: E402


def _polar(diagonal):
    r"""A patch placed where the polar stereographic coordinates put the
    ice, millions of metres from the origin."""
    return RectangleMesh(6, 6, 2.0e6, 2.1e6, originX=1.9e6, originY=2.0e6,
                         diagonal=diagonal)


def test_one_build_matches_itself_rebuilt():
    assert meshes_match(UnitSquareMesh(4, 4), UnitSquareMesh(4, 4))
    assert meshes_match(_polar("left"), _polar("left"))


def test_a_flipped_diagonal_is_another_build_with_the_same_counts():
    left, right = UnitSquareMesh(4, 4, diagonal="left"), UnitSquareMesh(4, 4, diagonal="right")
    assert left.cell_set.size == right.cell_set.size
    assert not meshes_match(left, right)
    assert not meshes_match(_polar("left"), _polar("right"))


def test_a_mesh_centred_on_the_pole_is_judged_on_the_size_of_its_terms():
    r"""Across the origin c_x c_y changes sign and its sum cancels, so the
    comparison is scaled by the |c|^2 sum, the size of the terms."""
    def pole(diagonal):
        return RectangleMesh(6, 6, 2.0e6, 2.0e6, originX=-1.0e6, originY=-1.0e6,
                             diagonal=diagonal)
    assert meshes_match(pole("left"), pole("left"))
    assert not meshes_match(pole("left"), pole("right"))


def test_different_counts_never_match():
    assert not meshes_match(UnitSquareMesh(4, 4), UnitSquareMesh(5, 4))
