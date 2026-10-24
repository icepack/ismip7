r"""A calving screen's legend states each member's parameter as the law its
checkpoint recorded, so a member name on another scale (``hfb015`` is
0.15 MPa, ``hfb086`` is 0.086 MPa) cannot mislabel it."""
import importlib.util
import os

import pytest

FIGURES = os.path.join(os.path.dirname(__file__), "..", "antarctica", "scripts", "figures")


@pytest.fixture(scope="module")
def figures():
    spec = importlib.util.spec_from_file_location(
        "calving_screen_figures", os.path.join(FIGURES, "calving_screen_figures.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("member, law, params, label", [
    ("hfb086", "hfb", dict(sigma_max=0.086), "resistive stress, 0.086 MPa"),
    ("hfb015", "hfb", dict(sigma_max=0.15), "resistive stress, 0.15 MPa"),
    ("vm130", "vonmises", dict(sigma_max_fl=0.13, sigma_max_gr=1.0), "von Mises, 0.13 MPa"),
    ("vm130g", "vonmises", dict(sigma_max_fl=0.13, sigma_max_gr=0.13),
     "von Mises, 0.13 MPa (grounded too)"),
    ("eig5e6", "eigencalving", dict(K=4.9e6), "eigencalving, K = 4.9e6"),
    ("thk235", "thickness", dict(hc=235), "thickness, 235 m"),
])
def test_the_label_is_the_recorded_parameter(figures, member, law, params, label):
    from icepack_tools import calving
    described = calving.make(law, marine_only=True, **params).describe()
    assert figures.member_label(member, figures.law_of(member), described) == label


def test_a_member_without_a_record_is_labelled_by_its_name(figures):
    assert figures.member_label("hfb086", "hfb") == "resistive stress (hfb086)"
