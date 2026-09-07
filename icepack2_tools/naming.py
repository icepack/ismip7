r"""MAP filename construction, shared by the inversion, the forward and the
preflight/launch gates.

A MAP is only valid for the (friction law, flow exponent, geometry space) it
was inverted under - the inversion absorbs whatever those get wrong into
theta/phi - so each of them is in the filename. The name is built in one place
because every copy of the rule is a chance for a gate to look for a file the
run will never write (or to bless a MAP the run will never load).

Pure Python: importable without Firedrake so the preflight stays fast.
"""

from icepack2_tools.runconfig import geometry_space, n_flow

_FRICTION_TAGS = {"regularized_coulomb": "_rc", "budd": "_budd"}


def map_n_tag():
    r"""Filename tag distinguishing MAPs inverted at different flow exponents
    so n=3 and n=4 MAPs coexist on disk. n=4 keeps the legacy untagged name
    (backward compatible with the `antarctica` MAPs); any other n gets
    `_n<N>` (e.g. `_n3`)."""
    n = n_flow()
    return "" if abs(n - 4.0) < 1e-9 else f"_n{int(round(n))}"


def map_geom_tag():
    r"""Filename tag for the geometry discretization a MAP was inverted under.

    CG1 keeps the legacy untagged name (every MAP on disk before Aug 2026);
    DG0 gets `_dg0`. They are NOT interchangeable: the inversion absorbs the
    front treatment into theta/phi, so a CG1 MAP driven by a DG0 forward has
    friction tuned to a calving front that the lumped lift made ~2x too thick.
    Separate names stop that from happening by accident.
    """
    return "" if geometry_space() == "cg1" else "_dg0"


def map_tag(friction, geometry=True):
    r"""Full MAP tag: friction law + flow exponent + geometry space."""
    return (_FRICTION_TAGS.get(friction, "") + map_n_tag()
            + (map_geom_tag() if geometry else ""))


def map_basename(friction, lc, geometry=True):
    r"""MAP filename (no directory) for this configuration."""
    return f"inversion_icepack2{map_tag(friction, geometry)}_{int(lc)}.h5"
