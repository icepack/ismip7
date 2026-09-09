#!/usr/bin/env python3
r"""Tiny distributed initialization test for transient solver modes.

This is not a physics or convergence qualification.  It exercises the actual
CG1-vector / DG0-symmetric-tensor / DG0-vector mixed layout, exact local
condensation, GAMG/MUMPS setup, and persistent transport-solver reuse without
loading Antarctic data.
"""

import argparse
import os
import sys

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)

import firedrake as fd  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402

from icepack2_tools.solverconfig import (  # noqa: E402
    diagnostic_solver_parameters,
    transport_solver_parameters,
)


def mixed_problem(mesh):
    velocity = fd.VectorFunctionSpace(mesh, "CG", 1)
    dg0 = fd.FiniteElement("DG", "triangle", 0)
    stress = fd.TensorFunctionSpace(mesh, dg0, symmetry=True)
    traction = fd.VectorFunctionSpace(mesh, dg0)
    mixed = velocity * stress * traction
    state = fd.Function(mixed)
    u, membrane, basal = fd.split(state)
    v, q, r = fd.TestFunctions(mixed)
    residual = (
        fd.inner(fd.grad(u), fd.grad(v))
        + fd.inner(u, v)
        + fd.inner(
            membrane - fd.sym(fd.grad(u)),
            q - fd.sym(fd.grad(v)),
        )
        + fd.inner(basal - u, r - v)
        - fd.inner(fd.as_vector((1.0, 1.0)), v)
    ) * fd.dx
    # Match simulation.py's structural-zero blocks for retained-first SCPC.
    zero = fd.Constant(0.0)
    residual += fd.derivative(
        zero * membrane[0, 0] * basal[0] * fd.dx, state
    )
    return fd.NonlinearVariationalProblem(residual, state), state


def test_mode(mesh, mode):
    os.environ["ISMIP7_DIAGNOSTIC_LINEAR_SOLVER"] = mode
    problem, state = mixed_problem(mesh)
    solver = fd.NonlinearVariationalSolver(
        problem,
        solver_parameters=diagnostic_solver_parameters(),
        options_prefix=f"ismip7_smoke_{mode}_",
    )
    solver.solve()
    reason = solver.snes.getConvergedReason()
    if reason <= 0:
        raise RuntimeError(f"{mode} diverged with SNES reason {reason}")
    PETSc.Sys.Print(
        f"PASS {mode}: snes_its={solver.snes.getIterationNumber()} "
        f"linear_its={solver.snes.getLinearSolveIterations()}"
    )
    return state


def test_persistent_transport(mesh):
    space = fd.FunctionSpace(mesh, "DG", 0)
    thickness = fd.Function(space).assign(1.0)
    old = fd.Function(space)
    source = fd.Function(space)
    dt = fd.Constant(0.1)
    trial = fd.TrialFunction(space)
    test = fd.TestFunction(space)
    form = ((trial - old) / dt * test + trial * test - source * test) * fd.dx
    problem = fd.LinearVariationalProblem(
        fd.lhs(form), fd.rhs(form), thickness
    )
    solver = fd.LinearVariationalSolver(
        problem,
        solver_parameters=transport_solver_parameters(),
        options_prefix="ismip7_transport_smoke_",
    )
    expected = 1.0
    for value, timestep in ((1.0, 0.1), (2.0, 0.05)):
        old.assign(thickness)
        source.assign(value)
        dt.assign(timestep)
        solver.solve()
        expected = (expected + timestep * value) / (1.0 + timestep)
        target = fd.Function(space).assign(expected)
        if fd.errornorm(target, thickness) > 1e-12:
            raise RuntimeError("persistent transport solver used stale coefficients")
    PETSc.Sys.Print(
        "PASS persistent transport solver: source and dt coefficient updates"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "modes",
        nargs="*",
        default=["scpc_gamg", "scpc_mumps", "full_mumps"],
    )
    args = parser.parse_args()
    mesh = fd.UnitSquareMesh(2, 2)
    for mode in args.modes:
        test_mode(mesh, mode)
    test_persistent_transport(mesh)


if __name__ == "__main__":
    main()
