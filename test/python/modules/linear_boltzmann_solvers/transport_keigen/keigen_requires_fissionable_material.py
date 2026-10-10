#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
k-eigenvalue solvers require fissionable material, checked on every rank.

A 1D slab (reflecting on both ends, 20 cells) is used with two material layouts:

1. No fissionable material: PowerIterationKEigenSolver and NonLinearKEigenSolver must reject the
   problem with a ValueError at Initialize() on every rank. Without the rule, the solve fails with
   zero fission production after work has been done.
2. Fissionable material only in the cells near zmax, which one rank owns in a 2-rank partition:
   the requirement is reduced across ranks, so both solvers must accept the problem on every rank
   and power iteration must return a finite, positive k.

Checks:
  KEIGEN_FISSION_REQUIREMENT_PASSED: number of cases (of 4) that behave as described.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
    )
    from pyopensn.xs import MultiGroupXS

REQUIREMENT = "requires fissionable material"


def make_problem(grid, xs_map):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
            }
        ],
        xs_map=xs_map,
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )


def make_solvers(problem):
    return [
        PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-8, max_iters=200),
        NonLinearKEigenSolver(problem=problem, nl_abs_tol=1.0e-8, nl_rel_tol=1.0e-8),
    ]


def initialize_outcome(solver):
    """Returns 'rejected' for a ValueError naming the requirement, 'accepted' if no error."""
    try:
        solver.Initialize()
    except ValueError as error:
        return "rejected" if REQUIREMENT in str(error) else "other"
    return "accepted"


def on_all_ranks(condition):
    return int(MPIAllReduce(int(condition), "min"))


if __name__ == "__main__":
    nodes = [i * 10.0 / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=8.0, zmax=10.0), 1, True
    )

    absorber = MultiGroupXS()
    absorber.CreateSimpleOneGroup(1.0, 0.5)
    fissile = MultiGroupXS()
    fissile.LoadFromOpenSn(
        os.path.join(os.path.dirname(__file__), "../../../../assets/xs/xs1g_delayed_sub_1p.cxs")
    )

    passed = 0

    # Case 1
    no_fission = make_problem(grid, [{"block_ids": [0, 1], "xs": absorber}])
    for solver in make_solvers(no_fission):
        passed += on_all_ranks(initialize_outcome(solver) == "rejected")

    # Case 2
    partial_fission = make_problem(
        grid,
        [{"block_ids": [0], "xs": absorber}, {"block_ids": [1], "xs": fissile}],
    )
    pi, nlke = make_solvers(partial_fission)
    pi_ok = initialize_outcome(pi) == "accepted"
    if pi_ok:
        pi.Execute()
        k = pi.GetEigenvalue()
        pi_ok = math.isfinite(k) and k > 0.0
    passed += on_all_ranks(pi_ok)
    passed += on_all_ranks(initialize_outcome(nlke) == "accepted")

    if rank == 0:
        print(f"KEIGEN_FISSION_REQUIREMENT_PASSED {passed}")
