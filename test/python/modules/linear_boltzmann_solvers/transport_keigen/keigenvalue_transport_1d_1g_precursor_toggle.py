#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D 1G delayed k-eigenvalue test.

Steady-state k should include delayed neutron production from the XS data
whether or not precursor storage is enabled.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver


def run_case(use_precursors):
    length = 100.0
    num_cells = 50
    dx = length / num_cells
    nodes = [i * dx for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/simple_fissile.xs")

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(
                    n_polar=32,
                    scattering_order=0,
                ),
                "inner_linear_method": "petsc_gmres",
                "l_max_its": 200,
                "l_abs_tol": 1.0e-10,
            }
        ],
        xs_map=[
            {
                "block_ids": [0],
                "xs": xs,
            }
        ],
        options={
            "use_precursors": use_precursors,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )

    solver = PowerIterationKEigenSolver(
        problem=problem,
        max_iters=2000,
        k_tol=1.0e-8,
    )
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue()


if __name__ == "__main__":
    expected_num_procs = 4
    if size != expected_num_procs:
        sys.exit(
            f"Incorrect number of processors. Expected {expected_num_procs} processors "
            f"but got {size}."
        )

    k_without_precursors = run_case(False)
    k_with_precursors = run_case(True)
    k_difference = abs(k_with_precursors - k_without_precursors)

    if rank == 0:
        print(f"Python use_precursors=False k-eigenvalue: {k_without_precursors}")
        print(f"Python use_precursors=True k-eigenvalue: {k_with_precursors}")
        print(f"Python precursor toggle k-difference: {k_difference}")

    if k_difference > 1.0e-10:
        raise RuntimeError(
            "Expected steady-state delayed k-eigenvalue to be independent of "
            f"use_precursors, difference is {k_difference}."
        )
