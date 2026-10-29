#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D 1G delayed k-eigenvalue with and without delayed neutrons.

use_precursors (alias delayed_neutrons) selects the fission model: with False, a k-eigenvalue solve
uses prompt neutrons only and gives the prompt eigenvalue k_p. The one-group data have
nu_p = 1.99 and nu_d = 0.01 with a single spectrum, so F_p = (nu_p / nu) F and k_p / k = 0.995
exactly.

Checks:
  k with use_precursors=True and False (0.9995433 and 0.9945456).
  Python k_p/k ratio error: |k_p / k - 0.995|.
  Python delayed_neutrons alias k-difference: |k(delayed_neutrons=False) - k(use_precursors=False)|.
  Python default options k-difference: |k(empty options) - k(use_precursors=True)|.
  Python both delayed-neutron options rejected: 1 when giving both names raises ValueError.
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


def run_case(delayed_options):
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
        # None gives an empty options block, all defaults
        options={}
        if delayed_options is None
        else {
            **delayed_options,
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

    k_prompt = run_case({"use_precursors": False})
    k_total = run_case({"use_precursors": True})
    k_alias = run_case({"delayed_neutrons": False})
    k_default = run_case(None)
    try:
        run_case({"use_precursors": True, "delayed_neutrons": True})
        both_rejected = 0
    except ValueError:
        both_rejected = 1

    if rank == 0:
        print(f"Python use_precursors=False k-eigenvalue: {k_prompt}")
        print(f"Python use_precursors=True k-eigenvalue: {k_total}")
        print(f"Python k_p/k ratio error: {abs(k_prompt / k_total - 0.995)}")
        print(f"Python delayed_neutrons alias k-difference: {abs(k_alias - k_prompt)}")
        print(f"Python default options k-difference: {abs(k_default - k_total)}")
        print(f"Python both delayed-neutron options rejected: {both_rejected}")
