#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D 1G keigenvalue test using power iteration
Test: Final k-eigenvalue: 0.9995433
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
    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Mesh variables
    L = 100.0      # Domain length
    n_cells = 50   # Number of cells
    dx = L / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Load cross-section data for a 1-group fissile material
    num_groups = 1
    xs_simple_fissile = MultiGroupXS()
    xs_simple_fissile.LoadFromOpenSn("simple_fissile.xs")

    # Angle information
    n_angles = 32  # Number of discrete angles
    scat_order = 0  # Scattering order

    # k-eigenvalue iteration parameters
    kes_max_iterations = 5000
    kes_tolerance = 1e-8

    # Source iteration parameters
    si_max_iterations = 500
    si_tolerance = 1e-8

    # Delayed neutrons
    use_precursors = True

    # Create and configure the discrete ordinates solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=n_angles),
                "inner_linear_method": "petsc_gmres",
                "l_max_its": si_max_iterations,
                "l_abs_tol": si_tolerance,
            }
        ],
        xs_map=[
            {
                "block_ids": [0],
                "xs": xs_simple_fissile,
            }
        ],
        options={
            "scattering_order": scat_order,
            "use_precursors": use_precursors,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
        }
    )
    k_solver = NonLinearKEigenSolver(
        lbs_problem=phys,
        nl_max_its=kes_max_iterations,
        nl_abs_tol=kes_tolerance,
    )
    k_solver.Initialize()
    k_solver.Execute()
