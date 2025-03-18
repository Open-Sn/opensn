#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 1D 1G KEigenvalue::Solver test using power iteration
# Test: Final k-eigenvalue: 0.9995433

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.xs import MultiGroupXS
    from pyopensn.solver import DiscreteOrdinatesSolver, NonLinearKEigen

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Mesh variables
    L = 100.0
    n_cells = 50

    # Transport angle information
    n_angles = 32
    scat_order = 0

    # k-eigenvalue iteration parameters
    kes_max_iterations = 5000
    kes_tolerance = 1e-8

    # Source iteration parameters
    si_max_iterations = 500
    si_tolerance = 1e-8

    # Delayed neutrons
    use_precursors = True

    # Setup mesh
    dx = L / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    # Cross sections
    xs_simple_fissile = MultiGroupXS()
    xs_simple_fissile.LoadFromOpenSn("simple_fissile.xs")

    # Setup Physics
    num_groups = 1
    phys = DiscreteOrdinatesSolver(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLProductQuadrature1DSlab(n_angles),
                "inner_linear_method": "petsc_gmres",
                "l_max_its": si_max_iterations,
                "l_abs_tol": si_tolerance,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_simple_fissile},
        ],
        options={
            "scattering_order": scat_order,
            "use_precursors": use_precursors,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
            "save_angular_flux": True,
        },
        sweep_type="CBC",
    )
    k_solver = NonLinearKEigen(
        lbs_problem=phys,
        nl_max_its=kes_max_iterations,
        nl_abs_tol=kes_tolerance,
    )
    k_solver.Initialize()
    k_solver.Execute()
