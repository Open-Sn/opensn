#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2G KEigenvalue::Solver test using NonLinearK
Test: Final k-eigenvalue: 0.5969127
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver

if __name__ == "__main__":

    with open("utils/qblock_mesh.py") as f:
        exec(f.read(), globals())
    with open("utils/qblock_materials.py") as f:
        exec(f.read(), globals())

    # Setup Physics
    pquad = GLCProductQuadrature2DXY(n_polar=8, n_azimuthal=16, scattering_order=2)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_max_its": 50,
                "gmres_restart_interval": 50,
                "l_abs_tol": 1.0e-10,
            },
        ],
        xs_map=xs_map,
        scattering_order=2,
        options={
            "boundary_conditions": [
                {"name": "xmin", "type": "reflecting"},
                {"name": "ymin", "type": "reflecting"},
            ],

            "use_precursors": False,

            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
        },
    )
    k_solver = NonLinearKEigenSolver(do_problem=phys)
    k_solver.Initialize()
    k_solver.Execute()
