#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2G KEigenvalue::Solver test using Power Iteration
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
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import SCDSAAcceleration

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
                "inner_linear_method": "petsc_richardson",
                "l_max_its": 2,
                "gmres_restart_interval": 50,
                "l_abs_tol": 1.0e-10,
            },
        ],
        xs_map=xs_map,
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
        ],
        options={
            "verbose_inner_iterations": True,
            "verbose_outer_iterations": True,
            "save_angular_flux": True,
        },
    )
    scdsa = SCDSAAcceleration(
        problem=phys,
        sdm="pwld"
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        acceleration=scdsa,
        k_tol=1.0e-8,
    )
    k_solver.Initialize()
    k_solver.Execute()
