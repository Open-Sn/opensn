#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D PWLD transport test with reflecting boundaries and multiple groupsets (CBC)
Test: Max-difference=0.0
"""

import os
import sys
import numpy as np
from mpi4py import MPI

if "opensn_console" not in globals():
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

if __name__ == "__main__":
    if size != 4:
        sys.exit(f"Incorrect number of processors. Expected 4 processors but got {size}.")

    nodes = [float(i) for i in range(5)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/diag_XS_64g_1mom_c0.99.xs")
    quadrature = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=0)
    solutions = []
    for use_gpus, intervals, save_angular_flux in (
        (False, ((0, 31), (32, 63)), False),
        (True, ((0, 63),), False),
        (True, ((0, 31), (32, 63)), False),
        (True, ((0, 31), (32, 63)), True),
    ):
        problem = DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=64,
            sweep_type="CBC",
            use_gpus=use_gpus,
            groupsets=[{
                "groups_from_to": interval,
                "angular_quadrature": quadrature,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 3,
            } for interval in intervals],
            xs_map=[{"block_ids": [0], "xs": xs}],
            boundary_conditions=[
                {"name": "xmin", "type": "isotropic", "group_strength": [1.0] * 64},
                {"name": "xmax", "type": "reflecting"},
                {"name": "ymax", "type": "reflecting"},
            ],
            options={"save_angular_flux": save_angular_flux},
        )
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()
        solver.Execute()
        solutions.append(np.array(problem.GetPhiOldLocal()))

    errors = [np.max(np.abs(phi - solutions[0])) for phi in solutions[1:]]
    for phi in solutions:
        group_fluxes = phi.reshape(-1, 64)
        errors.append(np.max(np.abs(group_fluxes - group_fluxes[:, :1])))
    difference = float(np.max(errors)) if np.all(np.isfinite(errors)) else float("inf")
    difference = MPI.COMM_WORLD.allreduce(difference, op=MPI.MAX)
    if rank == 0:
        print(f"Max-difference={difference:.12e}")
