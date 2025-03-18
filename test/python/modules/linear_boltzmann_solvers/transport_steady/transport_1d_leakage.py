#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 1D Transport leakage test
# Unit angular flux left boundary condition in a pure absorber with unit
# length and a unit absorption cross section. The analytic solution is:
# j^+ = \int_{0}^{1} \mu e^{-1/\mu} d\mu = 0.10969

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
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver

if __name__ == "__main__":
    # Check number of processors
    num_procs = 3
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    N = 100
    L = 1.0
    nodes = []
    for i in range(N + 1):
        nodes.append(i * L / N)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Add materials
    num_groups = 1
    sigma_t = 1.0

    xs1g = MultiGroupXS()
    xs1g.CreateSimpleOneGroup(sigma_t, 0.0)

    bsrc = [0.0 for _ in range(num_groups)]
    bsrc[0] = 1.0

    # Setup Physics
    pquad = GLProductQuadrature1DSlab(256)
    phys = DiscreteOrdinatesSolver(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs1g},
        ],
        options={
            "boundary_conditions": [
                {
                    "name": "zmin",
                    "type": "isotropic",
                    "group_strength": bsrc,
                },
            ],
            "scattering_order": 0,
            "save_angular_flux": True,
        },
    )

    ss_solver = SteadyStateSolver(lbs_problem=phys)
    # Solve the problem
    ss_solver.Initialize()
    ss_solver.Execute()

    # Compute the leakage
    leakage = phys.ComputeLeakage([])
    for k, v in leakage.items():
        print(f"{k}={v[0]:.5e}")
