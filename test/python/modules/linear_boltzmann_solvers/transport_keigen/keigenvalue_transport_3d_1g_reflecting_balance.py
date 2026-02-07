#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 1G infinite-medium k-eigen test with reflecting boundaries.
Analytic k_inf = nu_sigma_f / sigma_a = 1.0 for the chosen cross sections.
Balance table should close with k-normalized production.
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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver


if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Simple 3D orthogonal mesh
    n = 4
    length = 4.0
    nodes = [length * i / n for i in range(n + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    # 1-group fissile material with k_inf ~ 1.0
    xs = MultiGroupXS()
    xs.LoadFromOpenSn("simple_fissile.xs")

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=4, n_azimuthal=8, scattering_order=0
                ),
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs},
        ],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
        },
    )

    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        k_tol=1.0e-6,
        compute_balance=True,
    )
    k_solver.Initialize()
    k_solver.Execute()
