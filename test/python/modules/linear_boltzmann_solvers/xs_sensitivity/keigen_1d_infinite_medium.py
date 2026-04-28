#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Analytic 1D reflective infinite-medium k-eigenvalue check."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
else:
    rank = 0


if __name__ == "__main__":
    test_dir = os.path.dirname(__file__)
    sigma_t = 1.0
    sigma_s = 0.2
    production = 0.9
    analytic_k = production / (sigma_t - sigma_s)

    n_cells = 20
    nodes = [i / n_cells for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(test_dir, "xs_sensitivity_fissile_1g.xs"))

    quad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=0)
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-11,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )
    solver = PowerIterationKEigenSolver(problem=phys, max_iters=200, k_tol=1.0e-11)
    solver.Initialize()
    solver.Execute()
    k_eff = solver.GetEigenvalue()
    err = abs(k_eff - analytic_k)

    if rank == 0:
        print(f"XS_SENS_KEIGEN_K_ERR={err:.12e}")
