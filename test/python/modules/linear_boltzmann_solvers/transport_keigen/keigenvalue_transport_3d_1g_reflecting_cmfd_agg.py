#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 1G reflecting-boundary k-eigenvalue test using aggregated CMFD acceleration.

Test: Final k-eigenvalue: 1.0
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
    from pyopensn.solver import CMFDAcceleration

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    n = 4
    length = 4.0
    nodes = [length * i / n for i in range(n + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/simple_fissile.xs")

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
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        sweep_type=globals().get("sweep_type", "AAH"),
    )

    k_tolerance = 1.0e-6

    cmfd = CMFDAcceleration(
        problem=phys,
        aggregation_size=4,
        relaxation=1.0,
        l_abs_tol=1.0e-10,
        balance_residual_tolerance=10.0 * k_tolerance,
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        acceleration=cmfd,
        max_iters=400,
        k_tol=k_tolerance,
        compute_balance=True,
    )
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    sweeps = k_solver.GetNumSweeps()
    balance = k_solver.ComputeBalanceTable()
    if rank == 0:
        print(f"Python k-eigenvalue: {k}")
        print(f"Python sweeps: {sweeps}")
        print(f"Balance={balance['balance']:.6e}")

    if abs(k - 1.0) > 1.0e-4:
        raise RuntimeError(f"Expected k near 1.0, got {k}")
    if abs(balance["balance"]) > 5.0e-4:
        raise RuntimeError(f"Expected balance closure, got {balance['balance']}")
