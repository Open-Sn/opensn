#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2G k-eigenvalue test using power iteration with aggregated CMFD acceleration.

Test: Final k-eigenvalue: 0.5969127
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
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    N = 40
    L = 14.0
    xmin = 0.0
    dx = L / N
    nodes = [xmin + k * dx for k in range(N + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    vol1 = RPPLogicalVolume(
        xmin=-1000.0,
        xmax=10.0,
        ymin=-1000.0,
        ymax=10.0,
        infz=True,
    )
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    xss = {}
    xss["0"] = MultiGroupXS()
    xss["0"].LoadFromOpenSn("../../../../assets/xs/xs_water_g2.xs")
    xss["1"] = MultiGroupXS()
    xss["1"].LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")
    num_groups = xss["0"].num_groups

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
        xs_map=[
            {"block_ids": [0], "xs": xss["0"]},
            {"block_ids": [1], "xs": xss["1"]},
        ],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
        ],
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        sweep_type=globals().get("sweep_type", "AAH"),
    )

    k_tolerance = 1.0e-8

    cmfd = CMFDAcceleration(
        problem=phys,
        aggregation_size=8,
        allow_mixed_material_aggregation=globals().get(
            "allow_mixed_material_aggregation", False
        ),
        relaxation=1.0,
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
    # The eigenvector scale is arbitrary, so measure the balance residual relative to the
    # production rate.
    balance_table = k_solver.ComputeBalanceTable()
    relative_balance = balance_table["balance"] / balance_table["production_rate"]
    if rank == 0:
        print(f"Python k-eigenvalue: {k}")
        print(f"Python converged: {int(k_solver.IsConverged())}")
        print(f"Python CMFD skipped corrections: {cmfd.GetNumSkippedCorrections()}")
        print(f"Python sweeps: {sweeps}")
        print(f"Python relative balance: {relative_balance:.16e}")

    expected_k = 0.5969127
    if abs(k - expected_k) > 2.0e-6:
        raise RuntimeError(f"Expected k={expected_k}, got {k}")
    if abs(relative_balance) > 10.0 * k_tolerance:
        raise RuntimeError(f"Expected particle balance closure, got {relative_balance}")
