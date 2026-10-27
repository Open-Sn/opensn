#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CMFD on a large reflected core with a high dominance ratio.

A 2D two-group quarter core (fuel in [0, 200]^2 cm, water out to 220 cm, 2 cm cells, reflecting
xmin/ymin, vacuum xmax/ymax) has a dominance ratio close to one. CMFD with an identity coarse mesh
and one transport sweep per power iteration should converge in a number of power iterations that
does not grow with the core size: the coarse eigenvalue problem resolves the fundamental mode.
That requires the coarse power iteration to converge the coarse flux shape, not only the coarse
eigenvalue; stopping on the eigenvalue alone leaves the fundamental mode unaccelerated, and the
power iteration then needs several hundred iterations here.

Reference: unaccelerated power iteration with converged GMRES within-group solves, k_tol = 1e-12,
on the same problem (run this input with accel="none", 32 ranks): k = 0.9276389082.

Checks:
  CMFD_HD_K_ERR: |k / K_REFERENCE - 1|.
  CMFD_HD_ITERATIONS: power iterations; about 60 when the coarse problem is converged.
  CMFD_HD_SKIPPED: skipped CMFD corrections.
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
    from pyopensn.solver import (DiscreteOrdinatesProblem, PowerIterationKEigenSolver,
                                 CMFDAcceleration)
    from pyopensn.logvol import RPPLogicalVolume

K_REFERENCE = 0.9276389082

if __name__ == "__main__":
    accel = globals().get("accel", "cmfd")
    nodes = [2.0 * i for i in range(111)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=-1.0, xmax=200.0, ymin=-1.0, ymax=200.0, infz=True), 1, True)
    water = MultiGroupXS()
    water.LoadFromOpenSn("../../../../assets/xs/xs_water_g2.xs")
    fuel = MultiGroupXS()
    fuel.LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16,
                                                               scattering_order=0),
                "inner_linear_method": "classic_richardson" if accel == "cmfd" else "petsc_gmres",
                "l_abs_tol": 1.0e-12 if accel == "cmfd" else 1.0e-11,
                "l_max_its": 1 if accel == "cmfd" else 1000,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": water}, {"block_ids": [1], "xs": fuel}],
        boundary_conditions=[{"name": "xmin", "type": "reflecting"},
                             {"name": "ymin", "type": "reflecting"}],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )

    if accel == "cmfd":
        cmfd = CMFDAcceleration(problem=problem, coarse_mesh="identity",
                                balance_residual_tolerance=1.0e-7)
        solver = PowerIterationKEigenSolver(problem=problem, acceleration=cmfd, k_tol=1.0e-8,
                                            max_iters=1000)
    else:
        solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=20000)
    solver.Initialize()
    solver.Execute()

    k = solver.GetEigenvalue()
    if rank == 0:
        print(f"CMFD_HD_K {k:.10f}")
        if accel == "cmfd":
            print(f"CMFD_HD_K_ERR {abs(k / K_REFERENCE - 1.0):.6e}")
            print(f"CMFD_HD_ITERATIONS {solver.GetNumPowerIterations()}")
            print(f"CMFD_HD_SKIPPED {cmfd.GetNumSkippedCorrections()}")
