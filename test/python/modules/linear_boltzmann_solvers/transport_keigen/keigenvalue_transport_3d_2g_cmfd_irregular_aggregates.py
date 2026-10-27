#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CMFD on irregular coarse cells converges to a tight tolerance.

A 3D two-group fuel/water block (8^3 cells, reflecting xmin/ymin/zmin, vacuum elsewhere) is solved
with CMFD on rank-local aggregates of about 8 cells. Breadth-first aggregation on a partitioned
mesh gives irregular coarse cells, some with faces whose centroid lies to the side of the cell
centroid: the vector between them is nearly parallel to the face. The two-point diffusion coupling
of such a face must stay bounded (it is D A / |d| per side); a coupling that grows without
bound as d.n vanishes makes the CMFD iteration settle into a noise floor of about 1e-6 in k, so it
never meets k_tol = 1e-8.

Reference: unaccelerated power iteration on the same problem with converged GMRES within-group
solves and k_tol = 1e-10.

Checks:
  CMFD_IRREG_CONVERGED: 1 if the CMFD power iteration converged.
  CMFD_IRREG_K_ERR: |k_cmfd / k_reference - 1|.
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
    from pyopensn.solver import (DiscreteOrdinatesProblem, PowerIterationKEigenSolver,
                                 CMFDAcceleration)
    from pyopensn.logvol import RPPLogicalVolume


def solve(use_cmfd):
    nodes = [14.0 * i / 8 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=-1000.0, xmax=10.0, ymin=-1000.0, ymax=10.0, zmin=-1000.0,
                         zmax=10.0), 1, True)
    grid.SetOrthogonalBoundaries()
    water = MultiGroupXS()
    water.LoadFromOpenSn("../../../../assets/xs/xs_water_g2.xs")
    fuel = MultiGroupXS()
    fuel.LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": [0, 1],
                "angular_quadrature": GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8,
                                                                scattering_order=1),
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": water}, {"block_ids": [1], "xs": fuel}],
        boundary_conditions=[{"name": name, "type": "reflecting"}
                             for name in ("xmin", "ymin", "zmin")],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
        sweep_type=globals().get("sweep_type", "AAH"),
    )
    if use_cmfd:
        cmfd = CMFDAcceleration(problem=problem, aggregation_size=8,
                                balance_residual_tolerance=1.0e-7)
        solver = PowerIterationKEigenSolver(problem=problem, acceleration=cmfd, k_tol=1.0e-8,
                                            max_iters=400)
    else:
        solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-10, max_iters=5000)
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue(), solver.IsConverged()


if __name__ == "__main__":
    k_reference, _ = solve(False)
    k_cmfd, converged = solve(True)
    if rank == 0:
        print(f"CMFD_IRREG_K_REFERENCE {k_reference:.10f}")
        print(f"CMFD_IRREG_CONVERGED {int(converged)}")
        print(f"CMFD_IRREG_K_ERR {abs(k_cmfd / k_reference - 1.0):.6e}")
