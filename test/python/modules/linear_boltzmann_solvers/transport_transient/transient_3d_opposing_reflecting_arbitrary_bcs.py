#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D opposing reflecting and time-dependent arbitrary boundary-condition test.

Uses a small orthogonal mesh with:
  xmin/xmax: reflecting
  ymin/ymax: time-dependent arbitrary
  zmin/zmax: vacuum
"""

import os
import sys

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.math import AngularFluxTimeFunction
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    nodes = [i / 4 for i in range(5)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0, 1.0)

    quad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    def ymin_inflow(group, angle, time):
        omega = quad.omegas[angle]
        return 0.30 if omega.y > 0.0 and time <= 0.1 else 0.0

    def ymax_inflow(group, angle, time):
        omega = quad.omegas[angle]
        return 0.10 if omega.y < 0.0 and time <= 0.1 else 0.0

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        time_dependent=True,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 300,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {
                "name": "ymin",
                "type": "arbitrary",
                "time_function": AngularFluxTimeFunction(ymin_inflow),
            },
            {
                "name": "ymax",
                "type": "arbitrary",
                "time_function": AngularFluxTimeFunction(ymax_inflow),
            },
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )

    solver = TransientSolver(problem=problem, dt=0.1, theta=1.0, initial_state="zero")
    solver.Initialize()
    solver.SetTimeStep(0.1)
    solver.Advance()
    solver.Advance()

    phi = list(problem.GetPhiNewLocal())
    max_phi = comm.allreduce(max([abs(v) for v in phi] + [0.0]), op=MPI.MAX)

    if rank == 0:
        print(f"TD_OPPOSING_REFLECTING_ARBITRARY_3D_MAX_PHI {max_phi:.12e}")
        print("TD_OPPOSING_REFLECTING_ARBITRARY_3D_PASS 1")
