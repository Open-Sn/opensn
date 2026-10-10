#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Verify that steady-state restart moments survive a change in MPI partitioning.

The same heterogeneous absorber is solved on an x-striped and a y-striped four-rank KBA
partition.  A restart written by the x partition is loaded on the y partition and compared with
an independently converged y-partition solution.  Tight algebraic tolerances make the remaining
partition-dependent iteration error negligible compared with the acceptance tolerance.
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier

    def global_max(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier

    def global_max(value):
        return MPIAllReduce(value, "max")


def make_grid(partitioner):
    nodes = [i / 8.0 for i in range(9)]
    grid = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes], partitioner=partitioner
    ).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=0.5, ymin=0.0, ymax=1.0, infz=True), 1, True
    )
    return grid


def make_problem(grid, restart_read="", restart_write=""):
    xs_left = MultiGroupXS()
    xs_left.CreateSimpleOneGroup(1.3, 0.2)
    xs_right = MultiGroupXS()
    xs_right.CreateSimpleOneGroup(0.7, 0.4)
    quadrature = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=8, scattering_order=0)
    options = {
        "use_precursors": False,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }
    if restart_read:
        options["read_restart_path"] = restart_read
    if restart_write:
        options.update(
            {
                "restart_writes_enabled": True,
                "write_delayed_psi_to_restart": False,
                "write_restart_path": restart_write,
            }
        )
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 300,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs_right}, {"block_ids": [1], "xs": xs_left}],
        volumetric_sources=[VolumetricSource(block_ids=[1], group_strength=[1.0])],
        boundary_conditions=[
            {"name": "xmin", "type": "vacuum"},
            {"name": "xmax", "type": "vacuum"},
            {"name": "ymin", "type": "vacuum"},
            {"name": "ymax", "type": "vacuum"},
        ],
        options=options,
    )


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return np.array(problem.GetPhiNewLocal(), copy=True)


if __name__ == "__main__":
    if size != 4:
        sys.exit(f"Incorrect number of processors. Expected 4 but got {size}.")

    restart_stem = "restart_repartition/restart"
    x_partition = KBAGraphPartitioner(
        nx=4, ny=1, nz=1, xcuts=[0.25, 0.5, 0.75]
    )
    y_partition = KBAGraphPartitioner(
        nx=1, ny=4, nz=1, ycuts=[0.25, 0.5, 0.75]
    )

    solve(make_problem(make_grid(x_partition), restart_write=restart_stem))
    reference = solve(make_problem(make_grid(y_partition)))

    restarted_problem = make_problem(make_grid(y_partition), restart_read=restart_stem)
    restarted_solver = SteadyStateSourceSolver(problem=restarted_problem)
    restarted_solver.Initialize()
    restarted = np.array(restarted_problem.GetPhiNewLocal(), copy=True)

    local_difference = float(np.max(np.abs(reference - restarted)))
    difference = global_max(local_difference)
    if rank == 0:
        print(f"RESTART_REPARTITION_MAX_ABS_DIFF={difference:.16e}")
        print(f"RESTART_REPARTITION_PASS={int(difference < 1.0e-10)}")

    barrier()
    restart_file = f"{restart_stem}{rank}.restart.h5"
    if os.path.exists(restart_file):
        os.remove(restart_file)
    barrier()
    if rank == 0:
        os.rmdir(os.path.dirname(restart_stem))
