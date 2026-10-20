#!/usr/bin/env python3

"""Splitting AAH sweep messages must not change the solution.

The mesh is a strip of 2 x 1169 cells split along x between two ranks. With polar angle
aggregation each angle set holds the 8 polar angles of one azimuth, so every sweep message
across the interface carries 1169 faces x 2 PWLD nodes x 1 group x 8 angles = 18704 doubles.
max_mpi_message_size = 1025 bytes (any value of at least 1024 is accepted) splits each such
message into blocks of at most 128 doubles; the default of 32768 bytes into blocks of at most
4096 doubles.

Claim: both solutions agree to roundoff, since the blocks carry the same doubles (only the
order in which angle sets complete may change). The pure absorber is solved exactly by one
sweep, so the particle balance (production minus absorption and outflow) of the split run is
a roundoff-level conservation check.

Checks:
  AAHMessageSplitRelativeDifference: largest of the relative differences of the scalar-flux
    volume sum and maximum between the two message sizes.
  AAHMessageSplitRelativeBalance: |balance| / production of the split run.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def volume_value(field_function, logical_volume, operation):
    interpolation = FieldFunctionInterpolationVolume()
    interpolation.SetOperationType(operation)
    interpolation.SetLogicalVolume(logical_volume)
    interpolation.AddFieldFunction(field_function)
    interpolation.Execute()
    return interpolation.GetValue()


def solve(grid, xs, whole_domain, max_mpi_message_size):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=16, n_azimuthal=4, scattering_order=0
                ),
                "angle_aggregation_type": "polar",
                "inner_linear_method": "classic_richardson",
                "l_max_its": 2,
                "l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(logical_volume=whole_domain, group_strength=[1.0])],
        options={"max_mpi_message_size": max_mpi_message_size},
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    phi = problem.GetScalarFluxFieldFunction()[0]
    return (
        volume_value(phi, whole_domain, "sum"),
        volume_value(phi, whole_domain, "max"),
        solver.ComputeBalanceTable(),
    )


if __name__ == "__main__":
    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    num_cells_y = 1169
    xnodes = [0.0, 0.5, 1.0]
    ynodes = [0.01 * j for j in range(num_cells_y + 1)]
    grid = OrthogonalMeshGenerator(
        node_sets=[xnodes, ynodes],
        partitioner=KBAGraphPartitioner(nx=2, ny=1, nz=1, xcuts=[0.5]),
    ).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)
    whole_domain = RPPLogicalVolume(infx=True, infy=True, infz=True)

    sum_default, max_default, _ = solve(grid, xs, whole_domain, 32768)
    sum_split, max_split, balance = solve(grid, xs, whole_domain, 1025)

    difference = max(
        abs(sum_split - sum_default) / abs(sum_default),
        abs(max_split - max_default) / abs(max_default),
    )
    relative_balance = abs(balance["balance"]) / balance["production_rate"]
    if rank == 0:
        print(f"AAHMessageSplitRelativeDifference={difference:.8e}")
        print(f"AAHMessageSplitRelativeBalance={relative_balance:.8e}")
        sys.stdout.flush()
