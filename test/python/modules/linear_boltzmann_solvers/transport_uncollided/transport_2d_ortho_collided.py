#!/usr/bin/env python3
"""Parallel collided solve using a serially generated uncollided flux file."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.xs import MultiGroupXS


if __name__ == "__main__":
    nodes = [float(i - 10) * 5.0 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.1, c=0.5)

    quadrature = GLCProductQuadrature2DXY(
        n_polar=2,
        n_azimuthal=8,
        scattering_order=1,
    )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs},
        ],
        uncollided_flux=os.environ.get("OPENSN_UNCOLLIDED_FILE", "uncollided_flux.h5"),
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    first_solution = list(problem.GetPhiNewLocal())
    solver.Execute()
    second_solution = list(problem.GetPhiNewLocal())
    max_repeat_difference = max(
        abs(first - second) for first, second in zip(first_solution, second_solution)
    )
    balance = solver.ComputeBalanceTable()
    if rank == 0:
        print(f"CombinedBalance={balance['balance']:.6e}")
        print(f"RepeatedSolveMaxDifference={max_repeat_difference:.6e}")
    if abs(balance["balance"]) > 1.0e-5:
        raise RuntimeError(f"Combined balance failed: {balance['balance']}")
    if max_repeat_difference > 1.0e-5:
        raise RuntimeError(
            f"Repeated uncollided solve changed the solution: {max_repeat_difference}"
        )
