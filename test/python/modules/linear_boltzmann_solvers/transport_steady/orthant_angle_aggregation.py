#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Quadrant/octant aggregation on KBA and ParMETIS brick meshes.

1. A uniform one-group pure absorber with a unit isotropic source and
   reflecting boundaries has the exact scalar flux phi = q / sigma_a = 1. The
   2D and 3D KBA problems check quadrant and octant aggregation against it.

2. The infinite-medium fixed point is insensitive to sweep ordering, so a
   heterogeneous problem (scattering background, localized absorbing source
   region, mixed vacuum and reflecting boundaries) is solved on 2D and 3D
   meshes partitioned by the default ParMETIS partitioner. These partitions
   are not planar and produce cyclic inter-partition dependencies. The
   quadrant/octant solutions, with AAH and CBC sweeps, must match the
   single-angle AAH solution on the same discretization to within the
   iterative tolerance, and must conserve particles (production equals
   absorption plus net leakage).

3. Invalid configurations are rejected during setup: octant aggregation on a
   2D mesh, and Lebedev quadratures, whose rules all include the coordinate
   axis directions, with quadrant (2D) and octant (3D) aggregation.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size

    def MPIAllReduce(value, op="max"):
        operations = {"max": MPI.MAX}
        return MPI.COMM_WORLD.allreduce(value, op=operations[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import (
        GLCProductQuadrature2DXY,
        GLCProductQuadrature3DXYZ,
        LebedevQuadrature2DXY,
        LebedevQuadrature3DXYZ,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def solve(grid, dimension, aggregation, xs_map, source_blocks, boundary_conditions,
          sweep_type="AAH", quadrature=None):
    if quadrature is None and dimension == 2:
        quadrature = GLCProductQuadrature2DXY(
            n_polar=4, n_azimuthal=16, scattering_order=0
        )
    elif quadrature is None:
        quadrature = GLCProductQuadrature3DXYZ(
            n_polar=4, n_azimuthal=8, scattering_order=0
        )

    problem = DiscreteOrdinatesProblem(
        name=f"orthant_{dimension}d_{aggregation}_{sweep_type}",
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "angle_aggregation_type": aggregation,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 500,
            }
        ],
        xs_map=xs_map,
        volumetric_sources=[VolumetricSource(block_ids=source_blocks, group_strength=[1.0])],
        boundary_conditions=boundary_conditions,
        sweep_type=sweep_type,
        options={"verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal()), solver


def reflecting(dimension):
    names = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")[: 2 * dimension]
    return [{"name": name, "type": "reflecting"} for name in names]


def max_global_difference(left, right):
    local = max((abs(a - b) for a, b in zip(left, right)), default=0.0)
    return MPIAllReduce(local, "max")


def relative_balance(solver):
    table = solver.ComputeBalanceTable()
    return abs(table["balance"]) / table["production_rate"]


def aggregation_rejected(grid, dimension, aggregation, xs, expected, quadrature=None):
    try:
        solve(grid, dimension, aggregation, [{"block_ids": [0], "xs": xs}], [0],
              reflecting(dimension), quadrature=quadrature)
    except ValueError as error:
        return expected in str(error)
    return False


if __name__ == "__main__":
    if size != 4:
        sys.exit(f"Incorrect number of processors. Expected 4 processors but got {size}.")

    absorber = MultiGroupXS()
    absorber.CreateSimpleOneGroup(1.0, 0.0)
    absorber_map = [{"block_ids": [0], "xs": absorber}]

    # Analytic infinite medium on KBA partitions
    nodes = [0.0, 0.5, 1.0]
    partitioner_2d = KBAGraphPartitioner(
        nx=2, ny=2, nz=1, xcuts=[0.5], ycuts=[0.5]
    )
    grid_2d = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes], partitioner=partitioner_2d
    ).Execute()
    grid_2d.SetUniformBlockID(0)

    partitioner_3d = KBAGraphPartitioner(
        nx=2, ny=2, nz=1, xcuts=[0.5], ycuts=[0.5]
    )
    grid_3d = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes, nodes], partitioner=partitioner_3d
    ).Execute()
    grid_3d.SetUniformBlockID(0)

    quadrant_2d, _ = solve(grid_2d, 2, "quadrant", absorber_map, [0], reflecting(2))
    octant_3d, _ = solve(grid_3d, 3, "octant", absorber_map, [0], reflecting(3))
    exact_2d = max_global_difference(quadrant_2d, [1.0] * len(quadrant_2d))
    exact_3d = max_global_difference(octant_3d, [1.0] * len(octant_3d))

    # Heterogeneous problem on ParMETIS partitions
    background = MultiGroupXS()
    background.CreateSimpleOneGroup(1.0, 0.7)
    source_region = MultiGroupXS()
    source_region.CreateSimpleOneGroup(5.0, 0.2)
    heterogeneous_map = [
        {"block_ids": [0], "xs": background},
        {"block_ids": [1], "xs": source_region},
    ]

    heterogeneous = {}
    for dimension, num_cells, aggregation in ((2, 20, "quadrant"), (3, 8, "octant")):
        fine_nodes = [i / num_cells for i in range(num_cells + 1)]
        grid = OrthogonalMeshGenerator(node_sets=[fine_nodes] * dimension).Execute()
        grid.SetUniformBlockID(0)
        if dimension == 2:
            region = RPPLogicalVolume(xmin=0.2, xmax=0.55, ymin=0.3, ymax=0.7, infz=True)
            boundary_conditions = [
                {"name": "xmax", "type": "reflecting"},
                {"name": "ymax", "type": "reflecting"},
            ]
        else:
            region = RPPLogicalVolume(
                xmin=0.2, xmax=0.55, ymin=0.3, ymax=0.7, zmin=0.1, zmax=0.6
            )
            boundary_conditions = [
                {"name": "xmax", "type": "reflecting"},
                {"name": "zmin", "type": "reflecting"},
            ]
        grid.SetBlockIDFromLogicalVolume(region, 1, True)

        reference, _ = solve(
            grid, dimension, "single", heterogeneous_map, [1], boundary_conditions
        )
        for sweep_type in ("AAH", "CBC"):
            phi, solver = solve(
                grid, dimension, aggregation, heterogeneous_map, [1], boundary_conditions,
                sweep_type,
            )
            heterogeneous[(dimension, sweep_type)] = (
                max_global_difference(reference, phi),
                relative_balance(solver),
            )

    axis_direction = "is zero, so the direction is not inside any one"
    rejections = [
        aggregation_rejected(grid_2d, 2, "octant", absorber, "requires a 3D Cartesian mesh"),
        aggregation_rejected(
            grid_2d, 2, "quadrant", absorber, axis_direction,
            LebedevQuadrature2DXY(quadrature_order=7, scattering_order=0),
        ),
        aggregation_rejected(
            grid_3d, 3, "octant", absorber, axis_direction,
            LebedevQuadrature3DXYZ(quadrature_order=7, scattering_order=0),
        ),
    ]

    if rank == 0:
        print(f"ORTHANT_2D_EXACT_ERROR {exact_2d:.6e}")
        print(f"ORTHANT_3D_EXACT_ERROR {exact_3d:.6e}")
        for (dimension, sweep_type), (difference, balance) in heterogeneous.items():
            prefix = f"ORTHANT_PARMETIS_{dimension}D_{sweep_type}"
            print(f"{prefix}_SINGLE_DIFFERENCE {difference:.6e}")
            print(f"{prefix}_RELATIVE_BALANCE {balance:.6e}")
        print(f"ORTHANT_CONFIGURATION_REJECTIONS {sum(rejections)}")
