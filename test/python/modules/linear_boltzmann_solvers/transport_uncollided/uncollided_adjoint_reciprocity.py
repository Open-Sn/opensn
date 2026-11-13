#!/usr/bin/env python3
"""Verify the adjoint uncollided flux and its first-collision source against reciprocity.

A fast point source at A and a thermal point detector at B lie 2 cm apart near the center of an
8 cm cube, more than 3 mean free paths from every boundary. In the two-group cross sections the
only scattering is 0 -> 1 downscatter (sigma_s = 0.5, sigma_t = 1.0 and 0.8), so the forward
thermal flux at B is the single-scatter integral of the two point kernels. Its infinite-medium
value is evaluated independently in prolate spheroidal coordinates
(uncollided_unstructured_utils.single_scatter_point_response). The vacuum boundaries change it by
less than 1e-4: moving them 2 cm farther out changed both solutions by 7e-5.

The forward problem uses a forward uncollided flux from A and reads the thermal flux at B. The
adjoint problem uses an adjoint uncollided flux from a thermal detector at B, whose first-collision
source must use the transposed scattering matrix to reach the fast group, and reads the fast
adjoint flux at A. With the forward matrix the adjoint first-collision source is zero and the
adjoint flux at A vanishes.

The tolerances come from the angular discretization of the collided flux, which dominates. With
16 x 32 directions the forward and adjoint errors are 0.5% and 0.8%. With 12 x 24 they are 1.6%
and 0.4%, and with 8 x 16 they grow to 10% and 7%. Refining the mesh at fixed quadrature
increases them, as ray effects sharpen.

Checks:
  UncollidedReciprocityForwardError, UncollidedReciprocityAdjointError: relative errors with
    respect to the reference (below 0.03).
  UncollidedReciprocityDifference: |adjoint / forward - 1| (below 0.01).
"""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
        UncollidedSolver,
    )
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file
single_scatter_point_response = uncollided_utils.single_scatter_point_response

SOURCE = (3.03, 4.11, 3.93)
DETECTOR = (5.03, 4.11, 3.93)


def near_source(point):
    return RPPLogicalVolume(
        xmin=point[0] - 0.6,
        xmax=point[0] + 0.6,
        ymin=point[1] - 0.6,
        ymax=point[1] + 0.6,
        zmin=point[2] - 0.6,
        zmax=point[2] + 0.6,
    )


def solve(grid, xs, adjoint):
    # The forward source emits fast neutrons at A. The adjoint source is the thermal detector at B.
    location, strength = (DETECTOR, [0.0, 1.0]) if adjoint else (SOURCE, [1.0, 0.0])
    file_name = f"uncollided_adjoint_reciprocity_{'adjoint' if adjoint else 'forward'}.h5"
    remove_file(file_name)

    generator = UncollidedProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[{"groups_from_to": [0, 1]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(location), strength=strength)],
        near_source=[near_source(location)],
        scattering_order=0,
        options={"adjoint": adjoint},
    )
    generation = UncollidedSolver(problem=generator, file_name=file_name)
    generation.Initialize()
    generation.Execute()

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": [0, 1],
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=16, n_azimuthal=32, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 50,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=file_name,
        options={"adjoint": adjoint},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    remove_file(file_name)

    # Forward: thermal flux at the detector. Adjoint: fast importance at the source.
    probe, group = (SOURCE, 0) if adjoint else (DETECTOR, 1)
    scalar_flux = problem.GetScalarFluxFieldFunction()[group]
    return point_value(scalar_flux, probe, FieldFunctionInterpolationPoint, Vector3)


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    nodes = [0.5 * i for i in range(17)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "../../../../assets/xs/simple_2g_downscatter_td.cxs",
            )
        )
    )

    forward = solve(grid, xs, adjoint=False)
    adjoint = solve(grid, xs, adjoint=True)
    distance = sum((a - b) ** 2 for a, b in zip(SOURCE, DETECTOR)) ** 0.5
    reference = single_scatter_point_response(0.5, 1.0, 0.8, distance)

    forward_error = relative_error(forward, reference)
    adjoint_error = relative_error(adjoint, reference)
    difference = relative_error(adjoint, forward)

    if rank == 0:
        print(f"UncollidedReciprocityReference={reference:.8e}")
        print(f"UncollidedReciprocityForward={forward:.8e}")
        print(f"UncollidedReciprocityAdjoint={adjoint:.8e}")
        print(f"UncollidedReciprocityForwardError={forward_error:.8e}")
        print(f"UncollidedReciprocityAdjointError={adjoint_error:.8e}")
        print(f"UncollidedReciprocityDifference={difference:.8e}")
        sys.stdout.flush()
