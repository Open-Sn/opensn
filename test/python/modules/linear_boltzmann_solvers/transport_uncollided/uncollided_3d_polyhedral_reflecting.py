#!/usr/bin/env python3
"""Reflected uncollided transport on a true polyhedral pMesh.

The uniform-medium solution is the sum of the physical point source and its
image across x=0. This exercises reflected traversal through cells with up to
14 faces and polygonal faces with up to eight sides. The reflecting patch has
boundary ID zero, so an internal neighbor with global cell ID zero cannot be
mistaken for a boundary hit.
"""

import importlib
import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
utils = importlib.import_module("uncollided_unstructured_utils")


def analytic_flux(source, point, sigma_t):
    image = (-source[0], source[1], source[2])
    value = 0.0
    for location in (source, image):
        radius = math.dist(location, point)
        value += math.exp(-sigma_t * radius) / (4.0 * math.pi * radius**2)
    return value


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    mesh_path = utils.mesh_path("kobayashi/kobayashi_dog_leg_pmesh")
    grid = FromFileMeshGenerator(filename=mesh_path).Execute()
    grid.SetOrthogonalBoundaries()

    volumes = grid.ComputeVolumePerBlockID()
    volume_error = max(
        abs(volumes[block_id] - reference) / reference
        for block_id, reference in ((1, 1000.0), (2, 15000.0), (3, 344000.0))
    )

    sigma_t = 0.1
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)
    source = tuple(5.0 - 5.0 / math.sqrt(3.0) for _ in range(3))
    output_file = "uncollided_3d_polyhedral_reflecting.h5"
    utils.remove_file(output_file)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [1, 2, 3], "xs": xs}],
        point_sources=[PointSource(location=list(source), strength=[1.0])],
        near_source=[
            RPPLogicalVolume(
                xmin=0.0,
                xmax=10.0,
                ymin=0.0,
                ymax=10.0,
                zmin=0.0,
                zmax=10.0,
            )
        ],
        boundary_conditions=[{"name": "xmin", "type": "reflecting"}],
        exact_near_source_leakage=True,
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=output_file, progress_interval=0)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    points = [
        (17.3, 16.7, 12.9),
        (45.1, 44.2, 21.3),
        (12.4, 80.3, 42.1),
    ]
    errors = []
    for point in points:
        value = utils.point_value(
            scalar_flux,
            point,
            FieldFunctionInterpolationPoint,
            Vector3,
        )
        errors.append(utils.relative_error(value, analytic_flux(source, point, sigma_t)))

    utils.remove_file(output_file)
    if rank == 0:
        print(f"UncollidedPolyhedralVolumeError={volume_error:.8e}")
        print(f"UncollidedPolyhedralReflectingMaxError={max(errors):.8e}")
