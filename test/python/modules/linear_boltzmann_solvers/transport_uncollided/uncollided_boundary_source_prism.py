#!/usr/bin/env python3
"""Point source on a vacuum boundary face of a prism mesh.

The source sits on the top face z = H of a homogeneous box of extruded triangles. Half of its
strength leaves the domain directly; the other half is emitted into the domain, where the
uncollided flux is S exp(-sigma_t r) / (4 pi r^2). The whole mesh is ray traced.

The boundary plane z = H contains the source, so the faces in it are parallel to the flux and
carry no current. H = 0.7 is not exactly representable, and the centroid of a triangular face in
that plane, the average of three equal coordinates, differs from H in the last bit. The test
checks that this round-off does not turn those faces into inflow or outflow faces:

- the outflow matches the escape probability of the box, computed independently by integrating
  exp(-sigma_t R) over the inward directions, plus the half that leaves directly;
- the flux matches the point kernel at sample points, including next to the boundary plane.
"""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
extruded_mesh = uncollided_utils.extruded_mesh
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file
uncollided_3d = uncollided_utils.uncollided_3d
uncollided_escape_3d_box = uncollided_utils.uncollided_escape_3d_box
volume_integral = uncollided_utils.volume_integral


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    height = 0.7
    grid = extruded_mesh(
        ExtruderMeshGenerator,
        FromFileMeshGenerator,
        "triangle_mesh2x2.obj",
        height=height,
        num_layers=7,
    )
    grid.SetUniformBlockID(0)

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    source = (0.0441723, -0.0330477, height)
    strength = 1.0
    whole_domain = RPPLogicalVolume(
        xmin=-1.01, xmax=1.01, ymin=-1.01, ymax=1.01, zmin=-0.01, zmax=height + 0.01
    )

    file_name = "uncollided_boundary_source_prism.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(source), strength=[strength])],
        near_source=[whole_domain],
        exact_near_source_leakage=True,
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    scalar_integral = volume_integral(scalar_flux, whole_domain, FieldFunctionInterpolationVolume)
    computed_outflow = strength - sigma_t * scalar_integral
    # The escape integral skips the face that contains the source, so it covers the inward half
    # of the directions.
    reference_outflow = 0.5 * strength + strength * uncollided_escape_3d_box(
        sigma_t, source, [(-1.0, 1.0), (-1.0, 1.0), (0.0, height)]
    )
    outflow_error = relative_error(computed_outflow, reference_outflow)

    sample_points = [
        (0.46, 0.11, 0.37),
        (-0.39, 0.43, 0.12),
        (0.22, -0.58, 0.51),
        (0.61, 0.33, 0.67),
        (-0.52, -0.47, 0.68),
    ]
    max_error = 0.0
    for point in sample_points:
        value = point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3)
        exact = uncollided_3d(strength, sigma_t, source, point)
        max_error = max(max_error, relative_error(value, exact))

    if rank == 0:
        print(f"UncollidedBoundarySourcePrismOutflowRelativeError={outflow_error:.8e}")
        print(f"UncollidedBoundarySourcePrismMaxRelativeError={max_error:.8e}")
        sys.stdout.flush()

    remove_file(file_name)
    if outflow_error > 1.0e-3:
        raise RuntimeError(f"Boundary-source outflow error is too large: {outflow_error}")
    if max_error > 0.05:
        raise RuntimeError(f"Boundary-source scalar-flux error is too large: {max_error}")
