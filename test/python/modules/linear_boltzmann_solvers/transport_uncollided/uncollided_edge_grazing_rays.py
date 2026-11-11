#!/usr/bin/env python3
"""Uncollided ray tracing through rays that pass exactly through mesh edges.

The source sits at the centroid of a regular tetrahedron whose vertices are four corners
of a cube centered on the source. Gmsh meshes the surrounding cells of the 2 x 2 x 2 box
nearly mirror symmetric, so rays from some quadrature points to the source pass within
roundoff of an edge. Before the fix, such a ray was handed to a cell that it only touches
and the ray tracer stopped with "failed to advance". The mesh,
uncollided_edge_graze_cube.msh, was made with Gmsh 4.15.2 from an OpenCASCADE box with
those four points embedded, a maximum element size of 0.5, and an element size of 0.1 at
the points.

The test checks that the near-source ray tracing, with and without the exact near-source
leakage, completes, and that the outflow agrees with the exact escape fraction of the box.
"""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
remove_file = uncollided_utils.remove_file
uncollided_escape_3d_box = uncollided_utils.uncollided_escape_3d_box
volume_integral = uncollided_utils.volume_integral


def outflow_error(exact_near_source_leakage):
    grid = FromFileMeshGenerator(filename=mesh_path("uncollided_edge_graze_cube.msh")).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 0.8
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    # The whole box is ray traced, so every cell sends rays to the source.
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        zmin=-1.01,
        zmax=1.01,
    )
    file_name = f"uncollided_edge_grazing_rays_{int(exact_near_source_leakage)}.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=[0.0, 0.0, 0.0], strength=[1.0])],
        near_source=[whole_domain],
        exact_near_source_leakage=exact_near_source_leakage,
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name, progress_interval=0)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    integral = volume_integral(scalar_flux, whole_domain, FieldFunctionInterpolationVolume)
    computed_outflow = 1.0 - sigma_t * integral
    reference_outflow = uncollided_escape_3d_box(sigma_t, (0.0, 0.0, 0.0), [(-1.0, 1.0)] * 3)
    remove_file(file_name)
    return abs(computed_outflow - reference_outflow) / reference_outflow


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    quadrature_error = outflow_error(False)
    exact_error = outflow_error(True)
    if rank == 0:
        print(f"UncollidedEdgeGrazingOutflowError={quadrature_error:.8e}")
        print(f"UncollidedEdgeGrazingExactOutflowError={exact_error:.8e}")
