#!/usr/bin/env python3
"""Generate an unstructured uncollided field for MPI equivalence tests."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
remove_file = uncollided_utils.remove_file


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.8, c=0.55)
    # A small region around the point source, not the whole domain: see the
    # comment in uncollided_2d_multigroup_analytic.py. The source itself is
    # centered on its containing mesh cell (rather than an arbitrary
    # coordinate) so it doesn't sit flush against one of that cell's faces;
    # see UncollidedProblem's runtime warning for this.
    near_source_region = RPPLogicalVolume(
        xmin=0.0441723 - 0.08,
        xmax=0.0441723 + 0.08,
        ymin=-0.0330477 - 0.08,
        ymax=-0.0330477 + 0.08,
        infz=True,
    )
    file_name = "uncollided_parallel_equivalence.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=[0.0441723, -0.0330477, 0.0], strength=[1.0])],
        near_source=[near_source_region],
        scattering_order=1,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()
    if rank == 0:
        print("UncollidedParallelFileGenerated=1")
