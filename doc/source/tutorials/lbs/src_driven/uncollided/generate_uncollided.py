#!/usr/bin/env python3
"""Generate uncollided moments for a simple 3D first-collision tutorial."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../../"))
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
tutorial_common = importlib.import_module("uncollided_common")


if __name__ == "__main__":
    if size != 1:
        raise RuntimeError(
            "This tutorial generates the uncollided file in serial. "
            "Run the collided solve separately with MPI if desired."
        )

    tutorial_common.remove_file(tutorial_common.UNCOLLIDED_FILE)

    grid = tutorial_common.make_mesh(FromFileMeshGenerator)
    xs = tutorial_common.make_xs(MultiGroupXS)
    whole_domain = tutorial_common.make_whole_domain(RPPLogicalVolume)

    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[
            PointSource(location=list(tutorial_common.SOURCE_LOCATION), strength=[1.0])
        ],
        near_source=[whole_domain],
        scattering_order=0,
    )
    solver = UncollidedSolver(
        problem=problem,
        file_name=str(tutorial_common.UNCOLLIDED_FILE),
        progress_interval=25,
    )
    solver.Initialize()
    solver.Execute()

    if rank == 0:
        print(f"TutorialUncollidedFile={tutorial_common.UNCOLLIDED_FILE}")
        print("TutorialUncollidedFileGenerated=1")
