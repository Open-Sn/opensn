#!/usr/bin/env python3
"""Verify that the uncollided generator accepts only point sources."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
mesh_path = importlib.import_module("uncollided_unstructured_utils").mesh_path


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)
    rejected = False
    try:
        UncollidedProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{"groups_from_to": [0, 0]}],
            xs_map=[{"block_ids": [0], "xs": xs}],
            volumetric_sources=[
                VolumetricSource(block_ids=[0], group_strength=[1.0])
            ],
        )
    except Exception as error:
        rejected = "only point sources are supported" in str(error)

    if rank == 0:
        print(f"UncollidedVolumetricSourceRejected={int(rejected)}")
        sys.stdout.flush()

    if not rejected:
        raise RuntimeError("Uncollided generation accepted a volumetric source")
