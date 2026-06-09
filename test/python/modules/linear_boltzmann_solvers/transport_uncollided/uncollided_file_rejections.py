#!/usr/bin/env python3
"""Verify uncollided-file compatibility checks."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
remove_file = uncollided_utils.remove_file


def make_uncollided_file(grid, xs, scattering_order, file_name):
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=[0.13, -0.22, 0.0], strength=[1.0])],
        near_source=[whole_domain],
        scattering_order=scattering_order,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()


def make_collided_problem(grid, xs, scattering_order, file_name):
    quadrature = GLCProductQuadrature2DXY(
        n_polar=2,
        n_azimuthal=12,
        scattering_order=scattering_order,
    )
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=file_name,
    )


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)

    xs_reference = MultiGroupXS()
    xs_reference.CreateSimpleOneGroup(sigma_t=0.7, c=0.0)

    xs_mismatch = MultiGroupXS()
    xs_mismatch.CreateSimpleOneGroup(sigma_t=0.9, c=0.0)

    moment_file = "uncollided_file_reject_moment.h5"
    xs_file = "uncollided_file_reject_xs.h5"

    remove_file(moment_file)
    remove_file(xs_file)
    make_uncollided_file(grid, xs_reference, 0, moment_file)
    make_uncollided_file(grid, xs_reference, 0, xs_file)

    moment_order_rejected = False
    try:
        make_collided_problem(grid, xs_reference, 1, moment_file)
    except Exception as error:
        moment_order_rejected = (
            "moment order is lower than the collided scattering order" in str(error)
        )

    xs_mismatch_rejected = False
    try:
        make_collided_problem(grid, xs_mismatch, 0, xs_file)
    except Exception as error:
        xs_mismatch_rejected = "total cross-section mismatch" in str(error)

    remove_file(moment_file)
    remove_file(xs_file)

    if rank == 0:
        print(f"UncollidedMomentOrderRejected={int(moment_order_rejected)}")
        print(f"UncollidedXSMismatchRejected={int(xs_mismatch_rejected)}")
        sys.stdout.flush()

    if not moment_order_rejected:
        raise RuntimeError("Uncollided file with insufficient moment order was accepted.")
    if not xs_mismatch_rejected:
        raise RuntimeError("Uncollided file with mismatched total cross section was accepted.")
