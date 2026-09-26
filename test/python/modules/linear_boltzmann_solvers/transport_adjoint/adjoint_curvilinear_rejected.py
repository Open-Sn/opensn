#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Verify that curvilinear problems reject adjoint mode."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem
    from pyopensn.xs import MultiGroupXS


def make_problem(adjoint=False):
    grid = OrthogonalMeshGenerator(
        node_sets=[[0.0, 1.0], [0.0, 1.0]], coord_sys="cylindrical"
    ).Execute()
    grid.SetUniformBlockID(0)

    quadrature = GLCProductQuadrature2DRZ(
        n_polar=2, n_azimuthal=4, scattering_order=0
    )
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)

    return DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "angle_aggregation_type": "azimuthal",
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        options={"adjoint": adjoint},
    )


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Incorrect number of processors. Expected 1 processor but got {size}.")

    constructor_rejected = False
    try:
        make_problem(adjoint=True)
    except ValueError as error:
        constructor_rejected = "only for Cartesian geometry" in str(error)

    problem = make_problem()
    setter_rejected = False
    try:
        problem.SetAdjoint(True)
    except ValueError as error:
        setter_rejected = "only for Cartesian geometry" in str(error)

    if not constructor_rejected or not setter_rejected:
        sys.exit("Curvilinear adjoint mode was not rejected through every supported entry path.")

    if rank == 0:
        print("ADJOINT_CURVILINEAR_REJECTED=1")
