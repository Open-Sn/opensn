#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verify that Galerkin quadrature operators reject adjoint mode and cross-section sensitivities.

Transposing the cross sections yields the exact discrete adjoint of the scattering operator
M2D Sigma D2M only when G = W M2D^T diag(w) M2D commutes with Sigma. That holds for the standard
operators but not in general for Galerkin operators, so adjoint mode must be rejected for them.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import CrossSectionSensitivityPostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem
    from pyopensn.xs import MultiGroupXS


def make_problem(operator_method, adjoint=False):
    grid = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0], [0.0, 1.0]]).Execute()
    grid.SetUniformBlockID(0)

    quadrature = GLCProductQuadrature2DXY(
        n_polar=2, n_azimuthal=4, scattering_order=1, operator_method=operator_method
    )
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.5)

    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "angle_aggregation_type": "single",
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        options={"adjoint": adjoint},
    )


def is_galerkin_rejection(error):
    return "not supported with Galerkin quadrature operators" in str(error)


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Incorrect number of processors. Expected 1 processor but got {size}.")

    # Standard operators remain valid in adjoint mode.
    make_problem("standard", adjoint=True)
    make_problem("standard").SetAdjoint(True)

    for method in ["galerkin_one", "galerkin_three"]:
        constructor_rejected = False
        try:
            make_problem(method, adjoint=True)
        except ValueError as error:
            constructor_rejected = is_galerkin_rejection(error)

        problem = make_problem(method)
        setter_rejected = False
        try:
            problem.SetAdjoint(True)
        except ValueError as error:
            setter_rejected = is_galerkin_rejection(error)

        sensitivity_rejected = False
        try:
            CrossSectionSensitivityPostprocessor(problem=problem, sensitivity_type="sigma_t")
        except ValueError as error:
            sensitivity_rejected = is_galerkin_rejection(error)

        if not (constructor_rejected and setter_rejected and sensitivity_rejected):
            sys.exit(f"Galerkin operator '{method}' was not rejected through every entry path.")

    if rank == 0:
        print("ADJOINT_GALERKIN_REJECTED=1")
