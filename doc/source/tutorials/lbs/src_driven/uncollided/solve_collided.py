#!/usr/bin/env python3
"""Solve the collided component for a simple 3D first-collision tutorial."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../../"))
    )
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
tutorial_common = importlib.import_module("uncollided_common")


def point_value(field_function, point):
    interpolation = FieldFunctionInterpolationPoint()
    interpolation.SetPointOfInterest(Vector3(*point))
    interpolation.AddFieldFunction(field_function)
    interpolation.Execute()
    return interpolation.GetPointValue()


def volume_integral(field_function, logical_volume):
    interpolation = FieldFunctionInterpolationVolume()
    interpolation.SetOperationType("sum")
    interpolation.SetLogicalVolume(logical_volume)
    interpolation.AddFieldFunction(field_function)
    interpolation.Execute()
    return interpolation.GetValue()


if __name__ == "__main__":
    if not tutorial_common.UNCOLLIDED_FILE.exists():
        raise RuntimeError(
            f"Missing {tutorial_common.UNCOLLIDED_FILE}. Run generate_uncollided.py first."
        )

    grid = tutorial_common.make_mesh(FromFileMeshGenerator)
    xs = tutorial_common.make_xs(MultiGroupXS)
    whole_domain = tutorial_common.make_whole_domain(RPPLogicalVolume)

    quadrature = GLCProductQuadrature3DXYZ(
        n_polar=2,
        n_azimuthal=16,
        scattering_order=0,
    )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "gmres_restart_interval": 30,
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=str(tutorial_common.UNCOLLIDED_FILE),
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    point_flux = point_value(scalar_flux, tutorial_common.SAMPLE_POINT)
    total_integral = volume_integral(scalar_flux, whole_domain)
    balance = solver.ComputeBalanceTable()["balance"]

    if rank == 0:
        print(f"TutorialTotalFluxAtPoint={point_flux:.12e}")
        print(f"TutorialTotalFluxIntegral={total_integral:.12e}")
        print(f"TutorialBalance={balance:.12e}")
