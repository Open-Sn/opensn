#!/usr/bin/env python3
"""First-collision source and total-flux test on an extruded prism mesh."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
        UncollidedSolver,
    )
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
extruded_mesh = uncollided_utils.extruded_mesh
remove_file = uncollided_utils.remove_file
volume_integral = uncollided_utils.volume_integral


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = extruded_mesh(ExtruderMeshGenerator, FromFileMeshGenerator, "triangle_mesh2x2.obj")
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.8, c=0.55)
    source = PointSource(location=[0.0441723, -0.0330477, 1.11], strength=[1.0])
    near_source = RPPLogicalVolume(
        xmin=-0.22,
        xmax=0.22,
        ymin=-0.22,
        ymax=0.22,
        zmin=0.89,
        zmax=1.33,
    )
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        zmin=-0.01,
        zmax=2.01,
    )

    file_name = "uncollided_collided_3d_prism.h5"
    remove_file(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[source],
        near_source=[near_source],
        scattering_order=1,
    )
    uncollided_solver = UncollidedSolver(problem=uncollided, file_name=file_name)
    uncollided_solver.Initialize()
    uncollided_solver.Execute()
    uncollided_integral = volume_integral(
        uncollided.GetScalarFluxFieldFunction()[0],
        whole_domain,
        FieldFunctionInterpolationVolume,
    )

    quadrature = GLCProductQuadrature3DXYZ(
        n_polar=4,
        n_azimuthal=16,
        scattering_order=1,
    )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-9,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=file_name,
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    total_integral = volume_integral(
        problem.GetScalarFluxFieldFunction()[0],
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
    balance = solver.ComputeBalanceTable()

    collided_integral = total_integral - uncollided_integral
    collided_fraction = collided_integral / total_integral
    if rank == 0:
        print(f"UnstructuredCombinedBalance={balance['balance']:.8e}")
        print(f"UnstructuredCollidedFraction={collided_fraction:.8e}")

    remove_file(file_name)
    if abs(balance["balance"]) > 2.0e-5:
        raise RuntimeError(f"Combined balance failed: {balance['balance']}")
    if collided_fraction < 0.15:
        raise RuntimeError(
            f"First-collision source produced too little collided flux: {collided_fraction}"
        )
