#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPU steady-state test: update densities at runtime and verify response changes.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def make_problem():
    nodes = [i / 24.0 for i in range(25)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)
    quad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)
    src = VolumetricSource(block_ids=[0], group_strength=[1.0])

    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": quad,
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 200,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        density={"default_density": 1.0},
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
            "max_ags_iterations": 1,
        },
        use_gpus=True,
    )


def compute_max_scalar_flux(problem):
    fflist = problem.GetScalarFluxFieldFunction()
    vol = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("max")
    ffi.SetLogicalVolume(vol)
    ffi.AddFieldFunction(fflist[0])
    ffi.Initialize()
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    if size != 3:
        sys.exit(f"Incorrect number of processors. Expected 3 but got {size}.")

    problem = make_problem()
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    phi_rho1 = compute_max_scalar_flux(problem)

    problem.SetDensities(by_block=[{"block_ids": [0], "density": 2.0}])
    solver.Execute()
    phi_rho2 = compute_max_scalar_flux(problem)

    ratio = phi_rho2 / phi_rho1 if phi_rho1 != 0.0 else 0.0
    pass_flag = 1 if (abs(ratio - 0.5) < 1.0e-3 and phi_rho2 < phi_rho1) else 0

    if rank == 0:
        print(f"GPU_DENSITY_RATIO={ratio:.8e}")
        print(f"GPU_DENSITY_UPDATE_PASS {pass_flag}")
