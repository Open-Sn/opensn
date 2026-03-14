#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPU steady-state test: update XS map at runtime and verify response changes.
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


def make_problem(xs):
    nodes = [i / 24.0 for i in range(25)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

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

    xs_lo = MultiGroupXS()
    xs_lo.CreateSimpleOneGroup(1.0, 0.0)
    xs_hi = MultiGroupXS()
    xs_hi.CreateSimpleOneGroup(2.0, 0.0)

    problem = make_problem(xs_lo)
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    phi_lo = compute_max_scalar_flux(problem)

    problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_hi}])
    solver.Execute()
    phi_hi = compute_max_scalar_flux(problem)

    ratio = phi_hi / phi_lo if phi_lo != 0.0 else 0.0
    pass_flag = 1 if (abs(ratio - 0.5) < 1.0e-3 and phi_hi < phi_lo) else 0

    if rank == 0:
        print(f"GPU_XSMAP_RATIO={ratio:.8e}")
        print(f"GPU_XSMAP_SWAP_PASS {pass_flag}")
