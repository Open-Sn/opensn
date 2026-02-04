#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_unstructured_reflective_uniform_1g.py: Unstructured 1-group with vacuum boundaries
Expected: MAX=0.9065113759285 from 3D Cartesian reference problem
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume


if __name__ == "__main__":
    if size != 4:
        sys.exit(f"Incorrect number of processors. Expected 4 but got {size}.")

    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/rz_rect_unstructured.msh",
        coord_sys="cylindrical",
    )
    grid = meshgen.Execute()

    vol = RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.0, ymax=2.0, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol, 0, True)

    xs = MultiGroupXS()
    sigma_t = 5.0
    c = 0.8
    xs.CreateSimpleOneGroup(sigma_t, c)

    sigma_a = sigma_t * (1.0 - c)
    src = VolumetricSource(block_ids=[0], group_strength=[sigma_a])

    quad = GLCProductQuadrature2DRZ(n_polar=6, n_azimuthal=12, scattering_order=0)
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 100,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "rmin", "type": "reflecting"},
            {"name": "rmax", "type": "vacuum"},
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    ff = problem.GetScalarFieldFunctionList(only_scalar_flux=False)[0][0]

    ffi_max = FieldFunctionInterpolationVolume()
    ffi_max.SetOperationType("max")
    ffi_max.SetLogicalVolume(vol)
    ffi_max.AddFieldFunction(ff)
    ffi_max.Initialize()
    ffi_max.Execute()
    phi_max = ffi_max.GetValue()

    if rank == 0:
        print(f"MAX {phi_max:.12e}")
