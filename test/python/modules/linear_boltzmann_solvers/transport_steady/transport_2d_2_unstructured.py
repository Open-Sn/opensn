#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D PWLD unstructured mesh transport test with vacuum and incident-isotropic boundary conditions
Test: Max-value=0.51187 and 1.42458e-03
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":

    # Check number of processors
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Setup mesh using an unstructured triangle mesh file
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/TriangleMesh2x2Cuts.obj",
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            nz=1,
            xcuts=[0.0],
            ycuts=[0.0]
        )
    )
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetUniformBlockID(0)

    # Energy groups and cross-section data
    num_groups = 168
    xs_3_170 = MultiGroupXS()
    xs_3_170.LoadFromOpenSn("xs_168g.xs")

    # Sources
    strength = [0.0 for _ in range(num_groups)]
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)
    mg_src2 = VolumetricSource(block_ids=[2], group_strength=strength)

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=8, n_azimuthal=32, scattering_order=1)

    # Set up the boundary source
    bsrc = [0.0 for _ in range(num_groups)]
    bsrc[0] = 1.0

    # Create and configure the solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, 62),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
            {
                "groups_from_to": (63, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_3_170},
        ],
        scattering_order=1,
        options={
            "boundary_conditions": [
                {
                    "name": "xmin",
                    "type": "isotropic",
                    "group_strength": bsrc,
                },
            ],
            "max_ags_iterations": 1,
            "volumetric_sources": [mg_src1, mg_src2],
        }
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Retrieve field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Volume integration on the first field function
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5f}")

    # Volume integration on the tenth field function
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[9][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
