#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D PWLD transport test with vacuum and Incident-isotropic bundary conditions
Test: Max-value=0.50758 and 2.52527e-04
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
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    # Check number of processors
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/SquareMesh2x2QuadsBlock.obj",
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            nz=1,
            xcuts=[0.0],
            ycuts=[0.0],
        )
    )
    grid = meshgen.Execute()

    # Cross-section data
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)
    num_groups = 168
    xs_3_170 = MultiGroupXS()
    xs_3_170.LoadFromOpenSn("xs_168g.xs")

    # Volumetric sources
    strength = []
    for g in range(num_groups):
        strength.append(0.0)
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)
    mg_src2 = VolumetricSource(block_ids=[2], group_strength=strength)

    # Boundary sources
    bsrc = []
    for g in range(num_groups):
        bsrc.append(0.0)
    bsrc[0] = 1.0 / 4.0 / math.pi

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(2, 8)

    # Create solver
    phys = DiscreteOrdinatesSolver(
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
            {
                "block_ids": [0],
                "xs": xs_3_170
            }
        ],
        options={
            "boundary_conditions": [
                {"name": "xmin", "type": "isotropic", "group_strength": bsrc},
            ],
            "volumetric_sources": [mg_src1, mg_src2],
            "scattering_order": 1,
            "max_ags_iterations": 1,
            "save_angular_flux": True
        },
        sweep_type="CBC"
    )
    ss_solver = SteadyStateSolver(lbs_solver=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Volume integrations
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

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[159][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
