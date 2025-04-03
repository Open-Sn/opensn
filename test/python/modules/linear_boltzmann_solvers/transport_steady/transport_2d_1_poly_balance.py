#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 2D Transport test. Pure scatterer. Balance

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    # Check number of processors
    num_procs = 1
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    nodes = []
    N = 20
    L = 5
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)
    vol1 = RPPLogicalVolume(xmin=-1000.0, xmax=0.0, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    num_groups = 1
    xs_1g = MultiGroupXS()
    xs_1g.CreateSimpleOneGroup(1.0, 1.0)

    strength = [0.0 for _ in range(num_groups)]
    mg_src0 = VolumetricSource(block_ids=[0], group_strength=strength)
    strength[0] = 1.0
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)

    # Setup Physics
    fac = 1
    pquad = GLCProductQuadrature2DXY(6 * fac, 16 * fac)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_1g},
        ],
        options={
            "scattering_order": 0,
            "volumetric_sources": [mg_src0, mg_src1],
        },
    )

    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    phys.ComputeBalance()

    # Get field functions
    fflist = phys.GetScalarFieldFunctionList()

    # Volume integrations
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()

    if rank == 0:
        print(f"Max-value1={maxval:.5f}")

    # Volume integrations
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
