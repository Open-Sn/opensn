#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D PWLD transport test with vacuum and incident-isotropic boundary conditions
Test: Max-value=0.49903 and 7.18243e-4
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn import can_support_gpus
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationLine, FieldFunctionInterpolationVolume
    from pyopensn.math import Vector3
    from pyopensn.logvol import RPPLogicalVolume
    if can_support_gpus:
        from pyopensn.device import set_device, get_device_count

if __name__ == "__main__":

    # Check number of processors
    num_procs = 3
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Check for GPU support
    if not can_support_gpus:
        sys.exit("OpenSn was built without GPU support.")

    # Divide MPI ranks to multiple GPUs
    set_device(rank % get_device_count())

    # Setup mesh
    nodes = []
    N = 100
    L = 30.0
    xmin = 0.0
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()

    # Cross-section data
    num_groups = 168
    grid.SetUniformBlockID(0)
    xs_3_170 = MultiGroupXS()
    xs_3_170.LoadFromOpenSn("xs_168g.xs")

    # Volumetric sources
    strength = []
    for g in range(num_groups):
        strength.append(0.0)
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    # Boundary sources
    bsrc = []
    for g in range(num_groups):
        bsrc.append(0.0)
    bsrc[0] = 1.0

    # Angular quadrature
    pquad = GLProductQuadrature1DSlab(n_polar=80, scattering_order=5)

    # Create solver
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
            {
                "block_ids": [0],
                "xs": xs_3_170
            }
        ],
        scattering_order=5,
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "zmin", "type": "isotropic", "group_strength": bsrc},
        ],
        options={
            "max_ags_iterations": 1
        },
        use_gpus=True
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Line plot
    cline = FieldFunctionInterpolationLine()
    cline.SetInitialPoint(Vector3(0.0, 0.0, 0.0001 + xmin))
    cline.SetFinalPoint(Vector3(0.0, 0.0, 29.999 + xmin))
    cline.SetNumberOfPoints(50)
    cline.AddFieldFunction(fflist[164][0])
    cline.Initialize()
    cline.Execute()

    # Volume integrations
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
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

    ffi2 = FieldFunctionInterpolationVolume()
    curffi = ffi2
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[159][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
