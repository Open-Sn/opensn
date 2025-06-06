#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D PWLD unstructure mesh Transport test with vacuum boundary conditions
Test: Max-value1=6.55387e+00 Max-value2=1.02940e+00
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner, MeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/Sphere.case",
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            nz=1,
            xcuts=[0.0],
            ycuts=[0.0],
        )
    )
    grid = meshgen.Execute()

    # Define energy groups and load cross-section data
    num_groups = 5
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")

    # Create volumetric source
    strength = [0.0 for _ in range(num_groups)]
    mg_src0 = VolumetricSource(block_ids=[1], group_strength=strength)
    strength[0] = 1.0
    mg_src1 = VolumetricSource(block_ids=[0], group_strength=strength)

    # Setup the angular quadrature
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8)

    # Create and configure the discrete ordinates solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        options={
            # "restart_writes_enabled": True,
            # "write_delayed_psi_to_restart": True,
            # "write_restart_path": "transport_3d_5_cycles_2_restart/transport_3d_5_cycles_2",
            "read_restart_path": "transport_3d_5_cycles_2_restart/transport_3d_5_cycles_2",
            "scattering_order": 0,
            "volumetric_sources": [mg_src0, mg_src1],
        }
    )
    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Volume integration
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
        print(f"Max-value1={maxval:.5e}")

    # Volume integration
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[1][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
