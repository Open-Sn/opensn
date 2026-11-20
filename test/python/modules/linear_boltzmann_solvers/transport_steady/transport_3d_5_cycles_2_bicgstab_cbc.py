#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
Three-dimensional PWLD transport on an unstructured mesh with vacuum boundaries.

Compare the scalar-flux maxima in groups 0 and 1 with the AAH input
transport_3d_5_cycles_2_bicgstab.py.
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
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/sphere.case",
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            nz=1,
            xcuts=[0.0],
            ycuts=[0.0],
        )
    )
    grid = meshgen.Execute()

    num_groups = 5
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("../../../../assets/xs/xs_graphite_pure.xs")

    strength = [0.0 for _ in range(num_groups)]
    mg_src0 = VolumetricSource(block_ids=[1], group_strength=strength)
    strength[0] = 1.0
    mg_src1 = VolumetricSource(block_ids=[0], group_strength=strength)

    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_bicgstab",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 300,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        volumetric_sources=[mg_src0, mg_src1],
        sweep_type="CBC",
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    fflist = phys.GetScalarFluxFieldFunction(only_scalar_flux=False)

    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0][0])
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[1][0])
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
