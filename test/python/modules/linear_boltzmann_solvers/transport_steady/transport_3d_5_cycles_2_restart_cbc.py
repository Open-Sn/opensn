#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
Three-dimensional PWLD transport on an unstructured mesh with vacuum boundaries.

Compare the scalar-flux maxima in groups 0 and 1 with the AAH input
transport_3d_5_cycles_2_restart.py.

Load the checkpoint from transport_3d_5_cycles_2_cbc.py on the same four-rank
partition. The initial scalar-flux DOFs must agree exactly. The continued
flux must agree with the uninterrupted flux to 1e-9 relative error, with a
1e-12 linear tolerance in both solves.
"""

import os
import sys
import math

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    def global_max(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner, MeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
else:
    def global_max(value):
        return MPIAllReduce(value, "max")

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
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        volumetric_sources=[mg_src0, mg_src1],
        options={
            "read_restart_path": "cbc_cycle_checkpoint",
        },
        sweep_type="CBC",
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    reference = np.load(f"cbc_cycle_reference_{rank}.npy", allow_pickle=False)
    scale = global_max(float(np.max(np.abs(reference))))
    ss_solver.Initialize()
    restored = np.array(phys.GetPhiNewLocal(), copy=True)
    load_error = global_max(float(np.max(np.abs(restored - reference))))
    ss_solver.Execute()
    resumed = np.array(phys.GetPhiNewLocal(), copy=True)
    resume_error = global_max(float(np.max(np.abs(resumed - reference)))) / scale
    if rank == 0:
        print(f"CBC_RESTART_LOAD_MAX_ABS_DIFF={load_error:.16e}")
        print(f"CBC_RESTART_RESUME_MAX_REL_DIFF={resume_error:.16e}")

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
