#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D Transport test with point source FWD
SDM: PWLD
Test: QoI Value=2.90386e-05
      Inner Product=2.90543e-05
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import PointSource, VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.response import ResponseEvaluator
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
else:
    barrier = MPIBarrier

if __name__ == "__main__":

    # Check number of processors
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    N = 60
    L = 5.0
    ds = L / N
    nodes = [i * ds for i in range(N + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
        grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Define logical volumes and assign block IDs
    # Logical volume for block 1 (infinite in x, y from 0 to 0.8*L)
    vol1a = RPPLogicalVolume(infx=True, ymin=0.0, ymax=0.8 * L, infz=True)
    grid.SetBlockIDFromLogical(vol1a, 1, True)

    # Logical volume for block 0 (localized in x around 2.5)
    vol0 = RPPLogicalVolume(xmin=2.5 - 0.166666, xmax=2.5 + 0.166666, infy=True, infz=True)
    grid.SetBlockIDFromLogical(vol0, 0, True)

    # Logical volume for block 1 (second region: x centered about 2.5 and y from 0.9*L to L)
    vol1b = RPPLogicalVolume(xmin=2.5 - 1, xmax=2.5 + 1, ymin=0.9 * L, ymax=L, infz=True)
    grid.SetBlockIDFromLogical(vol1b, 1, True)

    # Add cross sections to materials
    xs_1g1 = MultiGroupXS()
    xs_1g1.CreateSimpleOneGroup(0.01, 0.01)

    xs_1g2 = MultiGroupXS()
    xs_1g2.CreateSimpleOneGroup(0.1 * 20, 0.8)

    # Create point source
    loc = [1.25 - 0.5 * ds, 1.5 * ds, 0.0]
    pt_src = PointSource(location=loc, strength=[1.0])

    # Create a 2D angular quadrature with 12 polar and 192 azimuthal angles.
    pquad = GLCProductQuadrature2DXY(12, 192)

    # Setup physics and solver
    phys = DiscreteOrdinatesSolver(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 500,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_1g1},
            {"block_ids": [1], "xs": xs_1g2},
        ],
        options={
            "scattering_order": 0,
            "point_sources": [pt_src],
        },
    )

    # Forward solve
    ss_solver = SteadyStateSolver(lbs_solver=phys)
        ss_solver.Initialize()
        ss_solver.Execute()

    # Get field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)
    ff_m0 = fflist[0][0]

    # Define QoI region and compute forward QoI
    qoi_vol = RPPLogicalVolume(xmin=0.5, xmax=0.8333, ymin=4.16666, ymax=4.33333, infz=True)

    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("sum")  # Using a sum operation (corresponding to OP_SUM)
    ffi.SetLogicalVolume(qoi_vol)
    ffi.AddFieldFunction(ff_m0)
    ffi.Initialize()
    ffi.Execute()
    fwd_qoi = ffi.GetValue()

    # Create adjoint source and switch to adjoint mode
    adj_src = VolumetricSource(logical_volume=qoi_vol, group_strength=[1.0])
    phys.SetOptions(adjoint=True, volumetric_sources=[adj_src])

    # Adjoint solve and write flux moments
        ss_solver.Execute()
    phys.WriteFluxMoments("adjoint_2d_2")

    # Create response evaluator and evaluate response
    evaluator = ResponseEvaluator(lbs_solver=phys)
    evaluator.SetOptions(
        buffers=[{'name': 'buff', 'file_prefixes': {'flux_moments': 'adjoint_2d_2'}}],
        sources={'point': [pt_src]}
    )
    adj_qoi = evaluator.EvaluateResponse("buff")

    # Print results
    if rank == 0:
        print(f"QoI Value={fwd_qoi:.5e}")
        print(f"Inner Product={adj_qoi:.5e}")

    # Cleanup
    barrier()
    if rank == 0:
        os.system("rm adjoint_2d_2*")
