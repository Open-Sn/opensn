#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D PWLD adjoint transport test with localized material source
Test: QoI Value=1.38399e-05
      Inner Product=1.38405e-05
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
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
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

    # Create mesh
    N = 60
    L = 5.0
    ds = L / N
    nodes = []
    for i in range(N + 1):
        nodes.append(i * ds)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Define a logical volume that spans in x (infinite) and in y from 0 to 0.8*L.
    vol1a = RPPLogicalVolume(infx=True, ymin=0.0, ymax=0.8 * L, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol1a, 1, True)

    # Define a volume for block 0: localized in x around 2.5.
    vol0 = RPPLogicalVolume(xmin=2.5 - 0.166666, xmax=2.5 + 0.166666, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    # Define a volume for block 2: localized in x and in y from 0 to 2*0.166666.
    vol2 = RPPLogicalVolume(
        xmin=2.5 - 0.166666,
        xmax=2.5 + 0.166666,
        ymin=0.0,
        ymax=2 * 0.166666,
        infz=True,
    )
    grid.SetBlockIDFromLogicalVolume(vol2, 2, True)

    # Define a second volume for block 1: localized in x around 2.5 and in y from 0.9*L to L.
    vol1b = RPPLogicalVolume(
        xmin=2.5 - 1,  # equivalent to (-1 + 2.5)
        xmax=2.5 + 1,  # equivalent to (1 + 2.5)
        ymin=0.9 * L,
        ymax=L,
        infz=True,
    )
    grid.SetBlockIDFromLogicalVolume(vol1b, 1, True)

    # Add cross sections to materials
    xs_1g1 = MultiGroupXS()
    xs_1g1.CreateSimpleOneGroup(0.01, 0.01)

    xs_1g2 = MultiGroupXS()
    xs_1g2.CreateSimpleOneGroup(0.1 * 20, 0.8)

    xs_1g3 = MultiGroupXS()
    xs_1g3.CreateSimpleOneGroup(0.3 * 20, 0.0)

    # Create sources
    # Forward source: only active in block 2 with a group strength of 3.0.
    fwd_src = VolumetricSource(block_ids=[2], group_strength=[3.0])

    # Create a 2D angular quadrature with 12 polar and 192 azimuthal angles.
    pquad = GLCProductQuadrature2DXY(n_polar=12, n_azimuthal=192, scattering_order=0)

    # Setup solver
    phys = DiscreteOrdinatesProblem(
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
            {"block_ids": [2], "xs": xs_1g3},
        ],
        scattering_order=0,
        volumetric_sources=[fwd_src],
    )

    # Forward solve
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Define QoI region and compute forward QoI
    qoi_vol = RPPLogicalVolume(xmin=0.5, xmax=0.8333, ymin=4.16666, ymax=4.33333, infz=True)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("sum")
    ffi.SetLogicalVolume(qoi_vol)
    ffi.AddFieldFunction(fflist[0][0])
    ffi.Initialize()
    ffi.Execute()
    fwd_qoi = ffi.GetValue()

    # Create adjoint source and switch to adjoint mode
    adj_src = VolumetricSource(logical_volume=qoi_vol, group_strength=[1.0])
    adjoint_options = {
        "adjoint": True,
        "volumetric_sources": [adj_src],
    }
    phys.SetOptions(
        adjoint=True
    )
    phys.SetVolumetricSources(volumetric_sources=[adj_src])

    # Adjoint solve and write flux moments
    ss_solver.Execute()
    phys.WriteFluxMoments("adjoint_2d_1")

    # Create response evaluator and evaluate response
    evaluator = ResponseEvaluator(problem=phys)
    evaluator.SetOptions(
        buffers=[{'name': 'buff', 'file_prefixes': {'flux_moments': 'adjoint_2d_1'}}],
        sources={'material': [{'block_id': 2, 'strength': [3.0]}]}
    )
    adj_qoi = evaluator.EvaluateResponse("buff")

    # Print results
    if rank == 0:
        print(f"QoI Value={fwd_qoi:.5e}")
        print(f"Inner Product={adj_qoi:.5e}")

    # Cleanup
    barrier()
    if rank == 0:
        os.system("rm adjoint_2d_1*")
