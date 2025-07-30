#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D PWLD adjoint transport test with point source Multigroup FWD
Test:
  QoI Value[0]= 1.12687e-06
  QoI Value[1]= 2.95934e-06
  QoI Value[2]= 3.92975e-06
  QoI Value[3]= 4.18474e-06
  QoI Value[4]= 3.89649e-06
  QoI Value[5]= 3.30482e-06
  QoI Value[6]= 1.54506e-06
  QoI Value[7]= 6.74868e-07
  QoI Value[8]= 3.06178e-07
  QoI Value[9]= 2.07284e-07
  sum(QoI Value)= 2.21354e-05
  Inner Product=3.30607e-06
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import PointSource, VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.response import ResponseEvaluator
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.math import VectorSpatialFunction
else:
    barrier = MPIBarrier

if __name__ == "__main__":

    # Check number of processors
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Create mesh
    N = 60
    L = 5.0
    ds = L / N
    nodes = [i * ds for i in range(N + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Set block IDs using logical volumes
    # Volume for block 1: infinite in x; y from 0 to 0.8*L.
    vol1a = RPPLogicalVolume(infx=True, ymin=0.0, ymax=0.8 * L, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol1a, 1, True)

    # Volume for block 0: localized in x around 2.5.
    vol0 = RPPLogicalVolume(xmin=2.5 - 0.166666, xmax=2.5 + 0.166666, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    # Volume for block 1 (second region): x centered about 2.5, y from 0.9*L to L.
    vol1b = RPPLogicalVolume(xmin=2.5 - 1, xmax=2.5 + 1, ymin=0.9 * L, ymax=L, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol1b, 1, True)

    # Add cross sections to materials
    num_groups = 10
    xs_1 = MultiGroupXS()
    xs_1.LoadFromOpenSn("response_2d_3_mat1.xs")

    xs_2 = MultiGroupXS()
    xs_2.LoadFromOpenSn("response_2d_3_mat2.xs")

    # Create point source (multigroup)
    src = [0.0] * num_groups
    src[0] = 1.0

    loc = [1.25 - 0.5 * ds, 1.5 * ds, 0.0]
    pt_src = PointSource(location=loc, strength=src)

    # Create a 2D angular quadrature with 4 polar and 48 azimuthal angles.
    pquad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=48, scattering_order=0)

    # Setup physics
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 500,
                "gmres_restart_interval": 100,
            }
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_1},
            {"block_ids": [1], "xs": xs_2},
        ],
        scattering_order=0,
        options={
            "point_sources": [pt_src],
        },
    )

    # Forward solve
    ss_solver = SteadyStateSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Define QoI region
    qoi_vol = RPPLogicalVolume(xmin=0.5, xmax=0.8333, ymin=4.16666, ymax=4.33333, infz=True)

    # Compute QoI for each group and sum them
    fwd_qois = []
    fwd_qoi_sum = 0.0
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)
    for g in range(num_groups):
        ff = fflist[g][0]
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType("sum")  # OP_SUM operation
        ffi.SetLogicalVolume(qoi_vol)
        ffi.AddFieldFunction(ff)
        ffi.Initialize()
        ffi.Execute()
        value = ffi.GetValue()
        fwd_qois.append(value)
        fwd_qoi_sum += value

    # Create adjoint source using a response function
    def ResponseFunction(xyz, mat_id):
        response = [0.0] * num_groups
        response[5] = 1.0
        return response
    response_func = VectorSpatialFunction(ResponseFunction)

    # Create the adjoint volumetric source with the response function over the QoI region.
    adjoint_source = VolumetricSource(
        logical_volume=qoi_vol,
        func=response_func
    )

    # Switch to adjoint mode
    phys.SetOptions(
        adjoint=True,
        volumetric_sources=[adjoint_source]
    )

    # Adjoint solve and write flux moments
    ss_solver.Execute()
    phys.WriteFluxMoments("adjoint_2d_3")

    # Create response evaluator and evaluate response
    evaluator = ResponseEvaluator(problem=phys)
    evaluator.SetOptions(
        buffers=[{"name": "buff", "file_prefixes": {"flux_moments": "adjoint_2d_3"}}],
        sources={"point": [pt_src]}
    )
    response_val = evaluator.EvaluateResponse("buff")

    # Print results
    if rank == 0:
        for g in range(num_groups):
            print(f"QoI Value[{g}]= {fwd_qois[g]:.5e}")
        print(f"sum(QoI Values)= {fwd_qoi_sum:.5e}")
        print(f"Inner Product= {response_val:.5e}")

    # Cleanup
    barrier()
    if rank == 0:
        os.system("rm adjoint_2d_3*")
