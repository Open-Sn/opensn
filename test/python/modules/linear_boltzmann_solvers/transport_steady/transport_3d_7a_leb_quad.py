#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D Transport test with 3D Ortho Mesh - Lebedev Quadrature order 3 Half Scattering Ratio.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import LebedevQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    nodes = []
    n_cells = 10
    length = 2.0
    xmin = - length / 2
    dx = length / n_cells
    for i in range(n_cells + 1):
        nodes.append(xmin + i * dx)

    meshgen = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes, nodes],
    )

    grid = meshgen.Execute()

    grid.SetUniformBlockID(0)

    xs_mat = MultiGroupXS()
    xs_mat.CreateSimpleOneGroup(sigma_t=1., c=0.5)

    num_groups = 1
    strength = []
    for g in range(num_groups):
        strength.append(1.0)
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    pquad = LebedevQuadrature3DXYZ(quadrature_order=3, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30
            }
        ],
        volumetric_sources=[mg_src],
        xs_map=[
            {
                "block_ids": [0],
                "xs": xs_mat
            }
        ]
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    logvol = RPPLogicalVolume(
        xmin=0., xmax=length / 2,
        ymin=0., ymax=length / 2,
        zmin=0., zmax=length / 2
    )

    logvol_whole_domain = RPPLogicalVolume(infx=True, infy=True, infz=True)

    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("avg")
    ffi.SetLogicalVolume(logvol)
    ffi.AddFieldFunction(fflist[0][0])
    ffi.Initialize()
    ffi.Execute()
    val = ffi.GetValue()
    if rank == 0:
        print(f"Max-value1={val:.5e}")

    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("max")
    ffi.SetLogicalVolume(logvol_whole_domain)
    ffi.AddFieldFunction(fflist[0][0])
    ffi.Initialize()
    ffi.Execute()
    val_whole = ffi.GetValue()
    if rank == 0:
        print(f"Max-value2={val_whole:.5e}")
