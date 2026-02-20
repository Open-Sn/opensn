#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_2g_scatter_vacuum.py: 2-group scatter with vacuum boundaries
Expected: G0_MAX=0.8429208738494, G1_MAX=0.3655486193072 from 3D Cartesian calculation
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume


if __name__ == "__main__":
    # Structured cylindrical mesh on 0<=r<=1, 0<=z<=2.
    nr = 50
    nz = 100
    r_nodes = [i / nr for i in range(nr + 1)]
    z_nodes = [2.0 * i / nz for i in range(nz + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[r_nodes, z_nodes], coord_sys="cylindrical")
    grid = meshgen.Execute()

    vol = RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.0, ymax=2.0, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol, 0, True)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("transport_2d_cyl_2g_2.xs")

    src = VolumetricSource(block_ids=[0], group_strength=[1.5, 0.8])

    quad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": quad,
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-7,
                "l_max_its": 5000,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "rmin", "type": "reflecting"},
            {"name": "rmax", "type": "vacuum"},
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    fflist = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)

    ffi_g0 = FieldFunctionInterpolationVolume()
    ffi_g0.SetOperationType("max")
    ffi_g0.SetLogicalVolume(vol)
    ffi_g0.AddFieldFunction(fflist[0][0])
    ffi_g0.Initialize()
    ffi_g0.Execute()
    phi_g0 = ffi_g0.GetValue()

    ffi_g1 = FieldFunctionInterpolationVolume()
    ffi_g1.SetOperationType("max")
    ffi_g1.SetLogicalVolume(vol)
    ffi_g1.AddFieldFunction(fflist[1][0])
    ffi_g1.Initialize()
    ffi_g1.Execute()
    phi_g1 = ffi_g1.GetValue()

    if rank == 0:
        print(f"G0_MAX {phi_g0:.12e}")
        print(f"G1_MAX {phi_g1:.12e}")
