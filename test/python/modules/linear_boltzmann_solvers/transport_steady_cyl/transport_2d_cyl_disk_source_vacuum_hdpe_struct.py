#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_disk_source_vacuum_hdpe_struct.py: Disk with HDPE OpenMC cross sections
Expected: PHI_MAX=1.224546944901e-01 from 3D Cartesian reference problem
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
    # Structured cylindrical mesh on 0<=r<=0.2, 0<=z<=0.1.
    nr = 100
    nz = 50
    r_nodes = [0.2 * i / nr for i in range(nr + 1)]
    z_nodes = [0.1 * i / nz for i in range(nz + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[r_nodes, z_nodes], coord_sys="cylindrical")
    grid = meshgen.Execute()

    vol = RPPLogicalVolume(xmin=0.0, xmax=0.2, ymin=0.0, ymax=0.1, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol, 0, True)

    xs = MultiGroupXS()
    xs.LoadFromOpenMC("HDPE.h5", "set1", 294.0)

    num_groups = 172
    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 1.0
    src = VolumetricSource(block_ids=[0], group_strength=strength)

    quad = GLCProductQuadrature2DRZ(n_polar=8, n_azimuthal=16, scattering_order=0)
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": quad,
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
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

    fflist = problem.GetScalarFieldFunctionList(only_scalar_flux=False)
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("max")
    ffi.SetLogicalVolume(vol)
    ffi.AddFieldFunction(fflist[0][0])
    ffi.Initialize()
    ffi.Execute()
    phi_max = ffi.GetValue()

    if rank == 0:
        print(f"PHI_MAX {phi_max:.12e}")
