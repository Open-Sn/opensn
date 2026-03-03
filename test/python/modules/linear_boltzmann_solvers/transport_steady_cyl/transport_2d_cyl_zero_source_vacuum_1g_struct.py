#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_zero_source_vacuum_1g_struct.py: Zero source, vacuum boundaries
Expected: PHI_MAX=0.0 from analytic zero-source solution
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
    nr = 40
    nz = 80
    r_nodes = [i * (1.0 / nr) for i in range(nr + 1)]
    z_nodes = [i * (2.0 / nz) for i in range(nz + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[r_nodes, z_nodes], coord_sys="cylindrical")
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)
    src = VolumetricSource(block_ids=[0], group_strength=[0.0])

    quad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)
    bcs = [
        {"name": "rmin", "type": "vacuum"},
        {"name": "rmax", "type": "vacuum"},
        {"name": "zmin", "type": "vacuum"},
        {"name": "zmax", "type": "vacuum"},
    ]

    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 50,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=bcs,
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    vol = RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.0, ymax=2.0, infz=True)
    fflist = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)

    ffi_max = FieldFunctionInterpolationVolume()
    ffi_max.SetOperationType("max")
    ffi_max.SetLogicalVolume(vol)
    ffi_max.AddFieldFunction(fflist[0][0])
    ffi_max.Initialize()
    ffi_max.Execute()
    phi_max = ffi_max.GetValue()

    if rank == 0:
        print(f"PHI_MAX {phi_max:.12e}")
