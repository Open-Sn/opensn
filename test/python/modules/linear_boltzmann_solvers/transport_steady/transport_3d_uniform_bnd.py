#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SDM: PWLD

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased, FieldFunctionInterpolationVolume
    from pyopensn.logvol import SphereLogicalVolume

if __name__ == "__main__":

    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/sphere_hex.e",
    )
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    grid.SetUniformBoundaryID("outside")

    # Material
    num_groups = 64
    xs_diag = MultiGroupXS()
    xs_diag.LoadFromOpenSn("diag_XS_64g_1mom_c0.99.xs")

    # Boundary source
    bsrc = []
    for g in range(num_groups):
        bsrc.append(0.0)
    bsrc[0] = 1.0

    # Quadrature
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=16, scattering_order=0)

    # Set up solver
    gs1 = [0, num_groups - 1]
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": gs1,
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_diag},
        ],
        scattering_order=0,
        boundary_conditions=[
            {"name": "outside", "type": "isotropic", "group_strength": bsrc},
        ],
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    vol0 = SphereLogicalVolume(r=1.1, x=0., y=0., z=0.)
    fflist = phys.GetScalarFieldFunctionList()

    ffi1 = FieldFunctionInterpolationVolume()
    ffi1.SetOperationType("sum")
    ffi1.SetLogicalVolume(vol0)
    ffi1.AddFieldFunction(fflist[0])
    ffi1.Initialize()
    ffi1.Execute()
    maxval = ffi1.GetValue()
    if rank == 0:
        print(f"Sum={maxval:.5f}")

    FieldFunctionGridBased.ExportMultipleToPVTU([fflist[0]], "test")
