#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SDM: PWLD

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume, FieldFunctionGridBased

if __name__ == "__main__":

    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/2blk_v2.msh",
    )
    grid = meshgen.Execute()

    # Material
    num_groups = 64
    xs_diag = MultiGroupXS()
    xs_diag.LoadFromOpenSn("diag_XS_64g_1mom_c0.99.xs")

    bsrc = []
    for g in range(num_groups):
        bsrc.append(0.0)
    bsrc[0] = 1.0

    # Quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16, scattering_order=1)

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
                "l_max_its": 300,
            },
        ],
        xs_map=[
            {"block_ids": [100, 101], "xs": xs_diag},
        ],
        boundary_conditions=[
            {"name": "left", "type": "isotropic", "group_strength": bsrc},
        ],
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    fflist = phys.GetScalarFieldFunctionList()

    vol0 = RPPLogicalVolume(xmin=-0.1, xmax=0.26, ymin=-0.1, ymax=1.1, infz=True)
    ffi1 = FieldFunctionInterpolationVolume()
    ffi1.SetOperationType("sum")
    ffi1.SetLogicalVolume(vol0)
    ffi1.AddFieldFunction(fflist[0])
    ffi1.Initialize()
    ffi1.Execute()
    maxval = ffi1.GetValue()
    if rank == 0:
        print(f"Sum={maxval:.5f}")
