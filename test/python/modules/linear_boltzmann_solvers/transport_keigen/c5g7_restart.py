#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# - Final k-eigenvalue    :         1.1925596 (265)

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver

if __name__ == "__main__":

    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="c5g7/mesh/2D_c5g7_coarse.msh",
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()

    # Create cross sections
    xss = []
    for m in range(0, 6 + 1):
        xss.append(MultiGroupXS())

    xss[0].LoadFromOpenSn("c5g7/materials/XS_water.xs")
    xss[1].LoadFromOpenSn("c5g7/materials/XS_UO2.xs")
    xss[2].LoadFromOpenSn("c5g7/materials/XS_7pMOX.xs")
    xss[3].LoadFromOpenSn("c5g7/materials/XS_guide_tube.xs")
    xss[4].LoadFromOpenSn("c5g7/materials/XS_4_3pMOX.xs")
    xss[5].LoadFromOpenSn("c5g7/materials/XS_8_7pMOX.xs")
    xss[6].LoadFromOpenSn("c5g7/materials/XS_fission_chamber.xs")

    num_groups = xss[0].num_groups

    # Create materials
    xs_map = []
    for m in range(0, 6 + 1):
        xs_map.append({"block_ids": [m], "xs": xss[m]})

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=1)

    # Solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_max_its": 5,
                "l_abs_tol": 1.0e-10,
                "angle_aggregation_type": "polar",
                "angle_aggregation_num_subsets": 1,
            },
        ],
        xs_map=xs_map,
        scattering_order=1,
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
        ],
        options={
            "verbose_outer_iterations": True,
            "verbose_inner_iterations": True,
            "power_field_function_on": True,
            "power_default_kappa": 1.0,
            "power_normalization": 1.0,
            "save_angular_flux": True,
            "read_restart_path": "c5g7_restart/c5g7",
            # "restart_writes_enabled": True,
            # "write_delayed_psi_to_restart": True,
            # "write_restart_time_interval": 60,
            # "write_restart_path": "c5g7_restart/c5g7",
        },
        sweep_type="CBC",
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        k_tol=1.0e-8,
    )
    k_solver.Initialize()
    k_solver.Execute()
