#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# - Final k-eigenvalue    :         1.1925596 (265)

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import PowerIterationKEigenSCDSASolver, PowerIterationKEigenSMMSolver, \
        NonLinearKEigenSolver

if __name__ == "__main__":

    # Mesh
    mesh_types = ['coarse', 'fine']
    if mesh_type not in mesh_types:
        if rank == 0:
            print("Error: 'mesh_type' must be one of: 'coarse', 'fine'")
        sys.exit(1)

    k_methods = ['pi', 'pi_scdsa', 'pi_scdsa_pwlc', 'pi_smm', 'pi_smm_pwld', 'jfnk']
    if k_method not in k_methods:
        if rank == 0:
            print("Error: 'k_method' must be one of: "
                  "'pi', 'pi_scdsa', 'pi_scdsa_pwlc','pi_smm', 'pi_smm_pwld', 'jfnk'")
        sys.exit(1)

    if rank == 0:
        print(f"Running C5G7 with mesh_type = {mesh_type} and k_method = {k_method}\n")

    if mesh_type == "coarse":
        mesh_file = "mesh/2D_c5g7_coarse.msh"
    else:
        mesh_file = "mesh/2D_c5g7_refined.msh"
    meshgen = FromFileMeshGenerator(filename=mesh_file)
    grid = meshgen.Execute()

    # Create cross sections
    xss = []
    for m in range(7):
        xss.append(MultiGroupXS())

    xss[0].LoadFromOpenSn("materials/XS_water.xs")
    xss[1].LoadFromOpenSn("materials/XS_UO2.xs")
    xss[2].LoadFromOpenSn("materials/XS_7pMOX.xs")
    xss[3].LoadFromOpenSn("materials/XS_guide_tube.xs")
    xss[4].LoadFromOpenSn("materials/XS_4_3pMOX.xs")
    xss[5].LoadFromOpenSn("materials/XS_8_7pMOX.xs")
    xss[6].LoadFromOpenSn("materials/XS_fission_chamber.xs")

    num_groups = xss[0].num_groups
    print("Num groups: ", num_groups)

    # Create materials
    xs_map = []
    for m in range(0, len(xss)):
        xs_map.append({"block_ids": [m], "xs": xss[m]})

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=1)

    # Solver
    if "scdsa" in k_method or "smm" in k_method:
        inner_linear_method = "classic_richardson"
        l_max_its = 2
    else:
        inner_linear_method = "petsc_gmres"
        l_max_its = 5

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": inner_linear_method,
                "l_max_its": l_max_its,
                "l_abs_tol": 1.0e-10,
                "angle_aggregation_type": "polar",
                "angle_aggregation_num_subsets": 1,
            },
        ],
        xs_map=xs_map,
        options={
            "boundary_conditions": [
                {"name": "xmin", "type": "reflecting"},
                {"name": "ymin", "type": "reflecting"},
            ],
            "scattering_order": 1,
            "verbose_outer_iterations": True,
            "verbose_inner_iterations": True,
            "power_field_function_on": True,
            "power_default_kappa": 1.0,
            "power_normalization": 1.0,
        },
        sweep_type="AAH",
    )

    # Execute Solver
    if k_method == "pi":
        k_solver = PowerIterationKEigenSolver(
            lbs_problem=phys,
            k_tol=1.0e-8,
        )
    elif k_method == "pi_scdsa":
        k_solver = PowerIterationKEigenSCDSASolver(
            lbs_problem=phys,
            diff_accel_sdm="pwld",
            accel_pi_verbose=True,
            k_tol=1.0e-8,
            accel_pi_k_tol=1.0e-8,
            accel_pi_max_its=50,
        )
    elif k_method == "pi_scdsa_pwlc":
        k_solver = PowerIterationKEigenSCDSASolver(
            lbs_problem=phys,
            diff_accel_sdm="pwlc",
            accel_pi_verbose=True,
            k_tol=1.0e-8,
            accel_pi_k_tol=1.0e-8,
            accel_pi_max_its=50,
        )
    elif k_method == "pi_smm":
        k_solver = PowerIterationKEigenSMMSolver(
            lbs_problem=phys,
            accel_pi_verbose=True,
            k_tol=1.0e-8,
            accel_pi_k_tol=1.0e-8,
            accel_pi_max_its=30,
            diff_sdm="pwlc",
        )
    elif k_method == "pi_smm_pwld":
        k_solver = PowerIterationKEigenSMMSolver(
            lbs_problem=phys,
            accel_pi_verbose=True,
            k_tol=1.0e-8,
            accel_pi_k_tol=1.0e-8,
            accel_pi_max_its=30,
            diff_sdm="pwld",
        )
    elif k_method == "jfnk":
        k_solver = NonLinearKEigenSolver(
            lbs_problem=phys,
            nl_max_its=50,
            nl_abs_tol=1.0e-9,
            nl_rel_tol=1.0e-9,
            l_abs_tol=1.0e-8,
            l_rel_tol=1.0e-8,
        )

    k_solver.Initialize()
    k_solver.Execute()
