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
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import NonLinearKEigenSolver
    from pyopensn.solver import CMFDAcceleration, SCDSAAcceleration, SMMAcceleration

def get_bool_option(name, default):
    if name in globals():
        return bool(globals()[name])
    return os.environ.get(name.upper(), "1" if default else "0") == "1"

def get_option(name, default):
    if name in globals():
        return globals()[name]
    return os.environ.get(name.upper(), default)

if __name__ == "__main__":
    asset_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../assets"))

    # Mesh
    mesh_types = ['coarse', 'fine']
    if mesh_type not in mesh_types:
        if rank == 0:
            print("Error: 'mesh_type' must be one of: 'coarse', 'fine'")
        sys.exit(1)

    k_methods = ['pi', 'pi_cmfd', 'pi_scdsa', 'pi_scdsa_pwlc', 'pi_smm', 'pi_smm_pwld', 'jfnk']
    if k_method not in k_methods:
        if rank == 0:
            print("Error: 'k_method' must be one of: "
                  "'pi', 'pi_cmfd', 'pi_scdsa', 'pi_scdsa_pwlc', 'pi_smm', "
                  "'pi_smm_pwld', 'jfnk'")
        sys.exit(1)

    if rank == 0:
        print(f"Running C5G7 with mesh_type = {mesh_type} and k_method = {k_method}\n")

    if mesh_type == "coarse":
        mesh_file = os.path.join(asset_dir, "mesh/c5g7/2d_c5g7_coarse.msh")
    else:
        mesh_file = os.path.join(asset_dir, "mesh/c5g7/2d_c5g7_refined.msh")
    meshgen = FromFileMeshGenerator(filename=mesh_file)
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()

    # Create cross sections
    xss = []
    for m in range(7):
        xss.append(MultiGroupXS())

    xss[0].LoadFromOpenSn(os.path.join(asset_dir, "xs/c5g7/XS_water.xs"))
    xss[1].LoadFromOpenSn(os.path.join(asset_dir, "xs/c5g7/XS_UO2.xs"))
    xss[2].LoadFromOpenSn(os.path.join(asset_dir, "xs/c5g7/XS_7pMOX.xs"))
    xss[3].LoadFromOpenSn(os.path.join(asset_dir, "xs/c5g7/XS_guide_tube.xs"))
    xss[4].LoadFromOpenSn(os.path.join(asset_dir, "xs/c5g7/XS_4_3pMOX.xs"))
    xss[5].LoadFromOpenSn(os.path.join(asset_dir, "xs/c5g7/XS_8_7pMOX.xs"))
    xss[6].LoadFromOpenSn(os.path.join(asset_dir, "xs/c5g7/XS_fission_chamber.xs"))

    num_groups = xss[0].num_groups

    # Create materials
    xs_map = []
    for m in range(0, len(xss)):
        xs_map.append({"block_ids": [m], "xs": xss[m]})

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=0)

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
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
        ],
        options={
            "verbose_outer_iterations": True,
            "verbose_inner_iterations": True,
            "save_angular_flux": "smm" in k_method,
            "power_default_kappa": 1.0,
        },
        sweep_type="AAH",
    )

    # Execute Solver
    if k_method == "pi":
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            k_tol=1.0e-8,
        )
    elif k_method == "pi_cmfd":
        cmfd = CMFDAcceleration(
            problem=phys,
            coarse_mesh=get_option("cmfd_coarse_mesh", "local_aggregation"),
            aggregation_size=int(get_option("cmfd_aggregation_size", 16)),
            update_scheme=get_bool_option("cmfd_update_scheme", True),
            update_wgs_max_its=int(get_option("cmfd_update_wgs_max_its", 4)),
            update_wgs_abs_tol=float(get_option("cmfd_update_wgs_abs_tol", 1.0e-4)),
            relaxation=float(get_option("cmfd_relaxation", 1.0)),
            adaptive_relaxation=get_bool_option("cmfd_adaptive_relaxation", True),
            adaptive_relaxation_min=float(get_option("cmfd_adaptive_relaxation_min", 0.25)),
            adaptive_relaxation_max=float(get_option("cmfd_adaptive_relaxation_max", 1.0)),
            adaptive_relaxation_growth=float(get_option("cmfd_adaptive_relaxation_growth", 1.25)),
            adaptive_relaxation_reduction=float(
                get_option("cmfd_adaptive_relaxation_reduction", 0.75)
            ),
            adaptive_relaxation_accept_fraction=float(
                get_option("cmfd_adaptive_relaxation_accept_fraction", 0.5)
            ),
            adaptive_relaxation_successes_to_grow=int(
                get_option("cmfd_adaptive_relaxation_successes_to_grow", 2)
            ),
            inactive_iterations=int(get_option("cmfd_inactive_iterations", 1)),
            l_abs_tol=float(get_option("cmfd_l_abs_tol", 1.0e-10)),
            max_iters=int(get_option("cmfd_l_max_its", 100)),
            pi_max_its=int(get_option("cmfd_pi_max_its", 5)),
            pi_k_tol=float(get_option("cmfd_pi_k_tol", 1.0e-10)),
            coarse_solver_policy=get_option("cmfd_coarse_solver_policy", "auto"),
            correction_max_attempts=int(get_option("cmfd_correction_max_attempts", 10)),
            correction_min_damping=float(get_option("cmfd_correction_min_damping", 1.0e-4)),
            negative_flux_tolerance=float(get_option("cmfd_negative_flux_tolerance", 1.0e-6)),
            verbose=get_bool_option("cmfd_verbose", False),
        )
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=cmfd,
            k_tol=1.0e-8,
        )
    elif k_method == "pi_scdsa":
        scdsa = SCDSAAcceleration(
            problem=phys,
            sdm="pwld",
            verbose=True,
            pi_k_tol=1.0e-8,
            pi_max_its=50
        )
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=scdsa,
            k_tol=1.0e-8
        )
    elif k_method == "pi_scdsa_pwlc":
        scdsa = SCDSAAcceleration(
            problem=phys,
            sdm="pwlc",
            verbose=True,
            pi_k_tol=1.0e-8,
            pi_max_its=50
        )
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=scdsa,
            k_tol=1.0e-8,
        )
    elif k_method == "pi_smm":
        smm = SMMAcceleration(
            problem=phys,
            verbose=True,
            pi_k_tol=1.0e-8,
            pi_max_its=30,
            sdm="pwlc"
        )
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=smm,
            k_tol=1.0e-8
        )
    elif k_method == "pi_smm_pwld":
        smm = SMMAcceleration(
            problem=phys,
            verbose=True,
            pi_k_tol=1.0e-8,
            pi_max_its=30,
            sdm="pwld"
        )
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=smm,
            k_tol=1.0e-8,
        )
    elif k_method == "jfnk":
        k_solver = NonLinearKEigenSolver(
            problem=phys,
            nl_max_its=50,
            nl_abs_tol=1.0e-9,
            nl_rel_tol=1.0e-9,
            l_abs_tol=1.0e-8,
            l_rel_tol=1.0e-8,
        )

    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    if rank == 0:
        print(f"Python method: {k_method}")
        print(f"Python k-eigenvalue: {k}")
        if hasattr(k_solver, "GetNumSweeps"):
            print(f"Python sweeps: {k_solver.GetNumSweeps()}")
