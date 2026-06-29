#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-rank 3D 2G heterogeneous extruded-unstructured benchmark.

Environment overrides:
  OPENSN_UNSTRUCT_BASE      2D OBJ base mesh name, default triangle_mesh2x2_cuts.obj
  OPENSN_UNSTRUCT_SCALE     x-y mesh scale, default 7.0
  OPENSN_UNSTRUCT_LAYERS    extrusion sub-layers, default 12
  OPENSN_UNSTRUCT_POLAR     polar angles, default 6
  OPENSN_UNSTRUCT_AZIMUTH   azimuthal angles, default 12
  OPENSN_UNSTRUCT_WGS_ATOL   within-group solver absolute tolerance, default 1.0e-8
  OPENSN_UNSTRUCT_WGS_MAXITS within-group solver maximum iterations, default 100
  OPENSN_UNSTRUCT_SOLVER    eigensolver/accelerator, pi, nlke, scdsa, or smm, default pi
  OPENSN_UNSTRUCT_SCATTERING_ORDER scattering order, default 1
  OPENSN_UNSTRUCT_NL_ATOL   NLKE nonlinear absolute tolerance, default 1.0e-8
  OPENSN_UNSTRUCT_NL_MAXITS NLKE nonlinear maximum iterations, default 400
  OPENSN_UNSTRUCT_NL_L_ATOL NLKE linear absolute tolerance, default OPENSN_UNSTRUCT_WGS_ATOL
  OPENSN_UNSTRUCT_NL_L_MAXITS NLKE linear maximum iterations, default 400
  OPENSN_UNSTRUCT_NL_INITIAL_PI NLKE initial power iterations, default 0
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
        SCDSAAcceleration,
        SMMAcceleration,
    )
    from pyopensn.logvol import RPPLogicalVolume


def get_string_option(name, default):
    if name in globals():
        return str(globals()[name])
    return os.environ.get(name.upper(), default)


if __name__ == "__main__":

    num_procs = 8
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    base_mesh = os.environ.get("OPENSN_UNSTRUCT_BASE", "triangle_mesh2x2_cuts.obj")
    xy_scale = float(os.environ.get("OPENSN_UNSTRUCT_SCALE", "7.0"))
    n_layers = int(os.environ.get("OPENSN_UNSTRUCT_LAYERS", "12"))
    n_polar = int(os.environ.get("OPENSN_UNSTRUCT_POLAR", "6"))
    n_azimuthal = int(os.environ.get("OPENSN_UNSTRUCT_AZIMUTH", "12"))
    scattering_order = int(os.environ.get("OPENSN_UNSTRUCT_SCATTERING_ORDER", "1"))
    wgs_atol = float(os.environ.get("OPENSN_UNSTRUCT_WGS_ATOL", "1.0e-8"))
    wgs_maxits = int(os.environ.get("OPENSN_UNSTRUCT_WGS_MAXITS", "100"))
    solver_name = get_string_option("solver_name", os.environ.get("OPENSN_UNSTRUCT_SOLVER", "pi"))
    solver_name = solver_name.lower()
    if solver_name in ("smm", "pi_smm") and "OPENSN_UNSTRUCT_SCATTERING_ORDER" not in os.environ:
        scattering_order = 0
    nl_atol = float(os.environ.get("OPENSN_UNSTRUCT_NL_ATOL", "1.0e-8"))
    nl_maxits = int(os.environ.get("OPENSN_UNSTRUCT_NL_MAXITS", "400"))
    nl_l_atol = float(os.environ.get("OPENSN_UNSTRUCT_NL_L_ATOL", str(wgs_atol)))
    nl_l_maxits = int(os.environ.get("OPENSN_UNSTRUCT_NL_L_MAXITS", "400"))
    nl_initial_pi = int(os.environ.get("OPENSN_UNSTRUCT_NL_INITIAL_PI", "0"))
    accel_sdm = os.environ.get("OPENSN_UNSTRUCT_ACCEL_SDM", "pwld")

    zmax = 14.0
    meshgen = ExtruderMeshGenerator(
        inputs=[
            FromFileMeshGenerator(
                filename=f"../../../../assets/mesh/{base_mesh}",
                scale=xy_scale,
            )
        ],
        layers=[
            {"z": zmax, "n": n_layers},
        ],
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            nz=2,
            xcuts=[0.0],
            ycuts=[0.0],
            zcuts=[0.5 * zmax],
        ),
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    fuel = RPPLogicalVolume(
        xmin=-1000.0,
        xmax=3.0,
        ymin=-1000.0,
        ymax=3.0,
        zmin=-1000.0,
        zmax=10.0,
    )
    grid.SetBlockIDFromLogicalVolume(fuel, 1, True)

    xss = {}
    xss["0"] = MultiGroupXS()
    xss["0"].LoadFromOpenSn("../../../../assets/xs/xs_water_g2.xs")
    xss["1"] = MultiGroupXS()
    xss["1"].LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")
    num_groups = xss["0"].num_groups

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, num_groups - 1],
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=n_polar, n_azimuthal=n_azimuthal, scattering_order=scattering_order
                ),
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": wgs_atol,
                "l_max_its": wgs_maxits,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xss["0"]},
            {"block_ids": [1], "xs": xss["1"]},
        ],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "save_angular_flux": solver_name in ("smm", "pi_smm"),
        },
    )

    if solver_name in ("pi", "power_iteration", "poweriteration"):
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            max_iters=400,
            k_tol=1.0e-7,
        )
        solver_label = "PowerIteration"
    elif solver_name in ("nlke", "nonlinear", "nonlinear_keigen"):
        k_solver = NonLinearKEigenSolver(
            problem=phys,
            nl_max_its=nl_maxits,
            nl_abs_tol=nl_atol,
            l_abs_tol=nl_l_atol,
            l_max_its=nl_l_maxits,
            num_initial_power_iterations=nl_initial_pi,
        )
        solver_label = "NonLinearKEigen"
    elif solver_name in ("scdsa", "pi_scdsa"):
        scdsa = SCDSAAcceleration(
            problem=phys,
            sdm=accel_sdm,
        )
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=scdsa,
            max_iters=400,
            k_tol=1.0e-7,
        )
        solver_label = "PowerIterationSCDSA"
    elif solver_name in ("smm", "pi_smm"):
        smm = SMMAcceleration(
            problem=phys,
            sdm=accel_sdm,
            pi_k_tol=1.0e-8,
            pi_max_its=30,
        )
        k_solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=smm,
            max_iters=400,
            k_tol=1.0e-7,
        )
        solver_label = "PowerIterationSMM"
    else:
        raise ValueError(f"Unknown OPENSN_UNSTRUCT_SOLVER '{solver_name}'.")

    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    if rank == 0:
        print(f"Python solver: {solver_label}")
        print(f"Python scattering order: {scattering_order}")
        print(f"Python k-eigenvalue: {k}")
        if hasattr(k_solver, "GetNumSweeps"):
            print(f"Python sweeps: {k_solver.GetNumSweeps()}")

    uses_default_benchmark = (
        base_mesh == "triangle_mesh2x2_cuts.obj"
        and abs(xy_scale - 7.0) < 1.0e-12
        and n_layers == 12
        and n_polar == 6
        and n_azimuthal == 12
        and scattering_order == 1
    )
    if uses_default_benchmark:
        expected_k = 1.0441134
        if abs(k - expected_k) > 5.0e-5:
            raise RuntimeError(f"Expected k near {expected_k}, got {k}")
