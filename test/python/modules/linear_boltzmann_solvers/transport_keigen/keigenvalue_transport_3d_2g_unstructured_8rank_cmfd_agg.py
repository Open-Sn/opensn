#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-rank 3D 2G heterogeneous extruded-unstructured benchmark using rank-local aggregated CMFD.

Environment overrides:
  OPENSN_UNSTRUCT_BASE      2D OBJ base mesh name, default triangle_mesh2x2_cuts.obj
  OPENSN_UNSTRUCT_SCALE     x-y mesh scale, default 7.0
  OPENSN_UNSTRUCT_LAYERS    extrusion sub-layers, default 12
  OPENSN_UNSTRUCT_POLAR     polar angles, default 6
  OPENSN_UNSTRUCT_AZIMUTH   azimuthal angles, default 12
  OPENSN_UNSTRUCT_WGS_ATOL   within-group solver absolute tolerance, default 1.0e-8
  OPENSN_UNSTRUCT_WGS_MAXITS within-group solver maximum iterations, default 100
  OPENSN_UNSTRUCT_CMFD_MESH CMFD coarse mesh, default local_aggregation
  OPENSN_UNSTRUCT_AGG       target fine cells per coarse cell, default 16
  OPENSN_UNSTRUCT_CMFD_VERBOSE enables CMFD timing diagnostics when set to 1
  OPENSN_UNSTRUCT_CMFD_ATOL  CMFD linear absolute tolerance, default 1.0e-10
  OPENSN_UNSTRUCT_CMFD_MAXITS CMFD linear maximum iterations, default 100
  OPENSN_UNSTRUCT_CMFD_UPDATE_SCHEME enables CMFD-owned loose transport updates, default 1
  OPENSN_UNSTRUCT_CMFD_UPDATE_WGS_ATOL update-scheme WGS tolerance, default 1.0e-4
  OPENSN_UNSTRUCT_CMFD_UPDATE_WGS_MAXITS update-scheme WGS max iterations, default 5
  OPENSN_UNSTRUCT_CMFD_PI_MAXITS CMFD inner PI maximum iterations, default 1
  OPENSN_UNSTRUCT_CMFD_PI_KTOL CMFD inner PI k tolerance, default 1.0e-10
  OPENSN_UNSTRUCT_CMFD_COARSE_SOLVER CMFD coarse solver policy, default auto
  OPENSN_UNSTRUCT_CMFD_PETSC_OPTIONS PETSc options for the CMFD KSP when requested
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
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration
    from pyopensn.logvol import RPPLogicalVolume

def get_bool_option(name, default):
    if name in globals():
        return bool(globals()[name])
    return os.environ.get(name.upper(), "1" if default else "0") == "1"

if __name__ == "__main__":

    num_procs = 8
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    base_mesh = os.environ.get("OPENSN_UNSTRUCT_BASE", "triangle_mesh2x2_cuts.obj")
    xy_scale = float(os.environ.get("OPENSN_UNSTRUCT_SCALE", "7.0"))
    n_layers = int(os.environ.get("OPENSN_UNSTRUCT_LAYERS", "12"))
    n_polar = int(os.environ.get("OPENSN_UNSTRUCT_POLAR", "6"))
    n_azimuthal = int(os.environ.get("OPENSN_UNSTRUCT_AZIMUTH", "12"))
    wgs_atol = float(os.environ.get("OPENSN_UNSTRUCT_WGS_ATOL", "1.0e-8"))
    wgs_maxits = int(os.environ.get("OPENSN_UNSTRUCT_WGS_MAXITS", "100"))
    cmfd_mesh = os.environ.get("OPENSN_UNSTRUCT_CMFD_MESH", "local_aggregation")
    aggregation_size = int(os.environ.get("OPENSN_UNSTRUCT_AGG", "16"))
    cmfd_verbose = get_bool_option("cmfd_verbose", False) or (
        os.environ.get("OPENSN_UNSTRUCT_CMFD_VERBOSE", "0") == "1"
    )
    cmfd_update_scheme = os.environ.get("OPENSN_UNSTRUCT_CMFD_UPDATE_SCHEME", "1") == "1"
    cmfd_update_wgs_atol = float(os.environ.get("OPENSN_UNSTRUCT_CMFD_UPDATE_WGS_ATOL", "1.0e-4"))
    cmfd_update_wgs_maxits = int(os.environ.get("OPENSN_UNSTRUCT_CMFD_UPDATE_WGS_MAXITS", "5"))
    cmfd_atol = float(os.environ.get("OPENSN_UNSTRUCT_CMFD_ATOL", "1.0e-10"))
    cmfd_maxits = int(os.environ.get("OPENSN_UNSTRUCT_CMFD_MAXITS", "100"))
    cmfd_pi_maxits = int(os.environ.get("OPENSN_UNSTRUCT_CMFD_PI_MAXITS", "1"))
    cmfd_pi_ktol = float(os.environ.get("OPENSN_UNSTRUCT_CMFD_PI_KTOL", "1.0e-10"))
    cmfd_coarse_solver = os.environ.get("OPENSN_UNSTRUCT_CMFD_COARSE_SOLVER", "auto")
    cmfd_correction_max_attempts = int(
        os.environ.get("OPENSN_UNSTRUCT_CMFD_CORRECTION_MAX_ATTEMPTS", "10")
    )
    cmfd_correction_min_damping = float(
        os.environ.get("OPENSN_UNSTRUCT_CMFD_CORRECTION_MIN_DAMPING", "1.0e-4")
    )
    cmfd_negative_flux_tolerance = float(
        os.environ.get("OPENSN_UNSTRUCT_CMFD_NEGATIVE_FLUX_TOLERANCE", "1.0e-6")
    )
    cmfd_petsc_options = os.environ.get(
        "OPENSN_UNSTRUCT_CMFD_PETSC_OPTIONS",
        "",
    )

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
                    n_polar=n_polar, n_azimuthal=n_azimuthal, scattering_order=1
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
        },
    )

    cmfd = CMFDAcceleration(
        problem=phys,
        coarse_mesh=cmfd_mesh,
        aggregation_size=aggregation_size,
        update_scheme=cmfd_update_scheme,
        update_wgs_abs_tol=cmfd_update_wgs_atol,
        update_wgs_max_its=cmfd_update_wgs_maxits,
        relaxation=1.0,
        l_abs_tol=cmfd_atol,
        max_iters=cmfd_maxits,
        pi_max_its=cmfd_pi_maxits,
        pi_k_tol=cmfd_pi_ktol,
        coarse_solver_policy=cmfd_coarse_solver,
        correction_max_attempts=cmfd_correction_max_attempts,
        correction_min_damping=cmfd_correction_min_damping,
        negative_flux_tolerance=cmfd_negative_flux_tolerance,
        verbose=cmfd_verbose,
        petsc_options=cmfd_petsc_options,
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        acceleration=cmfd,
        max_iters=400,
        k_tol=1.0e-7,
    )
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    sweeps = k_solver.GetNumSweeps()
    if rank == 0:
        print(f"Python k-eigenvalue: {k}")
        print(f"Python sweeps: {sweeps}")

    uses_default_benchmark = (
        base_mesh == "triangle_mesh2x2_cuts.obj"
        and abs(xy_scale - 7.0) < 1.0e-12
        and n_layers == 12
        and n_polar == 6
        and n_azimuthal == 12
        and cmfd_mesh == "local_aggregation"
        and aggregation_size == 16
    )
    if uses_default_benchmark:
        expected_k = 1.0441134
        if abs(k - expected_k) > 5.0e-5:
            raise RuntimeError(f"Expected k near {expected_k}, got {k}")
        if sweeps > 150:
            raise RuntimeError(f"Expected CMFD benchmark to use no more than 150 sweeps, got {sweeps}")
