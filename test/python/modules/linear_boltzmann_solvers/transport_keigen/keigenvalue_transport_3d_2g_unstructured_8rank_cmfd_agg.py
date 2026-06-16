#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-rank 3D 2G heterogeneous extruded-unstructured benchmark using rank-local aggregated CMFD.
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

if __name__ == "__main__":

    num_procs = 8
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    base_mesh = "triangle_mesh2x2_cuts.obj"
    xy_scale = 7.0
    n_layers = 12
    n_polar = 6
    n_azimuthal = 12
    wgs_atol = 1.0e-8
    wgs_maxits = 100
    cmfd_mesh = "local_aggregation"
    aggregation_size = 16
    cmfd_verbose = bool(globals().get("cmfd_verbose", False))
    cmfd_update_wgs_atol = 1.0e-4
    cmfd_update_wgs_maxits = 5

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
        sweep_type=globals().get("sweep_type", "AAH"),
    )

    k_tolerance = 1.0e-7

    cmfd = CMFDAcceleration(
        problem=phys,
        coarse_mesh=cmfd_mesh,
        aggregation_size=aggregation_size,
        update_wgs_abs_tol=cmfd_update_wgs_atol,
        update_wgs_max_its=cmfd_update_wgs_maxits,
        relaxation=1.0,
        balance_residual_tolerance=10.0 * k_tolerance,
        verbose=cmfd_verbose,
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        acceleration=cmfd,
        max_iters=400,
        k_tol=k_tolerance,
    )
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    sweeps = k_solver.GetNumSweeps()
    if rank == 0:
        print(f"Python k-eigenvalue: {k}")
        print(f"Python sweeps: {sweeps}")

    expected_k = 1.0441134
    if abs(k - expected_k) > 5.0e-5:
        raise RuntimeError(f"Expected k near {expected_k}, got {k}")
    if sweeps > 220:
        raise RuntimeError(
            f"Expected CMFD benchmark to use no more than 220 sweeps, got {sweeps}"
        )
