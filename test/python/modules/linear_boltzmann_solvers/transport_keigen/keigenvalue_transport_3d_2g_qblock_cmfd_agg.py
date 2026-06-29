#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 2G heterogeneous qblock k-eigenvalue benchmark using aggregated CMFD acceleration.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration
    from pyopensn.logvol import RPPLogicalVolume


def get_option(name, default):
    return globals()[name] if name in globals() else default


def get_bool_option(name, default):
    value = get_option(name, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("true", "1", "yes", "on")
    return bool(value)


if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    N = 8
    L = 14.0
    dx = L / N
    nodes = [i * dx for i in range(N + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    fuel = RPPLogicalVolume(
        xmin=-1000.0,
        xmax=10.0,
        ymin=-1000.0,
        ymax=10.0,
        zmin=-1000.0,
        zmax=10.0,
    )
    grid.SetBlockIDFromLogicalVolume(fuel, 1, True)
    grid.SetOrthogonalBoundaries()

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
                    n_polar=4, n_azimuthal=8, scattering_order=1
                ),
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 100,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xss["0"]},
            {"block_ids": [1], "xs": xss["1"]},
        ],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
        ],
        options={
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        sweep_type=get_option("sweep_type", "AAH"),
    )

    k_tolerance = float(get_option("k_tolerance", 1.0e-6))

    cmfd = CMFDAcceleration(
        problem=phys,
        aggregation_size=8,
        group_aggregation_size=int(get_option("cmfd_group_aggregation_size", 1)),
        relaxation=1.0,
        balance_residual_tolerance=float(
            get_option("cmfd_balance_residual_tolerance", 10.0 * k_tolerance)
        ),
        verbose=get_bool_option("cmfd_verbose", False),
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

    expected_k = 0.5028015
    if abs(k - expected_k) > 5.0e-5:
        raise RuntimeError(f"Expected k near {expected_k}, got {k}")
