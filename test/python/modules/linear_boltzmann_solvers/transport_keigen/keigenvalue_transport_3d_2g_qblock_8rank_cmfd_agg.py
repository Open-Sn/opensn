#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-rank 3D 2G heterogeneous qblock regression comparing PI with aggregated CMFD.
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


def get_bool_option(name, default):
    if name in globals():
        return bool(globals()[name])
    return default


def get_option(name, default):
    if name in globals():
        value = globals()[name]
        if isinstance(default, bool):
            if isinstance(value, str):
                return value.lower() in ("1", "true", "yes", "on")
            return bool(value)
        if isinstance(default, int) and not isinstance(default, bool):
            return int(value)
        if isinstance(default, float):
            return float(value)
        return value
    return default


if __name__ == "__main__":

    num_procs = 8
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    N = get_option("qblock_cells_per_axis", 12)
    n_polar = get_option("qblock_n_polar", 4)
    n_azimuthal = get_option("qblock_n_azimuthal", 8)
    aggregation_size = get_option("cmfd_aggregation_size", 24)
    cmfd_verbose = get_bool_option("cmfd_verbose", False)
    k_tolerance = get_option("k_tolerance", 1.0e-5)

    xss = {}
    xss["0"] = MultiGroupXS()
    xss["0"].LoadFromOpenSn("../../../../assets/xs/xs_water_g2.xs")
    xss["1"] = MultiGroupXS()
    xss["1"].LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")
    num_groups = xss["0"].num_groups

    def make_problem():
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

        return DiscreteOrdinatesProblem(
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
                    "l_abs_tol": get_option("wgs_atol", 1.0e-6),
                    "l_max_its": get_option("wgs_maxits", 60),
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

    def run_solver(use_cmfd):
        phys = make_problem()
        acceleration = None
        if use_cmfd:
            acceleration = CMFDAcceleration(
                problem=phys,
                aggregation_size=aggregation_size,
                relaxation=1.0,
                balance_residual_tolerance=10.0 * k_tolerance,
                verbose=cmfd_verbose,
            )
        solver_options = {
            "problem": phys,
            "max_iters": get_option("solver_max_iters", 200),
            "k_tol": k_tolerance,
        }
        if acceleration is not None:
            solver_options["acceleration"] = acceleration
        solver = PowerIterationKEigenSolver(**solver_options)
        solver.Initialize()
        solver.Execute()
        return solver.GetEigenvalue(), solver.GetNumSweeps()

    pi_k, pi_sweeps = run_solver(use_cmfd=False)
    cmfd_k, cmfd_sweeps = run_solver(use_cmfd=True)
    k_difference = abs(cmfd_k - pi_k)
    cmfd_sweep_ratio = cmfd_sweeps / pi_sweeps

    if rank == 0:
        print(f"Python PI k-eigenvalue: {pi_k}")
        print(f"Python PI sweeps: {pi_sweeps}")
        print(f"Python CMFD k-eigenvalue: {cmfd_k}")
        print(f"Python CMFD sweeps: {cmfd_sweeps}")
        print(f"Python CMFD-PI k-difference: {k_difference}")
        print(f"Python CMFD sweep ratio: {cmfd_sweep_ratio}")

    expected_k = get_option("expected_k", 0.5024441)
    if abs(pi_k - expected_k) > 5.0e-5:
        raise RuntimeError(f"Expected PI k near {expected_k}, got {pi_k}")
    if abs(cmfd_k - expected_k) > 5.0e-5:
        raise RuntimeError(f"Expected CMFD k near {expected_k}, got {cmfd_k}")
    if k_difference > k_tolerance:
        raise RuntimeError(f"Expected PI and CMFD k to agree, difference is {k_difference}")
    if cmfd_sweeps >= 0.75 * pi_sweeps:
        raise RuntimeError(
            f"Expected CMFD to reduce sweeps by at least 25%; "
            f"PI sweeps={pi_sweeps}, CMFD sweeps={cmfd_sweeps}"
        )
