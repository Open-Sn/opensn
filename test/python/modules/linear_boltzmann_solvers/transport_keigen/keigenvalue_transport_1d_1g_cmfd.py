#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D 1G k-eigenvalue smoke test using power iteration with CMFD acceleration.

Test: Final k-eigenvalue: 0.9945456
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration


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

    # Load cross-section data for a 1-group fissile material
    num_groups = 1
    xs_simple_fissile = MultiGroupXS()
    xs_simple_fissile.LoadFromOpenSn("../../../../assets/xs/simple_fissile.xs")
    xs_scaled_fissile = MultiGroupXS()
    xs_scaled_fissile.LoadFromOpenSn("../../../../assets/xs/simple_fissile.xs")
    xs_scaled_fissile.Scale(1.5)
    heterogeneous = get_option("heterogeneous", False)

    # Angle information
    n_angles = 32
    scat_order = 0
    pi_solver_max_iters = get_option("pi_solver_max_iters", 10000)
    cmfd_solver_max_iters = get_option("cmfd_solver_max_iters", 500)
    k_tolerance = 1.0e-10
    k_agreement_tolerance = 1.0e-7

    def make_problem():
        L = 100.0
        n_cells = 50
        dx = L / n_cells
        nodes = [i * dx for i in range(n_cells + 1)]
        meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
        grid = meshgen.Execute()
        grid.SetUniformBlockID(0)
        if heterogeneous:
            grid.SetBlockIDFromLogicalVolume(
                RPPLogicalVolume(zmin=50.0, zmax=100.0, infx=True, infy=True),
                1,
                True,
            )

        xs_map = [{"block_ids": [0], "xs": xs_simple_fissile}]
        if heterogeneous:
            xs_map.append({"block_ids": [1], "xs": xs_scaled_fissile})

        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=num_groups,
            groupsets=[
                {
                    "groups_from_to": (0, num_groups - 1),
                    "angular_quadrature": GLProductQuadrature1DSlab(
                        n_polar=n_angles,
                        scattering_order=scat_order
                    ),
                    "inner_linear_method": "petsc_richardson",
                    "l_max_its": 50,
                    "l_abs_tol": 1.0e-10,
                }
            ],
            xs_map=xs_map,
            options={
                "use_precursors": False,
                "verbose_inner_iterations": False,
                "verbose_outer_iterations": False,
            },
            sweep_type=get_option("sweep_type", "AAH"),
        )

    def run_power_iteration():
        phys = make_problem()
        solver = PowerIterationKEigenSolver(
            problem=phys,
            max_iters=pi_solver_max_iters,
            k_tol=k_tolerance,
        )
        solver.Initialize()
        solver.Execute()
        return solver.GetEigenvalue(), solver.GetNumPowerIterations(), solver.GetNumSweeps()

    def run_cmfd():
        phys = make_problem()
        coarse_mesh_type = get_option("cmfd_coarse_mesh", "local_aggregation")
        # In the heterogeneous case, the coarse cell [40, 55] deliberately straddles the material
        # interface at z=50, and on four ranks also a rank boundary. This exercises spatial
        # homogenization rather than merely grouping material-pure fine cells.
        if heterogeneous:
            overlay_edges = [0.0, 10.0, 20.0, 30.0, 40.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0]
        else:
            overlay_edges = [10.0 * i for i in range(11)]
        if coarse_mesh_type == "logical_volumes":
            coarse_mesh = {
                "type": "logical_volumes",
                "volumes": [
                    RPPLogicalVolume(
                        zmin=overlay_edges[i],
                        zmax=overlay_edges[i + 1],
                        infx=True,
                        infy=True,
                    )
                    for i in range(len(overlay_edges) - 1)
                ],
            }
        elif coarse_mesh_type == "cartesian":
            coarse_mesh = {"type": "cartesian", "z_edges": overlay_edges}
        elif coarse_mesh_type == "reactor_lattice":
            coarse_mesh = {"type": "reactor_lattice", "origin": [-1.0, -1.0], "pitch": [2.0, 2.0],
                           "shape": [1, 1], "axial_edges": overlay_edges}
        else:
            coarse_mesh = {"type": coarse_mesh_type,
                           "size": get_option("cmfd_aggregation_size", 2)}
        cmfd_kwargs = dict(
            problem=phys,
            coarse_mesh=coarse_mesh,
            verbose=get_option("cmfd_verbose", False),
            petsc_options="-CMFDAccelerationksp_type gmres -CMFDAccelerationpc_type jacobi",
        )
        cmfd = CMFDAcceleration(**cmfd_kwargs)
        solver = PowerIterationKEigenSolver(
            problem=phys,
            acceleration=cmfd,
            max_iters=cmfd_solver_max_iters,
            k_tol=k_tolerance,
        )
        solver.Initialize()
        solver.Execute()
        cmfd_status["converged"] = solver.IsConverged()
        cmfd_status["skipped"] = cmfd.GetNumSkippedCorrections()
        cmfd_status["coarse_cells"] = cmfd.GetNumCoarseCells()
        return solver.GetEigenvalue(), solver.GetNumPowerIterations(), solver.GetNumSweeps()

    cmfd_status = {}
    pi_k, pi_iterations, pi_sweeps = run_power_iteration()
    cmfd_k, cmfd_iterations, cmfd_sweeps = run_cmfd()
    k_difference = abs(cmfd_k - pi_k)
    cmfd_speedup = pi_sweeps / cmfd_sweeps

    if rank == 0:
        print(f"Python PI k-eigenvalue: {pi_k}")
        print(f"Python PI iterations: {pi_iterations}")
        print(f"Python PI sweeps: {pi_sweeps}")
        print(f"Python CMFD k-eigenvalue: {cmfd_k}")
        print(f"Python CMFD iterations: {cmfd_iterations}")
        print(f"Python CMFD sweeps: {cmfd_sweeps}")
        print(f"Python CMFD-PI k-difference: {k_difference}")
        print(f"Python CMFD sweep speedup: {cmfd_speedup}")
        print(f"Python CMFD converged: {int(cmfd_status['converged'])}")
        print(f"Python CMFD skipped corrections: {cmfd_status['skipped']}")
        print(f"Python CMFD coarse cells: {cmfd_status['coarse_cells']}")

    if not get_option("cmfd_limiter_stress", False):
        # The heterogeneous reference is the converged unaccelerated transport result from the
        # same discretization. CMFD changes only the iterative path.
        expected_k = 0.99470748 if heterogeneous else 0.9945456
        if abs(pi_k - expected_k) > 5.0e-6:
            raise RuntimeError(f"Expected PI k={expected_k}, got {pi_k}")
        if abs(cmfd_k - expected_k) > 5.0e-6:
            raise RuntimeError(f"Expected CMFD k={expected_k}, got {cmfd_k}")
        if k_difference > k_agreement_tolerance:
            raise RuntimeError(f"Expected PI and CMFD k to agree, difference is {k_difference}")
        if not heterogeneous and cmfd_sweeps >= 0.25 * pi_sweeps:
            raise RuntimeError(
                f"Expected CMFD to reduce sweeps by at least 4x; "
                f"PI sweeps={pi_sweeps}, CMFD sweeps={cmfd_sweeps}"
            )
