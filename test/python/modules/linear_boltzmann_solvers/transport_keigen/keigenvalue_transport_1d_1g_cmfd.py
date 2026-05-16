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
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration

def get_option(name, default):
    return globals()[name] if name in globals() else default

if __name__ == "__main__":

    # Mesh variables
    L = 100.0
    n_cells = 50
    dx = L / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # Load cross-section data for a 1-group fissile material
    num_groups = 1
    xs_simple_fissile = MultiGroupXS()
    xs_simple_fissile.LoadFromOpenSn("../../../../assets/xs/simple_fissile.xs")

    # Angle information
    n_angles = 32
    scat_order = 0

    phys = DiscreteOrdinatesProblem(
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
        xs_map=[
            {
                "block_ids": [0],
                "xs": xs_simple_fissile,
            }
        ],
        options={
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        }
    )
    cmfd = CMFDAcceleration(
        problem=phys,
        coarse_mesh="local_aggregation",
        aggregation_size=2,
        relaxation=get_option("cmfd_relaxation", 0.1),
        adaptive_relaxation=False,
        inactive_iterations=0,
        update_scheme=False,
        l_abs_tol=1.0e-10,
        max_iters=100,
        pi_max_its=50,
        correction_max_attempts=get_option("cmfd_correction_max_attempts", 10),
        correction_min_damping=get_option("cmfd_correction_min_damping", 1.0e-4),
        negative_flux_tolerance=get_option("cmfd_negative_flux_tolerance", 1.0e-6),
        verbose=get_option("cmfd_verbose", False),
        petsc_options="-CMFDAccelerationksp_type gmres -CMFDAccelerationpc_type jacobi",
    )
    k_solver = PowerIterationKEigenSolver(
        problem=phys,
        acceleration=cmfd,
        max_iters=400,
        k_tol=1.0e-8,
    )
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()
    sweeps = k_solver.GetNumSweeps()
    if rank == 0:
        print(f"Python k-eigenvalue: {k}")
        print(f"Python sweeps: {sweeps}")

    if not get_option("cmfd_limiter_stress", False):
        expected_k = 0.9945456
        if abs(k - expected_k) > 5.0e-6:
            raise RuntimeError(f"Expected k={expected_k}, got {k}")
