#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RZ k-eigenvalue test with vacuum boundaries and gmsh orthogonal mesh.

Cylinder with r=0.2 m and h=0.1 m, axisymmetric (rmin reflecting), vacuum
boundaries at rmax/zmin/zmax.
K_EIGEN = 0.0900968 is a historical baseline recorded in c51829ef2 as a 3D
transport result; its generating input and refinement data are not available here.
The acceptance checks also require AAH/CBC parity and global particle conservation,
P / k = A + L, for both solves. Diffusion is not an accuracy reference for this
pure absorber with optical radius 0.2 and optical height 0.1.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, PowerIterationKEigenSolver


if __name__ == "__main__":
    mesh_file = "../../../../assets/mesh/rz_disk_struct_quad.msh"
    meshgen = FromFileMeshGenerator(filename=mesh_file, coord_sys="cylindrical")
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    rmin = 0.0
    rmax = 0.2
    zmin = 0.0
    zmax = 0.1

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/simple_fissile_1g.xs")

    quad = GLCProductQuadrature2DRZ(n_polar=16, n_azimuthal=32, scattering_order=0)
    bcs = [
        {"name": "rmin", "type": "reflecting"},
        {"name": "rmax", "type": "vacuum"},
        {"name": "zmin", "type": "vacuum"},
        {"name": "zmax", "type": "vacuum"},
    ]

    problem_params = dict(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_max_its": 200,
                "l_abs_tol": 1.0e-10,
                "gmres_restart_interval": 100,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=bcs,
        options={
            "use_precursors": False,
            "verbose_inner_iterations": False,
        },
    )

    phys = DiscreteOrdinatesCurvilinearProblem(sweep_type="CBC", **problem_params)
    k_tolerance = 1.0e-8
    k_solver = PowerIterationKEigenSolver(problem=phys, k_tol=k_tolerance)
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()

    aah_problem = DiscreteOrdinatesCurvilinearProblem(sweep_type="AAH", **problem_params)
    aah_solver = PowerIterationKEigenSolver(problem=aah_problem, k_tol=k_tolerance)
    aah_solver.Initialize()
    aah_solver.Execute()
    k_aah = aah_solver.GetEigenvalue()
    # Each solve estimates its relative eigenvalue error below k_tolerance.
    # Add the two stopping-error budgets; the inner residual tolerance is 1e-10.
    parity_tol = 2.0 * k_tolerance
    parity_error = abs(k - k_aah) / abs(k_aah)
    parity_ok = parity_error <= parity_tol

    # The eigenvalue solver scales fission production by 1/k in this table.
    # Vacuum outer boundaries and the zero-area axis give zero physical inflow.
    # Use a relative conservation residual, independent of eigenvector normalization.
    balance_error = max(
        abs(table["balance"]) / abs(table["production_rate"])
        for table in (k_solver.ComputeBalanceTable(), aah_solver.ComputeBalanceTable())
    )
    # Budget: two k stopping tolerances + two inner tolerances + reduction roundoff
    # (2e-8 + 2e-10 + 1e-12 < 2.1e-8). Spatial/angular truncation does not relax balance.

    if rank == 0:
        print(f"K_EIGEN {k:.12e}")
        print(f"K_AAH {k_aah:.12e}")
        print(f"K_PARITY_ERROR {parity_error:.12e}")
        print(f"K_PARITY_OK {int(parity_ok)}")
        print(f"K_BALANCE_ERROR {balance_error:.12e}")
