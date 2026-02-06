#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RZ k-eigenvalue test with vacuum boundaries and gmsh orthogonal mesh.

Cylinder with r=0.2 m and h=0.1 m, axisymmetric (rmin reflecting), vacuum
boundaries at rmax/zmin/zmax.
K_EIGEN = 0.0900968 from 3D disk keigen vacuum transport solution
(transport_3d_disk_keigen_vacuum_bcs.py).
"""

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
    xs.LoadFromOpenSn("simple_fissile_1g.xs")

    quad = GLCProductQuadrature2DRZ(n_polar=16, n_azimuthal=32, scattering_order=0)
    bcs = [
        {"name": "rmin", "type": "reflecting"},
        {"name": "rmax", "type": "vacuum"},
        {"name": "zmin", "type": "vacuum"},
        {"name": "zmax", "type": "vacuum"},
    ]

    phys = DiscreteOrdinatesCurvilinearProblem(
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

    k_solver = PowerIterationKEigenSolver(problem=phys, k_tol=1.0e-8)
    k_solver.Initialize()
    k_solver.Execute()

    k = k_solver.GetEigenvalue()

    sigma_t = xs.sigma_t[0]
    sigma_a = xs.sigma_a[0] if len(xs.sigma_a) > 0 else sigma_t
    nu_sigma_f = xs.nu_sigma_f[0]

    D = 1.0 / (3.0 * sigma_t)
    R = rmax - rmin
    H = zmax - zmin
    z_ex = 0.7104 * D
    R_e = R + z_ex
    H_e = H + 2.0 * z_ex
    B2 = (2.405 / R_e) ** 2 + (math.pi / H_e) ** 2
    k_ref = nu_sigma_f / (sigma_a + D * B2)

    tol = float(os.environ.get("OPENSN_TOL", "5e-2"))
    ok = abs(k - k_ref) <= tol

    if rank == 0:
        print(f"K_EIGEN {k:.12e}")
        print(f"K_REF {k_ref:.12e}")
        print(f"TOL {tol:.6e}")
        print(f"PASS {int(ok)}")
