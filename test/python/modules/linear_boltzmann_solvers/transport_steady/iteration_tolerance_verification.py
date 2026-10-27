#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Source iteration and the pointwise across-groupset iteration must meet tight tolerances.

Both iterations stop when the pointwise change of the flux, scaled by (1 - rho) for the estimated
convergence rate rho, is below the tolerance. No fixed floor may loosen a tolerance the user asked
for: only changes at the round-off level of a relative change pass regardless.

  A. Classic Richardson, one group, sigma_t = 1, c = 0.9, q = 1, l_abs_tol = 1e-12, in a 3D
     infinite medium (reflecting on every face): phi = q / (sigma_t (1 - c)) = 10.
  B. Three groupsets (one group each) with down- and upscattering in a 2D infinite medium, with
     the pointwise across-groupset convergence check and ags_tolerance = 1e-12:
     phi = (diag(sigma_t) - S)^-1 q.

Checks:
  ITER_RICHARDSON_MAX_ERR, ITER_AGS_POINTWISE_MAX_ERR: max |phi / phi_exact - 1|.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

SIGMA_T = [1.0, 1.0, 1.0]
TRANSFER = [[0.5, 0.2, 0.1], [0.2, 0.5, 0.2], [0.1, 0.2, 0.6]]  # TRANSFER[g_to][g_from]
Q = [1.0, 0.0, 0.0]


def get_option(name, default):
    return globals()[name] if name in globals() else default


def exact_three_group_flux():
    n = len(Q)
    a = [[(SIGMA_T[i] if i == j else 0.0) - TRANSFER[i][j] for j in range(n)] + [Q[i]]
         for i in range(n)]
    for k in range(n):
        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            for j in range(k, n + 1):
                a[i][j] -= factor * a[k][j]
    phi = [0.0] * n
    for i in reversed(range(n)):
        phi[i] = (a[i][n] - sum(a[i][j] * phi[j] for j in range(i + 1, n))) / a[i][i]
    return phi


def max_relative_error(problem, exact):
    error = 0.0
    for g, phi_exact in enumerate(exact):
        for operation in ("max", "min"):
            interp = FieldFunctionInterpolationVolume()
            interp.SetOperationType(operation)
            interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
            interp.AddFieldFunction(problem.GetScalarFluxFieldFunction()[g])
            interp.Execute()
            error = max(error, abs(interp.GetValue() / phi_exact - 1.0))
    return error


def solve(grid, num_groups, groupsets, xs, q, faces, options):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=groupsets,
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=q)],
        boundary_conditions=[{"name": name, "type": "reflecting"} for name in faces],
        options=options,
        sweep_type=get_option("sweep_type", "AAH"),
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return problem


if __name__ == "__main__":
    # A. Classic Richardson
    nodes = [i / 6 for i in range(7)]
    grid_3d = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid_3d.SetUniformBlockID(0)
    one_group = MultiGroupXS()
    one_group.CreateSimpleOneGroup(1.0, 0.9)
    problem_a = solve(
        grid_3d, 1,
        [{"groups_from_to": (0, 0),
          "angular_quadrature": GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4,
                                                          scattering_order=0),
          "inner_linear_method": "classic_richardson", "l_abs_tol": 1.0e-12,
          "l_max_its": 2000}],
        one_group, [1.0], ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"),
        {"verbose_inner_iterations": False})
    err_a = max_relative_error(problem_a, [10.0])

    # B. Pointwise across-groupset iteration
    xs_file = "iteration_tolerance_verification.xs"
    if rank == 0:
        with open(xs_file, "w") as f:
            f.write("NUM_GROUPS 3\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g, value in enumerate(SIGMA_T):
                f.write(f"{g} {value}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(3):
                for g_from in range(3):
                    f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {TRANSFER[g_to][g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    three_group = MultiGroupXS()
    three_group.LoadFromOpenSn(xs_file)
    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)
    grid_2d = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid_2d.SetUniformBlockID(0)
    quad_2d = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)
    problem_b = solve(
        grid_2d, 3,
        [{"groups_from_to": (g, g), "angular_quadrature": quad_2d,
          "inner_linear_method": "petsc_gmres", "l_abs_tol": 1.0e-13, "l_max_its": 200}
         for g in range(3)],
        three_group, Q, ("xmin", "xmax", "ymin", "ymax"),
        {"verbose_inner_iterations": False, "max_ags_iterations": 2000,
         "ags_tolerance": 1.0e-12, "ags_convergence_check": "pointwise"})
    err_b = max_relative_error(problem_b, exact_three_group_flux())

    if rank == 0:
        print(f"ITER_RICHARDSON_MAX_ERR {err_a:.6e}")
        print(f"ITER_AGS_POINTWISE_MAX_ERR {err_b:.6e}")
