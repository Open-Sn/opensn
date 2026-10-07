#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D spherical discrete-ordinates verification against exact solutions.

Solid sphere of radius R = 1 unless noted (the center needs no boundary condition).
  A. Void (sigma_t = 0) with isotropic incident angular flux X = 2 on the surface: psi = X,
     so phi = 2. In the center cell the final direction mu = +1 has no inflow (zero area) and
     no absorption; its flux must come from the angular recursion.
  B. Hollow void sphere, r in [0.5, 1], with isotropic incident flux X = 2 on both surfaces:
     phi = 2 (the inner surface receives its boundary flux).
  C. Uniform medium (sigma_t = 1, c = 0.5, q = 1) with isotropic incident flux equal to the
     infinite-medium flux q / sigma_a = 2, and with a reflecting surface: phi = 2 exactly, which
     requires the angular redistribution to preserve an isotropic flux. Also checked with SINGLE
     angle aggregation and with linearly anisotropic scattering (P1). The balance-table production
     must equal q 4 pi R^3 / 3.
  D. Pure absorber (sigma_t = 1, R = 1, q = 1, vacuum surface): the escape fraction converges to
     the first-flight escape probability of a sphere of optical radius tau,
       P_esc = 3 / (8 tau^3) [2 tau^2 - 1 + (1 + 2 tau) exp(-2 tau)],
     as the number of directions increases (S8, S16, S32 on 40 cells).
  E. Two-group problem with downscattering: one groupset and one groupset per group agree.
  F. k-eigenvalue of a sphere with a reflecting surface: k = nu sigma_f / sigma_a.
  G. A point source at radius r represents a shell: the production rate equals its strength.

Checks: SPH_* values printed below.
"""

import math
import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource, PointSource
    from pyopensn.aquad import GLProductQuadrature1DSpherical
    from pyopensn.solver import (DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver,
                                 PowerIterationKEigenSolver)


def sphere_mesh(r_min, r_max, n_cells):
    nodes = [r_min + (r_max - r_min) * i / n_cells for i in range(n_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes], coord_sys="spherical").Execute()
    grid.SetUniformBlockID(0)
    return grid


def make_problem(grid, xs, num_groups, bcs, sources=None, n_polar=16, scattering_order=0,
                 agg="azimuthal", group_ranges=None, point_sources=None):
    quad = GLProductQuadrature1DSpherical(n_polar=n_polar, scattering_order=scattering_order)
    return DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": r,
                "angular_quadrature": quad,
                "angle_aggregation_type": agg,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 300,
            }
            for r in (group_ranges or [(0, num_groups - 1)])
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=sources or [],
        point_sources=point_sources or [],
        boundary_conditions=bcs,
        options={"verbose_inner_iterations": False, "max_ags_iterations": 100,
                 "ags_tolerance": 1.0e-12},
    )


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return np.array(problem.GetPhiNewLocal()), solver.ComputeBalanceTable()


def max_dev(phi, exact):
    local = float(np.max(np.abs(phi - exact))) if phi.size else 0.0
    return MPIAllReduce(local, "max") / exact


def one_group(sigma_t, c):
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, c)
    return xs


def write_xs(filename, text):
    if rank == 0:
        with open(filename, "w") as f:
            f.write(text)
    MPIBarrier()


if __name__ == "__main__":
    expected_ranks = int(globals().get("rank_count", 2))
    if size != expected_ranks:
        sys.exit(f"Incorrect number of processors. Expected {expected_ranks} but got {size}.")

    def iso(name, value):
        return {"name": name, "type": "isotropic", "group_strength": [value]}

    # A. Void sphere
    phi, _ = solve(make_problem(sphere_mesh(0.0, 1.0, 20), one_group(0.0, 0.0), 1,
                                [iso("zmax", 2.0)]))
    void_dev = max_dev(phi, 2.0)

    # B. Hollow void sphere
    phi, _ = solve(make_problem(sphere_mesh(0.5, 1.0, 20), one_group(0.0, 0.0), 1,
                                [iso("zmin", 2.0), iso("zmax", 2.0)]))
    hollow_dev = max_dev(phi, 2.0)

    # C. Uniform medium
    medium = one_group(1.0, 0.5)
    src = [VolumetricSource(block_ids=[0], group_strength=[1.0])]
    phi, table = solve(make_problem(sphere_mesh(0.0, 1.0, 20), medium, 1, [iso("zmax", 2.0)],
                                    src))
    uniform_dev = max_dev(phi, 2.0)
    production_err = abs(table["production_rate"] / (4.0 * math.pi / 3.0) - 1.0)
    phi, table = solve(make_problem(sphere_mesh(0.0, 1.0, 20), medium, 1,
                                    [{"name": "zmax", "type": "reflecting"}], src))
    reflecting_dev = max_dev(phi, 2.0)
    reflecting_balance = abs(table["balance"]) / table["production_rate"]
    phi, _ = solve(make_problem(sphere_mesh(0.0, 1.0, 20), medium, 1, [iso("zmax", 2.0)], src,
                                agg="single"))
    single_dev = max_dev(phi, 2.0)
    write_xs("sph_p1.xs", "NUM_GROUPS 1\nNUM_MOMENTS 2\n\nSIGMA_T_BEGIN\n0 1.0\nSIGMA_T_END\n\n"
             "TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.5\nM_GFROM_GTO_VAL 1 0 0 0.2\n"
             "TRANSFER_MOMENTS_END\n")
    p1 = MultiGroupXS()
    p1.LoadFromOpenSn("sph_p1.xs")
    phi, _ = solve(make_problem(sphere_mesh(0.0, 1.0, 20), p1, 1, [iso("zmax", 2.0)], src,
                                scattering_order=1))
    p1_dev = max_dev(phi[0::2], 2.0)

    # D. Escape probability of a pure-absorber sphere
    tau = 1.0
    p_esc = 3.0 / (8.0 * tau**3) * (2.0 * tau**2 - 1.0 + (1.0 + 2.0 * tau) * math.exp(-2.0 * tau))
    esc_err = {}
    for n_polar in (8, 16, 32):
        _, table = solve(make_problem(sphere_mesh(0.0, 1.0, 40), one_group(1.0, 0.0), 1,
                                      [{"name": "zmax", "type": "vacuum"}], src, n_polar=n_polar))
        esc_err[n_polar] = abs(table["outflow_rate"] / table["production_rate"] / p_esc - 1.0)
    esc_order = math.log(esc_err[8] / esc_err[32]) / math.log(4.0)

    # G. A point source at radius r represents a spherical shell; production equals its strength
    _, table = solve(make_problem(sphere_mesh(0.0, 1.0, 20), one_group(1.0, 0.0), 1,
                                  [{"name": "zmax", "type": "vacuum"}],
                                  point_sources=[PointSource(location=[0.0, 0.0, 0.425],
                                                             strength=[3.0])]))
    point_err = abs(table["production_rate"] / 3.0 - 1.0)

    # E. One groupset vs one groupset per group
    write_xs("sph_2g.xs", "NUM_GROUPS 2\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n0 1.0\n1 1.5\n"
             "SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.4\n"
             "M_GFROM_GTO_VAL 0 0 1 0.3\nM_GFROM_GTO_VAL 0 1 1 0.8\nTRANSFER_MOMENTS_END\n")
    two_group = MultiGroupXS()
    two_group.LoadFromOpenSn("sph_2g.xs")
    src2 = [VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])]
    phi_a, _ = solve(make_problem(sphere_mesh(0.0, 1.0, 20), two_group, 2, [], src2))
    phi_b, _ = solve(make_problem(sphere_mesh(0.0, 1.0, 20), two_group, 2, [], src2,
                                  group_ranges=[(0, 0), (1, 1)]))
    gs_diff = MPIAllReduce(float(np.max(np.abs(phi_a - phi_b) / np.abs(phi_a))), "max")

    # F. k-eigenvalue with a reflecting surface: k = nu sigma_f / sigma_a
    write_xs("sph_fissile.xs", "NUM_GROUPS 1\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n0 1.0\n"
             "SIGMA_T_END\n\nSIGMA_F_BEGIN\n0 0.1\nSIGMA_F_END\n\nNU_BEGIN\n0 2.5\nNU_END\n\n"
             "CHI_BEGIN\n0 1.0\nCHI_END\n\nTRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.7\n"
             "TRANSFER_MOMENTS_END\n")
    fissile = MultiGroupXS()
    fissile.LoadFromOpenSn("sph_fissile.xs")
    k_problem = make_problem(sphere_mesh(0.0, 1.0, 20), fissile, 1,
                             [{"name": "zmax", "type": "reflecting"}])
    k_solver = PowerIterationKEigenSolver(problem=k_problem, k_tol=1.0e-12)
    k_solver.Initialize()
    k_solver.Execute()
    k_err = abs(k_solver.GetEigenvalue() / (2.5 * 0.1 / 0.3) - 1.0)

    MPIBarrier()
    if rank == 0:
        for f in ("sph_p1.xs", "sph_2g.xs", "sph_fissile.xs"):
            os.remove(f)
        print(f"SPH_VOID_DEV {void_dev:.6e}")
        print(f"SPH_HOLLOW_DEV {hollow_dev:.6e}")
        print(f"SPH_UNIFORM_DEV {uniform_dev:.6e}")
        print(f"SPH_PRODUCTION_ERR {production_err:.6e}")
        print(f"SPH_REFLECTING_DEV {reflecting_dev:.6e}")
        print(f"SPH_REFLECTING_BALANCE {reflecting_balance:.6e}")
        print(f"SPH_SINGLE_AGG_DEV {single_dev:.6e}")
        print(f"SPH_P1_DEV {p1_dev:.6e}")
        print(f"SPH_ESCAPE_ERR_S8 {esc_err[8]:.6e}")
        print(f"SPH_ESCAPE_ERR_S32 {esc_err[32]:.6e}")
        print(f"SPH_ESCAPE_ORDER {esc_order:.4f}")
        print(f"SPH_POINT_PROD_ERR {point_err:.6e}")
        print(f"SPH_GROUPSET_DIFF {gs_diff:.6e}")
        print(f"SPH_KINF_ERR {k_err:.6e}")
