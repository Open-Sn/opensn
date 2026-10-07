#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fixed-source verification against analytic and semi-analytic references.

A. Isotropic incident flux on a 1D pure-absorber slab (sigma_t = 1, length L = 2, psi_in = 1 on
   zmin, vacuum on zmax). With OpenSn's quadrature weights (sum = 1) the discrete-ordinates
   solution gives, per unit area,
     inflow     = sum_{mu_n > 0} w_n mu_n psi_in
     outflow    = sum_{mu_n > 0} w_n mu_n psi_in exp(-sigma_t L / mu_n)
     absorption = inflow - outflow.
   PWLD is compared on 8 and 32 cells; the reported orders measure spatial convergence.
B. Infinite medium (reflecting on every face) in 1D, 2D, and 3D: a two-group material with down-
   and upscattering and a uniform source q has phi = (diag(sigma_t) - S)^-1 q everywhere.
   B2 repeats the 2D case with one groupset per group, so upscattering is resolved by the
   across-groupset iteration.
C. An arbitrary boundary that returns psi_in for every direction reproduces case A exactly.
D. A point source of strength S in a 3D pure-absorber cube with vacuum boundaries: the balance
   table must report production = S and absorption + outflow = S.

Checks:
  FSV_SLAB_ABS_ERR, FSV_SLAB_OUT_ERR: relative errors of absorption and outflow on 32 cells.
  FSV_SLAB_ABS_ORDER, FSV_SLAB_OUT_ORDER: observed orders between 8 and 32 cells.
  FSV_SLAB_IN_ERR: relative inflow error (exact, independent of the mesh).
  FSV_INF_MAX_ERR_1D/2D/3D: max |phi / phi_exact - 1| over all nodes and groups.
  FSV_INF_MAX_ERR_2D_AGS: the same for the 2D case with one groupset per group.
  FSV_ARBITRARY_DIFF: max relative difference between cases C and A.
  FSV_POINT_PROD_ERR, FSV_POINT_BALANCE: |production / S - 1| and |balance| / S.
"""

import math
import os
import sys

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
    from pyopensn.math import AngularFluxFunction
    from pyopensn.aquad import (GLProductQuadrature1DSlab, GLCProductQuadrature2DXY,
                                GLCProductQuadrature3DXYZ)
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


def get_option(name, default):
    return globals()[name] if name in globals() else default


def get_bool(value):
    return value if isinstance(value, bool) else str(value).lower() in ("1", "true", "yes")


def solve(grid, quad, xs_map, num_groups, bcs, vol_srcs=None, pt_srcs=None,
          group_ranges=None):
    group_ranges = group_ranges or [(0, num_groups - 1)]
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": group_range,
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
                "gmres_restart_interval": 100,
            }
            for group_range in group_ranges
        ],
        xs_map=xs_map,
        volumetric_sources=vol_srcs or [],
        point_sources=pt_srcs or [],
        boundary_conditions=bcs,
        options={"verbose_inner_iterations": False, "max_ags_iterations": 500,
                 "ags_tolerance": 1.0e-12},
        sweep_type=get_option("sweep_type", "AAH"),
        use_gpus=get_bool(get_option("use_gpus", False)),
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return problem, solver.ComputeBalanceTable()


def mesh(dim, n, length):
    nodes = [length * i / n for i in range(n + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes] * dim).Execute()
    grid.SetUniformBlockID(0)
    return grid


if __name__ == "__main__":
    if size != 2:
        sys.exit(f"Incorrect number of processors. Expected 2 but got {size}.")

    absorber = MultiGroupXS()
    absorber.CreateSimpleOneGroup(1.0, 0.0)

    # A. Isotropic incident flux on a slab
    length = 2.0
    slab_quad = GLProductQuadrature1DSlab(n_polar=8, scattering_order=0)
    incoming = [(w, om.z) for w, om in zip(slab_quad.weights, slab_quad.omegas) if om.z > 0.0]
    inflow_exact = sum(w * mu for w, mu in incoming)
    outflow_exact = sum(w * mu * math.exp(-length / mu) for w, mu in incoming)
    absorption_exact = inflow_exact - outflow_exact

    slab_bcs = [{"name": "zmin", "type": "isotropic", "group_strength": [1.0]},
                {"name": "zmax", "type": "vacuum"}]
    errors = {}
    for n in (8, 32):
        _, table = solve(mesh(1, n, length), slab_quad, [{"block_ids": [0], "xs": absorber}], 1,
                         slab_bcs)
        errors[n] = (abs(table["absorption_rate"] / absorption_exact - 1.0),
                     abs(table["outflow_rate"] / outflow_exact - 1.0),
                     abs(table["inflow_rate"] / inflow_exact - 1.0))
        if n == 32:
            table_a = table
    abs_order = math.log(errors[8][0] / errors[32][0]) / math.log(4.0)
    out_order = math.log(errors[8][1] / errors[32][1]) / math.log(4.0)

    # C. Constant arbitrary boundary must reproduce A
    arbitrary_bcs = [{"name": "zmin", "type": "arbitrary",
                      "function": AngularFluxFunction(lambda g, n: 1.0)},
                     {"name": "zmax", "type": "vacuum"}]
    _, table_c = solve(mesh(1, 32, length), slab_quad, [{"block_ids": [0], "xs": absorber}], 1,
                       arbitrary_bcs)
    arbitrary_diff = max(abs(table_c[k] / table_a[k] - 1.0)
                         for k in ("absorption_rate", "outflow_rate", "inflow_rate"))

    # B. Two-group infinite medium in 1D, 2D, and 3D
    sigma_t = [1.0, 2.0]
    transfer = [[0.5, 0.1], [0.4, 1.5]]  # transfer[g_to][g_from]
    q = [1.0, 0.5]
    a11, a12 = sigma_t[0] - transfer[0][0], -transfer[0][1]
    a21, a22 = -transfer[1][0], sigma_t[1] - transfer[1][1]
    det = a11 * a22 - a12 * a21
    phi_exact = [(a22 * q[0] - a12 * q[1]) / det, (a11 * q[1] - a21 * q[0]) / det]

    xs_file = "fsv_two_group.xs"
    if rank == 0:
        with open(xs_file, "w") as f:
            f.write("NUM_GROUPS 2\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g in range(2):
                f.write(f"{g} {sigma_t[g]}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(2):
                for g_from in range(2):
                    f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {transfer[g_to][g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    two_group = MultiGroupXS()
    two_group.LoadFromOpenSn(xs_file)

    faces = {1: ["zmin", "zmax"], 2: ["xmin", "xmax", "ymin", "ymax"],
             3: ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]}
    quads = {1: GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
             2: GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0),
             3: GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)}
    inf_err = {}
    for dim in (1, 2, 3):
        problem, _ = solve(mesh(dim, 4 if dim < 3 else 2, 1.0), quads[dim],
                           [{"block_ids": [0], "xs": two_group}], 2,
                           [{"name": f, "type": "reflecting"} for f in faces[dim]],
                           vol_srcs=[VolumetricSource(block_ids=[0], group_strength=q)])
        phi = list(problem.GetPhiNewLocal())
        local = max((abs(phi[i] / phi_exact[i % 2] - 1.0) for i in range(len(phi))),
                    default=0.0)
        inf_err[dim] = MPIAllReduce(local, "max")

    # B2. The same 2D medium with one groupset per group: upscattering couples the groupsets
    # through the across-groupset (AGS) iteration.
    problem, _ = solve(mesh(2, 4, 1.0), quads[2], [{"block_ids": [0], "xs": two_group}], 2,
                       [{"name": f, "type": "reflecting"} for f in faces[2]],
                       vol_srcs=[VolumetricSource(block_ids=[0], group_strength=q)],
                       group_ranges=[(0, 0), (1, 1)])
    phi = list(problem.GetPhiNewLocal())
    local = max((abs(phi[i] / phi_exact[i % 2] - 1.0) for i in range(len(phi))), default=0.0)
    inf_err_ags = MPIAllReduce(local, "max")

    # D. Point source balance in a 3D pure absorber
    strength = 2.5
    _, table_d = solve(mesh(3, 4, 2.0), GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8,
                                                                  scattering_order=0),
                       [{"block_ids": [0], "xs": absorber}], 1,
                       [{"name": f, "type": "vacuum"} for f in faces[3]],
                       pt_srcs=[PointSource(location=[0.9, 1.1, 1.05], strength=[strength])])
    point_prod_err = abs(table_d["production_rate"] / strength - 1.0)
    point_balance = abs(table_d["balance"]) / strength

    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)
        print(f"FSV_SLAB_IN_ERR {errors[32][2]:.6e}")
        print(f"FSV_SLAB_ABS_ERR {errors[32][0]:.6e}")
        print(f"FSV_SLAB_OUT_ERR {errors[32][1]:.6e}")
        print(f"FSV_SLAB_ABS_ORDER {abs_order:.4f}")
        print(f"FSV_SLAB_OUT_ORDER {out_order:.4f}")
        print(f"FSV_ARBITRARY_DIFF {arbitrary_diff:.6e}")
        for dim in (1, 2, 3):
            print(f"FSV_INF_MAX_ERR_{dim}D {inf_err[dim]:.6e}")
        print(f"FSV_INF_MAX_ERR_2D_AGS {inf_err_ags:.6e}")
        print(f"FSV_POINT_PROD_ERR {point_prod_err:.6e}")
        print(f"FSV_POINT_BALANCE {point_balance:.6e}")
