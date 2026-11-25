#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""Heterogeneous reflecting cube with the exact solution phi = (1, 100) and k = 1.

Both materials have the same absorption in the two groups, nu-sigma_f equal to that
absorption, chi = (1, 100) / 101, and only self-scattering. A spatially constant,
isotropic angular flux with group ratio 100 therefore solves the equation of every
material and every reflecting interface exactly, so the higher flux moments are zero
and k = 1. The materials have different totals and absorptions, and the mixed-material
aggregates homogenize both of them. The solve starts from the uniform scalar flux,
whose spectrum is far from the solution, uses one sweep per CMFD update, and is
compared with the analytic field after normalization by its maximum scalar flux.

The tolerances follow from the iteration tolerances. The power iteration stops at a
relative k change of 1e-12 and a fission-source change of 1e-10, and the coarse solves are
converged to a fraction of the transport change, so the field error is below 1e-6 and the k
error below 1e-9 with a wide margin. The partition puts mixed aggregates on every rank at 1,
2, and 4 ranks.

Optional argument: energy_bin, the number of transport groups per CMFD group (1 or 2).
"""
import os
import sys

import numpy as np

if "opensn_console" not in globals():
    import argparse
    import ast

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--py", action="append", default=[])
    arguments, _ = parser.parse_known_args()
    for assignment in arguments.py:
        key, value = assignment.split("=", 1)
        globals()[key] = ast.literal_eval(value)

    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    MPIBarrier = MPI.COMM_WORLD.Barrier

    def MPIAllReduce(value, operation):
        operations = {"max": MPI.MAX, "min": MPI.MIN, "sum": MPI.SUM}
        return MPI.COMM_WORLD.allreduce(value, op=operations[operation])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import (
        OrthogonalMeshGenerator, DistributedMeshGenerator, ExtruderMeshGenerator,
        KBAGraphPartitioner)
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import (
        DiscreteOrdinatesProblem, PowerIterationKEigenSolver, CMFDAcceleration)


def write_xs(filename, totals, absorption):
    with open(filename, "w") as stream:
        stream.write("NUM_GROUPS 2\nNUM_MOMENTS 1\nSIGMA_T_BEGIN\n")
        for g, value in enumerate(totals):
            stream.write(f"{g} {value}\n")
        stream.write("SIGMA_T_END\nNU_SIGMA_F_BEGIN\n")
        for g in range(2):
            stream.write(f"{g} {absorption}\n")
        stream.write("NU_SIGMA_F_END\nNU_BEGIN\n0 2.0\n1 2.0\nNU_END\n")
        stream.write("CHI_BEGIN\n0 0.009900990099009901\n1 0.990099009900990099\nCHI_END\n")
        stream.write("TRANSFER_MOMENTS_BEGIN\n")
        for g, value in enumerate(totals):
            stream.write(f"M_GFROM_GTO_VAL 0 {g} {g} {value - absorption}\n")
        stream.write("TRANSFER_MOMENTS_END\n")


if rank == 0:
    write_xs("cmfd_cube_material_a.xs", [1.0, 2.0], 0.2)
    write_xs("cmfd_cube_material_b.xs", [2.0, 3.0], 0.4)
MPIBarrier()
nx, nz = (1, 1) if size == 1 else ((2, 1) if size == 2 else (2, 2))
grid = DistributedMeshGenerator(
    inputs=[OrthogonalMeshGenerator(node_sets=[[0.0, 0.5, 1.0, 1.5, 2.0], [0.0, 1.0]]),
            ExtruderMeshGenerator(layers=[dict(z=1.0, n=1), dict(z=2.0, n=1)])],
    partitioner=KBAGraphPartitioner(
        nx=nx, ny=1, nz=nz, xcuts=[1.0] if nx == 2 else [],
        zcuts=[1.0] if nz == 2 else [])).Execute()
grid.SetOrthogonalBoundaries()
grid.SetBlockIDFromFunction(lambda point, old: int(point.x * 2.0) % 2)
a, b = MultiGroupXS(), MultiGroupXS()
a.LoadFromOpenSn("cmfd_cube_material_a.xs")
b.LoadFromOpenSn("cmfd_cube_material_b.xs")
problem = DiscreteOrdinatesProblem(
    mesh=grid, num_groups=2,
    groupsets=[dict(
        groups_from_to=(0, 1), angular_quadrature=GLCProductQuadrature3DXYZ(
            n_polar=2, n_azimuthal=4, scattering_order=1),
        inner_linear_method="classic_richardson", l_max_its=1, l_abs_tol=1e-12)],
    xs_map=[dict(block_ids=[0], xs=a), dict(block_ids=[1], xs=b)],
    boundary_conditions=[dict(name=name, type="reflecting")
                         for name in ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]],
    options=dict(verbose_inner_iterations=False, verbose_outer_iterations=False), sweep_type="AAH")
prefix = "mixed_material_"
acceleration = CMFDAcceleration(
    problem=problem, name=prefix, coarse_mesh={"type": "local_aggregation", "size": 4},
    coarse_groups=2 // globals().get("energy_bin", 1))
solver = PowerIterationKEigenSolver(
    problem=problem, acceleration=acceleration, max_iters=500, k_tol=1e-12,
    fission_source_tol=1e-10)
solver.Initialize()
solver.Execute()
converged = solver.IsConverged()
phi = np.asarray(problem.GetPhiNewLocal()).reshape(-1, 4, 2)
scale = MPIAllReduce(float(np.max(np.abs(phi[:, 0, 0]))), "max")
phi = phi / scale
expected = np.zeros((4, 2))
expected[0] = [1.0, 100.0]
field_error = MPIAllReduce(float(np.max(np.abs(phi - expected))), "max")
k_error = abs(solver.GetEigenvalue() - 1.0)
applied = acceleration.GetNumAppliedCorrections()
if rank == 0:
    print(f"CMFD_CUBE_CONVERGED {int(converged)}")
    print(f"CMFD_CUBE_APPLIED_CORRECTIONS_POSITIVE {int(applied > 0)}")
    print(f"CMFD_CUBE_K_ERROR {k_error:.12e}")
    print(f"CMFD_CUBE_FIELD_ERROR {field_error:.12e}")
