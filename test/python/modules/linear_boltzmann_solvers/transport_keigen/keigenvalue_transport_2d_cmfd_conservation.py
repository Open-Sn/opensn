#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""A frozen CMFD operator must conserve internal leakage for every trial vector.

The model has two unit-area cells with reflecting boundaries, and CMFD condenses its two
fine energy groups to one. Each material has the same absorption in both groups, so every
column sum of the CMFD loss matrix must equal the absorption of that column's cell, whatever
the fine spectrum after the sweep. This is an exact neutron-balance identity, not a computed
gold value. The coarse eigensolve starts unshifted, and this small coarse problem converges
before a Wielandt shift is applied, so every printed matrix is the loss matrix itself. Its
column sums must equal the fuel and water absorption. The two materials have contrasting
totals, so a face coupling that uses a flux-weighted diffusion coefficient on one side and an
unweighted one on the other breaks these identities by more than 1e-3. The two-rank variant
puts the interface on an MPI boundary.
"""
import os
import re
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    MPIBarrier = MPI.COMM_WORLD.Barrier

    def MPIAllReduce(value, operation):
        assert operation == "max"
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import (
        DiscreteOrdinatesProblem, PowerIterationKEigenSolver, CMFDAcceleration)
    from pyopensn.logvol import RPPLogicalVolume


def write_xs(filename, totals, absorption, fissile):
    with open(filename, "w") as stream:
        stream.write("NUM_GROUPS 2\nNUM_MOMENTS 1\nSIGMA_T_BEGIN\n")
        for g, total in enumerate(totals):
            stream.write(f"{g} {total}\n")
        stream.write("SIGMA_T_END\n")
        if fissile:
            stream.write("NU_SIGMA_F_BEGIN\n0 0.5\n1 0.5\nNU_SIGMA_F_END\n")
            stream.write("NU_BEGIN\n0 2.5\n1 2.5\nNU_END\n")
            stream.write("CHI_BEGIN\n0 1.0\n1 0.0\nCHI_END\n")
        stream.write("TRANSFER_MOMENTS_BEGIN\n")
        for g, total in enumerate(totals):
            stream.write(f"M_GFROM_GTO_VAL 0 {g} {g} {total - absorption}\n")
        stream.write("TRANSFER_MOMENTS_END\n")


if __name__ == "__main__":
    if rank == 0:
        write_xs("cmfd_conservation_fuel.xs", [1.0, 10.0], 0.4, True)
        write_xs("cmfd_conservation_water.xs", [2.0, 20.0], 0.02, False)
    MPIBarrier()
    maximum_error = 0.0
    for closure in ["partial"]:
        grid = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0, 2.0], [0.0, 1.0]]).Execute()
        grid.SetUniformBlockID(0)
        grid.SetBlockIDFromLogicalVolume(
            RPPLogicalVolume(xmin=-1.0, xmax=1.0, infy=True, infz=True), 1, True)
        fuel, water = MultiGroupXS(), MultiGroupXS()
        fuel.LoadFromOpenSn("cmfd_conservation_fuel.xs")
        water.LoadFromOpenSn("cmfd_conservation_water.xs")
        problem = DiscreteOrdinatesProblem(
            mesh=grid, num_groups=2,
            groupsets=[dict(
                groups_from_to=(0, 1),
                angular_quadrature=GLCProductQuadrature2DXY(
                    n_polar=4, n_azimuthal=8, scattering_order=0),
                inner_linear_method="classic_richardson", l_max_its=1, l_abs_tol=1e-12)],
            xs_map=[dict(block_ids=[1], xs=fuel), dict(block_ids=[0], xs=water)],
            boundary_conditions=[dict(name=name, type="reflecting")
                                 for name in ["xmin", "xmax", "ymin", "ymax"]],
            options=dict(verbose_inner_iterations=False, verbose_outer_iterations=False),
            sweep_type="AAH")
        prefix = f"conservation_{closure}_"
        matrix_path = f"cmfd_conservation_matrix_{closure}.txt"
        acceleration = CMFDAcceleration(
            problem=problem, name=prefix, coarse_mesh="identity", coarse_groups=1,
            petsc_options=f"-{prefix}ksp_view_mat ascii:{matrix_path}")
        solver = PowerIterationKEigenSolver(problem=problem, acceleration=acceleration,
                                            max_iters=1, k_tol=1e-10)
        solver.Initialize()
        solver.Execute()
        MPIBarrier()
        error = 0.0
        if rank == 0:
            with open(matrix_path) as stream:
                matrices = stream.read().split("Mat Object:")[1:]
            assert len(matrices) > 0, "Require at least one printed CMFD matrix"
            for matrix in matrices:
                column_sums = [0.0, 0.0]
                rows = re.findall(r"row (\d+):([^\n]+)", matrix)
                assert len(rows) == 2, "Require the complete two-cell matrix"
                for row, entries in rows:
                    assert int(row) in [0, 1]
                    for col, value in re.findall(r"\((\d+), ([\d.eE+-]+)\)", entries):
                        column_sums[int(col)] += float(value)
                error = max(error, abs(column_sums[0] - 0.4), abs(column_sums[1] - 0.02))
            # PETSc's default ASCII view rounds to six significant digits. The tolerance allows
            # that output rounding, not a physical balance error. The original defect exceeds 1e-3.
        error = MPIAllReduce(error, "max")
        assert error < 5e-6, (closure, error)
        maximum_error = max(maximum_error, error)
        MPIBarrier()
    if rank == 0:
        print(f"CMFD_CONSERVATION_COLUMN_SUM_ERROR {maximum_error:.12e}")
