#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""CBC no-lag tetrahedral transport and forced normal/delayed fragmentation.

A uniform pure absorber with unit source and matching isotropic inflow has the
exact solution phi_g = g+1, by sigma_a phi_g = q_g. The source is isotropic;
streaming is zero and the prescribed inflow preserves its constant value.
48 groups times three face nodes times eight bytes exceed the legal
1024-byte packet limit. A 48+1 groupset
split tests group offsets and both streams. No gold
value comes from OpenSn. Iteration tolerance 1e-12 is tighter than the 1e-9
relative acceptance limit. The existing tetrahedral geometry is used in its
native coordinates, with no mesh scale assumption.
"""
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank, size = MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size

    def MPIAllReduce(value, op="max"):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

if __name__ == "__main__":
    allow_cycles = globals().get("allow_cycles", False)
    grid = FromFileMeshGenerator(
        filename="../../../../assets/mesh/cube3.2.msh",
        partitioner=KBAGraphPartitioner(nx=size, xcuts=[0.032 * p / size for p in range(1, size)]),
    ).Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()
    # Diagonal 49-group pure absorption, written in the runner's isolated working directory.
    path = "cbc_equilibrium.xs"
    if rank == 0:
        with open(path, "w") as output:
            output.write("NUM_GROUPS 49\nNUM_MOMENTS 1\nSIGMA_T_BEGIN\n")
            for g in range(49):
                output.write(f"{g} 1.0\n")
            output.write("SIGMA_T_END\n")
    MPIBarrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    quad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=8, scattering_order=0)
    problem = DiscreteOrdinatesProblem(
        mesh=grid, num_groups=49,
        groupsets=[{"groups_from_to": interval, "angular_quadrature": quad,
                    "allow_cycles": allow_cycles,
                    "inner_linear_method": "petsc_gmres", "l_abs_tol": 1e-12,
                    "l_max_its": 200} for interval in [(0, 47), (48, 48)]],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(
            block_ids=[0], group_strength=[float(g) for g in range(1, 50)])],
        boundary_conditions=[{"name": name, "type": "isotropic",
                              "group_strength": [float(g) for g in range(1, 50)]}
                             for name in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")],
        sweep_type="CBC",
        options={"verbose_inner_iterations": False, "deterministic_sweeps": True,
                 "max_mpi_message_size": 1024},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    values = list(problem.GetPhiNewLocal())
    error = MPIAllReduce(max(
        (abs(phi / (i % 49 + 1) - 1.0) for i, phi in enumerate(values)), default=0.0))
    balance = solver.ComputeBalanceTable()
    relative_balance = abs(balance["balance"]) / balance["production_rate"]
    if rank == 0:
        print(f"CBC_EQUILIBRIUM_ERROR={error:.16e}")
        print(f"CBC_EQUILIBRIUM_BALANCE={relative_balance:.16e}")
