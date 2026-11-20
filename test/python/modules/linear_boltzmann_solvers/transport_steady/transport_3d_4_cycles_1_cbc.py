#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
Three-dimensional PWLD transport with vacuum and isotropic incident boundaries.

Particle conservation independently requires Q = A_eff + J_out - J_in.
A_eff integrates (sigma_t[g] - sum_h S0[h,g]) * phi[g] over volume,
with h restricted to the 21 transported groups. Scattering out of that range
counts as absorption. Internal currents cancel in the discrete PWLD balance.
CBC_CYCLE_RELATIVE_BALANCE is |Q + J_in - A_eff - J_out| / (Q + J_in).
The limit allows ten times the linear tolerance plus 1e-10 for summation.

The scalar-flux maxima in groups 0 and 19 also check agreement with the AAH
input transport_3d_4_cycles_1.py.
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    meshgen = ExtruderMeshGenerator(
        inputs=[
            FromFileMeshGenerator(
                filename="../../../../assets/mesh/square2x2_partition_cyclic3.obj"
            )
        ],
        layers=[{"z": 0.4, "n": 2},
                {"z": 0.8, "n": 2},
                {"z": 1.2, "n": 2},
                {"z": 1.6, "n": 2},
                ],  # layers
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            xcuts=[0.0],
            ycuts=[0.0], ),
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()

    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)
    vol1 = RPPLogicalVolume(xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    num_groups = 21
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("../../../../assets/xs/xs_graphite_pure.xs")

    strength = [0.0 for _ in range(num_groups)]
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)
    mg_src2 = VolumetricSource(block_ids=[2], group_strength=strength)

    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=1)

    bsrc = [0.0 for _ in range(num_groups)]
    bsrc[0] = 1.0

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 20],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        volumetric_sources=[mg_src1, mg_src2],
        boundary_conditions=[
            {"name": "zmax", "type": "isotropic", "group_strength": bsrc},
        ],
        sweep_type="CBC",
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    fflist = phys.GetScalarFluxFieldFunction()

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0])
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[19])
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")

    # Effective absorption includes removal into energy groups that are not transported.
    removal = list(xs_graphite.sigma_t)[:num_groups]
    with open("../../../../assets/xs/xs_graphite_pure.xs") as xs_file:
        for line in xs_file:
            values = line.split()
            if len(values) == 5 and values[0] == "M_GFROM_GTO_VAL":
                moment, g_from, g_to = map(int, values[1:4])
                if moment == 0 and g_from < num_groups and g_to < num_groups:
                    removal[g_from] -= float(values[4])
    absorption = []
    for g, field in enumerate(phys.GetScalarFluxFieldFunction()):
        integral = FieldFunctionInterpolationVolume()
        integral.SetOperationType("sum")
        integral.SetLogicalVolume(vol0)
        integral.AddFieldFunction(field)
        integral.Execute()
        absorption.append(removal[g] * integral.GetValue())
    balance = ss_solver.ComputeBalanceTable()
    gain = balance["production_rate"] + balance["inflow_rate"]
    residual = math.fsum((gain, -math.fsum(absorption), -balance["outflow_rate"]))
    relative_balance = abs(residual) / gain
    if rank == 0:
        print(f"CBC_CYCLE_RELATIVE_BALANCE={relative_balance:.16e}")
