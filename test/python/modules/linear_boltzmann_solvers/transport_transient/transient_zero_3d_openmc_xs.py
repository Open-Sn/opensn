#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Fixed-source time-dependent transport in a homogeneous cube using OpenMC-generated
# macroscopic, multigroup cross sections. The time-dependent solution is advanced to
# steady state and compared with the OpenSn steady state solution and the OpenMC
# steady state solution.
#
# OpenSn time-dependent solution: 50.7940
# OpenSn steady-state solution: 51.057722
# OpenMC steady-state solution: 50.96678

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    N = 10
    L = 10.0
    xmin = -L / 2.0
    dx = L / N
    nodes = [xmin + i * dx for i in range(N + 1)]

    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_water = MultiGroupXS()
    xs_water.LoadFromOpenMC(
        os.path.join(os.path.dirname(__file__), "xs_water.h5"), "set1", 294
    )
    num_groups = xs_water.num_groups

    strength = [0.0 for _ in range(num_groups)]
    strength[3] = 12.285
    src1 = VolumetricSource(block_ids=[0], group_strength=strength)

    pquad = GLCProductQuadrature3DXYZ(n_polar=8, n_azimuthal=16, scattering_order=1)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, num_groups - 1],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_water},
        ],
        volumetric_sources=[src1],
    )

    solver = TransientSolver(problem=phys, dt=0.01, theta=0.5, stop_time=0.1, initial_state="zero")
    solver.Initialize()
    solver.Execute()

    fflist = phys.GetScalarFluxFieldFunction()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType("max")
    field_interp.SetLogicalVolume(monitor_volume)
    field_interp.AddFieldFunction(fflist[3])
    field_interp.Initialize()
    field_interp.Execute()
    flux_max = field_interp.GetValue()

    if rank == 0:
        print(f"Max phi(0.1s) = {flux_max:.6f}")
