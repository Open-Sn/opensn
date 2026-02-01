#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Leakage sanity: 1D vacuum boundaries with source removal.

Initialize from a steady-state source solve, then remove the source and
advance one transient step. The scalar flux should decrease with leakage.
LEAKAGE_DECAY_PASS is 1 if phi decreases and remains non-negative.
"""

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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def max_phi(phys):
    fflist = phys.GetScalarFieldFunctionList()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType("max")
    field_interp.SetLogicalVolume(monitor_volume)
    field_interp.AddFieldFunction(fflist[0])
    field_interp.Initialize()
    field_interp.Execute()
    return field_interp.GetValue()


if __name__ == "__main__":
    dx = 2.0 / 40
    nodes = [i * dx for i in range(40 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 1.0
    Q = 2.0
    dt = 0.05

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0, 1.0)

    source = VolumetricSource(block_ids=[0], group_strength=[Q], start_time=0.0, end_time=10.0)

    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        boundary_conditions=[
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    steady = SteadyStateSourceSolver(problem=phys)
    steady.Initialize()
    steady.Execute()

    phi0 = max_phi(phys)

    phys.SetVolumetricSources(clear_volumetric_sources=True)

    solver = TransientSolver(problem=phys, dt=dt, theta=1.0, stop_time=dt, initial_state="existing")
    solver.Initialize()
    solver.Execute()

    phi1 = max_phi(phys)
    pass_flag = 1 if (phi1 >= 0.0 and phi1 < phi0) else 0

    if rank == 0:
        print(f"LEAKAGE_DECAY_PASS {pass_flag}")
