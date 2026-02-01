#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Zero-init transient: 1D pure absorber with a constant volumetric source.

A single Backward Euler step from zero should match the analytic
phi = Q*dt/(1 + sigma_t*dt) for a homogeneous infinite medium.
ZERO_INIT_PASS is 1 if the numeric and analytic values agree within 2%.
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
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientSolver
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
    dx = 1.0 / 10
    nodes = [i * dx for i in range(10 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 1.0
    Q = 1.5
    dt = 0.1

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0, 1.0)

    source = VolumetricSource(block_ids=[0], group_strength=[Q], start_time=0.0, end_time=1.0)

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
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = TransientSolver(problem=phys, dt=dt, theta=1.0, stop_time=dt, initial_state="zero")
    solver.Initialize()
    solver.Execute()

    phi_num = max_phi(phys)
    phi_analytic = Q * dt / (1.0 + sigma_t * dt)
    rel_err = abs(phi_num - phi_analytic) / phi_analytic
    pass_flag = 1 if rel_err < 0.02 else 0

    if rank == 0:
        print(f"ZERO_INIT_PASS {pass_flag}")
