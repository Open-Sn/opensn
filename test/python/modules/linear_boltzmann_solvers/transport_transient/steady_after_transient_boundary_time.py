#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A problem setter keeps time-dependent boundaries evaluated at the problem time.

A 1 cm pure-absorber slab has no source other than an isotropic boundary on zmin. A transient runs
to t = 2 s, then SetSteadyStateMode() rebuilds the problem's runtime objects and a steady-state
source solve follows at the same problem time.

1. The boundary is active only on [0, 1] s. At t = 2 s it is off, the problem has no source, and
   the exact steady solution is phi = 0. INACTIVE_WINDOW_PHI_MAX is max |phi|.
2. The boundary is active on [0, 10] s. At t = 2 s it is on, so phi must be positive (a control
   that the boundary path is exercised). ACTIVE_WINDOW_PHI_POSITIVE is 1 if min phi > 0.
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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def phi_extreme(phys, operation):
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType(operation)
    field_interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    field_interp.AddFieldFunction(phys.GetScalarFluxFieldFunction()[0])
    field_interp.Execute()
    return field_interp.GetValue()


def steady_after_transient(grid, xs, end_time):
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        time_dependent=True,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 100,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {
                "name": "zmin",
                "type": "isotropic",
                "group_strength": [1.0],
                "start_time": 0.0,
                "end_time": end_time,
            }
        ],
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    transient = TransientSolver(problem=phys, dt=0.5, theta=1.0, stop_time=2.0,
                                initial_state="zero")
    transient.Initialize()
    transient.Execute()

    phys.SetSteadyStateMode()
    steady = SteadyStateSourceSolver(problem=phys)
    steady.Initialize()
    steady.Execute()
    return phys


if __name__ == "__main__":
    grid = OrthogonalMeshGenerator(node_sets=[[i * 0.1 for i in range(11)]]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0, 1.0)

    inactive = steady_after_transient(grid, xs, end_time=1.0)
    inactive_phi_max = max(abs(phi_extreme(inactive, "max")), abs(phi_extreme(inactive, "min")))

    active = steady_after_transient(grid, xs, end_time=10.0)
    active_positive = phi_extreme(active, "min") > 0.0

    if rank == 0:
        print(f"INACTIVE_WINDOW_PHI_MAX {inactive_phi_max:.6e}")
        print(f"ACTIVE_WINDOW_PHI_POSITIVE {int(active_positive)}")
