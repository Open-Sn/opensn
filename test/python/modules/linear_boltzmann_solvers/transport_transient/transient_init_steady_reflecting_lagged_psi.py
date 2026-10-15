#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Entering transient mode keeps the lagged angular flux of opposing reflecting boundaries.

A 40 cm slab with reflecting boundaries on both ends is an infinite medium. For a pure absorber
with a uniform source, the exact (and discrete PWLD) solution is the constant phi = Q / sigma_a.
One end of the slab is a sweep cycle, so its incoming angular flux is lagged from the previous
sweep. Each solve does a single classic Richardson sweep (l_max_its = 1), so the lagged flux is
the only coupling between the two ends:

1. The first steady solve starts from a zero lagged flux and leaves a boundary layer. The slab is
   about 120 mean free paths thick for the most oblique direction, so the outgoing flux at the far
   end is the equilibrium value to round-off, and the second steady solve is exact.
2. SetTimeDependentMode() rebuilds the runtime objects. The lagged flux must survive it.
3. One backward-Euler step from the equilibrium state with the same source must stay at
   Q / sigma_a. If the lagged flux were reset to zero, the single sweep of that step would see a
   vacuum-like incoming flux at the lagged end and phi would drop there by O(1).

STEADY_DEVIATION and TRANSIENT_DEVIATION are max |phi - Q/sigma_a| / (Q/sigma_a) over the slab.
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


def phi_extreme(phys, operation):
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType(operation)
    field_interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    field_interp.AddFieldFunction(phys.GetScalarFluxFieldFunction()[0])
    field_interp.Execute()
    return field_interp.GetValue()


def max_deviation(phys, phi_exact):
    return max(abs(phi_extreme(phys, "max") - phi_exact),
               abs(phi_extreme(phys, "min") - phi_exact)) / phi_exact


if __name__ == "__main__":
    length = 40.0
    num_cells = 80
    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 1.0
    Q = 2.0
    dt = 1.0
    phi_exact = Q / sigma_t

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0, 1.0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-14,
                "l_max_its": 1,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[Q])],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )

    steady = SteadyStateSourceSolver(problem=phys)
    steady.Initialize()
    for _ in range(2):
        steady.Execute()
    steady_deviation = max_deviation(phys, phi_exact)

    phys.SetTimeDependentMode()

    transient = TransientSolver(
        problem=phys, dt=dt, theta=1.0, stop_time=dt, initial_state="existing"
    )
    transient.Initialize()
    transient.Execute()
    transient_deviation = max_deviation(phys, phi_exact)

    if rank == 0:
        print(f"STEADY_DEVIATION {steady_deviation:.6e}")
        print(f"TRANSIENT_DEVIATION {transient_deviation:.6e}")
