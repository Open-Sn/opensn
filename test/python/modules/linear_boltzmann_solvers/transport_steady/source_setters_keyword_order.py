#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SetVolumetricSources and SetPointSources clear the existing sources before adding the new ones,
whatever the keyword order.

Each call passes the new sources before the clear flag. A 10 cm slab with reflecting boundaries on
both ends is an infinite medium, and the material is a pure absorber (sigma_a = 1 cm^-1):

1. Volumetric: the original Q = 1 source is replaced by Q = 3, so phi = 3 / sigma_a everywhere
   (the exact solution, which PWLD reproduces). VOLUMETRIC_DEVIATION is
   max |phi - Q/sigma_a| / (Q/sigma_a).
2. Point: a S = 2 point source is replaced by S = 5. With no leakage, conservation gives an
   absorption rate sigma_a * integral(phi) = 5. POINT_BALANCE_ERROR is |absorption - S| / S.

If the clear were applied after the add, both problems would have no source and both metrics
would be 1.

Each setter also replaces its sources atomically: a replacement whose second source is invalid
(two group strengths in a one-group problem) must be rejected and leave the previous sources in
place. VOLUMETRIC_REJECTED and POINT_REJECTED are 1 if the call raised a ValueError;
VOLUMETRIC_ATOMIC_DEVIATION and POINT_ATOMIC_BALANCE_ERROR repeat the metrics above after the
rejected call (an incremental update would leave the valid first source installed: 1/3 and 0.4).
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
    from pyopensn.source import PointSource, VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def phi_operation(phys, operation):
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType(operation)
    field_interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    field_interp.AddFieldFunction(phys.GetScalarFluxFieldFunction()[0])
    field_interp.Execute()
    return field_interp.GetValue()


def rejected(action):
    try:
        action()
    except ValueError:
        return True
    return False


def volumetric_deviation(phys, phi_exact):
    return max(abs(phi_operation(phys, "max") - phi_exact),
               abs(phi_operation(phys, "min") - phi_exact)) / phi_exact


def solve(phys):
    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()
    solver.Execute()


if __name__ == "__main__":
    length = 10.0
    num_cells = 20
    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    sigma_a = 1.0
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_a, 0.0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )

    # Volumetric source replaced; the new sources come before the clear flag.
    Q = 3.0
    phys.SetVolumetricSources(
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[Q])],
        clear_volumetric_sources=True,
    )
    solve(phys)
    phi_exact = Q / sigma_a
    volumetric_dev = volumetric_deviation(phys, phi_exact)

    # A rejected replacement leaves the Q = 3 source in place.
    volumetric_rejected = rejected(lambda: phys.SetVolumetricSources(
        clear_volumetric_sources=True,
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[4.0]),
                            VolumetricSource(block_ids=[0], group_strength=[1.0, 1.0])],
    ))
    solve(phys)
    volumetric_atomic_dev = volumetric_deviation(phys, phi_exact)

    # Point source replaced, again with the new sources first. The source is inside a cell.
    S = 5.0
    phys.SetVolumetricSources(clear_volumetric_sources=True)
    phys.SetPointSources(point_sources=[PointSource(location=[0.0, 0.0, 2.25], strength=[2.0])])
    phys.SetPointSources(
        point_sources=[PointSource(location=[0.0, 0.0, 2.25], strength=[S])],
        clear_point_sources=True,
    )
    solve(phys)
    absorption = sigma_a * phi_operation(phys, "sum")
    point_balance_error = abs(absorption - S) / S

    # A rejected replacement leaves the S = 5 source in place.
    point_rejected = rejected(lambda: phys.SetPointSources(
        clear_point_sources=True,
        point_sources=[PointSource(location=[0.0, 0.0, 2.25], strength=[7.0]),
                       PointSource(location=[0.0, 0.0, 7.25], strength=[1.0, 1.0])],
    ))
    solve(phys)
    point_atomic_error = abs(sigma_a * phi_operation(phys, "sum") - S) / S

    if rank == 0:
        print(f"VOLUMETRIC_DEVIATION {volumetric_dev:.6e}")
        print(f"POINT_BALANCE_ERROR {point_balance_error:.6e}")
        print(f"VOLUMETRIC_REJECTED {int(volumetric_rejected)}")
        print(f"VOLUMETRIC_ATOMIC_DEVIATION {volumetric_atomic_dev:.6e}")
        print(f"POINT_REJECTED {int(point_rejected)}")
        print(f"POINT_ATOMIC_BALANCE_ERROR {point_atomic_error:.6e}")
