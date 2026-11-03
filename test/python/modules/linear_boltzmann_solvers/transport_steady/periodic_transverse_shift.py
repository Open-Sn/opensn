#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""Translation of a localized source translates the periodic flux in 2D/3D."""

import os
import sys

from mpi4py import MPI

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def solve(dim, shift):
    nodes = [0.0, 0.25, 0.5, 0.75, 1.0]
    mesh = OrthogonalMeshGenerator(node_sets=[nodes] * dim).Execute()
    mesh.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)
    if dim == 2:
        source_region = RPPLogicalVolume(
            xmin=shift, xmax=shift + 0.5, ymin=0.0, ymax=0.5, infz=True
        )
        quad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)
        names = ("xmin", "xmax", "ymin", "ymax")
    else:
        source_region = RPPLogicalVolume(
            xmin=0.0, xmax=0.5, ymin=0.0, ymax=0.5,
            zmin=shift, zmax=shift + 0.5
        )
        quad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)
        names = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
    source = VolumetricSource(logical_volume=source_region, group_strength=[1.0])
    problem = DiscreteOrdinatesProblem(
        mesh=mesh,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 300,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        boundary_conditions=[{"name": name, "type": "periodic"} for name in names],
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    ff = problem.GetScalarFluxFieldFunction()[0]
    interpolator = FieldFunctionInterpolationPoint()
    interpolator.AddFieldFunction(ff)
    samples = []
    for x in (0.125, 0.625):
        for y in (0.125, 0.625):
            for z in ((0.0,) if dim == 2 else (0.125, 0.625)):
                point = Vector3((x + shift) % 1.0 if dim == 2 else x,
                                y,
                                (z + shift) % 1.0 if dim == 3 else z)
                interpolator.SetPointOfInterest(point)
                interpolator.Execute()
                samples.append(interpolator.GetPointValue())
    return samples, solver.ComputeBalanceTable()["balance"]


if __name__ == "__main__":
    for dimension in (2, 3):
        baseline, baseline_balance = solve(dimension, 0.0)
        translated, translated_balance = solve(dimension, 0.5)
        error = max(abs(a - b) for a, b in zip(baseline, translated))
        if max(baseline) - min(baseline) < 0.01:
            raise RuntimeError("Localized source produced an unexpectedly uniform flux.")
        if MPI.COMM_WORLD.rank == 0:
            print(f"PeriodicTransverseError{dimension}D={error:.8e}")
            print(f"PeriodicTransverseBalance{dimension}D={baseline_balance:.8e}")
            print(f"PeriodicTransverseShiftedBalance{dimension}D={translated_balance:.8e}")
