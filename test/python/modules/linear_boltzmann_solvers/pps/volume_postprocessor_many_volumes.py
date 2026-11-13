#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VolumePostprocessor with many logical volumes of every kind on a 3D, two-group problem (3 ranks).

The volumes are a grid of RPPs, RPPs whose faces pass through cell centroids, infinite RPPs,
spheres (one with its surface through cell centroids), z-aligned and oblique RCCs, and a Boolean
volume. The postprocessor selects the cells of the volumes that have a bounding box through a
grid of cell centroids, and tests every cell for the others. Each average and integral is compared
with FieldFunctionInterpolationVolume, which tests every cell for every volume, and the integrals
over the RPP grid, which tiles the domain, must add up to the whole-domain integral.
"""

import os
import random
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import (
        BooleanLogicalVolume,
        RCCLogicalVolume,
        RPPLogicalVolume,
        SphereLogicalVolume,
    )
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import VolumePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def make_volumes(n, width):
    """Returns the logical volumes; the first m * m are the RPP grid that tiles the domain."""
    rng = random.Random(2026)
    volumes = []
    m = 5  # grid faces on cell faces (cells are 1/n wide)
    h = width / m
    for j in range(m):
        for i in range(m):
            volumes.append(
                RPPLogicalVolume(
                    xmin=i * h, xmax=(i + 1) * h, ymin=j * h, ymax=(j + 1) * h, zmin=0.0, zmax=width
                )
            )
    # Faces through cell centroids, which lie at (k + 0.5) / n.
    for k in range(0, int(width * n) - 3, 2):
        lo, hi = (k + 0.5) / n, (k + 3.5) / n
        volumes.append(RPPLogicalVolume(xmin=lo, xmax=hi, ymin=lo, ymax=width, zmin=lo, zmax=hi))
    volumes.append(RPPLogicalVolume(ymin=0.3, ymax=0.9, zmin=0.2, zmax=1.1, infx=True))
    volumes.append(RPPLogicalVolume(infx=True, infy=True, infz=True))
    for _ in range(12):
        x, y, z = (rng.uniform(0.0, width) for _ in range(3))
        volumes.append(SphereLogicalVolume(r=rng.uniform(0.2, 0.6 * width), x=x, y=y, z=z))
    c = 2.5 / n
    volumes.append(SphereLogicalVolume(r=2.0 / n, x=c, y=c, z=c))
    for _ in range(4):
        x0, y0 = rng.uniform(0.0, width), rng.uniform(0.0, width)
        volumes.append(
            RCCLogicalVolume(
                r=rng.uniform(0.2, 0.4 * width), x0=x0, y0=y0, z0=-0.5, vx=0.0, vy=0.0, vz=width
            )
        )
    for _ in range(6):
        x0, y0, z0 = (rng.uniform(0.0, width) for _ in range(3))
        vx, vy, vz = (rng.uniform(-width, width) for _ in range(3))
        volumes.append(
            RCCLogicalVolume(r=0.35 * width, x0=x0, y0=y0, z0=z0, vx=vx, vy=vy, vz=vz)
        )
    box = RPPLogicalVolume(
        xmin=0.1, xmax=width - 0.1, ymin=0.1, ymax=width - 0.1, zmin=0.1, zmax=width - 0.1
    )
    hole = SphereLogicalVolume(r=0.3 * width, x=width / 2, y=width / 2, z=width / 2)
    volumes.append(
        BooleanLogicalVolume(parts=[{"op": True, "lv": box}, {"op": False, "lv": hole}])
    )
    return volumes, m * m


if __name__ == "__main__":
    n, width = 5, 2.0
    nodes = [i / n for i in range(int(width * n) + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    source_volume = RPPLogicalVolume(xmin=0.4, xmax=1.0, ymin=0.6, ymax=1.4, zmin=0.8, zmax=1.6)
    grid.SetBlockIDFromLogicalVolume(source_volume, 1, True)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/transport_2d_cyl_2g.xs")
    num_groups = 2

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=2, n_azimuthal=4, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
            }
        ],
        xs_map=[{"block_ids": [0, 1], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[1], group_strength=[1.0, 0.5])],
        boundary_conditions=[
            {"name": b, "type": "vacuum"} for b in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
        ],
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    volumes, num_grid_volumes = make_volumes(n, width)
    flux = problem.GetScalarFluxFieldFunction()

    max_rel_diff = 0.0
    for value_type, op_type in (("avg", "avg"), ("integral", "sum")):
        pp = VolumePostprocessor(problem=problem, value_type=value_type, logical_volumes=volumes)
        pp.Execute()
        values = pp.GetValue()
        for i, volume in enumerate(volumes):
            for g in range(num_groups):
                ffi = FieldFunctionInterpolationVolume()
                ffi.SetOperationType(op_type)
                ffi.SetLogicalVolume(volume)
                ffi.AddFieldFunction(flux[g])
                ffi.Execute()
                reference = ffi.GetValue()
                rel_diff = abs(values[i][g] - reference) / abs(reference)
                max_rel_diff = max(max_rel_diff, rel_diff)
        if value_type == "integral":
            integrals = values

    whole = VolumePostprocessor(problem=problem, value_type="integral")
    whole.Execute()
    tiling_rel_err = max(
        abs(sum(integrals[i][g] for i in range(num_grid_volumes)) / whole.GetValue()[0][g] - 1.0)
        for g in range(num_groups)
    )

    if rank == 0:
        print(f"PP_NUM_VOLUMES={len(volumes)}")
        print(f"PP_MAX_REL_DIFF_VS_INTERPOLATION={max_rel_diff:.6e}")
        print(f"PP_GRID_TILING_REL_ERR={tiling_rel_err:.6e}")
