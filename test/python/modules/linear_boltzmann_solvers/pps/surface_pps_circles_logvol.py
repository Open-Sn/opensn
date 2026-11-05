#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2-material embedded circles transport test with fixed source.
Surface postprocessor acts on an arc of the outer circle (between y=0.25 and y=0.75),
restricted using RCCLogicalVolume over the uniform boundary 'outer'.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RCCLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.post import SurfacePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

if __name__ == "__main__":
    meshgen = FromFileMeshGenerator(filename="../../../../assets/mesh/circles.msh")
    grid = meshgen.Execute()
    grid.SetUniformBoundaryID("outer")

    # 2 materials: Block 1 is background, Block 2 is inner circle
    xs_bg = MultiGroupXS()
    xs_bg.CreateSimpleOneGroup(1.0, 0.9)

    xs_inner = MultiGroupXS()
    xs_inner.CreateSimpleOneGroup(2.0, 0.5)

    # Fixed source in inner circle (Block 2)
    src = VolumetricSource(block_ids=[2], group_strength=[10.0])

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)

    # Create solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
            }
        ],
        xs_map=[
            {"block_ids": [1], "xs": xs_bg},
            {"block_ids": [2], "xs": xs_inner},
        ],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "outer", "type": "vacuum"},
        ],
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Restrict surface postprocessor to the arc of the outer circle using RCCLogicalVolume
    lv_circle = RCCLogicalVolume(r=0.255, x0=0.0, y0=0.5, z0=-1.0, vx=0.0, vy=0.0, vz=2.0)

    pps_integral_net = SurfacePostprocessor(
        problem=phys,
        value_type="integral",
        current_type="net",
        boundaries=["outer"],
        logical_volumes=[lv_circle],
    )
    pps_integral_net.Execute()
    val_int_net = pps_integral_net.GetValue()[0][0]
    assert math.isclose(val_int_net, 0.16385257, rel_tol=1e-5)

    pps_integral_in = SurfacePostprocessor(
        problem=phys,
        value_type="integral",
        current_type="incoming",
        boundaries=["outer"],
        logical_volumes=[lv_circle],
    )
    pps_integral_in.Execute()
    assert math.isclose(pps_integral_in.GetValue()[0][0], 0., abs_tol=1e-15)

    pps_integral_out = SurfacePostprocessor(
        problem=phys,
        value_type="integral",
        current_type="outgoing",
        boundaries=["outer"],
        logical_volumes=[lv_circle],
    )
    pps_integral_out.Execute()
    val_int_out = pps_integral_out.GetValue()[0][0]
    assert math.isclose(val_int_out, 0.16385257, rel_tol=1e-5)

    pps_avg_net = SurfacePostprocessor(
        problem=phys,
        value_type="avg",
        current_type="net",
        boundaries=["outer"],
        logical_volumes=[lv_circle],
    )
    pps_avg_net.Execute()
    val_avg_net = pps_avg_net.GetValue()[0][0]
    assert math.isclose(val_avg_net, 0.31296494, rel_tol=1e-5)

    pps_max_net = SurfacePostprocessor(
        problem=phys,
        value_type="max",
        current_type="net",
        boundaries=["outer"],
        logical_volumes=[lv_circle],
    )
    pps_max_net.Execute()
    val_max_net = pps_max_net.GetValue()[0][0]
    assert math.isclose(val_max_net, 0.01543783, rel_tol=1e-5)

    pps_min_net = SurfacePostprocessor(
        problem=phys,
        value_type="min",
        current_type="net",
        boundaries=["outer"],
        logical_volumes=[lv_circle],
    )
    pps_min_net.Execute()
    val_min_net = pps_min_net.GetValue()[0][0]
    assert math.isclose(val_min_net, 0.00243616, rel_tol=1e-5)
