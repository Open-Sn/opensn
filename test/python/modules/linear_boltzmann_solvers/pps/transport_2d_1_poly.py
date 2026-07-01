#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D PWLD transport test with vacuum and Incident-isotropic bundary conditions
Test: Max-value=0.50758 and 2.52527e-04
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
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.post import SurfacePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

if __name__ == "__main__":
    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/square_mesh2x2_quads_block.obj",
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            nz=1,
            xcuts=[0.0],
            ycuts=[0.0],
        ),
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()

    # Cross-section data
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)
    num_groups = 168
    xs_3_170 = MultiGroupXS()
    xs_3_170.LoadFromOpenSn("../../../../assets/xs/xs_168g.xs")

    # Volumetric sources
    strength = []
    for g in range(num_groups):
        strength.append(0.0)
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)
    mg_src2 = VolumetricSource(block_ids=[2], group_strength=strength)

    # Boundary sources
    bsrc = []
    for g in range(num_groups):
        bsrc.append(0.0)
    bsrc[0] = 1.0

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=8, scattering_order=1)

    # Create solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, 62),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
            {
                "groups_from_to": (63, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_3_170}],
        volumetric_sources=[mg_src1, mg_src2],
        boundary_conditions=[
            {"name": "xmin", "type": "isotropic", "group_strength": bsrc},
        ],
        options={"max_ags_iterations": 1},
        sweep_type="CBC",
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # PPS over xmin boundary
    pps_xmin_min = SurfacePostprocessor(
        problem=phys,
        value_type="min",
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_min.Execute()

    min = pps_xmin_min.GetValue()[0]
    assert math.isclose(min[0], -0.3341564, rel_tol=1e-6)
    assert math.isclose(min[1], 2.09176149e-4, rel_tol=1e-6)

    pps_xmin_max = SurfacePostprocessor(
        problem=phys,
        value_type="max",
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_max.Execute()

    max = pps_xmin_max.GetValue()[0]
    assert math.isclose(max[0], -0.33314547, rel_tol=1e-6)
    assert math.isclose(max[1], 1.3370162e-3, rel_tol=1e-6)

    pps_xmin_int = SurfacePostprocessor(
        problem=phys,
        value_type="integral",
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_int.Execute()

    int = pps_xmin_int.GetValue()[0]
    assert math.isclose(int[0], -10.67881732, rel_tol=1e-6)
    assert math.isclose(int[1], 2.723323766e-2, rel_tol=1e-6)

    # PPS over xmin boundary for a single groups
    pps_xmin_g1_min = SurfacePostprocessor(
        problem=phys,
        value_type="min",
        group=4,
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_g1_min.Execute()

    min = pps_xmin_g1_min.GetValue()[0]
    assert math.isclose(min[0], 4.49308817e-4, rel_tol=1e-6)

    pps_xmin_g1_max = SurfacePostprocessor(
        problem=phys,
        value_type="max",
        group=4,
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_g1_max.Execute()

    max = pps_xmin_g1_max.GetValue()[0]
    assert math.isclose(max[0], 6.39729658e-4, rel_tol=1e-6)

    pps_xmin_g1_int = SurfacePostprocessor(
        problem=phys,
        value_type="integral",
        group=4,
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_g1_int.Execute()

    int = pps_xmin_g1_int.GetValue()[0]
    assert math.isclose(int[0], 1.81235559e-2, rel_tol=1e-6)

    # PPS over xmin boundary for a groupset 1
    pps_xmin_gs1_min = SurfacePostprocessor(
        problem=phys,
        value_type="min",
        groupset=1,
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_gs1_min.Execute()

    min = pps_xmin_gs1_min.GetValue()[0]
    assert math.isclose(min[0], 9.679619460480029e-05, rel_tol=1e-6)
    assert math.isclose(min[1], 8.67902372833718e-05, rel_tol=1e-6)

    pps_xmin_gs1_max = SurfacePostprocessor(
        problem=phys,
        value_type="max",
        groupset=1,
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_gs1_max.Execute()

    max = pps_xmin_gs1_max.GetValue()[0]
    assert math.isclose(max[0], 5.97642956e-4, rel_tol=1e-6)
    assert math.isclose(max[1], 5.38253292e-4, rel_tol=1e-6)

    pps_xmin_gs1_int = SurfacePostprocessor(
        problem=phys,
        value_type="integral",
        groupset=1,
        current_type="net",
        boundaries=["xmin"]
    )
    pps_xmin_gs1_int.Execute()

    int = pps_xmin_gs1_int.GetValue()[0]
    assert math.isclose(int[0], 1.32020085e-2, rel_tol=1e-6)
    assert math.isclose(int[1], 1.18703524e-2, rel_tol=1e-6)
