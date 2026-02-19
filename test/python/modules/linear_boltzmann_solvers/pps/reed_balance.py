#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Standard Reed 1D 1-group problem

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
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.post import VolumePostprocessor

if __name__ == "__main__":

    # Create Mesh
    widths = [2., 1., 2., 1., 2.]
    nrefs = [200, 200, 200, 200, 200]
    Nmat = len(widths)
    nodes = [0.]
    for imat in range(Nmat):
        dx = widths[imat] / nrefs[imat]
        for i in range(nrefs[imat]):
            nodes.append(nodes[-1] + dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()

    # Set block IDs
    z_min = 0.0
    z_max = widths[1]
    for imat in range(Nmat):
        z_max = z_min + widths[imat]
        if rank == 0:
            print("imat=", imat, ", zmin=", z_min, ", zmax=", z_max)
        lv = RPPLogicalVolume(infx=True, infy=True, zmin=z_min, zmax=z_max)
        grid.SetBlockIDFromLogicalVolume(lv, imat, True)
        z_min = z_max

    # Add cross sections to materials
    total = [50., 5., 0., 1., 1.]
    c = [0., 0., 0., 0.9, 0.9]
    xs_map = len(total) * [None]
    for imat in range(Nmat):
        xs_ = MultiGroupXS()
        xs_.CreateSimpleOneGroup(total[imat], c[imat])
        xs_map[imat] = {
            "block_ids": [imat], "xs": xs_,
        }

    # Create sources in 1st and 4th materials
    src0 = VolumetricSource(block_ids=[0], group_strength=[50.])
    src1 = VolumetricSource(block_ids=[3], group_strength=[1.])

    # Angular Quadrature
    gl_quad = GLProductQuadrature1DSlab(n_polar=128, scattering_order=0)

    # LBS block option
    num_groups = 1
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": gl_quad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-9,
                "l_max_its": 300,
                "gmres_restart_interval": 30,
            },
        ],
        xs_map=xs_map,
        volumetric_sources=[src0, src1],
        boundary_conditions=[
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"}
        ],
    )

    # Initialize and execute solver
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # PPS over whole domain
    pps_whole_min = VolumePostprocessor(problem=phys, value_type="min")
    pps_whole_min.Execute()

    pps_whole_max = VolumePostprocessor(problem=phys, value_type="max")
    pps_whole_max.Execute()

    pps_whole_avg = VolumePostprocessor(problem=phys, value_type="avg")
    pps_whole_avg.Execute()

    pps_whole_int = VolumePostprocessor(problem=phys, value_type="integral")
    pps_whole_int.Execute()

    if rank == 0:
        print(f"whole: min = {pps_whole_min.GetValue()[0][0]}")
        print(f"whole: max = {pps_whole_max.GetValue()[0][0]}")
        print(f"whole: avg = {pps_whole_avg.GetValue()[0][0]}")
        print(f"whole: int = {pps_whole_int.GetValue()[0][0]}")

    # PPS over a mesh block
    pps_blk4_min = VolumePostprocessor(
        problem=phys,
        value_type="min",
        block_ids=[4]
    )
    pps_blk4_min.Execute()

    pps_blk4_max = VolumePostprocessor(
        problem=phys,
        value_type="max",
        block_ids=[4]
    )
    pps_blk4_max.Execute()

    pps_blk4_avg = VolumePostprocessor(
        problem=phys,
        value_type="avg",
        block_ids=[4]
    )
    pps_blk4_avg.Execute()

    pps_blk4_int = VolumePostprocessor(
        problem=phys,
        value_type="integral",
        block_ids=[4]
    )
    pps_blk4_int.Execute()

    if rank == 0:
        print(f"blk4: min = {pps_blk4_min.GetValue()[0][0]}")
        print(f"blk4: max = {pps_blk4_max.GetValue()[0][0]}")
        print(f"blk4: avg = {pps_blk4_avg.GetValue()[0][0]}")
        print(f"blk4: int = {pps_blk4_int.GetValue()[0][0]}")

    # PPS over a logical volume
    lv1 = RPPLogicalVolume(infz=True, infy=True, zmin=4.0, zmax=7.)
    pps_lv1_min = VolumePostprocessor(
        problem=phys,
        value_type="min",
        logical_volumes=[lv1]
    )
    pps_lv1_min.Execute()

    pps_lv1_max = VolumePostprocessor(
        problem=phys,
        value_type="max",
        logical_volumes=[lv1]
    )
    pps_lv1_max.Execute()

    pps_lv1_avg = VolumePostprocessor(
        problem=phys,
        value_type="avg",
        logical_volumes=[lv1]
    )
    pps_lv1_avg.Execute()

    pps_lv1_int = VolumePostprocessor(
        problem=phys,
        value_type="integral",
        logical_volumes=[lv1]
    )
    pps_lv1_int.Execute()

    if rank == 0:
        print(f"lv1: min = {pps_lv1_min.GetValue()[0][0]}")
        print(f"lv1: max = {pps_lv1_max.GetValue()[0][0]}")
        print(f"lv1: avg = {pps_lv1_avg.GetValue()[0][0]}")
        print(f"lv1: int = {pps_lv1_int.GetValue()[0][0]}")

    # PPS over logical volume and mesh blocks
    lv2 = RPPLogicalVolume(infx=True, infy=True, zmin=6.1, zmax=6.5)
    lv3 = RPPLogicalVolume(infx=True, infy=True, zmin=6.9, zmax=7.2)

    pps_cmb_min = VolumePostprocessor(
        problem=phys,
        value_type="min",
        logical_volumes=[lv2, lv3],
        block_ids=[4]
    )
    pps_cmb_min.Execute()

    pps_cmb_max = VolumePostprocessor(
        problem=phys,
        value_type="max",
        logical_volumes=[lv2, lv3],
        block_ids=[4]
    )
    pps_cmb_max.Execute()

    pps_cmb_avg = VolumePostprocessor(
        problem=phys,
        value_type="avg",
        logical_volumes=[lv2, lv3],
        block_ids=[4]
    )
    pps_cmb_avg.Execute()

    pps_cmb_int = VolumePostprocessor(
        problem=phys,
        value_type="integral",
        logical_volumes=[lv2, lv3],
        block_ids=[4]
    )
    pps_cmb_int.Execute()

    if rank == 0:
        print(f"cmb: min[0] = {pps_cmb_min.GetValue()[0][0]}")
        print(f"cmb: min[1] = {pps_cmb_min.GetValue()[1][0]}")
        print(f"cmb: max[0] = {pps_cmb_max.GetValue()[0][0]}")
        print(f"cmb: max[1] = {pps_cmb_max.GetValue()[1][0]}")
        print(f"cmb: avg[0] = {pps_cmb_avg.GetValue()[0][0]}")
        print(f"cmb: avg[1] = {pps_cmb_avg.GetValue()[1][0]}")
        print(f"cmb: int[0] = {pps_cmb_int.GetValue()[0][0]}")
        print(f"cmb: int[1] = {pps_cmb_int.GetValue()[1][0]}")
