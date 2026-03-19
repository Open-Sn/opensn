#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SDM: PWLD

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.post import VolumePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

if __name__ == "__main__":
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/inclusions_gmsh_v4.msh",
    )
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()

    # Material
    num_groups = 8
    xs_diag = MultiGroupXS()
    xs_diag.LoadFromOpenSn("../../../../assets/xs/diag_XS_64g_1mom_c0.99.xs")

    # Source
    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 100.0
    strength[6] = 200.0
    mg_src = VolumetricSource(block_ids=[2], group_strength=strength)

    # Quadrature
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=16, scattering_order=1)

    # Set up solver
    gs1 = [0, 2]
    gs2 = [3, 6]
    gs3 = [7, num_groups - 1]
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": gs1,
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "gmres_restart_interval": 50,
                "l_abs_tol": 1.0e-6,
                "l_max_its": 100,
            },
            {
                "groups_from_to": gs2,
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "gmres_restart_interval": 50,
                "l_abs_tol": 1.0e-6,
                "l_max_its": 100,
            },
            {
                "groups_from_to": gs3,
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "gmres_restart_interval": 50,
                "l_abs_tol": 1.0e-6,
                "l_max_its": 100,
            },
        ],
        xs_map=[
            {"block_ids": [1, 2, 3], "xs": xs_diag},
        ],
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "max_ags_iterations": 1
        },
    )

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

    min = pps_whole_min.GetValue()[0]
    assert len(min) == num_groups
    assert math.isclose(min[0], 2505.23468, rel_tol=1e-6)
    assert math.isclose(min[6], 5010.469391, rel_tol=1e-6)

    max = pps_whole_max.GetValue()[0]
    assert len(max) == num_groups
    assert math.isclose(max[0], 2551.382249, rel_tol=1e-6)
    assert math.isclose(max[6], 5102.764544, rel_tol=1e-6)

    avg = pps_whole_avg.GetValue()[0]
    assert len(avg) == num_groups
    assert math.isclose(avg[0], 2525.183749, rel_tol=1e-6)
    assert math.isclose(avg[6], 5050.367541, rel_tol=1e-6)

    int = pps_whole_int.GetValue()[0]
    assert len(int) == num_groups
    assert math.isclose(int[0], 2525.183749, rel_tol=1e-6)
    assert math.isclose(int[6], 5050.367541, rel_tol=1e-6)

    # PPS over a mesh block
    pps_blk2_min = VolumePostprocessor(problem=phys, value_type="min", block_ids=[2])
    pps_blk2_min.Execute()

    pps_blk2_max = VolumePostprocessor(problem=phys, value_type="max", block_ids=[2])
    pps_blk2_max.Execute()

    pps_blk2_avg = VolumePostprocessor(problem=phys, value_type="avg", block_ids=[2])
    pps_blk2_avg.Execute()

    pps_blk2_int = VolumePostprocessor(problem=phys, value_type="integral", block_ids=[2])
    pps_blk2_int.Execute()

    min = pps_blk2_min.GetValue()[0]
    assert len(min) == num_groups
    assert math.isclose(min[0], 2512.2628697, rel_tol=1e-6)
    assert math.isclose(min[6], 5024.5257663, rel_tol=1e-6)

    max = pps_blk2_max.GetValue()[0]
    assert len(max) == num_groups
    assert math.isclose(max[0], 2551.3822487, rel_tol=1e-6)
    assert math.isclose(max[6], 5102.7645446, rel_tol=1e-6)

    avg = pps_blk2_avg.GetValue()[0]
    assert len(avg) == num_groups
    assert math.isclose(avg[0], 2535.7135624, rel_tol=1e-6)
    assert math.isclose(avg[6], 5071.4271681, rel_tol=1e-6)

    int = pps_blk2_int.GetValue()[0]
    assert len(int) == num_groups
    assert math.isclose(int[0], 640.3142799, rel_tol=1e-6)
    assert math.isclose(int[6], 1280.6285707, rel_tol=1e-6)

    # PPS over a logical volume
    lv1 = RPPLogicalVolume(xmin=0.5, xmax=1.0, ymin=0.0, ymax=0.25, zmin=0.7, zmax=1.0)
    pps_lv1_min = VolumePostprocessor(problem=phys, value_type="min", logical_volumes=[lv1])
    pps_lv1_min.Execute()

    pps_lv1_max = VolumePostprocessor(problem=phys, value_type="max", logical_volumes=[lv1])
    pps_lv1_max.Execute()

    pps_lv1_avg = VolumePostprocessor(problem=phys, value_type="avg", logical_volumes=[lv1])
    pps_lv1_avg.Execute()

    pps_lv1_int = VolumePostprocessor(problem=phys, value_type="integral", logical_volumes=[lv1])
    pps_lv1_int.Execute()

    min = pps_lv1_min.GetValue()[0]
    assert len(min) == num_groups
    assert math.isclose(min[0], 2510.0970706, rel_tol=1e-6)
    assert math.isclose(min[6], 5020.1941757, rel_tol=1e-6)

    max = pps_lv1_max.GetValue()[0]
    assert len(max) == num_groups
    assert math.isclose(max[0], 2529.948345, rel_tol=1e-6)
    assert math.isclose(max[6], 5059.896743, rel_tol=1e-6)

    avg = pps_lv1_avg.GetValue()[0]
    assert len(avg) == num_groups
    assert math.isclose(avg[0], 2518.3936498, rel_tol=1e-6)
    assert math.isclose(avg[6], 5036.7873447, rel_tol=1e-6)

    int = pps_lv1_int.GetValue()[0]
    assert len(int) == num_groups
    assert math.isclose(int[0], 94.8310452, rel_tol=1e-6)
    assert math.isclose(int[6], 189.662092, rel_tol=1e-6)

    # PPS over logical volume and mesh blocks
    pps_cmb_min = VolumePostprocessor(
        problem=phys,
        value_type="min",
        logical_volumes=[lv1],
        block_ids=[2])
    pps_cmb_min.Execute()

    pps_cmb_max = VolumePostprocessor(
        problem=phys,
        value_type="max",
        logical_volumes=[lv1],
        block_ids=[2])
    pps_cmb_max.Execute()

    pps_cmb_avg = VolumePostprocessor(
        problem=phys,
        value_type="avg",
        logical_volumes=[lv1],
        block_ids=[2])
    pps_cmb_avg.Execute()

    pps_cmb_int = VolumePostprocessor(
        problem=phys,
        value_type="integral",
        logical_volumes=[lv1],
        block_ids=[2])
    pps_cmb_int.Execute()

    min = pps_cmb_min.GetValue()[0]
    assert len(min) == num_groups
    assert math.isclose(min[0], 2516.3938277, rel_tol=1e-6)
    assert math.isclose(min[6], 5032.7876996, rel_tol=1e-6)

    max = pps_cmb_max.GetValue()[0]
    assert len(max) == num_groups
    assert math.isclose(max[0], 2529.9483446, rel_tol=1e-6)
    assert math.isclose(max[6], 5059.8967431, rel_tol=1e-6)

    avg = pps_cmb_avg.GetValue()[0]
    assert len(avg) == num_groups
    assert math.isclose(avg[0], 2524.7516068, rel_tol=1e-6)
    assert math.isclose(avg[6], 5049.5032621, rel_tol=1e-6)

    int = pps_cmb_int.GetValue()[0]
    assert len(int) == num_groups
    assert math.isclose(int[0], 18.776296, rel_tol=1e-6)
    assert math.isclose(int[6], 37.5525924, rel_tol=1e-6)

    # PPS over a groupset
    pps_gs1_min = VolumePostprocessor(problem=phys, value_type="min", groupset=1)
    pps_gs1_min.Execute()

    pps_gs1_max = VolumePostprocessor(problem=phys, value_type="max", groupset=1)
    pps_gs1_max.Execute()

    pps_gs1_avg = VolumePostprocessor(problem=phys, value_type="avg", groupset=1)
    pps_gs1_avg.Execute()

    pps_gs1_int = VolumePostprocessor(problem=phys, value_type="integral", groupset=1)
    pps_gs1_int.Execute()

    min = pps_gs1_min.GetValue()[0]
    assert len(min) == 4
    assert math.isclose(min[3], 5010.469391, rel_tol=1e-6)

    max = pps_gs1_max.GetValue()[0]
    assert len(max) == 4
    assert math.isclose(max[3], 5102.764544, rel_tol=1e-6)

    avg = pps_gs1_avg.GetValue()[0]
    assert len(avg) == 4
    assert math.isclose(avg[3], 5050.367541, rel_tol=1e-6)

    int = pps_gs1_int.GetValue()[0]
    assert len(int) == 4
    assert math.isclose(int[3], 5050.367541, rel_tol=1e-6)
