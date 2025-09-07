#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D Transport test with Vacuum BCs and a material source writing source moments.
SDM: PWLD
Test: Max-value=1.08320e-01 and 0.0
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    meshgen = ExtruderMeshGenerator(
        inputs=[
            FromFileMeshGenerator(
                filename="../../../../assets/mesh/SquareMesh2x2Quads.obj"
            )
        ],
        layers=[{"z": 0.4, "n": 2},
                {"z": 0.8, "n": 2},
                {"z": 1.2, "n": 2},
                {"z": 1.6, "n": 2},
                ],  # layers
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            xcuts=[0.0],
            ycuts=[0.0], ),
    )
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)
    vol1 = RPPLogicalVolume(
        xmin=-0.5 / 8,
        xmax=0.5 / 8,
        ymin=-0.5 / 8,
        ymax=0.5 / 8,
        infz=True,
    )
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    # Cross sections
    num_groups = 21
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")

    # Source
    strength = [0.0 for _ in range(num_groups)]
    mg_src0 = VolumetricSource(block_ids=[0], group_strength=strength)
    strength[0] = 1.0
    mg_src1 = VolumetricSource(block_ids=[1], group_strength=strength)

    # Setup Physics
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=1)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 20],
                "angular_quadrature": pquad,
                # "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 2,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        scattering_order=1,
        volumetric_sources=[mg_src0, mg_src1],
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    phys.CreateAndWriteSourceMoments("Qmoms")
