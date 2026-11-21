#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Write the adjoint data that response_2d_repartition_read.py reads on another number of ranks.

The model is that of response_2d_boundary_source.py with a source region (block 2, x < 1,
4 < y < 6) of the background material. The adjoint problem, driven by the detector cross section
over the detector (block 1, x > 9, 4 < y < 6) with vacuum boundaries, is solved on 4 ranks, and its
flux moments and angular fluxes are written as sets of 4 files. response_2d_repartition_read.py
depends on this test and removes the files.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    nodes = [i * 0.5 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    detector_volume = RPPLogicalVolume(xmin=9.0, xmax=10.0, ymin=4.0, ymax=6.0, infz=True)
    source_volume = RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=4.0, ymax=6.0, infz=True)
    grid.SetBlockIDFromLogicalVolume(detector_volume, 1, True)
    grid.SetBlockIDFromLogicalVolume(source_volume, 2, True)

    xs_background = MultiGroupXS()
    xs_background.CreateSimpleOneGroup(0.3, 0.5)
    xs_detector = MultiGroupXS()
    xs_detector.CreateSimpleOneGroup(0.8, 0.0)
    detector_sigma = 0.8

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=4, n_azimuthal=16, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[
            {"block_ids": [0, 2], "xs": xs_background},
            {"block_ids": [1], "xs": xs_detector},
        ],
        boundary_conditions=[
            {"name": name, "type": "vacuum"} for name in ("xmin", "xmax", "ymin", "ymax")
        ],
        options={"save_angular_flux": True},
    )
    problem.SetAdjoint(True)
    problem.SetVolumetricSources(
        volumetric_sources=[
            VolumetricSource(logical_volume=detector_volume, group_strength=[detector_sigma])
        ]
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    problem.WriteFluxMoments("RepartitionAdjPhi_p")
    problem.WriteAngularFluxes("RepartitionAdjPsi_p")

    if rank == 0:
        print("REPARTITION_WRITE_DONE=1")
