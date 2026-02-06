#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_disk_leakage_reflecting_bcs.py: Disk with radial vacuum
boundary and reflecting z boundaries.
Expected: From 3D Cartesian disk problem
(transport_3d_disk_leakage_reflecting_bcs.py):
Top=-6.190966583965322e-15
Bottom=-5.231802747012948e-16
Side=97.24938933424164
Total=97.24938933424164
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver


if __name__ == "__main__":
    rmin = 0.0
    rmax = 0.2
    zmin = 0.0
    zmax = 0.1
    nr = 100
    nz = 50
    r_nodes = [rmin + i * ((rmax - rmin) / nr) for i in range(nr + 1)]
    z_nodes = [zmin + i * ((zmax - zmin) / nz) for i in range(nz + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[r_nodes, z_nodes], coord_sys="cylindrical")
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)

    vol_src = VolumetricSource(block_ids=[0], group_strength=[9754.5])

    quad = GLCProductQuadrature2DRZ(n_polar=16, n_azimuthal=32, scattering_order=0)
    bcs = [
        {"name": "rmin", "type": "reflecting"},  # axis symmetry
        {"name": "rmax", "type": "vacuum"},
        {"name": "zmin", "type": "reflecting"},
        {"name": "zmax", "type": "reflecting"},
    ]

    phys = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 50,
                "gmres_restart_interval": 10,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=bcs,
        volumetric_sources=[vol_src],
        options={"save_angular_flux": True},
    )

    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()
    solver.Execute()

    leakage = phys.ComputeLeakage(["rmax", "zmin", "zmax"])
    lkg_side = leakage["rmax"]
    lkg_bottom = leakage["zmin"]
    lkg_top = leakage["zmax"]
    lkg_total = lkg_side + lkg_top + lkg_bottom
    scale = 2.0 * math.pi

    top_val = float(lkg_top.item()) * scale
    bottom_val = float(lkg_bottom.item()) * scale
    side_val = float(lkg_side.item()) * scale
    total_val = float(lkg_total.item()) * scale

    if rank == 0:
        volume = math.pi * rmax * rmax * (zmax - zmin)
        total_source = 9754.5 * volume
        print(f"Top leakage={top_val}")
        print(f"Bottom leakage={bottom_val}")
        print(f"Side leakage={side_val}")
        print(f"Total leakage={total_val}")
        print(f"Total source={total_source}")
        if total_source > 0.0:
            print(f"Leakage/Source={total_val / total_source}")
