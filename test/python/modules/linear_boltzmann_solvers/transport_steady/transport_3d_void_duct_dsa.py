#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D one-group WGDSA test of a void dog-leg duct through a scattering block. SDM: PWLD.

A 6 cm square duct with three legs runs through a 60 cm cube of a material with sigma_t = 1/cm
and scattering ratio 0.99. A source fills the duct entrance, and every boundary is vacuum.
WGDSA uses the void diffusion coefficient V / S of the duct, its volume over the area of its
walls and open ends. It changes only the convergence of the iteration, so the converged
solution satisfies the particle balance.

Checks:
  WGS status: converged within the within-group iteration limit (parameter wgs_max_its).
  BALANCE_REL: |production - absorption - net outflow| / production.
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
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    wgs_max_its = globals().get("wgs_max_its", 1000)

    # 60 cm cube of 2 cm cells
    nodes = [2.0 * i for i in range(31)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    # Dog-leg duct of side 6 cm at z in [26, 32]: along x, then y, then x to the boundary
    z0, z1 = 26.0, 32.0
    legs = [
        RPPLogicalVolume(xmin=0.0, xmax=46.0, ymin=10.0, ymax=16.0, zmin=z0, zmax=z1),
        RPPLogicalVolume(xmin=40.0, xmax=46.0, ymin=10.0, ymax=46.0, zmin=z0, zmax=z1),
        RPPLogicalVolume(xmin=40.0, xmax=60.0, ymin=40.0, ymax=46.0, zmin=z0, zmax=z1),
    ]
    for leg in legs:
        grid.SetBlockIDFromLogicalVolume(leg, 1, True)

    wall = MultiGroupXS()
    wall.CreateSimpleOneGroup(1.0, 0.99)
    void = MultiGroupXS()
    void.CreateSimpleOneGroup(0.0, 0.0)

    source = VolumetricSource(
        logical_volume=RPPLogicalVolume(xmin=0.0, xmax=6.0, ymin=10.0, ymax=16.0, zmin=z0, zmax=z1),
        group_strength=[1.0],
    )

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=4, n_azimuthal=8, scattering_order=0
                ),
                "inner_linear_method": "petsc_fgmres",
                "gmres_restart_interval": 30,
                "l_abs_tol": 1.0e-8,
                "l_max_its": wgs_max_its,
                "apply_wgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-1,
                "wgdsa_l_max_its": 5,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": wall}, {"block_ids": [1], "xs": void}],
        volumetric_sources=[source],
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    table = solver.ComputeBalanceTable()
    if rank == 0:
        net_outflow = table["outflow_rate"] - table["inflow_rate"]
        residual = table["production_rate"] - table["absorption_rate"] - net_outflow
        print(f"BALANCE_REL={abs(residual) / table['production_rate']:.6e}")
