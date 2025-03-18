#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2G KEigenvalue Solver test using Power Iteration
Test: Final k-eigenvalue: 0.5969127
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
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigen
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    # Setup mesh
    N = 40
    L = 14.0
    xmin = 0.0
    dx = L / N
    nodes = [xmin + k * dx for k in range(N + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    vol1 = RPPLogicalVolume(
        xmin=-1000.0,
        xmax=10.0,
        ymin=-1000.0,
        ymax=10.0,
        infz=True,
    )
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    # Cross-section data
    xss = {}
    xss["0"] = MultiGroupXS()
    xss["0"].LoadFromOpenSn("xs_water_g2.xs")
    xss["1"] = MultiGroupXS()
    xss["1"].LoadFromOpenSn("xs_fuel_g2.xs")
    num_groups = xss["0"].num_groups

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(8, 16)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_max_its": 50,
                "gmres_restart_interval": 50,
                "l_abs_tol": 1.0e-10,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xss["0"]},
            {"block_ids": [1], "xs": xss["1"]},
        ],
        options={
            "boundary_conditions": [
                {"name": "xmin", "type": "reflecting"},
                {"name": "ymin", "type": "reflecting"},
            ],
            "scattering_order": 2,
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": True,
        }
    )
    k_solver = PowerIterationKEigen(lbs_problem=phys)
    k_solver.Initialize()
    k_solver.Execute()
