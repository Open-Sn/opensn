#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script sets up a 64-group PWLD transport problem that is used for strong
scaling studies with OpenSn.
"""

import os
import sys
import math

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    from pyopensn.mesh import DistributedMeshGenerator, FromFileMeshGenerator, MeshGenerator
    from pyopensn.xs import LoadFromOpenSn
    from pyopensn.aquad import CreateProductQuadrature
    import pyopensn.lbs as lbs
    from pyopensn.mat import SetProperty, ISOTROPIC_MG_SOURCE, FROM_ARRAY
    from pyopensn import mat, materials

Ng = 64          # Number of energy groups
Npolar = 14      # Number of polar angles
Nazimuthal = 32  # Number of azimuthal angles

# Mesh
meshgen = DistributedMeshGenerator(
    inputs=[
        FromFileMeshGenerator(filename="strong_scaling.msh")
    ]
)
grid = meshgen.Execute()

# Cross-section datta
xs_diag = MultiGroupXS()
xs_diag.LoadFromOpenSn("diag_XS_64g_1mom_c0.99.xs")

# Boundary conditions
bsrc = [0.0 for _ in range(Ng)]
bsrc[0] = 1.0 / (4.0 * math.pi)

# Angular quadrature
pquad = GLCProductQuadrature3DXYZ(Npolar, Nazimuthal)

# Solver
phys = DiscreteOrdinatesProblem(
    mesh=grid,
    num_groups=Ng,
    groupsets=[
        {
            "groups_from_to": (0, Ng - 1),
            "angular_quadrature": pquad,
            "angle_aggregation_type": "single",
            "angle_aggregation_num_subsets": 1,
            "inner_linear_method": "petsc_richardson",
            "l_abs_tol": 1.0e-6,
            "l_max_its": 9,
        },
    ],
    xs_map=[
        {"block_ids": [1], "xs": xs_diag},
    ],
    options={
        "boundary_conditions": [
            {"name": "xmin", "type": "isotropic", "group_strength": bsrc},
        ],
        "scattering_order": 0,
    }
)
ss_solver = SteadyStateSolver(lbs_problem=phys)
ss_solver.Initialize()
ss_solver.Execute()
