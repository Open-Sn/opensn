#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun, February 23

Copyright (c) 2025 quocdang1998
"""

import os
import sys
import math

sys.path.append(os.path.abspath("../../../"))

from pyopensn.mesh import OrthogonalMeshGenerator, KBAGraphPartitioner
from pyopensn.mat import AddMaterial
from pyopensn.xs import LoadFromOpenSn
from pyopensn.source import VolumetricSource
from pyopensn.aquad import ProductQuadrature
from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
from pyopensn.fieldfunc import FieldFunctionGridBased

if __name__ == "__main__":

    # compute node coordinates
    nodes = []
    n_cells = 10
    length = 2.0
    xmin = - length / 2
    dx = length / n_cells
    for i in range(n_cells + 1):
        nodes.append(xmin + i * dx)

    # generate the mesh
    meshgen = OrthogonalMeshGenerator(
        node_sets = [nodes, nodes],
        partitioner = KBAGraphPartitioner(
            nx = 2,
            ny = 2,
            xcuts = [0.0],
            ycuts = [0.0],
        )
    )
    grid = meshgen.Execute()

    # set material
    grid.SetUniformMaterialID(0)

    # add material
    materials = []
    materials.append(AddMaterial("Material_A"))

    # load cross sections
    xs_matA = LoadFromOpenSn("xs_1g_MatA.xs")
    materials[0].SetTransportXSections(xs_matA)

    # create volumetric source
    num_groups = 1
    strength = []
    for g in range(num_groups):
        strength.append(1.0)
    mg_src = VolumetricSource(block_ids = [0], group_strength = strength)

    # initialize quadrature
    nazimu = 1
    npolar = 2
    pquad = ProductQuadrature("gauss-legendre-chebyshev", nazimu, npolar)
    pquad.OptimizeForPolarSymmetry(4.0 * math.pi)

    # create linear Boltzmann solver (LBS)
    phys = DiscreteOrdinatesSolver(
        mesh = grid,
        num_groups = num_groups,
        groupsets = [
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 30
            }
        ],
        options = {
            "volumetric_sources": [mg_src],
        }
    )

    # initialize steady state solver and execute
    ss_solver = SteadyStateSolver(lbs_solver = phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # retrieve the flux (field function)
    fflist = phys.GetScalarFieldFunctionList()
    vtk_basename = "first_example"
    FieldFunctionGridBased.ExportMultipleToVTK(
        [fflist[0][0]],  # export only the field function of group 0 (first []), moment 0 (second [])
        vtk_basename
    )
