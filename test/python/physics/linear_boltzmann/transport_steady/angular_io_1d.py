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
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.logvol import RPPLogicalVolume

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
    gl_quad = GLProductQuadrature1DSlab(n_polar=128)

    # LBS block option
    num_groups = 1
    solver_dict = {}
    solver_dict["mesh"] = grid
    solver_dict["num_groups"] = num_groups
    solver_dict["groupsets"] = [
        {
            "groups_from_to": (0, num_groups - 1),
            "angular_quadrature": gl_quad,
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-9,
            "l_max_its": 300,
            "gmres_restart_interval": 30,
        },
    ]
    solver_dict["xs_map"] = xs_map
    solver_dict["options"] = {
        "scattering_order": 0,
        "spatial_discretization": "pwld",
        "boundary_conditions": [
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"}
        ],
        "volumetric_sources": [src0, src1],
        "save_angular_flux": True,
    }

    phys1 = DiscreteOrdinatesProblem(**solver_dict)

    # Initialize and execute solver
    ss_solver = SteadyStateSolver(lbs_problem=phys1)
    ss_solver.Initialize()
    ss_solver.Execute()

    leakage_left_1 = phys1.ComputeLeakage(["zmin"])["zmin"][0]
    leakage_right_1 = phys1.ComputeLeakage(["zmax"])["zmax"][0]
    phys1.WriteAngularFluxes("angular_io")

    phys2 = DiscreteOrdinatesProblem(**solver_dict)
    ss_solver_2 = SteadyStateSolver(lbs_problem=phys2)
    ss_solver_2.Initialize()
    phys2.ReadAngularFluxes("angular_io")

    leakage_left_2 = phys2.ComputeLeakage(["zmin"])["zmin"][0]
    leakage_right_2 = phys2.ComputeLeakage(["zmax"])["zmax"][0]

    leakage_left_diff = leakage_left_1 - leakage_left_2
    leakage_right_diff = leakage_right_1 - leakage_right_2

    if rank == 0:
        print(f"Leakage-Diff1={leakage_left_diff:.5e}")
    if rank == 0:
        print(f"Leakage-Diff2={leakage_right_diff:.5e}")

    os.remove(f"angular_io{rank}.h5")
