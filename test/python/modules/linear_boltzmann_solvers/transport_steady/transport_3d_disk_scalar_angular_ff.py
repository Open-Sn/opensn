#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simple unstructured disk mesh example that writes scalar and angular flux field functions.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased


if __name__ == "__main__":
    meshgen = FromFileMeshGenerator(filename="../../../../assets/mesh/disk.msh")
    grid = meshgen.Execute()

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)

    vol_src = VolumetricSource(block_ids=[1], group_strength=[1.0])

    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=16, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [1], "xs": xs}],
        volumetric_sources=[vol_src],
        boundary_conditions=[
            {"name": "Outer", "type": "vacuum"},
            {"name": "Top", "type": "reflecting"},
            {"name": "Bottom", "type": "reflecting"},
        ],
        options={"save_angular_flux": True},
    )

    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()
    solver.Execute()

    angular_ff_list = phys.GetAngularFieldFunctionList(groups=[0], angles=[0])
    FieldFunctionGridBased.ExportMultipleToPVTU(angular_ff_list, "angular_flux_disk")

    scalar_ff = phys.GetScalarFluxFieldFunction()[0]
    FieldFunctionGridBased.ExportMultipleToPVTU([scalar_ff], "scalar_flux_disk")

    if rank == 0:
        print("Wrote scalar_flux_disk.pvtu and angular_flux_disk.pvtu")
