#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adjoint responses on 3 ranks from adjoint data written on 4 ranks.

response_2d_repartition_write.py writes the adjoint flux moments and angular fluxes of its model as
sets of 4 files. This test solves the forward problems of the same model on 3 ranks, so the sets
hold more files than there are ranks: rank 0 reads two files, ranks 1 and 2 one each, and every
cell's values must reach the rank that owns it here. By forward-adjoint duality,

* for an isotropic boundary source on xmin, the forward detector response int_V sigma_d phi dV
  equals int_{xmin} int_{Omega . n < 0} |Omega . n| psi_b psi^dagger dOmega dA, which reads the
  adjoint angular fluxes;
* for a uniform volumetric source q in the source region, the forward detector response equals
  int_V q phi^dagger dV, which reads the adjoint flux moments.

Both identity errors are checked, then the files are removed.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


if __name__ == "__main__":
    num_procs = 3
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
    boundary_strength = 1.0
    source_strength = 1.0
    phi_prefix = "RepartitionAdjPhi_p"
    psi_prefix = "RepartitionAdjPsi_p"
    num_files = 4

    def solve_forward(boundary_source, volumetric_sources):
        bcs = [{"name": name, "type": "vacuum"} for name in ("xmax", "ymin", "ymax")]
        bcs.append(
            {"name": "xmin", "type": "isotropic", "group_strength": [boundary_strength]}
            if boundary_source
            else {"name": "xmin", "type": "vacuum"}
        )
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
            boundary_conditions=bcs,
            volumetric_sources=volumetric_sources,
        )
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()
        interpolator = FieldFunctionInterpolationVolume()
        interpolator.SetOperationType("sum")
        interpolator.SetLogicalVolume(detector_volume)
        interpolator.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
        interpolator.Execute()
        return problem, detector_sigma * interpolator.GetValue()

    # Boundary source: the adjoint response reads the 4-file angular-flux set.
    boundary_problem, boundary_forward_qoi = solve_forward(True, [])
    boundary_evaluator = ResponseEvaluator(problem=boundary_problem)
    boundary_evaluator.SetOptions(
        data=[
            {
                "name": "detector",
                "file_prefixes": {"flux_moments": phi_prefix, "angular_fluxes": psi_prefix},
            }
        ],
        sources={
            "boundary": [
                {"name": "xmin", "type": "isotropic", "group_strength": [boundary_strength]}
            ]
        },
    )
    boundary_adjoint_qoi = boundary_evaluator.EvaluateResponse("detector")

    # Volumetric source in block 2: the adjoint response reads the 4-file flux-moment set.
    material_problem, material_forward_qoi = solve_forward(
        False, [VolumetricSource(block_ids=[2], group_strength=[source_strength])]
    )
    material_evaluator = ResponseEvaluator(problem=material_problem)
    material_evaluator.SetOptions(
        data=[{"name": "detector", "file_prefixes": {"flux_moments": phi_prefix}}],
        sources={"material": [{"block_id": 2, "strength": [source_strength]}]},
    )
    material_adjoint_qoi = material_evaluator.EvaluateResponse("detector")

    if rank == 0:
        print(f"Boundary forward QoI={boundary_forward_qoi:.8e}")
        print(f"Boundary adjoint QoI={boundary_adjoint_qoi:.8e}")
        print(
            "Boundary identity error="
            f"{abs(boundary_adjoint_qoi / boundary_forward_qoi - 1.0):.3e}"
        )
        print(f"Material forward QoI={material_forward_qoi:.8e}")
        print(f"Material adjoint QoI={material_adjoint_qoi:.8e}")
        print(
            "Material identity error="
            f"{abs(material_adjoint_qoi / material_forward_qoi - 1.0):.3e}"
        )

    barrier()
    if rank == 0:
        for prefix in (phi_prefix, psi_prefix):
            for k in range(num_files):
                os.remove(f"{prefix}{k}.h5")
