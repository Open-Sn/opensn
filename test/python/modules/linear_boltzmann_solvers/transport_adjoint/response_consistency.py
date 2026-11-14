#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Response evaluation consistent with the forward solve, and validated adjoint data (2 ranks).

A 2D one-group problem is solved forward and in adjoint mode for the flux integrated over a QoI
region. By the discrete forward-adjoint identity, a source folded against the adjoint flux equals
the forward QoI of that source to the solver tolerance (1e-12):

- RESP_FUNC_REL_DIFF: a spatially varying (func) volumetric source, which the forward solve
  applies through the consistent mass matrix.
- RESP_WINDOW_REL_DIFF: the same source plus a source with start_time = 1, which a steady
  forward solve at t = 0 does not apply.
- RESP_FORWARD_DATA_REL_DIFF: the forward flux used as a data set, folded against the QoI
  source, which gives the forward QoI directly.

Data sets that do not match the problem are rejected on every rank (each flag is 1 only if every
rank raised the error):

- RESP_QUADRATURE_REJECTED: angular fluxes from a different quadrature with the same number of
  directions.
- RESP_PARTITION_REJECTED: flux moments written with a different mesh partition.
- RESP_MISSING_FILE_REJECTED: a flux-moment file missing on one rank only.
- RESP_INPUT_REJECTED: unknown data sets, duplicate boundary sources, nonexistent blocks, and
  nonfinite source strengths.
- RESP_SOLVER_ORIENTATION_REJECTED: an adjoint angular-flux file written before the adjoint
  solve, so still in solver orientation, whose directions the boundary term would mispair.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    global_sum = MPI.COMM_WORLD.allreduce
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import VectorSpatialFunction
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier
    global_sum = MPIAllReduce


def linear_source(p, num_groups):
    return [1.0 + 3.0 * p.x + 2.0 * p.y] * num_groups


def make_grid(cut_axis):
    nodes = [i / 8 for i in range(9)]
    partitioner = (KBAGraphPartitioner(nx=2, ny=1, xcuts=[0.5]) if cut_axis == "x"
                   else KBAGraphPartitioner(nx=1, ny=2, ycuts=[0.5]))
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes], partitioner=partitioner).Execute()
    grid.SetUniformBlockID(0)
    return grid


def make_problem(grid, xs, n_polar, n_azimuthal, sources):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=n_polar, n_azimuthal=n_azimuthal, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=sources,
        options={"save_angular_flux": True},
    )


def raises_on_all_ranks(action, error_type):
    try:
        action()
        raised = 0
    except error_type:
        raised = 1
    return int(global_sum(raised) == size)


if __name__ == "__main__":
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.5)
    qoi_volume = RPPLogicalVolume(xmin=0.5, xmax=1.0, ymin=0.5, ymax=1.0, infz=True)

    func_source = VolumetricSource(block_ids=[0], func=VectorSpatialFunction(linear_source))
    late_source = VolumetricSource(block_ids=[0], group_strength=[5.0], start_time=1.0)

    # Forward solve at t = 0: the late source is inactive.
    grid = make_grid("x")
    problem = make_problem(grid, xs, 4, 8, [func_source, late_source])
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType("sum")
    ffi.SetLogicalVolume(qoi_volume)
    ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
    ffi.Execute()
    forward_qoi = ffi.GetValue()
    problem.WriteFluxMoments("response_consistency_forward_phi")

    # Adjoint solve for the QoI. An angular flux written before the solve is in solver orientation.
    problem.SetAdjoint(True)
    qoi_source = VolumetricSource(logical_volume=qoi_volume, group_strength=[1.0])
    problem.SetVolumetricSources(volumetric_sources=[qoi_source])
    problem.WriteAngularFluxes("response_consistency_solver_psi")
    solver.Execute()
    problem.WriteFluxMoments("response_consistency_phi")
    problem.WriteAngularFluxes("response_consistency_psi")

    evaluator = ResponseEvaluator(problem=problem)
    evaluator.SetOptions(
        data=[{"name": "phi", "file_prefixes": {"flux_moments": "response_consistency_phi"}}],
        sources={"volumetric": [func_source]},
    )
    func_rel_diff = abs(evaluator.EvaluateResponse("phi") / forward_qoi - 1.0)
    evaluator.SetOptions(sources={"volumetric": [late_source]})
    window_rel_diff = abs(evaluator.EvaluateResponse("phi") / forward_qoi - 1.0)

    forward_data = ResponseEvaluator(problem=problem)
    forward_data.SetOptions(
        data=[{"name": "fwd",
               "file_prefixes": {"flux_moments": "response_consistency_forward_phi"}}],
        sources={"volumetric": [qoi_source]},
    )
    forward_data_rel_diff = abs(forward_data.EvaluateResponse("fwd") / forward_qoi - 1.0)

    solver_orientation_rejected = raises_on_all_ranks(
        lambda: ResponseEvaluator(problem=problem).SetOptions(
            data=[{"name": "psi",
                   "file_prefixes": {"angular_fluxes": "response_consistency_solver_psi"}}]
        ),
        RuntimeError,
    )

    # Same number of directions (16), different quadrature.
    other_quadrature = make_problem(grid, xs, 8, 4, [])
    quadrature_rejected = raises_on_all_ranks(
        lambda: ResponseEvaluator(problem=other_quadrature).SetOptions(
            data=[{"name": "psi",
                   "file_prefixes": {"angular_fluxes": "response_consistency_psi"}}]
        ),
        RuntimeError,
    )

    # Same mesh and rank count, different partition.
    other_partition = make_problem(make_grid("y"), xs, 4, 8, [])
    partition_rejected = raises_on_all_ranks(
        lambda: ResponseEvaluator(problem=other_partition).SetOptions(
            data=[{"name": "phi",
                   "file_prefixes": {"flux_moments": "response_consistency_phi"}}]
        ),
        RuntimeError,
    )

    # A file that exists on every rank but the last.
    if rank == size - 1:
        os.remove(f"response_consistency_phi{rank}.h5")
    barrier()
    missing_file_rejected = raises_on_all_ranks(
        lambda: ResponseEvaluator(problem=problem).SetOptions(
            data=[{"name": "phi",
                   "file_prefixes": {"flux_moments": "response_consistency_phi"}}]
        ),
        RuntimeError,
    )

    unknown_data_rejected = raises_on_all_ranks(
        lambda: evaluator.EvaluateResponse("no_such_data"), ValueError
    )
    duplicate_boundary_rejected = raises_on_all_ranks(
        lambda: ResponseEvaluator(problem=problem).SetOptions(
            sources={"boundary": [
                {"name": "xmin", "type": "isotropic", "group_strength": [1.0]},
                {"name": "xmin", "type": "isotropic", "group_strength": [2.0]},
            ]}
        ),
        ValueError,
    )
    missing_block_rejected = raises_on_all_ranks(
        lambda: ResponseEvaluator(problem=problem).SetOptions(
            sources={"material": [{"block_id": 99, "strength": [1.0]}]}
        ),
        ValueError,
    )
    nonfinite_source_rejected = raises_on_all_ranks(
        lambda: ResponseEvaluator(problem=problem).SetOptions(
            sources={"material": [{"block_id": 0, "strength": [float("nan")]}]}
        ),
        ValueError,
    )
    barrier()
    if rank == 0:
        print(f"RESP_FUNC_REL_DIFF={func_rel_diff:.6e}")
        print(f"RESP_WINDOW_REL_DIFF={window_rel_diff:.6e}")
        print(f"RESP_FORWARD_DATA_REL_DIFF={forward_data_rel_diff:.6e}")
        print(f"RESP_QUADRATURE_REJECTED={quadrature_rejected}")
        print(f"RESP_PARTITION_REJECTED={partition_rejected}")
        print(f"RESP_MISSING_FILE_REJECTED={missing_file_rejected}")
        input_rejected = (
            unknown_data_rejected
            * duplicate_boundary_rejected
            * missing_block_rejected
            * nonfinite_source_rejected
        )
        print(f"RESP_INPUT_REJECTED={input_rejected}")
        print(f"RESP_SOLVER_ORIENTATION_REJECTED={solver_orientation_rejected}")
        for name in os.listdir("."):
            if name.startswith("response_consistency_") and name.endswith(".h5"):
                os.remove(name)
