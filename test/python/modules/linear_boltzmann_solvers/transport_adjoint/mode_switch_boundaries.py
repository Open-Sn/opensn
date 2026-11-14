#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A mode change keeps vacuum and reflecting boundaries and clears the boundary sources and the source
moments.

The problem is a one-group 1D slab with sigma_t = 1 and scattering ratio c = 0.5.

Reflecting boundaries are kept. With both faces reflecting the slab is an infinite medium, so a
uniform unit source gives the flat flux phi = q / sigma_a = 2 exactly, forward and adjoint. The
forward problem is solved, switched to adjoint mode with only a new volumetric source, and then
switched back with only the forward source. Without the reflecting boundaries the slab leaks and
the flux falls well below 2.
  MODE_SWITCH_FWD_REFLECTING_ERR, MODE_SWITCH_ADJ_REFLECTING_ERR, MODE_SWITCH_FWD_AGAIN_ERR:
    the largest relative deviation of phi from 2 (0 to the solver tolerance).

Boundary sources are cleared. A forward problem driven only by an isotropic boundary source on zmin
and an arbitrary boundary source on zmax is switched to adjoint mode. With no adjoint source and
vacuum boundaries, the adjoint flux is exactly zero. A kept boundary source would drive it.
  MODE_SWITCH_FWD_BOUNDARY_DRIVEN: 1 when the forward flux is positive everywhere, which shows the
    boundary sources were active.
  MODE_SWITCH_ADJ_SOURCE_FREE_MAX: the largest |phi| of the source-free adjoint solve (0).
  MODE_SWITCH_ADJ_VACUUM_REL_DIFF: with a unit adjoint source added, the relative difference of the
    average flux from a fresh adjoint problem built with vacuum boundaries (0 to the solver
    tolerance). Reflecting boundaries would give 2 instead.

Source moments are cleared. The reflecting slab is driven by unit source moments read from a file
(use_source_moments=True), which gives phi = 2 again. After the switch to adjoint mode the solve
must fail until new moments are read. A kept forward source would silently give phi = 2. With the
moments read again, the one-group adjoint is the forward flux.
  MODE_SWITCH_SRC_MOMENTS_FWD_ERR, MODE_SWITCH_SRC_MOMENTS_ADJ_ERR: the largest relative deviation
    of phi from 2 before the switch and after reading the moments again (0 to the solver tolerance).
  MODE_SWITCH_SRC_MOMENTS_CLEARED: 1 when the adjoint solve without new moments raises ValueError.
With use_gpus=True the sweeps run on the GPU, so the mode switches also rebuild the device data.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import AngularFluxFunction
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


INFINITE_MEDIUM_PHI = 2.0
USE_GPUS = globals().get("use_gpus", False)


def phi_reduce(problem, operation):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType(operation)
    ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
    ffi.Execute()
    return ffi.GetValue()


def infinite_medium_err(problem):
    return max(
        abs(phi_reduce(problem, op) - INFINITE_MEDIUM_PHI) / INFINITE_MEDIUM_PHI
        for op in ("max", "min")
    )


if __name__ == "__main__":
    nodes = [i * 0.1 for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.5)
    quadrature = GLProductQuadrature1DSlab(n_polar=8, scattering_order=0)

    def unit_source():
        return VolumetricSource(block_ids=[0], group_strength=[1.0])

    def make_problem(boundaries, sources, adjoint=False, problem_xs=None, use_moments=False):
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[
                {
                    "groups_from_to": (0, 0),
                    "angular_quadrature": quadrature,
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-12,
                    "l_max_its": 500,
                }
            ],
            xs_map=[{"block_ids": [0], "xs": problem_xs or xs}],
            volumetric_sources=sources,
            boundary_conditions=boundaries,
            options={
                "adjoint": adjoint,
                "use_source_moments": use_moments,
                "verbose_inner_iterations": False,
                "verbose_outer_iterations": False,
            },
            **({"use_gpus": True} if USE_GPUS else {}),
        )

    def solve(problem):
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()

    # Reflecting boundaries are kept across both switches.
    reflecting = [{"name": "zmin", "type": "reflecting"}, {"name": "zmax", "type": "reflecting"}]
    problem = make_problem(reflecting, [unit_source()])
    solve(problem)
    fwd_err = infinite_medium_err(problem)
    problem.SetAdjoint(True)
    problem.SetVolumetricSources(volumetric_sources=[unit_source()])
    solve(problem)
    adj_err = infinite_medium_err(problem)
    problem.SetAdjoint(False)
    problem.SetVolumetricSources(volumetric_sources=[unit_source()])
    solve(problem)
    fwd_again_err = infinite_medium_err(problem)

    # Isotropic and arbitrary boundary sources are cleared and become vacuum.
    boundary_sources = [
        {"name": "zmin", "type": "isotropic", "group_strength": [1.0]},
        {"name": "zmax", "type": "arbitrary", "function": AngularFluxFunction(lambda g, n: 1.0)},
    ]
    problem = make_problem(boundary_sources, [])
    solve(problem)
    boundary_driven = 1 if phi_reduce(problem, "min") > 0.1 else 0
    problem.SetAdjoint(True)
    solve(problem)
    source_free_max = max(abs(phi_reduce(problem, op)) for op in ("max", "min"))
    problem.SetVolumetricSources(volumetric_sources=[unit_source()])
    solve(problem)
    switched_avg = phi_reduce(problem, "avg")

    vacuum = [{"name": "zmin", "type": "vacuum"}, {"name": "zmax", "type": "vacuum"}]
    fresh = make_problem(vacuum, [unit_source()], adjoint=True)
    solve(fresh)
    fresh_avg = phi_reduce(fresh, "avg")
    vacuum_rel_diff = abs(switched_avg - fresh_avg) / abs(fresh_avg)

    # Source moments are cleared. A pure absorber with a unit source gives phi = 1, whose flux
    # moments are read back as unit source moments.
    absorber = MultiGroupXS()
    absorber.CreateSimpleOneGroup(1.0, 0.0)
    helper = make_problem(reflecting, [unit_source()], problem_xs=absorber)
    solve(helper)
    moments_prefix = "mode_switch_boundaries_moments"
    helper.WriteFluxMoments(moments_prefix)
    problem = make_problem(reflecting, [], use_moments=True)
    problem.ReadSourceMoments(moments_prefix, False)
    solve(problem)
    moments_fwd_err = infinite_medium_err(problem)
    problem.SetAdjoint(True)
    try:
        solve(problem)
        moments_cleared = 0
    except ValueError as error:
        moments_cleared = int("requires source moments" in str(error))
    problem.ReadSourceMoments(moments_prefix, False)
    solve(problem)
    moments_adj_err = infinite_medium_err(problem)
    os.remove(f"{moments_prefix}{rank}.h5")

    if rank == 0:
        print(f"MODE_SWITCH_FWD_REFLECTING_ERR={fwd_err:.6e}")
        print(f"MODE_SWITCH_ADJ_REFLECTING_ERR={adj_err:.6e}")
        print(f"MODE_SWITCH_FWD_AGAIN_ERR={fwd_again_err:.6e}")
        print(f"MODE_SWITCH_FWD_BOUNDARY_DRIVEN={boundary_driven}")
        print(f"MODE_SWITCH_ADJ_SOURCE_FREE_MAX={source_free_max:.6e}")
        print(f"MODE_SWITCH_ADJ_VACUUM_REL_DIFF={vacuum_rel_diff:.6e}")
        print(f"MODE_SWITCH_SRC_MOMENTS_FWD_ERR={moments_fwd_err:.6e}")
        print(f"MODE_SWITCH_SRC_MOMENTS_CLEARED={moments_cleared}")
        print(f"MODE_SWITCH_SRC_MOMENTS_ADJ_ERR={moments_adj_err:.6e}")
