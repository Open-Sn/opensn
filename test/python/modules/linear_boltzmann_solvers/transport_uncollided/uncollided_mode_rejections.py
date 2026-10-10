#!/usr/bin/env python3
"""
Verify that a problem with an uncollided flux rejects time-dependent and adjoint modes, both at
construction and when the mode is changed afterwards, rejects replacement cross sections (the
first-collision source was computed from the original ones), rejects other sources (the
first-collision source replaces fixed sources, which would be silently ignored), both at
construction and when added afterwards, and that a rejected change leaves the problem unchanged.
Uncollided generation itself rejects a point source outside its near-source logical volume at
construction (it used to fail partway through UncollidedSolver.Execute()).

A rejected SetTimeDependentMode(), SetAdjoint(True), SetXSMap(), or SetVolumetricSources() must
leave the problem as it was, so a steady-state solve afterwards must reproduce the solve of an
identical fresh problem exactly.

Checks:
  UncollidedCtorTimeDependentRejected, UncollidedCtorAdjointRejected,
  UncollidedSetTimeDependentRejected, UncollidedSetAdjointRejected, UncollidedSetXSMapRejected,
  UncollidedCtorSourceRejected, UncollidedLaterSourceRejected: 1 if rejected with the
    uncollided-flux error.
  UncollidedNearSourceRejected: 1 if the misplaced near-source volume is rejected.
  UncollidedRejectedStateUnchanged: 1 if the problem is still steady-state and forward.
  UncollidedRejectedMaxDiff: max |phi - phi_fresh| after the rejected changes.
"""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
        UncollidedSolver,
    )
    from pyopensn.source import PointSource, VolumetricSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
remove_file = uncollided_utils.remove_file

STEADY_ONLY = "uncollided flux is only supported for steady-state"
FORWARD_ONLY = "uncollided flux is not supported for adjoint"


def make_uncollided_file(grid, xs, file_name):
    whole_domain = RPPLogicalVolume(xmin=-1.01, xmax=1.01, ymin=-1.01, ymax=1.01, infz=True)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=[0.13, -0.22, 0.0], strength=[1.0])],
        near_source=[whole_domain],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()


def near_source_outside_rejected(grid, xs):
    # The near-source volume does not contain the point source at (0.13, -0.22).
    away = RPPLogicalVolume(xmin=0.5, xmax=1.01, ymin=0.5, ymax=1.01, infz=True)
    return rejected(
        lambda: UncollidedProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{"groups_from_to": [0, 0]}],
            xs_map=[{"block_ids": [0], "xs": xs}],
            point_sources=[PointSource(location=[0.13, -0.22, 0.0], strength=[1.0])],
            near_source=[away],
            scattering_order=0,
        ),
        "lies outside its near-source logical volume",
    )


def make_problem(grid, xs, file_name, **kwargs):
    options = {"save_angular_flux": True, "verbose_inner_iterations": False}
    options.update(kwargs.pop("options", {}))
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=12, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=file_name,
        options=options,
        **kwargs,
    )


def rejected(action, expected):
    try:
        action()
    except ValueError as error:
        return expected in str(error)
    return False


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal())


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)

    # The velocity makes the cross sections valid for time-dependent mode, so the uncollided flux
    # is the only reason for the rejection.
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.7, c=0.3, velocity=1.0)

    file_name = "uncollided_mode_rejections.h5"
    remove_file(file_name)
    make_uncollided_file(grid, xs, file_name)

    ctor_td = rejected(lambda: make_problem(grid, xs, file_name, time_dependent=True), STEADY_ONLY)
    ctor_adj = rejected(
        lambda: make_problem(grid, xs, file_name, options={"adjoint": True}), FORWARD_ONLY
    )

    problem = make_problem(grid, xs, file_name)
    set_td = rejected(problem.SetTimeDependentMode, STEADY_ONLY)
    set_adj = rejected(lambda: problem.SetAdjoint(True), FORWARD_ONLY)
    other_xs = MultiGroupXS()
    other_xs.CreateSimpleOneGroup(sigma_t=0.7, c=0.3, velocity=1.0)
    set_xs = rejected(
        lambda: problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": other_xs}]),
        "cross sections cannot be replaced after loading an uncollided flux file",
    )
    unchanged = not problem.IsTimeDependent() and not problem.IsAdjoint()

    def volumetric_source():
        return VolumetricSource(block_ids=[0], group_strength=[10.0])

    ctor_source = rejected(
        lambda: make_problem(grid, xs, file_name, volumetric_sources=[volumetric_source()]),
        "An uncollided flux file replaces fixed sources",
    )
    later_source = rejected(
        lambda: problem.SetVolumetricSources(volumetric_sources=[volumetric_source()]),
        "An uncollided flux file replaces fixed sources",
    )

    near_source = near_source_outside_rejected(grid, xs)

    max_diff = float("inf")
    if unchanged:
        phi = solve(problem)
        phi_fresh = solve(make_problem(grid, xs, file_name))
        max_diff = max(abs(a - b) for a, b in zip(phi, phi_fresh))

    remove_file(file_name)

    if rank == 0:
        print(f"UncollidedCtorTimeDependentRejected={int(ctor_td)}")
        print(f"UncollidedCtorAdjointRejected={int(ctor_adj)}")
        print(f"UncollidedSetTimeDependentRejected={int(set_td)}")
        print(f"UncollidedSetAdjointRejected={int(set_adj)}")
        print(f"UncollidedSetXSMapRejected={int(set_xs)}")
        print(f"UncollidedCtorSourceRejected={int(ctor_source)}")
        print(f"UncollidedLaterSourceRejected={int(later_source)}")
        print(f"UncollidedNearSourceRejected={int(near_source)}")
        print(f"UncollidedRejectedStateUnchanged={int(unchanged)}")
        print(f"UncollidedRejectedMaxDiff={max_diff:.8e}")
        sys.stdout.flush()
