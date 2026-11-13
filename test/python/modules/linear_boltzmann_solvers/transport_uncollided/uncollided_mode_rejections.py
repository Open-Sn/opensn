#!/usr/bin/env python3
"""
Verify that a problem with an uncollided flux rejects time-dependent mode, and a change of the
transport mode, both at construction and when the mode is changed afterwards, rejects replacement
cross sections (the first-collision source was computed from the original ones), rejects other
sources (the first-collision source replaces fixed sources, which would be silently ignored), both
at construction and when added afterwards, and that a rejected change leaves the problem unchanged.
Uncollided generation itself rejects a point source outside its near-source logical volume at
construction (it used to fail partway through UncollidedSolver.Execute()).

The generator (UncollidedProblem) records its transport mode in the file. A forward file is
rejected by an adjoint problem and an adjoint file by a forward problem, because the point sources
of the two modes are different physical sources. An adjoint problem accepts an adjoint file.
SetAdjoint(True) on the generator clears its point sources, as for any problem, so generation then
rejects the missing point source.

Incoming-flux boundaries are fixed sources too: with an uncollided flux the sweep does not apply
boundary sources, so an isotropic boundary would be silently ignored. It must be rejected at
construction and when set afterwards.

A rejected SetTimeDependentMode(), SetAdjoint(True), SetXSMap(), or SetVolumetricSources() must
leave the problem as it was, so a steady-state solve afterwards must reproduce the solve of an
identical fresh problem exactly.

Checks:
  UncollidedCtorTimeDependentRejected, UncollidedSetTimeDependentRejected,
  UncollidedSetAdjointRejected, UncollidedSetXSMapRejected, UncollidedCtorSourceRejected,
  UncollidedLaterSourceRejected: 1 if rejected with the uncollided-flux error.
  UncollidedCtorAdjointRejected: 1 if an adjoint problem rejects the forward file.
  UncollidedNearSourceRejected: 1 if the misplaced near-source volume is rejected.
  UncollidedRejectedStateUnchanged: 1 if the problem is still steady-state and forward.
  UncollidedRejectedMaxDiff: max |phi - phi_fresh| after the rejected changes.
  UncollidedAdjointFileForwardRejected: 1 if a forward problem rejects an adjoint file.
  UncollidedAdjointFileAdjointAccepted: 1 if an adjoint problem accepts an adjoint file.
  UncollidedGeneratorSetAdjointClearsSources: 1 if SetAdjoint(True) on the generator succeeds and
    generation then rejects the missing point source.
  UncollidedCtorBoundaryRejected, UncollidedLaterBoundaryRejected: 1 if an isotropic boundary is
    rejected with the uncollided-flux error at construction and when set afterwards.
"""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator
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
extruded_mesh = uncollided_utils.extruded_mesh
remove_file = uncollided_utils.remove_file

STEADY_ONLY = "uncollided flux is only supported for steady-state"
MODE_CHANGE = "transport mode cannot be changed after loading an uncollided flux file"
FORWARD_FILE_IN_ADJOINT = "was generated in forward mode, but the problem is in adjoint mode"
ADJOINT_FILE_IN_FORWARD = "was generated in adjoint mode, but the problem is in forward mode"


def make_generator(grid, xs, options=None):
    whole_domain = RPPLogicalVolume(xmin=-1.01, xmax=1.01, ymin=-1.01, ymax=1.01, infz=True)
    return UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=[0.13, -0.22, 1.05], strength=[1.0])],
        near_source=[whole_domain],
        scattering_order=0,
        options=options or {},
    )


def run_generator(problem, file_name):
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()
    remove_file(file_name)


def accepted(action):
    try:
        action()
    except ValueError:
        return False
    return True


def make_uncollided_file(grid, xs, file_name):
    problem = make_generator(grid, xs)
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
            point_sources=[PointSource(location=[0.13, -0.22, 1.05], strength=[1.0])],
            near_source=[away],
            scattering_order=0,
        ),
        "lies outside its near-source logical volume",
    )


def make_problem(grid, xs, file_name, **kwargs):
    options = {"save_angular_flux": True, "verbose_inner_iterations": False}
    options.update(kwargs.pop("options", {}))
    kwargs["uncollided_flux"] = file_name
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=4, n_azimuthal=12, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
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

    grid = extruded_mesh(
        ExtruderMeshGenerator, FromFileMeshGenerator, "triangle_mesh2x2.obj", num_layers=4
    )
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    # The velocity makes the cross sections valid for time-dependent mode, so the uncollided flux
    # is the only reason for the rejection.
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.7, c=0.3, velocity=1.0)

    file_name = "uncollided_mode_rejections.h5"
    remove_file(file_name)
    make_uncollided_file(grid, xs, file_name)

    ctor_td = rejected(lambda: make_problem(grid, xs, file_name, time_dependent=True), STEADY_ONLY)
    ctor_adj = rejected(
        lambda: make_problem(grid, xs, file_name, options={"adjoint": True}),
        FORWARD_FILE_IN_ADJOINT,
    )

    problem = make_problem(grid, xs, file_name)
    set_td = rejected(problem.SetTimeDependentMode, STEADY_ONLY)
    set_adj = rejected(lambda: problem.SetAdjoint(True), MODE_CHANGE)
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

    # An adjoint file holds the uncollided flux of a detector, which only an adjoint problem uses.
    adjoint_file_name = "uncollided_mode_rejections_adjoint.h5"
    remove_file(adjoint_file_name)
    adjoint_generator = make_generator(grid, xs, {"adjoint": True})
    solver_adjoint = UncollidedSolver(problem=adjoint_generator, file_name=adjoint_file_name)
    solver_adjoint.Initialize()
    solver_adjoint.Execute()
    adjoint_file_forward = rejected(
        lambda: make_problem(grid, xs, adjoint_file_name), ADJOINT_FILE_IN_FORWARD
    )
    adjoint_file_adjoint = accepted(
        lambda: make_problem(grid, xs, adjoint_file_name, options={"adjoint": True})
    )
    remove_file(adjoint_file_name)

    # SetAdjoint clears the point sources, so generation has nothing to trace.
    generator = make_generator(grid, xs)
    generator.SetAdjoint(True)
    generator_cleared = generator.IsAdjoint() and rejected(
        lambda: run_generator(generator, "uncollided_mode_rejections_generator.h5"),
        "at least one point source is required",
    )

    max_diff = float("inf")
    if unchanged:
        phi = solve(problem)
        phi_fresh = solve(make_problem(grid, xs, file_name))
        max_diff = max(abs(a - b) for a, b in zip(phi, phi_fresh))

    # An incoming-flux boundary would be silently ignored with an uncollided flux.
    boundary = [{"name": "xmin", "type": "isotropic", "group_strength": [1.0]}]
    ctor_boundary = rejected(
        lambda: make_problem(grid, xs, file_name, boundary_conditions=boundary),
        "incoming-flux",
    )
    later_boundary = rejected(
        lambda: problem.SetBoundaryOptions(boundary_conditions=boundary), "incoming-flux"
    )

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
        print(f"UncollidedAdjointFileForwardRejected={int(adjoint_file_forward)}")
        print(f"UncollidedAdjointFileAdjointAccepted={int(adjoint_file_adjoint)}")
        print(f"UncollidedGeneratorSetAdjointClearsSources={int(generator_cleared)}")
        print(f"UncollidedCtorBoundaryRejected={int(ctor_boundary)}")
        print(f"UncollidedLaterBoundaryRejected={int(later_boundary)}")
        sys.stdout.flush()
