#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Problem configuration rules reject invalid input on every rank, and a rejected change leaves the
problem unchanged.

Each case builds an invalid problem and checks that construction raises a ValueError naming the
violated rule:

1. groupsets that omit the first group,
2. groupsets that have a gap,
3. groupsets that omit the last group,
4. cells without cross sections, on one rank only (the rule reduces across ranks, so every rank
   raises instead of only the rank that owns the cells),
5. cross sections with fewer groups than the problem,
6. restart writes without a restart path,
7. a groupset without an angular quadrature,
8. two violations at once (both are reported),
9. a spherical mesh and
10. an RZ mesh given to the Cartesian DiscreteOrdinatesProblem (only the curvilinear problem
   discretizes their angular-derivative terms),
11. azimuthal angle aggregation on an unstructured RZ mesh (which requires single aggregation),
12. a reflecting boundary that is not planar. Its faces lie on xmin and ymin and are owned by one
    rank, but the rule reduces across ranks so every rank raises.
13. a restart path whose directory is an existing file (only rank 0 creates the directory, and
    every rank must raise).

Finally, SetBoundaryOptions must reject a reflecting rmax boundary on an RZ problem (reflecting
axis), and a non-planar reflecting boundary on an XY problem, before the boundaries change: a
fixed-source solve afterwards must reproduce the solve of an identical fresh problem to round-off.

Checks:
  CONFIG_REJECTIONS_PASSED: number of cases (of 15) that behave as described.
  CONFIG_RZ_REJECTED_MAX_DIFF: max |phi - phi_fresh| after the rejected RZ boundary change.
  CONFIG_NONPLANAR_REJECTED_MAX_DIFF: max |phi - phi_fresh| after the rejected XY boundary change.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import (
        GLCProductQuadrature2DRZ,
        GLCProductQuadrature2DXY,
        GLProductQuadrature1DSlab,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator, OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesCurvilinearProblem,
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
    )
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def xy_problem(grid, xs, num_groups, groupsets, options=None):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=groupsets,
        xs_map=[{"block_ids": [0], "xs": xs}],
        options=options or {},
    )


def groupset(first, last, quadrature=True):
    block = {"groups_from_to": (first, last)}
    if quadrature:
        block["angular_quadrature"] = GLCProductQuadrature2DXY(
            n_polar=2, n_azimuthal=4, scattering_order=0
        )
    return block


def rejected(action, *expected):
    try:
        action()
    except ValueError as error:
        return all(text in str(error) for text in expected)
    return False


def rz_problem(grid, xs):
    return DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DRZ(
                    n_polar=2, n_azimuthal=4, scattering_order=0
                ),
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=[
            {"name": "rmin", "type": "reflecting"},
            {"name": "rmax", "type": "vacuum"},
        ],
        options={"verbose_inner_iterations": False},
    )


def xy_source_problem(grid, xs, boundary_conditions=None):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[dict(groupset(0, 0), inner_linear_method="petsc_gmres", l_abs_tol=1.0e-12)],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=boundary_conditions or [],
        options={"verbose_inner_iterations": False},
    )


def max_diff_from_fresh(problem, make_fresh):
    phi = solve(problem)
    phi_fresh = solve(make_fresh())
    local_diff = max((abs(a - b) for a, b in zip(phi, phi_fresh)), default=0.0)
    return MPIAllReduce(local_diff, "max")


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal())


if __name__ == "__main__":
    nodes = [i / 8 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs1 = MultiGroupXS()
    xs1.CreateSimpleOneGroup(1.0, 0.5)

    passed = 0
    passed += rejected(
        lambda: xy_problem(grid, xs1, 3, [groupset(1, 2)]),
        "first groupset must start at group 0",
    )
    passed += rejected(
        lambda: xy_problem(grid, xs1, 3, [groupset(0, 0), groupset(2, 2)]),
        "should start at group 1",
    )
    passed += rejected(
        lambda: xy_problem(grid, xs1, 3, [groupset(0, 1)]),
        "last groupset must end at group 2",
    )

    # Block 1 occupies only the cells near x = 1, which one rank owns.
    partial_grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    partial_grid.SetUniformBlockID(0)
    partial_grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.9, xmax=2.0, ymin=0.9, ymax=2.0, infz=True), 1, True
    )
    local_ok = rejected(
        lambda: xy_problem(partial_grid, xs1, 1, [groupset(0, 0)]), "invalid material id"
    )
    passed += int(MPIAllReduce(int(local_ok), "min"))

    passed += rejected(lambda: xy_problem(grid, xs1, 2, [groupset(0, 1)]), "fewer groups")
    passed += rejected(
        lambda: xy_problem(grid, xs1, 1, [groupset(0, 0)], {"restart_writes_enabled": True}),
        "non-empty `write_restart_path`",
    )
    passed += rejected(
        lambda: xy_problem(grid, xs1, 1, [groupset(0, 0, quadrature=False)]),
        "does not have an associated quadrature set",
    )
    passed += rejected(
        lambda: xy_problem(grid, xs1, 3, [groupset(0, 0), groupset(2, 2)],
                           {"restart_writes_enabled": True}),
        "invalid configuration",
        "should start at group 1",
        "non-empty `write_restart_path`",
    )

    rz_grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes], coord_sys="cylindrical").Execute()
    rz_grid.SetUniformBlockID(0)

    sphere_grid = OrthogonalMeshGenerator(node_sets=[nodes], coord_sys="spherical").Execute()
    sphere_grid.SetUniformBlockID(0)
    slab_groupset = {
        "groups_from_to": (0, 0),
        "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
    }
    passed += rejected(
        lambda: xy_problem(sphere_grid, xs1, 1, [slab_groupset]),
        "Spherical meshes require DiscreteOrdinatesCurvilinearProblem",
    )
    passed += rejected(
        lambda: xy_problem(rz_grid, xs1, 1, [groupset(0, 0)]),
        "Cylindrical meshes require DiscreteOrdinatesCurvilinearProblem",
    )

    unstructured_rz_grid = FromFileMeshGenerator(
        filename=os.path.join(os.path.dirname(__file__),
                              "../../../../assets/mesh/rz_annulus_single.msh"),
        coord_sys="cylindrical",
    ).Execute()
    unstructured_rz_grid.SetUniformBlockID(0)
    passed += rejected(
        lambda: DiscreteOrdinatesCurvilinearProblem(
            mesh=unstructured_rz_grid,
            num_groups=1,
            groupsets=[
                {
                    "groups_from_to": (0, 0),
                    "angular_quadrature": GLCProductQuadrature2DRZ(
                        n_polar=2, n_azimuthal=4, scattering_order=0
                    ),
                    "angle_aggregation_type": "azimuthal",
                }
            ],
            xs_map=[{"block_ids": [0], "xs": xs1}],
        ),
        'unstructured RZ meshes require angle_aggregation_type "single"',
    )

    # One boundary name on the faces near the corner (0, 0), on both xmin and ymin.
    corner_grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    corner_grid.SetUniformBlockID(0)
    corner_grid.SetBoundaryIDFromLogicalVolume(
        RPPLogicalVolume(xmin=-1.0, xmax=0.2, ymin=-1.0, ymax=0.2, infz=True), "corner", True
    )
    corner_reflecting = [{"name": "corner", "type": "reflecting"}]
    passed += rejected(
        lambda: xy_source_problem(corner_grid, xs1, corner_reflecting),
        'reflecting boundary "corner" is not planar',
    )

    # The restart directory would be an existing regular file.
    blocker = "problem_configuration_rejections_restart_blocker"
    if rank == 0:
        with open(blocker, "w") as blocker_file:
            blocker_file.write("not a directory\n")
    MPIAllReduce(0)
    passed += rejected(
        lambda: xy_problem(
            grid,
            xs1,
            1,
            [groupset(0, 0)],
            {"restart_writes_enabled": True, "write_restart_path": blocker + "/restart"},
        ),
        "could not create the restart directory",
    )
    MPIAllReduce(0)
    if rank == 0:
        os.remove(blocker)

    problem = rz_problem(rz_grid, xs1)
    passed += rejected(
        lambda: problem.SetBoundaryOptions(
            boundary_conditions=[{"name": "rmax", "type": "reflecting"}]
        ),
        "Reflecting boundary on rmax",
    )
    rz_max_diff = max_diff_from_fresh(problem, lambda: rz_problem(rz_grid, xs1))

    problem = xy_source_problem(corner_grid, xs1)
    passed += rejected(
        lambda: problem.SetBoundaryOptions(boundary_conditions=corner_reflecting),
        'reflecting boundary "corner" is not planar',
    )
    nonplanar_max_diff = max_diff_from_fresh(problem, lambda: xy_source_problem(corner_grid, xs1))

    if rank == 0:
        print(f"CONFIG_REJECTIONS_PASSED {passed}")
        print(f"CONFIG_RZ_REJECTED_MAX_DIFF {rz_max_diff:.6e}")
        print(f"CONFIG_NONPLANAR_REJECTED_MAX_DIFF {nonplanar_max_diff:.6e}")
