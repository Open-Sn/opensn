#!/usr/bin/env python3
"""Check the CMFD input contract and that skipped CMFD corrections reproduce unaccelerated PI."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration


if __name__ == "__main__":
    # Load cross-section data for a 1-group fissile material
    num_groups = 1
    xs_simple_fissile = MultiGroupXS()
    xs_simple_fissile.LoadFromOpenSn("../../../../assets/xs/simple_fissile.xs")

    # Angle information
    n_angles = 32
    scat_order = 0

    def make_problem(wgs_max_its=50, wgs_abs_tol=1.0e-10):
        L = 100.0
        n_cells = 50
        dx = L / n_cells
        nodes = [i * dx for i in range(n_cells + 1)]
        meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
        grid = meshgen.Execute()
        grid.SetUniformBlockID(0)

        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=num_groups,
            groupsets=[
                {
                    "groups_from_to": (0, num_groups - 1),
                    "angular_quadrature": GLProductQuadrature1DSlab(
                        n_polar=n_angles,
                        scattering_order=scat_order
                    ),
                    "inner_linear_method": "petsc_richardson",
                    "l_max_its": wgs_max_its,
                    "l_abs_tol": wgs_abs_tol,
                }
            ],
            xs_map=[
                {
                    "block_ids": [0],
                    "xs": xs_simple_fissile,
                }
            ],
            options={
                "use_precursors": False,
                "verbose_inner_iterations": False,
                "verbose_outer_iterations": False,
            },
            sweep_type="AAH",
        )

    def execute(phys, acceleration=None, iterations=1):
        kwargs = dict(problem=phys, max_iters=iterations, k_tol=1.0e-10)
        if acceleration is not None:
            kwargs["acceleration"] = acceleration
        solver = PowerIterationKEigenSolver(**kwargs)
        solver.Initialize()
        before = solver.GetNumSweeps()
        solver.Execute()
        return solver.GetEigenvalue(), solver.GetNumSweeps() - before

    def assert_same(actual, expected):
        assert abs(actual[0] - expected[0]) < 1.0e-12, (actual, expected)
        assert actual[1] == expected[1], (actual, expected)

    def assert_rejected(description, construct):
        try:
            construct()
        except (ValueError, RuntimeError, TypeError) as error:
            return str(error)
        raise AssertionError(f"CMFD accepted {description}")

    # Values that contradict the fixed method, and malformed coarse-mesh and coarse-group input,
    # are rejected when the accelerator is constructed.
    phys = make_problem()
    for description, options in (
        ("the net closure", dict(current_closure="net")),
        ("a relaxation other than 1", dict(relaxation=0.5)),
        ("inactive iterations", dict(inactive_iterations=1)),
        ("an unknown coarse-mesh type", dict(coarse_mesh="octree")),
        ("an unknown coarse-mesh key",
         dict(coarse_mesh={"type": "local_aggregation", "pitch": [1.0, 1.0]})),
        ("a lattice without its geometry", dict(coarse_mesh="reactor_lattice")),
        ("an aggregation size given twice",
         dict(coarse_mesh={"type": "local_aggregation", "size": 4}, aggregation_size=4)),
        ("zero coarse groups", dict(coarse_groups=0)),
        ("a coarse-group string other than auto", dict(coarse_groups="all")),
        ("both coarse_groups and group_aggregation_size",
         dict(coarse_groups=1, group_aggregation_size=1)),
    ):
        assert_rejected(description, lambda: CMFDAcceleration(problem=phys, **options))

    # The solver rejects sweeps_per_update without CMFD.
    for description, options in (
        ("sweeps_per_update without CMFD", dict(sweeps_per_update=2)),
    ):
        def construct_and_initialize():
            solver = PowerIterationKEigenSolver(problem=phys, **options)
            solver.Initialize()
        assert_rejected(description, construct_and_initialize)

    # Each deprecated name maps to its replacement: the same coarse mesh, coarse groups and
    # transport sweeps give the same iterates.
    def cmfd_result(cmfd_options, solver_options=None):
        phys = make_problem()
        cmfd = CMFDAcceleration(problem=phys, **cmfd_options)
        solver = PowerIterationKEigenSolver(problem=phys, acceleration=cmfd, max_iters=6,
                                            k_tol=1.0e-10, **(solver_options or {}))
        solver.Initialize()
        solver.Execute()
        return solver.GetEigenvalue(), solver.GetNumSweeps(), cmfd.GetNumCoarseCells()

    for deprecated, replacement, solver_options in (
        (dict(aggregation_size=5), dict(coarse_mesh={"type": "local_aggregation", "size": 5}),
         None),
        (dict(group_aggregation_size=1), dict(coarse_groups=1), None),
        (dict(update_wgs_max_its=2), dict(), dict(sweeps_per_update=2)),
    ):
        old = cmfd_result(deprecated)
        new = cmfd_result(replacement, solver_options)
        assert abs(old[0] - new[0]) < 1.0e-12 and old[1:] == new[1:], (deprecated, old, new)

    # A correction whose coarse linear solve fails is skipped, so the iterates are those of the
    # unaccelerated power iteration with the same number of sweeps per update.
    reference = execute(make_problem(1, 1.0e-12), iterations=3)
    for petsc_options in (
        "-CMFDAccelerationksp_type gmres -CMFDAccelerationksp_max_it 0 "
        "-CMFDAccelerationpc_type jacobi",
        # PETSc reports this one-iteration Jacobi solve as converged because it computes no norm.
        # Its positive but inaccurate solution passes the other correction safeties, so only the
        # true-residual check rejects it. Without that check the result differs from the
        # unaccelerated reference.
        "-CMFDAccelerationksp_type richardson -CMFDAccelerationksp_max_it 1 "
        "-CMFDAccelerationksp_norm_type none -CMFDAccelerationpc_type jacobi",
    ):
        phys = make_problem()
        cmfd = CMFDAcceleration(problem=phys, petsc_options=petsc_options)
        assert_same(execute(phys, cmfd, iterations=3), reference)
        assert cmfd.GetNumSkippedCorrections() == 3, cmfd.GetNumSkippedCorrections()

    # petsc_options configure only the accelerator that was given them. An accelerator with the
    # same name created afterwards, without options, applies its corrections.
    phys = make_problem()
    cmfd = CMFDAcceleration(problem=phys)
    execute(phys, cmfd, iterations=3)
    assert cmfd.GetNumSkippedCorrections() == 0, cmfd.GetNumSkippedCorrections()

    if rank == 0:
        print("CMFD transport state checks passed")
