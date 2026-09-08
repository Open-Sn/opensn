#!/usr/bin/env python3
"""Check CMFD inactive/rejected updates against PI and preserve WGS settings on reinitialization."""

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

    phys = make_problem()
    for option, value in (("relaxation", 0.0), ("pi_max_its", 0)):
        try:
            CMFDAcceleration(problem=phys, **{option: value})
        except (ValueError, RuntimeError) as error:
            assert option in str(error), str(error)
        else:
            raise AssertionError(f"CMFD accepted {option}=0, which can freeze k_eff")

    reference = execute(make_problem())
    phys = make_problem()
    cmfd = CMFDAcceleration(problem=phys, inactive_iterations=1)
    assert_same(execute(phys, cmfd), reference)
    # Reach the active phase, then initialize a new executor on the same problem/acceleration.
    execute(phys, cmfd, iterations=2)
    assert_same(execute(phys, cmfd), reference)

    reference = execute(make_problem(1, 1.0e-12), iterations=3)
    for rejection_options in (
        dict(max_iters=0, coarse_solver_policy="petsc_options",
             petsc_options="-CMFDAccelerationksp_type gmres -CMFDAccelerationpc_type jacobi"),
        dict(relaxation=100.0, correction_max_attempts=2, correction_min_damping=1.0),
    ):
        phys = make_problem()
        cmfd = CMFDAcceleration(problem=phys, **rejection_options)
        assert_same(execute(phys, cmfd, iterations=3), reference)

    if rank == 0:
        print("CMFD transport state checks passed")
