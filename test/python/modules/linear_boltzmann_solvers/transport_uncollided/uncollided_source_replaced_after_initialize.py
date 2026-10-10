#!/usr/bin/env python3
"""
UncollidedSolver.Execute() uses the point sources the problem has when it runs.

Replacing the point source between Initialize() and Execute() (same location, so the configuration
stays valid) must change the generated uncollided flux. Transport is linear in the source
strength, so the collided solve driven by a file generated after replacing a unit source with a
strength-2 source must be exactly twice the solve driven by a file generated with the unit source.
A generator that used a copy of the sources taken at Initialize() would give the same flux (a
relative error of 0.5).

Checks:
  UncollidedReplacedSourceRelErr: max |phi_replaced - 2 phi_unit| / max |2 phi_unit|.
  UncollidedUnitFluxNonzero: 1 if the unit-source flux is nonzero.
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
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
remove_file = uncollided_utils.remove_file

LOCATION = [0.13, -0.22, 0.0]


def make_generator(grid, xs, strength):
    whole_domain = RPPLogicalVolume(xmin=-1.01, xmax=1.01, ymin=-1.01, ymax=1.01, infz=True)
    return UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=LOCATION, strength=[strength])],
        near_source=[whole_domain],
        scattering_order=0,
    )


def collided_phi(grid, xs, file_name):
    problem = DiscreteOrdinatesProblem(
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
        options={"verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal())


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.7, c=0.3)

    unit_file = "uncollided_source_replaced_unit.h5"
    replaced_file = "uncollided_source_replaced_double.h5"
    remove_file(unit_file)
    remove_file(replaced_file)

    unit_solver = UncollidedSolver(problem=make_generator(grid, xs, 1.0), file_name=unit_file)
    unit_solver.Initialize()
    unit_solver.Execute()

    # Replace the unit source with a strength-2 source after Initialize().
    generator = make_generator(grid, xs, 1.0)
    replaced_solver = UncollidedSolver(problem=generator, file_name=replaced_file)
    replaced_solver.Initialize()
    generator.SetPointSources(
        clear_point_sources=True,
        point_sources=[PointSource(location=LOCATION, strength=[2.0])],
    )
    replaced_solver.Execute()

    phi_unit = collided_phi(grid, xs, unit_file)
    phi_replaced = collided_phi(grid, xs, replaced_file)
    scale = max(abs(2.0 * v) for v in phi_unit)
    unit_nonzero = scale > 0.0
    rel_err = (max(abs(a - 2.0 * b) for a, b in zip(phi_replaced, phi_unit)) / scale
               if unit_nonzero else float("inf"))

    remove_file(unit_file)
    remove_file(replaced_file)

    if rank == 0:
        print(f"UncollidedReplacedSourceRelErr={rel_err:.8e}")
        print(f"UncollidedUnitFluxNonzero={int(unit_nonzero)}")
        sys.stdout.flush()
