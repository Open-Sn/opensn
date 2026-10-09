#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
Tiny source-driven GPU transport input for compute-sanitizer memcheck.

compute-sanitizer memcheck adds a large fixed cost to every kernel launch, and
a sweep launches kernels per angle set and per dependency level. This input
keeps the mesh, angle count, and iteration count small while still exercising
the GPU sweep with multiple groups, P1 scattering, fission, a volumetric
source, and one of three boundary configurations selected by ``bc_case``:

- ``source``: vacuum boundaries plus an isotropic incident flux on xmin.
- ``reflecting``: a single reflecting boundary (xmin), vacuum elsewhere.
- ``opposing_reflecting``: reflecting xmin/xmax and ymin/ymax, which creates
  cyclic boundary dependencies between angle sets.

The inner solve is capped at two GMRES iterations, so the flux is not
converged. Cycle breaking is disabled for CBC so CPU and GPU use the same
algebraic iteration.
The same capped iteration sequence is run on the CPU and on the GPU
and the scalar-flux vectors are compared. With identical algebra the two agree
to round-off. This is a CPU/GPU equivalence check of the GPU sweep, not an
accuracy check; accuracy of the converged solution is covered by the other
transport regression tests.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    comm = MPI.COMM_WORLD

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX}

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


BOUNDARY_CASES = {
    "source": [
        {"name": "xmin", "type": "isotropic", "group_strength": [1.0, 0.5]},
    ],
    "reflecting": [
        {"name": "xmin", "type": "reflecting"},
    ],
    "opposing_reflecting": [
        {"name": "xmin", "type": "reflecting"},
        {"name": "xmax", "type": "reflecting"},
        {"name": "ymin", "type": "reflecting"},
        {"name": "ymax", "type": "reflecting"},
    ],
}


def get_option(name, default):
    return globals().get(name, default)


def solve(boundary_conditions, sweep_type, use_gpus):
    nodes = [i / 3.0 for i in range(4)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=4,
                    n_azimuthal=8,
                    scattering_order=1,
                ),
                "inner_linear_method": "petsc_gmres",
                # Unreachable on purpose: the solve always stops at l_max_its.
                "l_abs_tol": 1.0e-14,
                "l_max_its": 2,
                "allow_cycles": sweep_type != "CBC",
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])],
        boundary_conditions=boundary_conditions,
        sweep_type=sweep_type,
        use_gpus=use_gpus,
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return list(problem.GetPhiNewLocal())


def compare(reference, candidate):
    """Return the global max-abs and relative-L2 differences of two local flux vectors."""
    if len(reference) != len(candidate):
        raise RuntimeError("CPU and GPU scalar-flux vectors have different sizes.")

    local_max_abs = max((abs(a - b) for a, b in zip(reference, candidate)), default=0.0)
    local_diff_sq = sum((a - b) ** 2 for a, b in zip(reference, candidate))
    local_ref_sq = sum(value * value for value in reference)

    max_abs = MPIAllReduce(local_max_abs, "max")
    diff_norm = math.sqrt(MPIAllReduce(local_diff_sq, "sum"))
    ref_norm = math.sqrt(MPIAllReduce(local_ref_sq, "sum"))
    if not math.isfinite(ref_norm) or ref_norm <= 0.0:
        raise RuntimeError(f"CPU reference scalar flux has invalid norm {ref_norm}.")
    return max_abs, diff_norm / ref_norm


if __name__ == "__main__":
    if size != 2:
        sys.exit(f"Incorrect number of processors. Expected 2 processors but got {size}.")

    bc_case = get_option("bc_case", "source")
    sweep_type = get_option("sweep_type", "AAH")
    if bc_case not in BOUNDARY_CASES:
        sys.exit(f"Unknown bc_case '{bc_case}'. Expected one of {sorted(BOUNDARY_CASES)}.")

    boundary_conditions = BOUNDARY_CASES[bc_case]
    cpu_phi = solve(boundary_conditions, sweep_type, False)
    gpu_phi = solve(boundary_conditions, sweep_type, True)
    max_abs, relative_l2 = compare(cpu_phi, gpu_phi)

    if rank == 0:
        print(f"GPU_CPU_PHI_MAX_ABS_DIFF= {max_abs:.16e}")
        print(f"GPU_CPU_PHI_REL_L2_DIFF= {relative_l2:.16e}")
