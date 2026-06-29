#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Regression for runtime boundary-condition mutation.

Each case is first solved with a freshly constructed problem to establish the
reference answer. A single reusable problem is then mutated through the same
boundary-condition sequence and must reproduce the fresh answers.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    comm = MPI.COMM_WORLD

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN, "bor": MPI.BOR}

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op=_mpi_ops[op])

    def MPIBarrier():
        comm.Barrier()
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.math import AngularFluxFunction
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume


def require_procs(expected):
    if size != expected:
        sys.exit(f"Incorrect number of processors. Expected {expected} processors but got {size}.")


def build_mesh():
    nx = 10
    ny = 8
    lx = 2.5
    ly = 2.0
    nodes_x = [lx * i / nx for i in range(nx + 1)]
    nodes_y = [ly * i / ny for i in range(ny + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes_x, nodes_y]).Execute()
    grid.SetUniformBlockID(0)
    return grid


def build_xs():
    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/simple_upscatter.xs")
    return xs


def build_source():
    src_vol = RPPLogicalVolume(xmin=0.65, xmax=1.55, ymin=0.55, ymax=1.45, infz=True)
    return VolumetricSource(logical_volume=src_vol, group_strength=[1.0, 0.25, 0.05])


def build_quadrature():
    return GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=0)


def build_arbitrary_bc(pquad, num_groups):
    mask = np.zeros(len(pquad.omegas), dtype=bool)
    normalization = 0.0

    for n, (omega, weight) in enumerate(zip(pquad.omegas, pquad.weights)):
        if omega.x > 0.0 and omega.y > 0.0:
            mask[n] = True
            normalization += omega.x * weight * num_groups

    if normalization <= 0.0:
        raise RuntimeError("Arbitrary boundary setup found no incoming directions.")

    strength = 0.35 / normalization

    def incident_flux(group_index, angle_index):
        if mask[angle_index]:
            return strength * (1.0 + 0.15 * group_index)
        return 0.0

    return AngularFluxFunction(incident_flux)


def make_boundary_cases(arbitrary_bc):
    vacuum = [
        {"name": "xmin", "type": "vacuum"},
        {"name": "xmax", "type": "vacuum"},
        {"name": "ymin", "type": "vacuum"},
        {"name": "ymax", "type": "vacuum"},
    ]
    isotropic_drive = [
        {"name": "xmin", "type": "isotropic", "group_strength": [0.20, 0.05, 0.01]},
        {"name": "xmax", "type": "vacuum"},
        {"name": "ymin", "type": "vacuum"},
        {"name": "ymax", "type": "vacuum"},
    ]
    mixed_isotropic_vacuum_reflecting = [
        {"name": "xmin", "type": "isotropic", "group_strength": [0.12, 0.03, 0.02]},
        {"name": "xmax", "type": "vacuum"},
        {"name": "ymin", "type": "reflecting"},
        {"name": "ymax", "type": "vacuum"},
    ]
    all_reflecting = [
        {"name": "xmin", "type": "reflecting"},
        {"name": "xmax", "type": "reflecting"},
        {"name": "ymin", "type": "reflecting"},
        {"name": "ymax", "type": "reflecting"},
    ]
    return [
        {
            "name": "all_vacuum_explicit",
            "clear": True,
            "boundaries": vacuum,
        },
        {
            "name": "isotropic_drive",
            "clear": False,
            "boundaries": isotropic_drive,
        },
        {
            "name": "arbitrary_and_isotropic",
            "clear": False,
            "boundaries": [
                {"name": "xmin", "type": "arbitrary", "function": arbitrary_bc},
                {"name": "xmax", "type": "vacuum"},
                {"name": "ymin", "type": "vacuum"},
                {"name": "ymax", "type": "isotropic", "group_strength": [0.02, 0.03, 0.04]},
            ],
        },
        {
            "name": "reflecting_axes",
            "clear": False,
            "boundaries": [
                {"name": "xmin", "type": "reflecting"},
                {"name": "xmax", "type": "vacuum"},
                {"name": "ymin", "type": "reflecting"},
                {"name": "ymax", "type": "vacuum"},
            ],
        },
        {
            "name": "all_reflecting",
            "clear": False,
            "boundaries": all_reflecting,
        },
        {
            "name": "clear_then_partial_default_vacuum",
            "clear": True,
            "boundaries": [
                {"name": "xmin", "type": "isotropic", "group_strength": [0.10, 0.00, 0.00]},
            ],
        },
        {
            "name": "clear_to_default_vacuum",
            "clear": True,
            "boundaries": [],
        },
        {
            "name": "return_to_isotropic_drive",
            "clear": True,
            "boundaries": isotropic_drive,
        },
        {
            "name": "return_to_all_vacuum_explicit",
            "clear": True,
            "boundaries": vacuum,
        },
        {
            "name": "final_mixed_isotropic_vacuum_reflecting",
            "clear": True,
            "boundaries": mixed_isotropic_vacuum_reflecting,
        },
        {
            "name": "final_all_vacuum_after_mixed",
            "clear": True,
            "boundaries": vacuum,
        },
        {
            "name": "final_all_reflecting_after_vacuum",
            "clear": True,
            "boundaries": all_reflecting,
        },
    ]


def build_problem(grid, xs, source, pquad, boundaries, sweep_type):
    num_groups = xs.num_groups
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, num_groups - 1],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-9,
                "l_max_its": 300,
                "gmres_restart_interval": 50,
                "apply_wgdsa": True,
                "apply_tgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-7,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=boundaries,
        volumetric_sources=[source],
        options={
            "max_ags_iterations": 1,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        sweep_type=sweep_type,
    )


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    leakage = problem.ComputeLeakage(["xmin", "xmax", "ymin", "ymax"])
    return {
        "phi": list(problem.GetPhiNewLocal()),
        "leakage": {name: [float(v) for v in leakage[name]] for name in leakage},
    }


def compare_results(ref, got):
    phi_diff = 0.0
    for a, b in zip(ref["phi"], got["phi"]):
        phi_diff = max(phi_diff, abs(a - b))

    leakage_diff = 0.0
    for name in ref["leakage"]:
        for a, b in zip(ref["leakage"][name], got["leakage"][name]):
            leakage_diff = max(leakage_diff, abs(a - b))

    local = (phi_diff, leakage_diff)
    return tuple(MPIAllReduce(value, "max") for value in local)


if __name__ == "__main__":
    require_procs(4)

    grid = build_mesh()
    xs = build_xs()
    source = build_source()
    pquad = build_quadrature()
    boundary_cases = make_boundary_cases(build_arbitrary_bc(pquad, xs.num_groups))

    max_case_diff = 0.0
    for sweep_type in ["AAH", "CBC"]:
        references = []
        for case in boundary_cases:
            ref_problem = build_problem(grid, xs, source, pquad, case["boundaries"], sweep_type)
            references.append(solve(ref_problem))

        mutable_problem = build_problem(
            grid, xs, source, pquad, boundary_cases[0]["boundaries"], sweep_type
        )
        mutable_solver = SteadyStateSourceSolver(problem=mutable_problem)
        mutable_solver.Initialize()

        for i, case in enumerate(boundary_cases):
            if i > 0:
                mutable_problem.SetBoundaryOptions(
                    clear_boundary_conditions=case["clear"],
                    boundary_conditions=case["boundaries"],
                )

            mutable_problem.ZeroPhi()
            mutable_solver.Execute()
            mutated = {
                "phi": list(mutable_problem.GetPhiNewLocal()),
                "leakage": {
                    name: [float(v) for v in values]
                    for name, values in mutable_problem.ComputeLeakage(
                        ["xmin", "xmax", "ymin", "ymax"]
                    ).items()
                },
            }
            phi_diff, leakage_diff = compare_results(references[i], mutated)
            case_diff = max(phi_diff, leakage_diff)
            max_case_diff = max(max_case_diff, case_diff)

            if rank == 0:
                print(f"BOUNDARY_CASE_{sweep_type}_{i}_{case['name']}_MAX_DIFF={case_diff:.8e}")
                print(f"BOUNDARY_CASE_{sweep_type}_{i}_{case['name']}_PHI_DIFF={phi_diff:.8e}")
                print(
                    f"BOUNDARY_CASE_{sweep_type}_{i}_{case['name']}_LEAKAGE_DIFF="
                    f"{leakage_diff:.8e}"
                )

    match_tol = 5.0e-8
    match = 1 if max_case_diff <= match_tol else 0
    if rank == 0:
        print(f"BOUNDARY_MUTATION_MAX_DIFF={max_case_diff:.8e}")
        print(f"BOUNDARY_MUTATION_MATCH_TOL={match_tol:.8e}")
        print(f"BOUNDARY_MUTATION_MATCH={match}")
