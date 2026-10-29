#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Null transients from k-eigenvalue and fixed-source steady states with delayed neutrons.

A k-eigenvalue solve leaves the eigenvector phi and the precursors C_j = gamma_j nu_d sigma_f phi /
(k lambda_j). With the fissionable cross sections replaced by ones with nu_p and nu_d divided by k,
phi is a critical eigenstate and the precursors are at equilibrium, so a transient from it must
stay stationary for any theta: the steady fission source (prompt plus delayed) and the transient
one (prompt, plus decay of the precursors and the precursors produced in the step) agree.

keigen cases: a 2D quarter core, fuel (two groups with upscatter, two precursor families whose
delayed spectra differ from each other and from the prompt spectrum) in [0, 1.5]^2 and a
non-fissionable reflector to 3, reflecting on xmin and ymin and vacuum on xmax and ymax. Power
iteration (AAH with theta = 1, CBC with theta = 1/2, AAH with one groupset per group) and the
nonlinear k-eigenvalue solver (AAH, theta = 1/2), with delayed neutrons, and power iteration
(AAH, theta = 1/2) with use_precursors=False, prompt neutrons only in the k solve and the
transient.

openmc case: the slab of keigenvalue_transport_1d_openmcxs_delayed_slab.py (OpenMC MGXS fuel with
six delayed groups and water, vacuum ends, k = 0.358) with a source in the fuel. Its steady state,
with the precursors from the steady-state solver, is stationary for theta = 1. The fast inverse
velocity of these data is zero.

Check:
  NULL_KEIGEN_MAX_REL_DRIFT: max over cases, steps and nodes of |phi^n - phi^0| / max phi^0.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLProductQuadrature1DSlab
    from pyopensn.source import VolumetricSource
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
        SteadyStateSourceSolver,
        TransientSolver,
    )

FUEL = {
    "sigma_t": [1.0, 2.0],
    # (from-group, to-group, value), with upscatter
    "transfers": [(0, 0, 0.6), (0, 1, 0.25), (1, 1, 1.5), (1, 0, 0.02)],
    "sigma_f": [0.05, 0.25],
    "nu_prompt": [2.5, 2.4],
    "nu_delayed": [0.05, 0.04],
    "chi_prompt": [0.9, 0.1],
    "inv_velocity": [0.1, 1.0],
    # (emission spectrum, fractional yield, decay constant)
    "families": [([0.5, 0.5], 0.4, 0.0124), ([0.1, 0.9], 0.6, 0.305)],
}
REFLECTOR = {
    "sigma_t": [0.8, 2.5],
    "transfers": [(0, 0, 0.6), (0, 1, 0.18), (1, 1, 2.45), (1, 0, 0.001)],
    "inv_velocity": [0.1, 1.0],
}
NUM_STEPS = 4


def load_xs(medium, tag, nu_factor=1.0):
    """Writes the medium, with nu_p and nu_d times nu_factor, as an OpenSn file and loads it."""
    families = medium.get("families", [])
    lines = ["NUM_GROUPS 2", "NUM_MOMENTS 1", f"NUM_PRECURSORS {len(families)}"]

    def block(name, values):
        lines.append(f"{name}_BEGIN")
        lines.extend(f"{g} {v!r}" for g, v in enumerate(values))
        lines.append(f"{name}_END")

    absorption = list(medium["sigma_t"])
    for gf, _, v in medium["transfers"]:
        absorption[gf] -= v
    block("SIGMA_T", medium["sigma_t"])
    block("SIGMA_A", absorption)
    block("INV_VELOCITY", medium["inv_velocity"])
    if families:
        block("SIGMA_F", medium["sigma_f"])
        block("NU_PROMPT", [nu_factor * v for v in medium["nu_prompt"]])
        block("NU_DELAYED", [nu_factor * v for v in medium["nu_delayed"]])
        block("CHI_PROMPT", medium["chi_prompt"])
        lines.append("CHI_DELAYED_BEGIN")
        for j, (spectrum, _, _) in enumerate(families):
            lines.extend(f"G_PRECURSOR_VAL {g} {j} {v!r}" for g, v in enumerate(spectrum))
        lines.append("CHI_DELAYED_END")
        block("PRECURSOR_DECAY_CONSTANTS", [decay for _, _, decay in families])
        block("PRECURSOR_FRACTIONAL_YIELDS", [yld for _, yld, _ in families])
    lines.append("TRANSFER_MOMENTS_BEGIN")
    lines.extend(f"M_GFROM_GTO_VAL 0 {gf} {gt} {v!r}" for gf, gt, v in medium["transfers"])
    lines.append("TRANSFER_MOMENTS_END")

    path = f"transient_null_transient_from_keigen_{tag}_r{rank}.xs"
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)
    return xs


def max_drift(problem, solver, phi0):
    """Advances NUM_STEPS steps; returns max |phi^n - phi^0| / max phi^0."""
    scale = MPIAllReduce(float(np.max(np.abs(phi0))) if phi0.size else 0.0, "max")
    drift = 0.0
    for _ in range(NUM_STEPS):
        solver.Advance()
        phi = np.array(problem.GetPhiNewLocal())
        local = float(np.max(np.abs(phi - phi0))) if phi.size else 0.0
        drift = max(drift, MPIAllReduce(local, "max") / scale)
    return drift


def keigen_drift(eigen_solver, sweep_type, theta, groupsets, use_precursors=True):
    tag = f"{eigen_solver}_{sweep_type}_{theta}_{groupsets}_{use_precursors}"
    nodes = [0.25 * i for i in range(13)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(1)
    fuel_region = RPPLogicalVolume(xmin=0.0, xmax=1.5, ymin=0.0, ymax=1.5, infz=True)
    grid.SetBlockIDFromLogicalVolume(fuel_region, 0, True)
    reflector = load_xs(REFLECTOR, tag + "_reflector")

    group_ranges = [[0, 1]] if groupsets == 1 else [[0, 0], [1, 1]]
    options = {
        "save_angular_flux": True,
        "use_precursors": use_precursors,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }
    if groupsets > 1:
        options["ags_tolerance"] = 1.0e-13
        options["max_ags_iterations"] = 1000
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        sweep_type=sweep_type,
        groupsets=[
            {
                "groups_from_to": group_range,
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=8, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 500,
            }
            for group_range in group_ranges
        ],
        xs_map=[
            {"block_ids": [0], "xs": load_xs(FUEL, tag + "_fuel")},
            {"block_ids": [1], "xs": reflector},
        ],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "xmax", "type": "vacuum"},
            {"name": "ymax", "type": "vacuum"},
        ],
        options=options,
    )
    if eigen_solver == "pi":
        solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-13, max_iters=2000)
    else:
        solver = NonLinearKEigenSolver(
            problem=problem, nl_abs_tol=1.0e-13, nl_rel_tol=1.0e-13, nl_max_its=200
        )
    solver.Initialize()
    solver.Execute()
    k = solver.GetEigenvalue()
    phi0 = np.array(problem.GetPhiNewLocal(), copy=True)

    problem.SetTimeDependentMode()
    transient = TransientSolver(
        problem=problem, dt=0.05, stop_time=1.0e9, theta=theta, initial_state="existing",
        verbose=False,
    )
    transient.Initialize()
    critical_fuel = load_xs(FUEL, tag + "_critical", nu_factor=1.0 / k)
    problem.SetXSMap(
        xs_map=[{"block_ids": [0], "xs": critical_fuel}, {"block_ids": [1], "xs": reflector}]
    )
    return max_drift(problem, transient, phi0)


def openmc_drift():
    xs_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../transport_keigen/slab_delayed_mgxs.h5"
    )
    xs_fuel = MultiGroupXS()
    xs_fuel.LoadFromOpenMC(xs_file, "fuel", 294.0)
    xs_water = MultiGroupXS()
    xs_water.LoadFromOpenMC(xs_file, "water", 294.0)

    nodes = [-2.0 + 0.2 * i for i in range(61)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(1)
    fuel_region = RPPLogicalVolume(infx=True, infy=True, zmin=-2.0, zmax=2.0)
    grid.SetBlockIDFromLogicalVolume(fuel_region, 0, True)
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 500,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_fuel}, {"block_ids": [1], "xs": xs_water}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])],
        boundary_conditions=[
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    steady = SteadyStateSourceSolver(problem=problem)
    steady.Initialize()
    steady.Execute()
    phi0 = np.array(problem.GetPhiNewLocal(), copy=True)

    problem.SetTimeDependentMode()
    transient = TransientSolver(
        problem=problem, dt=0.05, stop_time=1.0e9, theta=1.0, initial_state="existing",
        verbose=False,
    )
    transient.Initialize()
    return max_drift(problem, transient, phi0)


if __name__ == "__main__":
    drift = 0.0
    for case in (("pi", "AAH", 1.0, 1), ("pi", "CBC", 0.5, 1), ("pi", "AAH", 1.0, 2),
                 ("nlke", "AAH", 0.5, 1), ("pi", "AAH", 0.5, 1, False)):
        case_drift = keigen_drift(*case)
        if rank == 0:
            print(f"NULL_KEIGEN_CASE {case} drift {case_drift:.6e}")
        drift = max(drift, case_drift)
    openmc = openmc_drift()
    if rank == 0:
        print(f"NULL_KEIGEN_CASE openmc drift {openmc:.6e}")
    drift = max(drift, openmc)

    if rank == 0:
        print(f"NULL_KEIGEN_MAX_REL_DRIFT {drift:.6e}")
