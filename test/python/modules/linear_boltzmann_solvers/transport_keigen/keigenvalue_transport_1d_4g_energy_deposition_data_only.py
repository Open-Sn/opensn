#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
ENERGY_DEPOSITION in OpenSn cross-section files is data only: it must not change the solve.

A 1D slab [0, 100] cm with 100 cells has fuel in [0, 50] cm (reflecting at z = 0) and a reflector
in [50, 100] cm (vacuum at z = 100). The synthetic four-group data have two neutron groups (0 fast,
1 thermal) and two photon groups (2, 3) fed by neutron-to-photon transfers, solved by power
iteration with one groupset for the neutrons and one for the photons. The same problem is solved
twice, from files without and with an ENERGY_DEPOSITION section (otherwise identical).

Checks:
  KHET_K_REL_DIFF, KHET_PHI_MAX_REL_DIFF: the relative difference of k and the largest difference
    of any flux moment, over the largest flux moment, between the two solves (zero).
  KHET_ITS_DIFF: the difference of their power iteration counts (zero).
  KHET_EDEP_FUEL_REL_ERR, KHET_EDEP_REFL_REL_ERR: the average of the "energy_deposition" field
    function over each region against sum_g ED_g times the region average of phi_g, with the ED_g
    of the file (agreement to roundoff).
  KHET_EDEP_ABSENT_MAX: the largest magnitude of that field function without the section (zero).
  KHET_DETERMINISTIC_SWEEPS: 1 when both solves use the deterministic_sweeps option, else 0.

On several ranks the default sweeps add the flux moments in message-arrival order, so two solves of
the same problem can differ in the last bits. The two-rank entry passes deterministic_sweeps=True,
which makes the comparison exact; on one rank the default sweeps are already reproducible.
"""

import math
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
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

# Set to True on the command line (deterministic_sweeps=True) for bitwise reproducible sweeps.
DETERMINISTIC_SWEEPS = globals().get("deterministic_sweeps", False)
if not isinstance(DETERMINISTIC_SWEEPS, bool):
    raise TypeError("deterministic_sweeps must be True or False")

NUM_GROUPS = 4

# Transfers (moment 0) as (from, to, value).
PHOTON_TRANSFERS = [(2, 2, 0.005), (2, 3, 0.01), (3, 3, 0.02)]
MATERIALS = {
    "fuel": {
        "sigma_t": [0.20, 0.80, 0.05, 0.10],
        "sigma_a": [0.015, 0.10, 0.035, 0.08],
        "transfers": [(0, 0, 0.17), (0, 1, 0.015), (1, 1, 0.70)] + PHOTON_TRANSFERS
        + [(0, 2, 0.002), (0, 3, 0.006), (1, 2, 0.02), (1, 3, 0.06)],
        "fission": {"sigma_f": [0.004, 0.08, 0.0, 0.0], "nu": [2.43, 2.43, 0.0, 0.0],
                    "chi": [1.0, 0.0, 0.0, 0.0]},
        "energy_deposition": [0.01, 0.005, 0.002, 0.001],
    },
    "reflector": {
        "sigma_t": [0.25, 1.20, 0.05, 0.10],
        "sigma_a": [0.002, 0.02, 0.035, 0.08],
        "transfers": [(0, 0, 0.22), (0, 1, 0.028), (1, 1, 1.18)] + PHOTON_TRANSFERS
        + [(1, 2, 0.004), (1, 3, 0.012)],
        "fission": None,
        "energy_deposition": [0.001, 0.002, 0.002, 0.001],
    },
}


def section(name, values):
    return ([f"{name}_BEGIN"] + [f"{g} {v!r}" for g, v in enumerate(values)]
            + [f"{name}_END"])


def xs_text(material, with_energy_deposition):
    data = MATERIALS[material]
    lines = [f"NUM_GROUPS {NUM_GROUPS}", "NUM_MOMENTS 1"]
    lines += section("SIGMA_T", data["sigma_t"]) + section("SIGMA_A", data["sigma_a"])
    lines.append("TRANSFER_MOMENTS_BEGIN")
    lines += [f"M_GFROM_GTO_VAL 0 {gf} {gt} {v!r}" for gf, gt, v in data["transfers"]]
    lines.append("TRANSFER_MOMENTS_END")
    if data["fission"] is not None:
        lines += section("SIGMA_F", data["fission"]["sigma_f"])
        lines += section("NU", data["fission"]["nu"]) + section("CHI", data["fission"]["chi"])
    if with_energy_deposition:
        lines += section("ENERGY_DEPOSITION", data["energy_deposition"])
    return "\n".join(lines) + "\n"


def load_xs(material, with_energy_deposition):
    # Each rank reads its own copy of the file.
    path = f"khet_{material}_{int(with_energy_deposition)}_r{rank}.xs"
    with open(path, "w") as f:
        f.write(xs_text(material, with_energy_deposition))
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)
    return xs


def volume_value(ff, op_type, logical_volume):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType(op_type)
    ffi.SetLogicalVolume(logical_volume)
    ffi.AddFieldFunction(ff)
    ffi.Execute()
    return ffi.GetValue()


def solve(with_energy_deposition):
    nodes = [float(i) for i in range(101)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(1)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=-1.0, zmax=50.0), 0, True)
    quad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=0, operator_method="standard")
    groupset = {
        "angular_quadrature": quad,
        "inner_linear_method": "petsc_gmres",
        "l_abs_tol": 1.0e-12,
        "l_max_its": 1000,
        "gmres_restart_interval": 30,
    }
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=NUM_GROUPS,
        groupsets=[dict(groupset, groups_from_to=(0, 1)), dict(groupset, groups_from_to=(2, 3))],
        xs_map=[
            {"block_ids": [0], "xs": load_xs("fuel", with_energy_deposition)},
            {"block_ids": [1], "xs": load_xs("reflector", with_energy_deposition)},
        ],
        boundary_conditions=[{"name": "zmin", "type": "reflecting"}],
        options={
            "use_precursors": False,
            "ags_convergence_check": "l2",
            "ags_tolerance": 1.0e-12,
            "max_ags_iterations": 100,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "deterministic_sweeps": DETERMINISTIC_SWEEPS,
        },
    )
    solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=1000,
                                        reset_phi0=True)
    solver.Initialize()
    solver.Execute()
    return problem, solver.GetEigenvalue(), solver.GetNumPowerIterations()


def energy_deposition_rel_err(problem, material, zmin, zmax):
    vol = RPPLogicalVolume(infx=True, infy=True, zmin=zmin, zmax=zmax)
    edep = problem.CreateFieldFunction(f"edep_{material}", "energy_deposition")
    value = volume_value(edep, "avg", vol)
    phi = problem.GetScalarFluxFieldFunction(only_scalar_flux=True)
    reference = sum(ed * volume_value(phi[g], "avg", vol)
                    for g, ed in enumerate(MATERIALS[material]["energy_deposition"]))
    return abs(value / reference - 1.0)


if __name__ == "__main__":
    problem_0, k_0, its_0 = solve(False)
    problem_1, k_1, its_1 = solve(True)

    phi_0 = list(problem_0.GetPhiNewLocal())
    phi_1 = list(problem_1.GetPhiNewLocal())
    local_diff = 0.0
    if len(phi_0) != len(phi_1):
        local_diff = math.inf
    else:
        for a, b in zip(phi_0, phi_1):
            # A NaN must fail the check, and max() would skip it.
            d = abs(a - b)
            local_diff = max(local_diff, d) if math.isfinite(d) else math.inf
    phi_diff = MPIAllReduce(local_diff, "max")
    phi_max = MPIAllReduce(max(abs(a) for a in phi_0), "max")

    fuel_err = energy_deposition_rel_err(problem_1, "fuel", -1.0, 50.0)
    refl_err = energy_deposition_rel_err(problem_1, "reflector", 50.0, 101.0)

    whole = RPPLogicalVolume(infx=True, infy=True, infz=True)
    absent = problem_0.CreateFieldFunction("edep_absent", "energy_deposition")
    absent_max = max(abs(volume_value(absent, "max", whole)),
                     abs(volume_value(absent, "min", whole)))

    if rank == 0:
        print(f"KHET_K_WITHOUT={k_0:.17e}")
        print(f"KHET_K_WITH={k_1:.17e}")
        print(f"KHET_ITS_WITHOUT={its_0}")
        print(f"KHET_ITS_WITH={its_1}")
        print(f"KHET_K_REL_DIFF={abs(k_1 / k_0 - 1.0):.17e}")
        print(f"KHET_ITS_DIFF={its_1 - its_0}")
        print(f"KHET_PHI_MAX_REL_DIFF={phi_diff / phi_max:.17e}")
        print(f"KHET_EDEP_FUEL_REL_ERR={fuel_err:.17e}")
        print(f"KHET_EDEP_REFL_REL_ERR={refl_err:.17e}")
        print(f"KHET_EDEP_ABSENT_MAX={absent_max:.17e}")
        print(f"KHET_DETERMINISTIC_SWEEPS={int(DETERMINISTIC_SWEEPS)}")
