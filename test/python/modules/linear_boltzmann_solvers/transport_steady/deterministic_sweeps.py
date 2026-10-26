#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""
Bitwise reproducibility of ``options={"deterministic_sweeps": True}``.

Without the option, the order in which a sweep adds the directions' contributions to the scalar
flux moments and outflows follows MPI message arrival (several CPU ranks) or GPU thread and stream
scheduling, so repeated solves differ in the last bits. With it, the moments and outflows are
summed from the stored angular flux in a fixed angle order.

Each case solves two independent copies of the same problem with the option and requires the
scalar flux on every rank, the boundary leakage, and (for ``case="keigen"``) the eigenvalue to
agree bitwise. It also solves the problem without the option and requires agreement to the
convergence level: the two modes differ only in the rounding of the same sums. The steady solves
are converged to an absolute GMRES tolerance of 1e-12, so their scalar flux and leakage agree to
about 1e-12 relative; the eigenvalues are converged to k_tol = 1e-10, and they agree to about that.

Agreement with the default mode does not show that both are correct, so each forward case also
requires the deterministic solution to conserve particles: production plus inflow equals
absorption plus outflow, relative to production plus inflow, with the CSDA deposition out of the
lowest group in the CSDA case and the inventory change over the step in the transient case. The
rebuilt moments give the absorption and the rebuilt outflows give the leakage, so the balance ties
the two rebuilt quantities together. It closes to the same convergence level. The adjoint case has
no such check: the balance table uses the forward absorption cross section, which does not
conserve the adjoint equation's group sum.

Cases: ``steady`` (3D fixed source), ``keigen`` (3D power iteration with CMFD), ``rz`` (RZ fixed
source), ``transient`` (two theta = 0.5 steps, so the rebuild uses theta psi^{n+1} + (1 - theta)
psi^n), ``csda`` (3D CSDA fixed source, which also rebuilds the energy-slope moments), and
``adjoint`` (3D adjoint fixed source). The option stores the angular flux without enabling
``save_angular_flux``, so the adjoint case also requires the stored angular flux to be reoriented
like that of ``save_angular_flux``: they must agree to the same round-off level.

Options: ``case``, ``sweep_type`` ("AAH" or "CBC"), ``use_gpus``, and ``fused_aao_sweep``.
"""

import hashlib
import math
import os
import struct
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    comm = MPI.COMM_WORLD

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DRZ, GLCProductQuadrature3DXYZ
    from pyopensn.mesh import OrthogonalMeshGenerator, KBAGraphPartitioner
    from pyopensn.solver import (
        CMFDAcceleration,
        DiscreteOrdinatesCurvilinearProblem,
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        SteadyStateSourceSolver,
        TransientSolver,
    )
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def get_option(name, default):
    return globals().get(name, default)


CASE = get_option("case", "steady")
SWEEP_TYPE = get_option("sweep_type", "AAH")
USE_GPUS = get_option("use_gpus", False)
FUSED = get_option("fused_aao_sweep", False)
TIME_STEP = 0.1


def load_xs():
    xs = MultiGroupXS()
    if CASE == "transient":
        # Two groups with velocities, as the time-dependent mode requires.
        xs.LoadFromOpenSn("../../../../assets/xs/simple_2g_downscatter_td.cxs")
    elif CASE == "csda":
        xs.LoadFromCEPXS(
            "../../../../assets/xs/cepxs_synthetic_csda_3g.bxslib", material_id=0, csda_format=True
        )
    else:
        xs.LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")
    return xs


def solve(deterministic, save_angular_flux=False):
    # A 5 cm cube (a 5 cm square in RZ): the eigenvalue iteration converges robustly on 1-4 ranks.
    nodes = [5.0 * i / 6.0 for i in range(7)]
    nx = 2 if size >= 2 else 1
    ny = size // nx
    partitioner = KBAGraphPartitioner(
        nx=nx, ny=ny, nz=1, xcuts=[0.5] if nx == 2 else [], ycuts=[0.5] if ny == 2 else []
    )
    if CASE == "rz":
        grid = OrthogonalMeshGenerator(
            node_sets=[nodes, nodes], partitioner=partitioner, coord_sys="cylindrical"
        ).Execute()
        quadrature = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)
        angle_aggregation = "azimuthal"
        # The axis is reflecting; a reflecting zmin exercises the delayed angular fluxes.
        boundary_conditions = [
            {"name": "rmin", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
        ]
    else:
        grid = OrthogonalMeshGenerator(
            node_sets=[nodes, nodes, nodes], partitioner=partitioner
        ).Execute()
        quadrature = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=1)
        angle_aggregation = "single"
        # A reflecting boundary exercises the delayed angular fluxes.
        boundary_conditions = [{"name": "xmin", "type": "reflecting"}]
    grid.SetUniformBlockID(0)

    num_groups = 3 if CASE == "csda" else 2
    inner = {"inner_linear_method": "petsc_gmres", "l_abs_tol": 1.0e-12, "l_max_its": 300}
    problem_args = dict(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": quadrature,
                "angle_aggregation_type": angle_aggregation,
                **inner,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": load_xs()}],
        boundary_conditions=boundary_conditions,
        options={
            "deterministic_sweeps": deterministic,
            "save_angular_flux": save_angular_flux or CASE == "transient",
            "csda_enabled": CASE == "csda",
            "adjoint": CASE == "adjoint",
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        sweep_type=SWEEP_TYPE,
    )
    if CASE != "keigen":
        problem_args["volumetric_sources"] = [
            VolumetricSource(block_ids=[0], group_strength=[1.0] + [0.0] * (num_groups - 1))
        ]
    if CASE == "transient":
        problem_args["time_dependent"] = True
    if CASE == "rz":
        problem = DiscreteOrdinatesCurvilinearProblem(**problem_args)
    else:
        problem_args["options"]["fused_aao_sweep"] = FUSED
        problem_args["use_gpus"] = USE_GPUS
        problem = DiscreteOrdinatesProblem(**problem_args)

    k_eff = None
    if CASE == "transient":
        solver = TransientSolver(problem=problem, dt=TIME_STEP, theta=0.5, initial_state="zero")
        solver.Initialize()
        solver.Advance()
        solver.Advance()
    elif CASE != "keigen":
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()
    else:
        cmfd = CMFDAcceleration(
            problem=problem, coarse_mesh="local_aggregation", aggregation_size=8
        )
        solver = PowerIterationKEigenSolver(
            problem=problem, acceleration=cmfd, k_tol=1.0e-10, max_iters=200
        )
        solver.Initialize()
        solver.Execute()
        k_eff = solver.GetEigenvalue()

    relative_balance = None if CASE == "adjoint" else relative_particle_balance(solver)

    phi = list(problem.GetPhiNewLocal())
    boundary_leakage = problem.ComputeLeakage([])
    leakage = [value for name in sorted(boundary_leakage) for value in boundary_leakage[name]]
    psi = [float(value) for value in problem.GetPsi()[0]] if CASE == "adjoint" else []
    return phi, leakage, k_eff, psi, relative_balance


def relative_particle_balance(solver):
    # Particle conservation of the rebuilt moments and outflows, relative to the particle gain.
    table = solver.ComputeBalanceTable()
    gain = table["production_rate"] + table["inflow_rate"]
    if CASE == "csda":
        # Includes the particles that slow down out of the lowest group; already relative.
        return table["csda_particle_balance"]
    if CASE == "transient":
        # Includes the storage term: inventory change minus dt times the time-centered net rate.
        return table["inventory_residual"] / (TIME_STEP * gain)
    return table["balance"] / gain


def digest(values):
    return hashlib.sha256(b"".join(struct.pack("<d", v) for v in values)).hexdigest()


def relative_l2(reference, candidate, distributed):
    diff_sq = sum((a - b) ** 2 for a, b in zip(reference, candidate))
    ref_sq = sum(a * a for a in reference)
    if distributed:
        diff_sq, ref_sq = MPIAllReduce(diff_sq, "sum"), MPIAllReduce(ref_sq, "sum")
    return math.sqrt(diff_sq) / max(math.sqrt(ref_sq), 1.0e-300)


if __name__ == "__main__":
    if CASE not in ("steady", "keigen", "rz", "transient", "csda", "adjoint"):
        sys.exit(
            f"Unknown case '{CASE}'. "
            "Expected 'steady', 'keigen', 'rz', 'transient', 'csda', or 'adjoint'."
        )

    first = solve(True)
    second = solve(True)
    default = solve(False)

    phi_match = MPIAllReduce(1 if digest(first[0]) == digest(second[0]) else 0, "min")
    leakage_match = 1 if digest(first[1]) == digest(second[1]) else 0
    phi_diff = relative_l2(first[0], default[0], distributed=True)
    leakage_diff = relative_l2(first[1], default[1], distributed=False)
    adjoint_psi_diff = None
    if CASE == "adjoint":
        saved = solve(False, save_angular_flux=True)
        adjoint_psi_diff = relative_l2(saved[3], first[3], distributed=True)

    if rank == 0:
        print(f"DETERMINISTIC_PHI_REPEAT_MATCH= {phi_match}")
        print(f"DETERMINISTIC_LEAKAGE_REPEAT_MATCH= {leakage_match}")
        print(f"DETERMINISTIC_VS_DEFAULT_PHI_REL_DIFF= {phi_diff:.6e}")
        print(f"DETERMINISTIC_VS_DEFAULT_LEAKAGE_REL_DIFF= {leakage_diff:.6e}")
        if CASE != "adjoint":
            print(f"DETERMINISTIC_RELATIVE_BALANCE= {first[4]:.6e}")
        if CASE == "keigen":
            k_match = 1 if first[2] == second[2] else 0
            print(f"DETERMINISTIC_K_REPEAT_MATCH= {k_match}")
            k_diff = abs(first[2] - default[2]) / default[2]
            print(f"DETERMINISTIC_VS_DEFAULT_K_REL_DIFF= {k_diff:.6e}")
        if CASE == "adjoint":
            print(f"DETERMINISTIC_VS_SAVED_ADJOINT_PSI_REL_DIFF= {adjoint_psi_diff:.6e}")
