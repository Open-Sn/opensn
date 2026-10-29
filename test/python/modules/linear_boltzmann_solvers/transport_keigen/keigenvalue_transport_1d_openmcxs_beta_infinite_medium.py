#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenMC MGXS delayed-neutron data in OpenMC's beta form, in an infinite medium.

beta_form_mgxs.h5 (written by make_beta_form_mgxs.py, whose DATA this test uses) gives a
two-group fuel with upscatter by nu-fission, chi and beta[j][g], the delayed fraction of each of
two delayed groups by fissioning group. All neutrons are born with chi, so the production matrix is
F = chi (nu sigma_f)^T with delayed neutrons and F_p = chi ((1 - sum_j beta_j) nu sigma_f)^T
without. A 1D slab with reflecting faces is an infinite medium, so with L = diag(sigma_t) - S
(S row: to-group) numpy gives exact references:

  k_inf = 0.8841026, the dominant eigenvalue of L^-1 F (use_precursors=True), and
  k_p = 0.8784718, that of L^-1 F_p (use_precursors=False), which depends on how the import splits
  nu-fission into prompt and delayed parts by group;
  the fixed-source flux (L - F)^-1 q, and a transient from that steady state, whose precursors the
  import derives from beta and decay-rate, stays stationary for theta = 1.

Checks:
  OPENMC_BETA_K_REL_ERR: |k / k_inf - 1|.
  OPENMC_BETA_K_PROMPT_REL_ERR: |k_p / k_p,ref - 1|.
  OPENMC_BETA_SOURCE_PHI_REL_ERR: max |phi - phi_ref| / max phi_ref.
  OPENMC_BETA_NULL_TRANSIENT_DRIFT: max over steps and nodes of |phi^n - phi^0| / max phi^0.
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
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.source import VolumetricSource
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        SteadyStateSourceSolver,
        TransientSolver,
    )

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from make_beta_form_mgxs import DATA  # noqa: E402

XS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "beta_form_mgxs.h5")


def operators():
    """L = diag(sigma_t) - S, F and F_p."""
    S = np.zeros((2, 2))
    for gf, gt, v in DATA["transfers"]:
        S[gt, gf] = v
    nu_sigma_f = np.array(DATA["nu"]) * np.array(DATA["fission"])
    prompt_fraction = 1.0 - np.array(DATA["beta"]).sum(axis=0)
    F = np.outer(DATA["chi"], nu_sigma_f)
    F_p = np.outer(DATA["chi"], prompt_fraction * nu_sigma_f)
    return np.diag(DATA["total"]) - S, F, F_p


def dominant(L, F):
    return float(np.max(np.linalg.eigvals(np.linalg.solve(L, F)).real))


def build_problem(use_precursors, source=None):
    xs = MultiGroupXS()
    xs.LoadFromOpenMC(XS_FILE, "fuel", 294.0)
    nodes = [0.25 * i for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    settings = {}
    if source is not None:
        settings["volumetric_sources"] = [VolumetricSource(block_ids=[0], group_strength=source)]
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": use_precursors,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
        **settings,
    )


def eigenvalue(use_precursors):
    problem = build_problem(use_precursors)
    solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-13, max_iters=1000)
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue()


def max_rel_err(phi, phi_ref):
    local = float(np.max(np.abs(phi - phi_ref))) if phi.size else 0.0
    return MPIAllReduce(local, "max") / float(np.max(np.abs(phi_ref)))


if __name__ == "__main__":
    L, F, F_p = operators()
    results = {
        "K_REL_ERR": abs(eigenvalue(True) / dominant(L, F) - 1.0),
        "K_PROMPT_REL_ERR": abs(eigenvalue(False) / dominant(L, F_p) - 1.0),
    }

    q = np.array([1.0, 0.5])
    problem = build_problem(True, source=q.tolist())
    steady = SteadyStateSourceSolver(problem=problem)
    steady.Initialize()
    steady.Execute()
    phi0 = np.array(problem.GetPhiNewLocal()).reshape(-1, 2)
    phi_ref = np.linalg.solve(L - F, q)
    results["SOURCE_PHI_REL_ERR"] = max_rel_err(phi0, np.broadcast_to(phi_ref, phi0.shape))

    problem.SetTimeDependentMode()
    transient = TransientSolver(
        problem=problem, dt=0.1, stop_time=1.0e9, theta=1.0, initial_state="existing",
        verbose=False,
    )
    transient.Initialize()
    drift = 0.0
    for _ in range(4):
        transient.Advance()
        phi = np.array(problem.GetPhiNewLocal()).reshape(-1, 2)
        drift = max(drift, max_rel_err(phi, phi0))
    results["NULL_TRANSIENT_DRIFT"] = drift

    if rank == 0:
        for key, value in results.items():
            print(f"OPENMC_BETA_{key}={value:.6e}")
