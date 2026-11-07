#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TGDSA removal of the alternating Jacobi error mode.

The groups of a groupset are iterated with Jacobi iteration in energy. With
converged within-group scattering its infinite-medium iteration matrix
C = (diag(sigma_t) - diag(S))^-1 offdiag(S) has the Perron eigenvalue rho and
can also have an eigenvalue near -rho, whose error mode changes sign every
iteration. The two-grid diffusion correction represents only the Perron mode.
The local energy correction of TGDSA removes the other modes.

Two infinite media (a 1D slab with reflecting boundaries) have a source in
group 0, so the exact solution is uniform, phi = (diag(sigma_t) - S)^-1 q.
  A. The C5G7 moderator (7 groups, one groupset) has rho = 0.8719 and the
     eigenvalue -0.8715.
  B. Two groups with strong upscattering have the eigenvalues +/-0.99775.
The reference is solved here from the cross sections with numpy.

Each case runs exactly 20 iterations of classic or PETSc Richardson with WGDSA
and TGDSA, so the result measures the acceleration, not the stopping test.
The error is then below 2e-7 in every case. Without the local correction it is
6e-4 to 8e-3: case B stalls (the -0.99775 mode decays by 0.998 per iteration)
and case A converges at about 0.87 per iteration.

Check:
  TGDSA_ALTERNATING_MAX_REL_ERR: max over cases, inner solvers and groups of
  |phi / phi_exact - 1| (min and max of phi over the domain).
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def read_sigma_t_and_s0(path):
    """Total cross sections and isotropic transfer matrix S[g_to, g_from] of an OpenSn xs file."""
    num_groups, sigma_t, S, section = 0, None, None, None
    with open(path) as f:
        for line in f:
            words = line.split()
            if not words:
                continue
            if words[0] == "NUM_GROUPS":
                num_groups = int(words[1])
                sigma_t, S = np.zeros(num_groups), np.zeros((num_groups, num_groups))
            elif words[0] in ("SIGMA_T_BEGIN", "TRANSFER_MOMENTS_BEGIN"):
                section = words[0]
            elif words[0].endswith("_END"):
                section = None
            elif section == "SIGMA_T_BEGIN":
                sigma_t[int(words[0])] = float(words[1])
            elif section == "TRANSFER_MOMENTS_BEGIN" and int(words[1]) == 0:
                S[int(words[3]), int(words[2])] = float(words[4])
    return sigma_t, S


def max_relative_error(problem, phi_exact):
    error = 0.0
    for g, exact in enumerate(phi_exact):
        for operation in ("max", "min"):
            interp = FieldFunctionInterpolationVolume()
            interp.SetOperationType(operation)
            interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
            interp.AddFieldFunction(problem.GetScalarFluxFieldFunction()[g])
            interp.Execute()
            value = interp.GetValue()
            error = max(error, abs(value / exact - 1.0) if value == value else 1.0e300)
    return error


if __name__ == "__main__":
    # B. Two groups with upscattering, written to a file in the working directory
    two_group_file = "tgdsa_alternating_mode_two_group.xs"
    if rank == 0:
        sigma_t2 = [3.0, 6.0]
        transfer = [[1.8, 3.0], [1.197, 2.994]]  # transfer[g_to][g_from]
        with open(two_group_file, "w") as f:
            f.write("NUM_GROUPS 2\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g in range(2):
                f.write(f"{g} {sigma_t2[g]}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(2):
                for g_from in range(2):
                    f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {transfer[g_to][g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()

    cases = (("c5g7_water", "../../../../assets/xs/c5g7/XS_water.xs"),
             ("two_group", two_group_file))

    nodes = [i * 10.0 / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    max_err = 0.0
    for name, path in cases:
        sigma_t, S = read_sigma_t_and_s0(path)
        num_groups = len(sigma_t)
        q = np.zeros(num_groups)
        q[0] = 1.0
        phi_exact = np.linalg.solve(np.diag(sigma_t) - S, q)

        xs = MultiGroupXS()
        xs.LoadFromOpenSn(path)
        for method in ("classic_richardson", "petsc_richardson"):
            problem = DiscreteOrdinatesProblem(
                mesh=grid,
                num_groups=num_groups,
                groupsets=[
                    {
                        "groups_from_to": (0, num_groups - 1),
                        "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8,
                                                                        scattering_order=0),
                        "inner_linear_method": method,
                        "l_abs_tol": 1.0e-14,
                        "l_max_its": 20,
                        "apply_wgdsa": True,
                        "wgdsa_l_abs_tol": 1.0e-12,
                        "apply_tgdsa": True,
                        "tgdsa_l_abs_tol": 1.0e-12,
                    }
                ],
                xs_map=[{"block_ids": [0], "xs": xs}],
                volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=list(q))],
                boundary_conditions=[{"name": "zmin", "type": "reflecting"},
                                     {"name": "zmax", "type": "reflecting"}],
                options={"verbose_inner_iterations": False},
            )
            solver = SteadyStateSourceSolver(problem=problem)
            solver.Initialize()
            solver.Execute()
            err = max_relative_error(problem, phi_exact)
            max_err = max(max_err, err)
            if rank == 0:
                print(f"CASE {name} {method} max_rel_err {err:.6e}")

    MPIBarrier()
    if rank == 0:
        os.remove(two_group_file)
        print(f"TGDSA_ALTERNATING_MAX_REL_ERR {max_err:.6e}")
