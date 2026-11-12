#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D test of an infinite medium of graphite (168 groups) with WGDSA in both groupsets and TGDSA in
the thermal groupset. SDM: PWLD. With energy_coupling=True, both groupsets couple the WGDSA
corrections in energy (wgdsa_energy_coupling).

One slab cell with reflecting faces and a source in group 0 is an infinite medium, so the scalar
flux is the solution of the multigroup balance phi = (diag(sigma_t) - S_0)^-1 q, with S_0 the
isotropic transfer matrix, which is evaluated here from the cross-section file independently of
the transport solve. The within-group tolerance (1e-10) is tighter than the across-groupset
tolerance, which is absolute (1e-6 on a flux norm of about 600), so the across-groupset iteration
can converge to it.

Checks:
  AGS_CONVERGED: 1 when the steady-state solve reports convergence.
  PHI_MAX_REL_ERR: max over nodes and groups of |phi / phi_exact - 1|.
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

if __name__ == "__main__":

    num_procs = 1

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    nodes = []
    N = 1
    L = 1e6
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    num_groups = 168
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("../../../../assets/xs/xs_graphite_pure.xs")
    xs_air = MultiGroupXS()
    xs_air.LoadFromOpenSn("../../../../assets/xs/xs_air50RH.xs")

    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 1.0
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    # Setup Physics
    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=1)

    # Optional WGDSA energy coupling (default False)
    dsa_options = {}
    if "energy_coupling" in globals():
        dsa_options["wgdsa_energy_coupling"] = energy_coupling

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 62],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 1000,
                "gmres_restart_interval": 30,
                "apply_wgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-2,
                **dsa_options,
            },
            {
                "groups_from_to": [63, num_groups - 1],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 1000,
                "gmres_restart_interval": 30,
                "apply_wgdsa": True,
                "apply_tgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-2,
                **dsa_options,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_graphite},
            {"block_ids": [1], "xs": xs_air},
        ],
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
    )

    # Initialize and Execute Solver
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Infinite-medium reference from the cross-section file: sigma_t and the P0 transfer matrix.
    sigma_t = np.zeros(num_groups)
    transfer = np.zeros((num_groups, num_groups))
    section = None
    with open("../../../../assets/xs/xs_graphite_pure.xs") as f:
        for line in f:
            words = line.split()
            if not words or words[0].startswith("#"):
                continue
            if words[0] in ("SIGMA_T_BEGIN", "TRANSFER_MOMENTS_BEGIN"):
                section = words[0]
            elif words[0].endswith("_END"):
                section = None
            elif section == "SIGMA_T_BEGIN":
                sigma_t[int(words[0])] = float(words[1])
            elif section == "TRANSFER_MOMENTS_BEGIN" and words[1] == "0":
                # M_GFROM_GTO_VAL moment g_from g_to value
                transfer[int(words[3]), int(words[2])] += float(words[4])
    q = np.array(strength)
    phi_exact = np.linalg.solve(np.diag(sigma_t) - transfer, q)

    # Local unknowns are ordered node, moment, group; scattering order 1 in 1D gives 2 moments.
    phi = np.array(phys.GetPhiNewLocal())
    num_moments = 2
    phi0 = phi.reshape(-1, num_moments, num_groups)[:, 0, :]
    max_rel_err = float(np.max(np.abs(phi0 / phi_exact - 1.0)))
    if rank == 0:
        print(f"AGS_CONVERGED={int(ss_solver.IsConverged())}")
        print(f"PHI_MAX_REL_ERR={max_rel_err:.6e}")
