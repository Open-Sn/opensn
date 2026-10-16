#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Source moments made from a flux-moment file.

ReadFluxMomentsAndMakeSourceMoments(file_base, single_file) reads flux moments phi_f from file
and sets the external source moments to the scattering and fission source they produce,
(S + F) phi_f. It leaves the problem's own flux vectors (phi_new and phi_old) unchanged.

Infinite medium (1D slab, both faces reflecting), one group, isotropic scattering, sigma_t = 1,
scattering ratio c = 0.25. With a uniform isotropic source q the flux is flat and exact for any
quadrature and mesh: phi = q / (sigma_t (1 - c)).

1. A problem with the volumetric source q = 1 has phi_1 = 4/3, written to file.
2. A second problem with no volumetric source and use_source_moments=True reads phi_1 from file.
   Its source moments are S phi_1 = c sigma_t phi_1 = 1/3, so its flux is phi_2 = 4/9. A source
   made from the second problem's own flux (zero before its solve) gives phi_2 = 0; using the
   file flux itself as the source gives phi_2 = 16/9.
3. A solved problem: a third problem with use_source_moments=True takes phi_1 itself as its source
   moments (ReadSourceMoments), so its flux is phi_3 = 16/9. It then reads phi_1 with
   ReadFluxMomentsAndMakeSourceMoments and is solved again, giving phi_4 = 4/9 as in 2. A source
   made from the problem's own solution phi_3 gives phi_4 = 16/27.
4. ReadFluxMoments loads phi_1 into one of the third problem's flux vectors, so that one holds
   phi_4 = 4/9 and the other phi_1 = 4/3, and the problem reads phi_1 once more.
5. Two groups, one groupset per group, so that the source made from the file has scattering and
   fission terms across groupsets: sigma_t = (1, 2), transfers 0->0 1/2, 0->1 1/4, 1->0 1/8,
   1->1 1, sigma_f = (1/16, 1/4), nu = 2, chi = (3/4, 1/4), and A = diag(sigma_t) - S - F with
   F = chi (nu sigma_f)^T. A problem with the source q = (1, 0) has phi_1 = A^-1 q =
   (224, 72)/55, written to file. A new problem made from that file as in 2 has the source
   moments (S + F) phi_1, so phi_2 = A^-1 (S + F) phi_1 = (56288, 27144)/3025.

Checks:
  FLUX_MOMENTS_PHI1_MAX_REL_ERR: max |phi_1 / (4/3) - 1| over all nodes.
  FLUX_MOMENTS_PHI2_MAX_REL_ERR: max |phi_2 / (4/9) - 1| over all nodes.
  FLUX_MOMENTS_PHI3_MAX_REL_ERR: max |phi_3 / (16/9) - 1| over all nodes.
  FLUX_MOMENTS_PHI4_MAX_REL_ERR: max |phi_4 / (4/9) - 1| over all nodes.
  FLUX_MOMENTS_PHI_NEW_UNCHANGED: the number of the three ReadFluxMomentsAndMakeSourceMoments calls
    (in 2, 3 and 4) after which phi_new is bit for bit the vector before the call, on every rank.
  FLUX_MOMENTS_PHI_OLD_UNCHANGED: the same for phi_old.
  FLUX_MOMENTS_2G_PHI1_MAX_REL_ERR, FLUX_MOMENTS_2G_PHI2_MAX_REL_ERR: phi_1 and phi_2 of 5,
    max |phi_g / phi_g,exact - 1| over all nodes and both groups.
"""

import os
import sys
import numpy as np

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
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


def build_problem(grid, xs, sources, options, groupsets=((0, 0),)):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=groupsets[-1][1] + 1,
        groupsets=[
            {
                "groups_from_to": groups,
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
            for groups in groupsets
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=sources,
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options=options,
    )


def max_rel_err(problem, expected):
    expected = np.atleast_1d(expected)
    phi = np.array(problem.GetPhiNewLocal()).reshape(-1, expected.size)
    local = float(np.max(np.abs(phi / expected - 1.0))) if phi.size > 0 else 0.0
    return MPIAllReduce(local, "max")


def unchanged_on_every_rank(before, after):
    local = 1 if before.tobytes() == after.tobytes() else 0
    return MPIAllReduce(local, "min")


def write_xs(file_name, sigma_t, transfer, sigma_f, nu, chi):
    """Write an OpenSn cross-section file with one moment; transfer[g][gp] is gp -> g."""
    lines = [f"NUM_GROUPS {len(sigma_t)}", "NUM_MOMENTS 1"]
    for name, values in (("SIGMA_T", sigma_t), ("SIGMA_F", sigma_f), ("NU", nu), ("CHI", chi)):
        lines += [f"{name}_BEGIN", *(f"{g} {v!r}" for g, v in enumerate(values)), f"{name}_END"]
    lines.append("TRANSFER_MOMENTS_BEGIN")
    for g, row in enumerate(transfer):
        lines += [f"M_GFROM_GTO_VAL 0 {gp} {g} {v!r}" for gp, v in enumerate(row)]
    lines.append("TRANSFER_MOMENTS_END")
    with open(file_name, "w") as f:
        f.write("\n".join(lines) + "\n")


def read_flux_moments_and_make_source_moments(problem, file_base):
    """Make the source moments from the file; return whether phi_new, phi_old are unchanged."""
    phi_new = np.array(problem.GetPhiNewLocal())
    phi_old = np.array(problem.GetPhiOldLocal())
    problem.ReadFluxMomentsAndMakeSourceMoments(file_base, False)
    return (
        unchanged_on_every_rank(phi_new, np.array(problem.GetPhiNewLocal())),
        unchanged_on_every_rank(phi_old, np.array(problem.GetPhiOldLocal())),
    )


if __name__ == "__main__":
    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    sigma_t = 1.0
    c = 0.25
    q = 1.0
    phi1_exact = q / (sigma_t * (1.0 - c))
    phi2_exact = c * phi1_exact / (1.0 - c)
    phi3_exact = phi1_exact / (sigma_t * (1.0 - c))

    nodes = [0.1 * i for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, c)

    options = {"verbose_inner_iterations": False, "verbose_outer_iterations": False}
    source_moments_options = {**options, "use_source_moments": True}

    # Problem 1: uniform source, flat flux phi_1, written to file.
    source = VolumetricSource(block_ids=[0], group_strength=[q])
    phys1 = build_problem(grid, xs, [source], dict(options))
    solver1 = SteadyStateSourceSolver(problem=phys1)
    solver1.Initialize()
    solver1.Execute()
    err1 = max_rel_err(phys1, phi1_exact)

    file_base = "flux_moments_to_source_moments_1d_"
    phys1.WriteFluxMoments(file_base)

    # Problem 2: source moments made from the file flux phi_1.
    phys2 = build_problem(grid, xs, [], dict(source_moments_options))
    solver2 = SteadyStateSourceSolver(problem=phys2)
    solver2.Initialize()
    unchanged = [read_flux_moments_and_make_source_moments(phys2, file_base)]
    solver2.Execute()
    err2 = max_rel_err(phys2, phi2_exact)

    # Problem 3: solved with phi_1 as its source moments, then given the source made from phi_1.
    phys3 = build_problem(grid, xs, [], dict(source_moments_options))
    solver3 = SteadyStateSourceSolver(problem=phys3)
    solver3.Initialize()
    phys3.ReadSourceMoments(file_base, False)
    solver3.Execute()
    err3 = max_rel_err(phys3, phi3_exact)
    unchanged.append(read_flux_moments_and_make_source_moments(phys3, file_base))
    solver3.Execute()
    err4 = max_rel_err(phys3, phi2_exact)

    # Problem 3 with phi_new and phi_old different.
    phys3.ReadFluxMoments(file_base, False)
    unchanged.append(read_flux_moments_and_make_source_moments(phys3, file_base))

    phi_new_unchanged = sum(new for new, _ in unchanged)
    phi_old_unchanged = sum(old for _, old in unchanged)

    # Two groups, one groupset each, with fission: source moments from the file flux phi_1.
    sigma_t_2g = [1.0, 2.0]
    transfer = [[0.5, 0.125], [0.25, 1.0]]
    sigma_f, nu, chi = [0.0625, 0.25], [2.0, 2.0], [0.75, 0.25]
    production = np.outer(chi, np.multiply(nu, sigma_f))
    a_2g = np.diag(sigma_t_2g) - np.array(transfer) - production
    phi1_2g_exact = np.linalg.solve(a_2g, [1.0, 0.0])
    phi2_2g_exact = np.linalg.solve(a_2g, (np.array(transfer) + production) @ phi1_2g_exact)

    xs_file = f"flux_moments_to_source_moments_1d_2g_{rank}.xs"
    write_xs(xs_file, sigma_t_2g, transfer, sigma_f, nu, chi)
    xs_2g = MultiGroupXS()
    xs_2g.LoadFromOpenSn(xs_file)
    os.remove(xs_file)

    groupsets_2g = ((0, 0), (1, 1))
    options_2g = {**options, "max_ags_iterations": 100, "ags_tolerance": 1.0e-10}
    source_2g = VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])
    phys5 = build_problem(grid, xs_2g, [source_2g], options_2g, groupsets_2g)
    solver5 = SteadyStateSourceSolver(problem=phys5)
    solver5.Initialize()
    solver5.Execute()
    err5 = max_rel_err(phys5, phi1_2g_exact)

    file_base_2g = "flux_moments_to_source_moments_1d_2g_"
    phys5.WriteFluxMoments(file_base_2g)

    phys6 = build_problem(grid, xs_2g, [], {**options_2g, "use_source_moments": True}, groupsets_2g)
    solver6 = SteadyStateSourceSolver(problem=phys6)
    solver6.Initialize()
    phys6.ReadFluxMomentsAndMakeSourceMoments(file_base_2g, False)
    solver6.Execute()
    err6 = max_rel_err(phys6, phi2_2g_exact)

    # Each rank wrote and read only its own files.
    for own_file in (f"{file_base}{rank}.h5", f"{file_base_2g}{rank}.h5"):
        if os.path.exists(own_file):
            os.remove(own_file)

    if rank == 0:
        print(f"FLUX_MOMENTS_PHI1_MAX_REL_ERR={err1:.6e}")
        print(f"FLUX_MOMENTS_PHI2_MAX_REL_ERR={err2:.6e}")
        print(f"FLUX_MOMENTS_PHI3_MAX_REL_ERR={err3:.6e}")
        print(f"FLUX_MOMENTS_PHI4_MAX_REL_ERR={err4:.6e}")
        print(f"FLUX_MOMENTS_PHI_NEW_UNCHANGED={phi_new_unchanged}")
        print(f"FLUX_MOMENTS_PHI_OLD_UNCHANGED={phi_old_unchanged}")
        print(f"FLUX_MOMENTS_2G_PHI1_MAX_REL_ERR={err5:.6e}")
        print(f"FLUX_MOMENTS_2G_PHI2_MAX_REL_ERR={err6:.6e}")
