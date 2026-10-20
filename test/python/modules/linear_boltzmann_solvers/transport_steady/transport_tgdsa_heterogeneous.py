#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TGDSA without WGDSA converges in a heterogeneous multigroup problem.

A 1D slab of 20 cm (40 cells, vacuum on both ends) holds a 12-group HDPE
moderator with a 2 cm He-3 layer at 14-16 cm, and a fission-spectrum source in
the HDPE. The data are the reduced HDPE and He-3 libraries of the
numerical-methods tutorials.

The groupset (all 12 groups) uses classic Richardson with TGDSA alone. The
TGDSA spectrum is that of Jacobi iteration with converged within-group
scattering, which is thermal in both materials. The spectrum of the single-sweep
iteration (all scattering lagged) is instead concentrated in the most strongly
self-scattering group of each material: the thermal groups in HDPE but groups
0-1 in He-3. With that spectrum, the collapsed correction is mapped to different
groups on the two sides of the interface and the iteration diverges (the
pointwise flux change stays above 1). The solve here needs about 125
iterations; 300 are allowed.

The reference is an independent solve of the same problem with unaccelerated
GMRES converged to 1e-12.

Check:
  TGDSA_HETEROGENEOUS_MAX_REL_DIFF: max |phi - phi_ref| / max |phi_ref| over all
  groups and nodes.
"""

import os
import sys
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


SOURCE_SPECTRUM = [
    8.70673882298e-1, 1.27603333571e-1, 1.69874391340e-3, 2.39751217396e-5,
    6.37406413603e-8, 8.61581497551e-10, 1.06916189569e-10, 2.94232444417e-10,
    5.81630205686e-11, 2.77856961803e-11, 4.45769378270e-12, 1.77083905175e-12,
]


def load_xs(name):
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "../../../../assets/xs", name))
    return xs


def solve(grid, hdpe, he3, method, tgdsa, max_its, tol):
    groupset = {
        "groups_from_to": (0, 11),
        "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
        "inner_linear_method": method,
        "l_abs_tol": tol,
        "l_max_its": max_its,
        "gmres_restart_interval": 200,
        "apply_tgdsa": tgdsa,
        "tgdsa_l_abs_tol": 1.0e-12,
    }
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=12,
        groupsets=[groupset],
        xs_map=[{"block_ids": [0], "xs": hdpe}, {"block_ids": [1], "xs": he3}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=SOURCE_SPECTRUM)],
        boundary_conditions=[
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={"verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return np.array(problem.GetPhiNewLocal())


if __name__ == "__main__":
    nodes = [i * 20.0 / 40 for i in range(41)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=14.0, zmax=16.0), 1, True
    )

    hdpe = load_xs("xs12g_hdpe_reduced.cxs")
    he3 = load_xs("xs12g_he3_reduced.cxs")

    phi_ref = solve(grid, hdpe, he3, "petsc_gmres", False, 2000, 1.0e-12)
    phi = solve(grid, hdpe, he3, "classic_richardson", True, 300, 1.0e-10)

    local_diff = float(np.max(np.abs(phi - phi_ref))) if phi.size else 0.0
    local_scale = float(np.max(np.abs(phi_ref))) if phi_ref.size else 0.0
    diff = MPIAllReduce(local_diff if np.isfinite(local_diff) else 1.0e300, "max")
    scale = MPIAllReduce(local_scale, "max")
    if rank == 0:
        print(f"TGDSA_HETEROGENEOUS_MAX_REL_DIFF {diff / scale:.6e}")
