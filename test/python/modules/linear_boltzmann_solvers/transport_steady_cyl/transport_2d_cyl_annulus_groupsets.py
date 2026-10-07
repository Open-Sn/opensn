#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RZ annulus inner-boundary inflow, angular ordering, and multiple groupsets.

A. Void annulus, r in [0.5, 1], z in [0, 1], with isotropic incident angular flux X = 2 on every
   boundary: psi = X, so phi = 2 to O(sigma_t L). The inner radial boundary is not the symmetry
   axis and must receive its boundary flux; without a reflecting axis, the angle sets of each
   polar level must still be swept in azimuthal order. Checked with AZIMUTHAL and SINGLE angle
   aggregation.
B. Two-group cylinder with downscattering and leakage: one groupset and one groupset per group
   give the same flux.

Checks:
  RZ_ANNULUS_DEV_AZIMUTHAL, RZ_ANNULUS_DEV_SINGLE: max |phi / 2 - 1|.
  RZ_GROUPSET_DIFF: max relative difference between the groupset structures.
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

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver


def solve(grid, xs, num_groups, bcs, agg, group_ranges, sources=None):
    quad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": r,
                "angular_quadrature": quad,
                "angle_aggregation_type": agg,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
            }
            for r in group_ranges
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=sources or [],
        boundary_conditions=bcs,
        options={"verbose_inner_iterations": False, "max_ags_iterations": 100,
                 "ags_tolerance": 1.0e-12},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return np.array(problem.GetPhiNewLocal())


if __name__ == "__main__":
    if size != 2:
        sys.exit(f"Incorrect number of processors. Expected 2 but got {size}.")

    void = MultiGroupXS()
    void.CreateSimpleOneGroup(1.0e-8, 0.0)
    annulus = OrthogonalMeshGenerator(
        node_sets=[[0.5 + 0.05 * i for i in range(11)], [0.1 * i for i in range(11)]],
        coord_sys="cylindrical").Execute()
    annulus.SetUniformBlockID(0)
    bcs = [{"name": n, "type": "isotropic", "group_strength": [2.0]}
           for n in ("rmin", "rmax", "zmin", "zmax")]
    dev = {}
    for agg in ("azimuthal", "single"):
        phi = solve(annulus, void, 1, bcs, agg, [(0, 0)])
        dev[agg] = MPIAllReduce(float(np.max(np.abs(phi / 2.0 - 1.0))), "max")

    xs_file = "rz_annulus_2g.xs"
    if rank == 0:
        with open(xs_file, "w") as f:
            f.write("NUM_GROUPS 2\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n0 1.0\n1 1.5\nSIGMA_T_END\n\n"
                    "TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.4\nM_GFROM_GTO_VAL 0 0 1 0.3\n"
                    "M_GFROM_GTO_VAL 0 1 1 0.8\nTRANSFER_MOMENTS_END\n")
    MPIBarrier()
    two_group = MultiGroupXS()
    two_group.LoadFromOpenSn(xs_file)
    cylinder = OrthogonalMeshGenerator(node_sets=[[0.1 * i for i in range(11)]] * 2,
                                       coord_sys="cylindrical").Execute()
    cylinder.SetUniformBlockID(0)
    src = [VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])]
    axis = [{"name": "rmin", "type": "reflecting"}]
    phi_a = solve(cylinder, two_group, 2, axis, "azimuthal", [(0, 1)], src)
    phi_b = solve(cylinder, two_group, 2, axis, "azimuthal", [(0, 0), (1, 1)], src)
    gs_diff = MPIAllReduce(float(np.max(np.abs(phi_a - phi_b) / np.abs(phi_a))), "max")

    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)
        print(f"RZ_ANNULUS_DEV_AZIMUTHAL {dev['azimuthal']:.6e}")
        print(f"RZ_ANNULUS_DEV_SINGLE {dev['single']:.6e}")
        print(f"RZ_GROUPSET_DIFF {gs_diff:.6e}")
