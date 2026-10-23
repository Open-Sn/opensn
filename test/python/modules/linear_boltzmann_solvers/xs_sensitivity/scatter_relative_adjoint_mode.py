#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Relative scattering sensitivities computed with the problem in adjoint mode.

Two groups, sigma_t = (1, 1), within-group scattering 0.5 in each group, downscattering 0.3 from
group 0 to group 1 and no upscattering, on a 1D slab of length 1 with both faces reflecting (an
infinite medium). The response R = int phi_1 of the source q = (1, 0) has the forward flux
phi = (2, 1.2) and the adjoint flux phi* = (1.2, 2) (source (0, 1)). The sensitivity of R to
sigma_s(g -> g') is int phi*_g' phi_g: 4 for 0 -> 1 and 1.44 for 1 -> 0. The relative sensitivity
is sigma_s(g -> g') times that: 1.2 for 0 -> 1 and 0 for 1 -> 0. Its coefficient is the forward
cross section whatever the mode; read from the transposed matrix it would give 0 and 0.432.

The sensitivities are computed after the adjoint solve, with the problem in adjoint mode, and again
after the problem is switched back to forward mode.

Checks (MODE = ADJ_MODE or FWD_MODE):
  XS_SENS_SCATTER_<MODE>_ABS_ERR: max over 0 -> 1 and 1 -> 0 of |S / S_exact - 1|.
  XS_SENS_SCATTER_<MODE>_RELATIVE_ERR: max of |S_rel(0 -> 1) / 1.2 - 1| and |S_rel(1 -> 0)|.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.post import CrossSectionSensitivityPostprocessor
else:
    barrier = MPIBarrier

XS_TEXT = """NUM_GROUPS 2
NUM_MOMENTS 1
SIGMA_T_BEGIN
0 1.0
1 1.0
SIGMA_T_END
SIGMA_A_BEGIN
0 0.2
1 0.5
SIGMA_A_END
TRANSFER_MOMENTS_BEGIN
M_GFROM_GTO_VAL 0 0 0 0.5
M_GFROM_GTO_VAL 0 0 1 0.3
M_GFROM_GTO_VAL 0 1 1 0.5
TRANSFER_MOMENTS_END
"""


def remove_rank_file(prefix):
    try:
        os.remove(prefix + str(rank) + ".h5")
    except FileNotFoundError:
        pass


def scatter_sensitivity(problem, from_group, to_group, relative):
    pp = CrossSectionSensitivityPostprocessor(
        problem=problem,
        sensitivity_type="scatter",
        moment=0,
        from_group=from_group,
        to_group=to_group,
        relative=relative,
        forward_flux_moments=fwd_phi_prefix,
        adjoint_flux_moments=adj_phi_prefix,
    )
    pp.Execute()
    return pp.GetValue()[0][0]


def errors(problem):
    abs_err = max(
        abs(scatter_sensitivity(problem, 0, 1, False) / 4.0 - 1.0),
        abs(scatter_sensitivity(problem, 1, 0, False) / 1.44 - 1.0),
    )
    relative_err = max(
        abs(scatter_sensitivity(problem, 0, 1, True) / 1.2 - 1.0),
        abs(scatter_sensitivity(problem, 1, 0, True)),
    )
    return abs_err, relative_err


if __name__ == "__main__":
    nodes = [0.1 * i for i in range(11)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs_path = f"xs_sens_scatter_adjoint_mode_r{rank}.xs"
    with open(xs_path, "w") as f:
        f.write(XS_TEXT)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_path)
    os.remove(xs_path)

    boundaries = [{"name": "zmin", "type": "reflecting"}, {"name": "zmax", "type": "reflecting"}]
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": [0, 1],
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0])],
        boundary_conditions=boundaries,
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    solver = SteadyStateSourceSolver(problem=phys)
    solver.Initialize()
    solver.Execute()
    fwd_phi_prefix = "xs_sens_scatter_adjoint_mode_fwd_phi_"
    phys.WriteFluxMoments(fwd_phi_prefix)

    phys.SetAdjoint(True)
    phys.SetBoundaryOptions(boundary_conditions=boundaries)
    phys.SetVolumetricSources(
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[0.0, 1.0])]
    )
    solver.Execute()
    adj_phi_prefix = "xs_sens_scatter_adjoint_mode_adj_phi_"
    phys.WriteFluxMoments(adj_phi_prefix)

    adj_mode_abs_err, adj_mode_relative_err = errors(phys)
    phys.SetAdjoint(False)
    fwd_mode_abs_err, fwd_mode_relative_err = errors(phys)

    if rank == 0:
        print(f"XS_SENS_SCATTER_ADJ_MODE_ABS_ERR={adj_mode_abs_err:.12e}")
        print(f"XS_SENS_SCATTER_ADJ_MODE_RELATIVE_ERR={adj_mode_relative_err:.12e}")
        print(f"XS_SENS_SCATTER_FWD_MODE_ABS_ERR={fwd_mode_abs_err:.12e}")
        print(f"XS_SENS_SCATTER_FWD_MODE_RELATIVE_ERR={fwd_mode_relative_err:.12e}")

    barrier()
    for prefix in [fwd_phi_prefix, adj_phi_prefix]:
        remove_rank_file(prefix)
