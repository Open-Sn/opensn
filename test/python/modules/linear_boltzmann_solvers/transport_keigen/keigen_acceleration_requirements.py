#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
k-eigenvalue acceleration requirements are checked collectively by the owning solver's
Initialize(), which raises a ValueError on every rank, and valid neighboring configurations are
accepted.

Each pair below is a configuration that violates a requirement and, where useful, the closest
configuration that satisfies it:

1. SCDSA with two groupsets (rejected).
2. SMM with P1 scattering and without stored angular fluxes (rejected; both errors are reported).
3. SMM with a zero-strength isotropic boundary (rejected: the SMM boundary closure covers only
   vacuum and reflecting boundaries), and
4. power iteration without acceleration on the same problem (accepted: a zero-strength isotropic
   boundary is not an external source).
5. SCDSA with an AAH sweep and l_max_its=1 (rejected: unstable with lagged angular fluxes), and
6. the same with l_max_its=2 (accepted).
7. An acceleration object created for a different problem than the solver's (rejected).
8. NonLinearKEigenSolver with WGDSA on two groupsets (rejected), and
9. the same groups in one groupset (accepted).
10. CMFD in 1D spherical geometry (rejected: the coarse balance mixes spherical face areas with
    slab-length volumes), and
11. CMFD in RZ geometry (accepted).
12. SCDSA with sdm="pwlc" in 1D spherical geometry (rejected: the correction diverges), and
13. SCDSA with sdm="pwld" in 1D spherical geometry (accepted).
14. Mixed-material aggregation requested with an identity CMFD mesh (rejected).

Checks:
  KEIGEN_ACCELERATION_REQUIREMENTS_PASSED: number of cases (of 14) that behave as described.
"""

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
    from pyopensn.aquad import (
        GLCProductQuadrature2DRZ,
        GLCProductQuadrature2DXY,
        GLProductQuadrature1DSpherical,
    )
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        CMFDAcceleration,
        DiscreteOrdinatesCurvilinearProblem,
        DiscreteOrdinatesProblem,
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
        SCDSAAcceleration,
        SMMAcceleration,
    )
    from pyopensn.xs import MultiGroupXS


def groupset(groups, quadrature, **extra):
    block = {
        "groups_from_to": groups,
        "angular_quadrature": quadrature,
        "inner_linear_method": "petsc_richardson",
        "l_max_its": 2,
        "l_abs_tol": 1.0e-10,
    }
    block.update(extra)
    return block


def xy_problem(grid, xs, num_groups, groupsets, boundary_conditions=None, options=None):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=groupsets,
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=boundary_conditions or [],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False,
                 **(options or {})},
    )


def curvilinear_problem(grid, xs, quadrature):
    return DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[groupset((0, 0), quadrature, angle_aggregation_type="azimuthal")],
        xs_map=[{"block_ids": [0], "xs": xs}],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )


def outcome(solver, *expected):
    """Returns 'rejected' for a ValueError naming every expected text, 'accepted' if no error."""
    try:
        solver.Initialize()
    except ValueError as error:
        return "rejected" if all(text in str(error) for text in expected) else "other"
    return "accepted"


def on_all_ranks(condition):
    return int(MPIAllReduce(int(condition), "min"))


def pi(problem, acceleration=None):
    kwargs = {"problem": problem, "k_tol": 1.0e-8, "max_iters": 50}
    if acceleration is not None:
        kwargs["acceleration"] = acceleration
    return PowerIterationKEigenSolver(**kwargs)


if __name__ == "__main__":
    xs_dir = os.path.join(os.path.dirname(__file__), "../../../../assets/xs")
    fissile = MultiGroupXS()
    fissile.LoadFromOpenSn(os.path.join(xs_dir, "simple_fissile_1g.xs"))

    xs2_file = "keigen_acceleration_requirements_2g.xs"
    if rank == 0:
        with open(xs2_file, "w") as f:
            f.write("NUM_GROUPS 2\nNUM_MOMENTS 1\n\n"
                    "SIGMA_T_BEGIN\n0 1.0\n1 1.0\nSIGMA_T_END\n\n"
                    "SIGMA_F_BEGIN\n0 0.05\n1 0.05\nSIGMA_F_END\n\n"
                    "NU_PROMPT_BEGIN\n0 2.4\n1 2.4\nNU_PROMPT_END\n\n"
                    "CHI_PROMPT_BEGIN\n0 1.0\n1 0.0\nCHI_PROMPT_END\n\n"
                    "TRANSFER_MOMENTS_BEGIN\n"
                    "M_GFROM_GTO_VAL 0 0 0 0.5\nM_GFROM_GTO_VAL 0 0 1 0.4\n"
                    "M_GFROM_GTO_VAL 0 1 1 0.9\n"
                    "TRANSFER_MOMENTS_END\n")
    MPIAllReduce(0)
    fissile2 = MultiGroupXS()
    fissile2.LoadFromOpenSn(xs2_file)

    nodes = [i * 10.0 / 8 for i in range(9)]
    xy_grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    xy_grid.SetUniformBlockID(0)
    rz_grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes], coord_sys="cylindrical").Execute()
    rz_grid.SetUniformBlockID(0)
    sph_grid = OrthogonalMeshGenerator(node_sets=[[i * 10.0 / 16 for i in range(17)]],
                                       coord_sys="spherical").Execute()
    sph_grid.SetUniformBlockID(0)

    def xy_quad(scattering_order=0):
        return GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4,
                                        scattering_order=scattering_order)

    def rz_quad():
        return GLCProductQuadrature2DRZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    def sph_quad():
        return GLProductQuadrature1DSpherical(n_polar=8, scattering_order=0)

    save_psi = {"save_angular_flux": True}
    passed = 0

    # Case 1
    problem = xy_problem(xy_grid, fissile2, 2,
                         [groupset((0, 0), xy_quad()), groupset((1, 1), xy_quad())])
    passed += on_all_ranks(
        outcome(pi(problem, SCDSAAcceleration(problem=problem)), "single groupset") == "rejected"
    )

    # Case 2
    problem = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad(scattering_order=1))])
    passed += on_all_ranks(
        outcome(pi(problem, SMMAcceleration(problem=problem)),
                "invalid configuration", "P0 scattering", "save_angular_flux") == "rejected"
    )

    # Cases 3 and 4
    isotropic_zero = [{"name": "xmin", "type": "isotropic", "group_strength": [0.0]}]
    problem = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad())], isotropic_zero,
                         save_psi)
    passed += on_all_ranks(
        outcome(pi(problem, SMMAcceleration(problem=problem)),
                "only vacuum and reflecting boundaries") == "rejected"
    )
    problem = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad())], isotropic_zero)
    passed += on_all_ranks(outcome(pi(problem)) == "accepted")

    # Cases 5 and 6
    problem = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad(), l_max_its=1)])
    passed += on_all_ranks(
        outcome(pi(problem, SCDSAAcceleration(problem=problem)),
                "one sweep per iteration") == "rejected"
    )
    problem = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad())])
    passed += on_all_ranks(
        outcome(pi(problem, SCDSAAcceleration(problem=problem))) == "accepted"
    )

    # Case 7
    other = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad())])
    problem = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad())])
    passed += on_all_ranks(
        outcome(pi(problem, SCDSAAcceleration(problem=other)), "different problem") == "rejected"
    )

    # Cases 8 and 9
    def nlke(problem):
        return NonLinearKEigenSolver(problem=problem, nl_abs_tol=1.0e-8, nl_rel_tol=1.0e-8)

    problem = xy_problem(xy_grid, fissile2, 2,
                         [groupset((0, 0), xy_quad(), apply_wgdsa=True),
                          groupset((1, 1), xy_quad(), apply_wgdsa=True)])
    passed += on_all_ranks(outcome(nlke(problem), "single-groupset problems") == "rejected")
    problem = xy_problem(xy_grid, fissile2, 2, [groupset((0, 1), xy_quad(), apply_wgdsa=True)])
    passed += on_all_ranks(outcome(nlke(problem)) == "accepted")

    # Cases 10 and 11
    problem = curvilinear_problem(sph_grid, fissile, sph_quad())
    passed += on_all_ranks(
        outcome(pi(problem, CMFDAcceleration(problem=problem)),
                "not supported in 1D spherical geometry") == "rejected"
    )
    problem = curvilinear_problem(rz_grid, fissile, rz_quad())
    passed += on_all_ranks(outcome(pi(problem, CMFDAcceleration(problem=problem))) == "accepted")

    # Cases 12 and 13
    problem = curvilinear_problem(sph_grid, fissile, sph_quad())
    passed += on_all_ranks(
        outcome(pi(problem, SCDSAAcceleration(problem=problem, sdm="pwlc")),
                'sdm="pwlc" is not supported in 1D spherical geometry') == "rejected"
    )
    problem = curvilinear_problem(sph_grid, fissile, sph_quad())
    passed += on_all_ranks(
        outcome(pi(problem, SCDSAAcceleration(problem=problem, sdm="pwld"))) == "accepted"
    )

    # Case 14
    problem = xy_problem(xy_grid, fissile, 1, [groupset((0, 0), xy_quad())])
    passed += on_all_ranks(
        outcome(
            pi(
                problem,
                CMFDAcceleration(
                    problem=problem,
                    coarse_mesh="identity",
                    allow_mixed_material_aggregation=True,
                ),
            ),
            "accepts allow_mixed_material_aggregation=true only with local_aggregation or "
            "global_aggregation",
        )
        == "rejected"
    )

    MPIAllReduce(0)
    if rank == 0:
        os.remove(xs2_file)
        print(f"KEIGEN_ACCELERATION_REQUIREMENTS_PASSED {passed}")
