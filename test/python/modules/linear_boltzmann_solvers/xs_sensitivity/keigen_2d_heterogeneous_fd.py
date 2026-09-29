"""2D heterogeneous k-eigenvalue sensitivities checked against central finite differences.

A one-group fissile corner region, a P1-scattering reflector, and an absorber strip with vacuum
boundaries give a strongly anisotropic flux. For each selected coefficient x, the adjoint-based
sensitivity dk/dx from CrossSectionSensitivityPostprocessor (with k-eigenvalue scaling) is
compared with the central difference (k(x + h) - k(x - h)) / (2 h) of the forward eigenvalue.

The adjoint is solved by sweeping the forward directions with transposed cross sections and
must be reoriented to the physical adjoint by both k-eigenvalue solvers. The sigma_t sensitivity
pairs forward and adjoint angular fluxes, and the P1 scatter sensitivity pairs odd flux moments,
so both require that reorientation. The central-difference truncation error is O(h^2) and the
eigenvalue tolerance is 1e-13, so the relative agreement is limited to about 1e-7.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import CrossSectionSensitivityPostprocessor
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
    )
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


PREFIX = "xs_sens_keigen_2d_"
BASE = {"nu_sigma_f": 0.35, "sigma_s0_refl": 0.9, "sigma_s1_refl": 0.3, "sigma_t_abs": 2.0}


def write_xs(path, sigma_t, sigma_s0, sigma_s1, sigma_f=0.0, nu_sigma_f=0.0):
    with open(path, "w") as f:
        f.write(f"NUM_GROUPS 1\nNUM_MOMENTS 2\n\nSIGMA_T_BEGIN\n0 {sigma_t!r}\nSIGMA_T_END\n\n")
        if sigma_f > 0.0:
            f.write(f"SIGMA_F_BEGIN\n0 {sigma_f!r}\nSIGMA_F_END\n\nCHI_BEGIN\n0 1.0\nCHI_END\n\n")
            f.write(
                f"PRODUCTION_MATRIX_BEGIN\nGPRIME_G_VAL 0 0 {nu_sigma_f!r}\n"
                "PRODUCTION_MATRIX_END\n\n"
            )
        f.write(
            f"TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 {sigma_s0!r}\n"
            f"M_GFROM_GTO_VAL 1 0 0 {sigma_s1!r}\nTRANSFER_MOMENTS_END\n"
        )


def load_xs(tag, params):
    files = [f"{PREFIX}{tag}_{name}.xs" for name in ("refl", "core", "abs")]
    if rank == 0:
        write_xs(files[0], 1.0, params["sigma_s0_refl"], params["sigma_s1_refl"])
        write_xs(files[1], 1.0, 0.6, 0.1, sigma_f=0.1, nu_sigma_f=params["nu_sigma_f"])
        write_xs(files[2], params["sigma_t_abs"], 0.5, 0.0)
    barrier()
    xs = []
    for path in files:
        x = MultiGroupXS()
        x.LoadFromOpenSn(path)
        xs.append(x)
    barrier()
    if rank == 0:
        for path in files:
            os.remove(path)
    return xs


def remove_rank_file(prefix):
    try:
        os.remove(prefix + str(rank) + ".h5")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    nodes = [i * 0.5 for i in range(17)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=4.0, ymin=0.0, ymax=4.0, infz=True), 1, True
    )
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=5.0, xmax=6.0, ymin=0.0, ymax=8.0, infz=True), 2, True
    )
    quad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=1)
    boundaries = [{"name": n, "type": "vacuum"} for n in ("xmin", "xmax", "ymin", "ymax")]

    def make_problem(tag, params):
        xs = load_xs(tag, params)
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[
                {
                    "groups_from_to": (0, 0),
                    "angular_quadrature": quad,
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-13,
                    "l_max_its": 500,
                }
            ],
            xs_map=[{"block_ids": [i], "xs": xs[i]} for i in range(3)],
            boundary_conditions=boundaries,
            options={
                "save_angular_flux": True,
                "verbose_inner_iterations": False,
                "verbose_outer_iterations": False,
            },
        )

    def power_iteration(problem):
        solver = PowerIterationKEigenSolver(problem=problem, max_iters=2000, k_tol=1.0e-13)
        solver.Initialize()
        solver.Execute()
        return solver.GetEigenvalue()

    def perturbed_k(key, delta):
        params = dict(BASE)
        params[key] += delta
        return power_iteration(make_problem(f"{key}_{'p' if delta > 0 else 'm'}", params))

    h = 1.0e-4

    def central_difference(key):
        return (perturbed_k(key, h) - perturbed_k(key, -h)) / (2.0 * h)

    fd = {key: central_difference(key) for key in ("sigma_t_abs", "nu_sigma_f", "sigma_s1_refl")}

    # Forward and adjoint solves at the base state.
    problem = make_problem("base", BASE)
    k_eff = power_iteration(problem)
    fwd_psi, fwd_phi = PREFIX + "fwd_psi_", PREFIX + "fwd_phi_"
    problem.WriteAngularFluxes(fwd_psi)
    problem.WriteFluxMoments(fwd_phi)

    problem.SetAdjoint(True)
    power_iteration(problem)
    adj_psi, adj_phi = PREFIX + "adj_psi_", PREFIX + "adj_phi_"
    problem.WriteAngularFluxes(adj_psi)
    problem.WriteFluxMoments(adj_phi)

    nl_solver = NonLinearKEigenSolver(
        problem=problem,
        nl_abs_tol=1.0e-13,
        nl_rel_tol=1.0e-13,
        nl_max_its=100,
        l_abs_tol=1.0e-14,
        l_max_its=100,
    )
    nl_solver.Initialize()
    nl_solver.Execute()
    nl_adj_psi, nl_adj_phi = PREFIX + "nl_adj_psi_", PREFIX + "nl_adj_phi_"
    problem.WriteAngularFluxes(nl_adj_psi)
    problem.WriteFluxMoments(nl_adj_phi)

    def sensitivity(adjoint_psi, adjoint_phi, **options):
        pp = CrossSectionSensitivityPostprocessor(
            problem=problem,
            forward_angular_fluxes=fwd_psi,
            adjoint_angular_fluxes=adjoint_psi,
            forward_flux_moments=fwd_phi,
            adjoint_flux_moments=adjoint_phi,
            **options,
        )
        pp.Execute()
        pp.ApplyKEigenvalueScaling(k_eff)
        return pp.GetValue()[0][0]

    sigma_t_options = {"sensitivity_type": "sigma_t", "block_ids": [2], "group": 0}
    results = {
        "SIGMA_T": (sensitivity(adj_psi, adj_phi, **sigma_t_options), fd["sigma_t_abs"]),
        "NLKE_SIGMA_T": (
            sensitivity(nl_adj_psi, nl_adj_phi, **sigma_t_options),
            fd["sigma_t_abs"],
        ),
        "PRODUCTION": (
            sensitivity(
                adj_psi, adj_phi, sensitivity_type="production", block_ids=[1], group=0
            ),
            fd["nu_sigma_f"],
        ),
        "SCATTER_P1": (
            sensitivity(
                adj_psi,
                adj_phi,
                sensitivity_type="scatter",
                block_ids=[0],
                moment=1,
                from_group=0,
                to_group=0,
            ),
            fd["sigma_s1_refl"],
        ),
    }

    if rank == 0:
        for name, (adjoint_value, fd_value) in results.items():
            rel_err = abs(adjoint_value - fd_value) / abs(fd_value)
            print(f"XS_SENS_KEIGEN_2D_{name}_REL_ERR={rel_err:.12e}")

    barrier()
    for prefix in (fwd_psi, fwd_phi, adj_psi, adj_phi, nl_adj_psi, nl_adj_phi):
        remove_rank_file(prefix)
