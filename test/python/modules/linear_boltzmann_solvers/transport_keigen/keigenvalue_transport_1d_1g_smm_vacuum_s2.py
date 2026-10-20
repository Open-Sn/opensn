#!/usr/bin/env python3

"""Power iteration with SMM must converge to the S2 eigenvalue of a slab with vacuum faces.

This is the slab of keigenvalue_transport_1d_1g_smm_vacuum.py, solved with the Gauss-Legendre
S2 set (mu = +-1/sqrt(3)) in place of S32. At a vacuum face the low-order SMM equation uses the
boundary factor sum(w |mu|) / sum(w) of the angular quadrature, 1/sqrt(3) = 0.5774 for S2,
where the Marshak condition of a diffusion solver has the coefficient 1/2. A vacuum term with
the coefficient 1/2 in place of the quadrature's factor moves the SMM eigenvalue by 1.1e-3 here;
for S32 the factor is 0.5004, within 4e-4 of 1/2. With mu^2 = 1/3 and |mu| equal to the
boundary factor, the second-moment and boundary corrections of SMM vanish for this set, so its
eigenvalue is that of the low-order diffusion discretization alone.

Problem: 1D slab of length 10 cm, vacuum at both faces, 50 cells, one group with
sigma_t = 1, sigma_s = 0.3, sigma_f = 0.35 and nu = 1.99 (nu sigma_f = 0.6965 /cm), GL S2.
The input writes these cross sections itself. They carry no delayed-neutron data.

Reference: on a slab the S2 equations are exactly a diffusion equation with
D = mu^2 / sigma_t = 1/(3 sigma_t) and, at a vacuum face, an outgoing current mu phi. Its
eigenvalue is k = nu sigma_f / (sigma_a + D B^2) with B tan(B L / 2) = sigma_t / mu, so
B = 0.28189208 /cm and k = 0.95872234. Power iteration without acceleration reproduces it to
5.2e-7 on this mesh.

Checks (relative eigenvalue differences):
  SMMVacuumS2PIError: power iteration against the S2 reference (tolerance 1e-5).
  SMMVacuumS2PWLDError, SMMVacuumS2PWLCError: power iteration with SMM (sdm "pwld", "pwlc")
    against power iteration, with the tolerance 5e-4 of the S32 input. Both SMM eigenvalues are
    those of the diffusion discretizations, 9.1e-6 from power iteration on this mesh.

The test runs on two ranks, where the SMM diffusion solver selects its iterative policy for
this small problem.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import SMMAcceleration
else:
    barrier = MPIBarrier


XS_FILE = "keigenvalue_transport_1d_1g_smm_vacuum_s2.xs"
XS_TEXT = """NUM_GROUPS 1
NUM_MOMENTS 1
NUM_PRECURSORS 0

SIGMA_T_BEGIN
0 1.0
SIGMA_T_END

SIGMA_F_BEGIN
0 0.35
SIGMA_F_END

NU_BEGIN
0 1.99
NU_END

CHI_BEGIN
0 1.0
CHI_END

TRANSFER_MOMENTS_BEGIN
M_GFROM_GTO_VAL 0 0 0 0.3
TRANSFER_MOMENTS_END
"""


def load_xs():
    """Write the one-group cross sections, load them on every rank and remove the file."""
    if rank == 0:
        with open(XS_FILE, "w", encoding="ascii") as f:
            f.write(XS_TEXT)
    barrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(XS_FILE)
    barrier()
    if rank == 0:
        os.remove(XS_FILE)
    return xs


def solve_k(grid, xs, sdm):
    """Power iteration, with SMM when sdm is given."""
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=2, scattering_order=0),
                "inner_linear_method": "classic_richardson",
                "l_max_its": 2 if sdm else 1000,
                "l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        options={
            "use_precursors": False,
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    if sdm:
        smm = SMMAcceleration(problem=problem, pi_k_tol=1.0e-12, pi_max_its=50, sdm=sdm)
        solver = PowerIterationKEigenSolver(problem=problem, acceleration=smm, k_tol=1.0e-12,
                                            max_iters=2000)
    else:
        solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=20000)
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue()


if __name__ == "__main__":

    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    length = 10.0
    num_cells = 50
    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = load_xs()

    k_reference = 0.95872234
    k_pi = solve_k(grid, xs, None)
    k_pwld = solve_k(grid, xs, "pwld")
    k_pwlc = solve_k(grid, xs, "pwlc")

    if rank == 0:
        print(f"SMMVacuumS2PIK={k_pi:.12f}")
        print(f"SMMVacuumS2PWLDK={k_pwld:.12f}")
        print(f"SMMVacuumS2PWLCK={k_pwlc:.12f}")
        print(f"SMMVacuumS2PIError={abs(k_pi - k_reference) / k_reference:.8e}")
        print(f"SMMVacuumS2PWLDError={abs(k_pwld - k_pi) / k_pi:.8e}")
        print(f"SMMVacuumS2PWLCError={abs(k_pwlc - k_pi) / k_pi:.8e}")
        sys.stdout.flush()
