#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D delayed transient with a ramped XS (monotonic sanity).

Test intent
- Similar to the prompt ramp, but with delayed neutrons enabled; ensures precursor coupling remains stable
  as XS changes stepwise over time.

Physics
- 1-group, 1 precursor. nu*Sigma_f ramps upward in time. The fission rate should grow monotonically for
  this case under reflecting BCs.

Gold values
- No fixed numeric gold; we require monotonic non-decreasing FR and finite values over the ramp.

What we check and why
- TRANSIENT_OK enforces finite response and non-decreasing FR, catching delayed-source instability or
  XS update timing issues.
"""

import math
import os


def xs_path(name):
    return os.path.join(os.path.dirname(__file__), name)


def build_mesh_3d(n, length):
    dx = length / n
    nodes = [i * dx for i in range(n + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    return grid


def read_block_value(file_path, block_begin, block_end):
    in_block = False
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == block_begin:
                in_block = True
                continue
            if line == block_end:
                in_block = False
                continue
            if in_block:
                parts = line.split()
                if len(parts) == 2 and parts[0] == "0":
                    return float(parts[1])
    return None


def load_xs_scalar_params(xs_file):
    nu_prompt = read_block_value(xs_file, "NU_PROMPT_BEGIN", "NU_PROMPT_END")
    nu_delayed = read_block_value(xs_file, "NU_DELAYED_BEGIN", "NU_DELAYED_END")
    lam = read_block_value(xs_file, "PRECURSOR_DECAY_CONSTANTS_BEGIN",
                           "PRECURSOR_DECAY_CONSTANTS_END")
    frac_yield = read_block_value(xs_file, "PRECURSOR_FRACTIONAL_YIELDS_BEGIN",
                                  "PRECURSOR_FRACTIONAL_YIELDS_END")
    return nu_prompt, nu_delayed, lam, frac_yield


if __name__ == "__main__":
    grid = build_mesh_3d(n=4, length=8.0)

    xs_list = []
    for i in range(5):
        xs = MultiGroupXS()
        xs.LoadFromOpenSn(xs_path(f"xs1g_delayed_ramp_{i}.cxs"))
        xs_list.append(xs)
    xs_crit = xs_list[0]
    xs_super = xs_list[-1]
    xs_scalar_file = xs_path("xs1g_delayed_ramp_0.cxs")

    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    num_groups = 1
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
                "gmres_restart_interval": 50,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = TransientKEigenSolver(problem=phys, max_iters=200, k_tol=1.0e-10)
    solver.Initialize()

    # XS-based kinetics parameters (1-group, 1-precursor)
    sigma_a = xs_super.sigma_a[0]
    nu_prompt, nu_delayed, lam, frac_yield = load_xs_scalar_params(xs_scalar_file)
    # Use beta from precursor fractional yield and nu_total from nu_prompt + nu_delayed
    beta = frac_yield
    v = 1.0 / xs_super.inv_velocity[0]

    def k_from_nu_sigma_f(nu_sigma_f):
        return nu_sigma_f / sigma_a

    def rho_from_k(k):
        return (k - 1.0) / k

    # Ramp parameters
    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    t_end = 0.2
    ramp_time = t_end

    def mix_factor(t):
        if t <= 0.0:
            return 0.0
        if t >= ramp_time:
            return 1.0
        return t / ramp_time

    def xs_index(t):
        f = mix_factor(t)
        return min(int(f * (len(xs_list) - 1) + 1.0e-12), len(xs_list) - 1)

    def nu_sigma_f_of_t(t):
        sigma_f = xs_list[xs_index(t)].sigma_f[0]
        nu_total = nu_prompt + nu_delayed
        return sigma_f * nu_total

    def rho_of_t(t):
        k = k_from_nu_sigma_f(nu_sigma_f_of_t(t))
        return rho_from_k(k)

    def Lambda_of_t(t):
        return 1.0 / (v * nu_sigma_f_of_t(t))

    nu_sigma_f_crit = nu_sigma_f_of_t(0.0)

    fr0 = phys.ComputeFissionProduction("new")

    rel_tol = 2.0e-2
    ok = True
    growth_ok = True
    last_ratio = None
    print("step time ratio_numeric ratio_analytic")
    step = 0
    while phys.GetTime() < t_end:
        step += 1
        t_from = phys.GetTime()

        # Update XS mix for current step (piecewise-constant over dt)
        idx = xs_index(t_from)
        phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_list[idx]}])

        solver.Step()
        fr_new = phys.ComputeFissionProduction("new")
        solver.Advance()
        t_to = phys.GetTime()

        ratio_num = fr_new / fr0
        ratio_ana = 1.0

        print(f"{step:4d} {t_to:10.4e} {ratio_num:12.6e} {ratio_ana:12.6e}")
        if (not math.isfinite(fr_new)) or (not math.isfinite(ratio_num)):
            ok = False
        elif last_ratio is not None and ratio_num < (last_ratio - 1.0e-6):
            growth_ok = False
        last_ratio = ratio_num

    ok = ok and growth_ok
    print(f"TRANSIENT_OK {1 if ok else 0}")
