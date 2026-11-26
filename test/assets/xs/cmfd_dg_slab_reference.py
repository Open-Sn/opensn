#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""Generate independent two-group DG-Sn eigenpairs for the CMFD slab tests.

No OpenSn code or exported operator is used. Each slab has linear elements of
width 0.5 cm with the exact mass matrix h/6 [[2,1],[1,2]], weak streaming,
upwind traces, specular reflection at z=0 and vacuum at the right end. Four
Gauss-Legendre directions use the normalized dmu/2 measure. The scattering
matrices below are indexed [destination, source]. Each case writes its
cross-section files and a reference file that holds k, the two outgoing group
currents, and then the [cell, node, group] scalar flux in cell order. The flux
is normalized so that the integrated nu-fission production is one.

Run this script in this directory to reproduce the checked-in assets.
"""
import hashlib
from pathlib import Path

import numpy as np
from scipy.linalg import eig

WIDTH = 0.5
MU, WEIGHTS = np.polynomial.legendre.leggauss(4)
WEIGHTS = WEIGHTS / 2.0
MASS = WIDTH * np.array([[2.0, 1.0], [1.0, 2.0]]) / 6.0

CASES = {
    # Alternating materials A and B in eight cells with full two-group scattering.
    "cmfd_material_dg_reference.txt": dict(
        cell_materials=[0, 1, 0, 1, 0, 1, 0, 1],
        xs_files=["cmfd_material_a_2g.xs", "cmfd_material_b_2g.xs"],
        totals=[[1.0, 2.0], [3.0, 4.0]],
        scatter=[[[0.6, 0.02], [0.1, 1.8]], [[2.5, 0.03], [0.2, 3.5]]],
        nu_sigma_f=[[0.2, 0.1], [0.4, 0.05]],
        chi=[[0.8, 0.2], [0.1, 0.9]]),
    # Four fuel cells next to two strong absorber cells. The thermal spectrum and the
    # vacuum boundary make the exact cell-average flux of the last cell negative.
    "cmfd_signed_slab_dg_reference.txt": dict(
        cell_materials=[0, 0, 0, 0, 1, 1],
        xs_files=["cmfd_signed_fuel_2g.xs", "cmfd_signed_absorber_2g.xs"],
        totals=[[1.0, 2.0], [10.0, 12.0]],
        scatter=[[[0.8, 0.0], [0.0, 1.8]], [[0.0, 0.0], [0.0, 0.0]]],
        nu_sigma_f=[[0.4, 0.4], [0.0, 0.0]],
        chi=[[1.0 / 101.0, 100.0 / 101.0], [0.0, 0.0]]),
}


def write_xs(path, totals, scatter, nu_sigma_f, chi):
    lines = ["NUM_GROUPS 2", "NUM_MOMENTS 1", "SIGMA_T_BEGIN"]
    lines += [f"{g} {value!r}" for g, value in enumerate(totals)]
    lines += ["SIGMA_T_END"]
    if any(value != 0.0 for value in nu_sigma_f):
        lines += ["NU_SIGMA_F_BEGIN"] + [f"{g} {v!r}" for g, v in enumerate(nu_sigma_f)]
        lines += ["NU_SIGMA_F_END", "NU_BEGIN", "0 2.0", "1 2.0", "NU_END", "CHI_BEGIN"]
        lines += [f"{g} {value!r}" for g, value in enumerate(chi)] + ["CHI_END"]
    lines += ["TRANSFER_MOMENTS_BEGIN"]
    for source in range(2):
        for destination in range(2):
            value = scatter[destination][source]
            if value != 0.0:
                lines.append(f"M_GFROM_GTO_VAL 0 {source} {destination} {value!r}")
    lines += ["TRANSFER_MOMENTS_END"]
    path.write_text("\n".join(lines) + "\n")


def solve(cell_materials, totals, scatter, nu_sigma_f, chi):
    cells, angles, groups = len(cell_materials), len(MU), 2
    totals, scatter = np.asarray(totals), np.asarray(scatter)
    nu_sigma_f, chi = np.asarray(nu_sigma_f), np.asarray(chi)
    size = cells * angles * groups * 2
    loss, fission = np.zeros((size, size)), np.zeros((size, size))

    def index(cell, angle, group, node):
        return ((cell * angles + angle) * groups + group) * 2 + node

    for c, material in enumerate(cell_materials):
        for m, direction in enumerate(MU):
            incoming = 0 if direction > 0 else 1
            neighbor = c - 1 if direction > 0 else c + 1
            for g in range(groups):
                block = direction * np.array([[-0.5, 0.5], [-0.5, 0.5]])
                block += totals[material, g] * MASS
                block[incoming, incoming] += abs(direction)
                for i in range(2):
                    for j in range(2):
                        loss[index(c, m, g, i), index(c, m, g, j)] += block[i, j]
                if 0 <= neighbor < cells:
                    loss[index(c, m, g, incoming), index(neighbor, m, g, 1 - incoming)] -= abs(
                        direction)
                elif c == 0 and direction > 0:
                    reflected = int(np.argmin(abs(MU + direction)))
                    loss[index(c, m, g, 0), index(c, reflected, g, 0)] -= direction
                for mp in range(angles):
                    for gp in range(groups):
                        for i in range(2):
                            for j in range(2):
                                coefficient = WEIGHTS[mp] * MASS[i, j]
                                loss[index(c, m, g, i), index(c, mp, gp, j)] -= (
                                    scatter[material, g, gp] * coefficient)
                                fission[index(c, m, g, i), index(c, mp, gp, j)] += (
                                    chi[material, g] * nu_sigma_f[material, gp] * coefficient)

    values, vectors = eig(fission, loss)
    eligible = np.flatnonzero((abs(values.imag) < 1e-10) & (values.real > 0))
    mode = int(eligible[np.argmax(values[eligible].real)])
    k = float(values[mode].real)
    assert k >= np.max(np.abs(np.delete(values, mode))) * (1.0 - 1e-10)
    psi = vectors[:, mode].real.reshape(cells, angles, groups, 2)
    phi = np.einsum("m,cmgn->cng", WEIGHTS, psi)
    production = sum(WIDTH * np.sum(phi[c].mean(axis=0) * nu_sigma_f[material])
                     for c, material in enumerate(cell_materials))
    psi /= production
    phi /= production
    residual = loss @ psi.ravel() - fission @ psi.ravel() / k
    relative_residual = np.linalg.norm(residual) / np.linalg.norm(fission @ psi.ravel() / k)
    assert relative_residual < 1e-12
    leakage = np.einsum("m,mg->g", WEIGHTS[MU > 0] * MU[MU > 0], psi[-1, MU > 0, :, 1])
    absorption = sum(
        WIDTH * np.sum(phi[c].mean(axis=0) * (totals[material] - scatter[material].sum(axis=0)))
        for c, material in enumerate(cell_materials))
    assert abs(absorption + leakage.sum() - 1.0 / k) < 1e-12
    return k, leakage, phi, relative_residual


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for filename, case in CASES.items():
        for material, xs_file in enumerate(case["xs_files"]):
            write_xs(here / xs_file, case["totals"][material], case["scatter"][material],
                     case["nu_sigma_f"][material], case["chi"][material])
        k, leakage, phi, relative_residual = solve(
            case["cell_materials"], case["totals"], case["scatter"], case["nu_sigma_f"],
            case["chi"])
        np.savetxt(here / filename, np.r_[k, leakage, phi.ravel()], fmt="%.17e",
                   header="k, outgoing group currents, then scalar flux [cell,node,group]; "
                          "integrated nu-fission production=1; "
                          "generator_sha256=" + digest)
        print(filename, "k =", k, "relative residual =", relative_residual)
