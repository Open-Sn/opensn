#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""Generate the synthetic CEPXS fixture used by the CSDA regressions.

This is artificial, nonphysical three-group data with energy bounds in MeV,
macroscopic total/transfer cross sections in cm^-1, stopping power in MeV/cm,
and responses in MeV/cm (energy) and cm^-1 (charge). The table is destination
major with eleven rows per destination group, encoded as little-endian float64
inside Fortran records with little-endian uint32 markers. It is original OpenSn
test data under MIT. Tests only read the checked-in file; regeneration is explicit:

    python3 test/assets/xs/generate_cepxs_synthetic.py /tmp/cepxs_synthetic_csda_3g.bxslib
"""

import argparse
import struct

G = 3
E_BOUNDS = [3.0, 2.0, 1.0, 0.0]
SIGMA_T = [1.0, 1.2, 1.5]
SIGMA_EDEP = [0.5, 0.35, 0.2]
SIGMA_CDEP = [0.05, 0.03, 0.01]
STOPPING = [0.4, 0.3, 0.2]
Q = [1.0, 0.0, 0.0]

# Transfer matrix in OpenSn convention: arrival g, departing gp.
S0 = [
    [0.20, 0.00, 0.00],
    [0.10, 0.15, 0.00],
    [0.03, 0.12, 0.10],
]


def write_fortran_record(f, payload):
    n = len(payload)
    f.write(struct.pack("<I", n))
    f.write(payload)
    f.write(struct.pack("<I", n))


def write_library(path):
    n_groups = G
    n_materials = 1
    n_entries = 11
    total_xs_row_1b = 8
    self_scatter_row_1b = 9
    n_moments = 1
    n_tables = n_materials * n_moments

    meta = struct.pack(
        "<8i",
        n_groups,
        n_materials,
        n_entries,
        total_xs_row_1b,
        self_scatter_row_1b,
        n_moments,
        0,
        n_tables,
    )

    table = [0.0] * (n_groups * n_entries)
    for g_to in range(n_groups):
        base = g_to * n_entries
        table[base + 1] = SIGMA_CDEP[g_to]  # row 2
        table[base + 2] = SIGMA_EDEP[g_to]  # row 3
        table[base + 4] = STOPPING[g_to]    # row 5
        table[base + 7] = SIGMA_T[g_to]     # row 8 total
        table[base + 8] = S0[g_to][g_to]    # row 9 self
        if g_to >= 1:
            table[base + 9] = S0[g_to][g_to - 1]   # row 10 from g-1
        if g_to >= 2:
            table[base + 10] = S0[g_to][g_to - 2]  # row 11 from g-2

    with open(path, "wb") as f:
        write_fortran_record(f, b"SYNTHETIC 3G CSDA TEST")
        write_fortran_record(f, meta)
        write_fortran_record(f, struct.pack(f"<{len(E_BOUNDS)}d", *E_BOUNDS))
        write_fortran_record(f, struct.pack(f"<{len(table)}d", *table))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    write_library(parser.parse_args().output)
