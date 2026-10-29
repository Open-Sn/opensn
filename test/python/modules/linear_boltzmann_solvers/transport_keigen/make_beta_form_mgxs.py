#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Writes beta_form_mgxs.h5, an OpenMC MGXS library in OpenMC's beta form for
keigenvalue_transport_1d_openmcxs_beta_infinite_medium.py.

One fissionable, isotropic two-group set "fuel" at 294 K with upscatter and two delayed groups:
"nu-fission" and "chi" with "beta" by delayed group and fissioning group ([delayed group][group]),
"decay-rate", and no prompt/delayed datasets. The data are in DATA, which the test also uses for
its references. Regenerate with `python3 make_beta_form_mgxs.py` (needs h5py).
"""

import os

DATA = {
    "total": [1.0, 2.0],
    # (from-group, to-group, value)
    "transfers": [(0, 0, 0.6), (0, 1, 0.25), (1, 1, 1.5), (1, 0, 0.02)],
    "fission": [0.04, 0.2],
    "nu": [2.5, 2.4],
    "chi": [0.9, 0.1],
    # beta[j][g]: delayed fraction of delayed group j for fission in group g
    "beta": [[0.002, 0.003], [0.004, 0.0035]],
    "decay-rate": [0.0124, 0.305],
    "inverse-velocity": [0.1, 1.0],
}


def main():
    import h5py
    import numpy as np

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "beta_form_mgxs.h5")
    num_groups = len(DATA["total"])
    with h5py.File(path, "w") as f:
        f.attrs["filetype"] = np.bytes_("mgxs")
        f.attrs["version"] = np.array([1, 0])
        f.attrs["energy_groups"] = np.int64(num_groups)
        f.attrs["delayed_groups"] = np.int64(len(DATA["beta"]))
        f.attrs["group structure"] = np.array([1.0e-5, 0.625, 2.0e7])

        fuel = f.create_group("fuel")
        fuel.attrs["fissionable"] = True
        fuel.attrs["order"] = np.int64(0)
        fuel.attrs["representation"] = np.bytes_("isotropic")
        fuel.attrs["scatter_format"] = np.bytes_("legendre")
        fuel.attrs["scatter_shape"] = np.bytes_("[G][G'][Order]")

        xs = fuel.create_group("294K")
        for name in ("total", "fission", "chi", "decay-rate", "inverse-velocity"):
            xs.create_dataset(name, data=np.array(DATA[name]))
        xs.create_dataset("nu-fission", data=np.array(DATA["nu"]) * np.array(DATA["fission"]))
        xs.create_dataset("beta", data=np.array(DATA["beta"]))

        # OpenMC scatter data: for each incoming group, the outgoing groups g_min..g_max (1-based)
        scatter = np.zeros((num_groups, num_groups))
        for gf, gt, v in DATA["transfers"]:
            scatter[gf, gt] = v
        scatter_data = xs.create_group("scatter_data")
        scatter_data.create_dataset("g_min", data=np.ones(num_groups, dtype=np.int32))
        scatter_data.create_dataset("g_max", data=np.full(num_groups, num_groups, dtype=np.int32))
        scatter_data.create_dataset("scatter_matrix", data=scatter.flatten())


if __name__ == "__main__":
    main()
