import os
import sys
import numpy as np
"""
# Multigroup Cross Sections

## Background

OpenSn is not supplied with cross-section libraries.
Users are expected to supply their own cross-section data.
One may use other open-source software to generate multigroup cross sections
(e.g., NJOY, Dragon, OpenMC).

Here, we load the cross sections from an OpenMC HDF5 file.
"""
if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.xs import MultiGroupXS

if __name__ == "__main__":
    # load cross sections
    xs_uo2 = MultiGroupXS()
    filepath = "../../../test/python/modules/linear_boltzmann_solvers/transport_keigen/uo2.h5"
    xs_uo2.LoadFromOpenMC(filepath, "set1", 294.0)

    # Retrieve properties
    ng = xs_uo2.num_groups
    sca_order = xs_uo2.scattering_order
    print("num_groups       = ", ng)
    print("scattering_order = ", sca_order)

    # note cross sections are read-only objects of type <memoryview>
    # retrieve a numpy array
    siga = np.array(xs_uo2.sigma_a)
    print("siga = ", siga)
    # retrieve as list
    sigt = list(xs_uo2.sigma_t)
    print("sigt = ", sigt)
