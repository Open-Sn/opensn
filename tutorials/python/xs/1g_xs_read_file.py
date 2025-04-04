import os
import sys
import numpy as np
"""
# Multigroup Cross Sections

OpenSn is not supplied with cross-section libraries.
Users are expected to supply their own cross-section data.
One may use other open-source software to generate multigroup cross sections
(e.g., NJOY, Dragon, OpenMC).

Here, we load the cross sections from a plain text file (OpenSn format).
"""
if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.xs import MultiGroupXS

if __name__ == "__main__":

    # load cross sections
    xs_matA = MultiGroupXS()
    xs_matA.LoadFromOpenSn("xs_1g.xs")

    # Retrieve properties
    ng = xs_matA.num_groups
    sca_order = xs_matA.scattering_order
    print("num_groups       = ", ng)
    print("scattering_order = ", sca_order)

    # note cross sections are read-only objects of type <memoryview>
    siga = xs_matA.sigma_a
    print("type = ", type(siga))
    # retrieve a numpy array
    siga = np.array(xs_matA.sigma_a)
    print("siga = ", siga)
    # retrieve as list
    sigt = list(xs_matA.sigma_t)
    print("sigt = ", sigt)
