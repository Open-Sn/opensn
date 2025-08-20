import sys
import os
import glob
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    # Append parent directory to locate the pyopensn modules
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator

if __name__ == "__main__":

    # Test Hex
    meshgen = FromFileMeshGenerator(
        filename="openfoam_meshes/Hex"
    )
    grid = meshgen.Execute()
    grid.ExportToPVTU("Hex")

    # Test Tets
    meshgen = FromFileMeshGenerator(
        filename="openfoam_meshes/Tets"
    )
    grid = meshgen.Execute()
    grid.ExportToPVTU("Tets")

    # Test Pyr
    meshgen = FromFileMeshGenerator(
        filename="openfoam_meshes/Pyr"
    )
    grid = meshgen.Execute()
    grid.ExportToPVTU("Pyr")

    # Test Wedge
    meshgen = FromFileMeshGenerator(
        filename="openfoam_meshes/Wedge"
    )
    grid = meshgen.Execute()
    grid.ExportToPVTU("Wedge")

for ending in ("*.pvtu", "*.vtu"):
    for f in glob.glob(ending):
        os.remove(f)
