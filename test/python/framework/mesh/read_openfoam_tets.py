import sys
import os

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    # Append parent directory to locate the pyopensn modules
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator

if __name__ == "__main__":

    # Test Tets
    meshgen = FromFileMeshGenerator(
        filename="openfoam_meshes/tets"
    )
    grid = meshgen.Execute()
