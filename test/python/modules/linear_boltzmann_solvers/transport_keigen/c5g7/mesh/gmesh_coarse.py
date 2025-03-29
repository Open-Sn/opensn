if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator

    meshgen = FromFileMeshGenerator(
          filename = "mesh/2D_c5g7_coarse.msh",
    )
    grid = meshgen.Execute()
