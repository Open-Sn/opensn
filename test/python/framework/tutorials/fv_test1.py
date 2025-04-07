import os

# Setup mesh
N = 100
L = 2.0
xmin = -L / 2
dx = L / N
nodes = [xmin + k * dx for k in range(N + 1)]
meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
grid = meshgen.Execute()
grid.SetUniformBlockID(0)

SimTest01_FV(grid)
MPIBarrier()
if rank == 0:
    os.system("rm CodeTut1_FV*")
