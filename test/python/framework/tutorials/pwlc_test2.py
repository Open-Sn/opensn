import os
import math

# Setup mesh
N = 10
L = 2.0
xmin = -L / 2
dx = L / N
nodes = [xmin + k * dx for k in range(N + 1)]
meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
grid = meshgen.Execute()
grid.SetUniformBlockID(0)


def MMS_phi(pt):
    return math.cos(math.pi * pt[0]) + math.cos(math.pi * pt[1])


def MMS_q(pt):
    return math.pi**2 * (math.cos(math.pi * pt[0]) + math.cos(math.pi * pt[1]))


SimTest04_PWLC(grid)
MPIBarrier()
if rank == 0:
    os.system("rm CodeTut4_PWLC*")
