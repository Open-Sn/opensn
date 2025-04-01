# Setup mesh
N = 30
L = 2.0
xmin = -L / 2
dx = L / N
nodes = [xmin + k * dx for k in range(N + 1)]
meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
grid = meshgen.Execute()

math_SDM_Test01_Continuous(grid, "PWLC", False)
