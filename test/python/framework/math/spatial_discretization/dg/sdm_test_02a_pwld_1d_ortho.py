N = 100
L = 2.0
xmin = -L / 2
dx = L / N
nodes = [xmin + k * dx for k in range(N + 1)]
meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
grid = meshgen.Execute()

math_SDM_Test02_Discontinuous(grid, "PWLD", False)
