meshgen = FromFileMeshGenerator(
    filename="../../../../../assets/mesh/triangle_mesh2x2.obj",
)
grid = meshgen.Execute()

math_SDM_Test02_Discontinuous(grid, "PWLD", False)
