meshgen = FromFileMeshGenerator(
    filename="../../../../../assets/mesh/gmsh_all_tets.vtu",
)
grid = meshgen.Execute()

math_SDM_Test02_Discontinuous(grid, "PWLD", False)
