meshgen = FromFileMeshGenerator(
    filename="../../../../../assets/mesh/gmsh_all_hexes.vtu",
)
grid = meshgen.Execute()

math_SDM_Test02_Discontinuous(grid, "PWLD", False)
