meshgen = FromFileMeshGenerator(
    filename="../../../../../assets/mesh/gmsh_all_tets.vtu",
)
grid = meshgen.Execute()

math_SDM_Test01_Continuous(grid, "PWLC", False)
