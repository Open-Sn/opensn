meshgen = FromFileMeshGenerator(
    filename="../../../../../assets/mesh/GMSH_AllTets.vtu",
)
grid = meshgen.Execute()

math_SDM_Test01_Continuous(grid, "PWLC", False)
