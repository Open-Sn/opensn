meshgen = FromFileMeshGenerator(
    filename="../../../../../assets/mesh/TriangleMesh2x2.obj",
)
grid = meshgen.Execute()

math_SDM_Test01_Continuous(grid, "PWLC", False)
