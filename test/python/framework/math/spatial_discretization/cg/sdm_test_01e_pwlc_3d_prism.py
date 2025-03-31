meshgen = ExtruderMeshGenerator(
    inputs=[
        FromFileMeshGenerator(
            filename="../../../../../assets/mesh/TriangleMesh2x2.obj",
        )
    ],
    layers=[
        {"z": 0.4, "n": 2},
        {"z": 0.8, "n": 2},
        {"z": 1.2, "n": 2},
        {"z": 1.6, "n": 2}
    ]
)
grid = meshgen.Execute()

math_SDM_Test01_Continuous(grid, "PWLC", False)
