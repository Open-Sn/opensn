# Setup mesh
N = 11
L = 11.0
xmin = -L / 2
dx = L / N
nodes = [xmin + k * dx for k in range(N + 1)]
meshgen1 = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
grid = meshgen1.Execute()
grid.SetUniformBlockID(0)
grid.SetOrthogonalBoundaries()

# Run a unit simulation test for ray tracing
SimTest93_RayTracing(grid)

num_groups = 1
xs1g = MultiGroupXS()
xs1g.CreateSimpleOneGroup(1.0, 0.0)

# Add point source
src = []
for g in range(1, num_groups + 1):
    src.append(0.0)
src[0] = 1.0
pt_src = PointSource(
    location=[0.0, 0.0, 0.0],
    strength=src,
)

# Create a 2D XY quadrature
pquad = GLCProductQuadrature2DXY(12 * 4 * 2, 12 * 2 * 4 * 4)

# Setup Physics
solver_name = "LBS"
phys1 = DiscreteOrdinatesSolver(
    name=solver_name,
    mesh=grid,
    num_groups=num_groups,
    groupsets=[
        {
            "groups_from_to": (0, num_groups - 1),
            "angular_quadrature": pquad,
            "inner_linear_method": "petsc_richardson",
            "l_abs_tol": 1.0e-6,
            "l_max_its": 0,
        },
    ],
    xs_map=[
        {"block_ids": [0], "xs": xs1g},
    ],
    options={
        "scattering_order": 0,
        "point_sources": [pt_src],
        "field_function_prefix": solver_name,
    },
)
ss_solver = SteadyStateSolver(lbs_problem=phys1)
ss_solver.Initialize()
ss_solver.Execute()

ff_m0 = phys1.GetScalarFieldFunctionList()

# Export the first scalar field function to VTK
FieldFunctionGridBased.ExportMultipleToVTK([ff_m0[0]], "SimTest_93_" + solver_name)

MPIBarrier()
if rank == 0:
    import os
    os.system("rm SimTest_93*")
