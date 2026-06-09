#!/usr/bin/env python3

def run(api, rank):
    GLCProductQuadrature2DXY = api["GLCProductQuadrature2DXY"]
    FieldFunctionInterpolationVolume = api["FieldFunctionInterpolationVolume"]
    RPPLogicalVolume = api["RPPLogicalVolume"]
    FromFileMeshGenerator = api["FromFileMeshGenerator"]
    DiscreteOrdinatesProblem = api["DiscreteOrdinatesProblem"]
    SteadyStateSourceSolver = api["SteadyStateSourceSolver"]
    MultiGroupXS = api["MultiGroupXS"]
    from uncollided_unstructured_utils import mesh_path, volume_integral

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.8, c=0.55)
    quadrature = GLCProductQuadrature2DXY(
        n_polar=2,
        n_azimuthal=24,
        scattering_order=1,
    )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-9,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux="uncollided_parallel_equivalence.h5",
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )
    integral = volume_integral(
        problem.GetScalarFluxFieldFunction()[0],
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
    balance = solver.ComputeBalanceTable()["balance"]
    if rank == 0:
        print(f"ParallelEquivalentIntegral={integral:.12e}")
        print(f"ParallelEquivalentBalance={balance:.12e}")
    if abs(balance) > 2.0e-5:
        raise RuntimeError(f"Combined balance failed: {balance}")
