-- SDM: PWLD

-- Groups
Ng = 64

-- Angles
Npolar = 7
Nazimuthal = 8

meshgen1 = mesh.DistributedMeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "1n_64rpn_256cr.msh",
    }),
  },
})
mesh.MeshGenerator.Execute(meshgen1)

-- Material
materials = {}
materials[0] = mat.AddMaterial("Test Material")
mat.SetProperty(materials[0], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "../diag_XS_64g_1mom_c0.99.xs")
src = {}
for g = 1, Ng do
  src[g] = 0.0
end
mat.SetProperty(materials[0], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

-- Boundary conditions
bsrc = {}
for g = 1, Ng do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0 / 4.0 / math.pi
lbs_options = {
  boundary_conditions = { { name = "xmin", type = "isotropic", group_strength = bsrc } },
  scattering_order = 0,
}

-- Quadrature
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, Npolar, Nazimuthal)

-- Set up solver
gs1 = { 0, Ng - 1 }
lbs_block = {
  num_groups = Ng,
  groupsets = {
    {
      groups_from_to = gs1,
      angular_quadrature_handle = pquad0,
      angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "krylov_richardson",
      l_abs_tol = 1.0e-6,
      l_max_its = 9,
      gmres_restart_interval = 10,
    },
  },
}
phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys1 })

-- Solve
solver.Initialize(ss_solver)
solver.Execute(ss_solver)
