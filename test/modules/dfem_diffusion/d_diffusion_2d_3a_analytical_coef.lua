-- Setup mesh
nodes = {}
N = 100
L = 2
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen1:Execute()

-- Set Material IDs
grid:SetUniformMaterialID(0)

function D_coef(i, pt)
  return 3.0 + pt.x + pt.y
end

function Q_ext(i, pt)
  return pt.x * pt.x
end

function Sigma_a(i, pt)
  return pt.x * pt.y * pt.y
end

-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = logvol.RPPLogicalVolume.Create({ xmin = 0.99999, xmax = 1000.0, infy = true, infz = true })
w_vol =
  logvol.RPPLogicalVolume.Create({ xmin = -1000.0, xmax = -0.99999, infy = true, infz = true })
n_vol = logvol.RPPLogicalVolume.Create({ ymin = 0.99999, ymax = 1000.0, infx = true, infz = true })
s_vol =
  logvol.RPPLogicalVolume.Create({ ymin = -1000.0, ymax = -0.99999, infx = true, infz = true })

e_bndry = "0"
w_bndry = "1"
n_bndry = "2"
s_bndry = "3"

grid:SetBoundaryIDFromLogicalVolume(e_vol, e_bndry, true)
grid:SetBoundaryIDFromLogicalVolume(w_vol, w_bndry, true)
grid:SetBoundaryIDFromLogicalVolume(n_vol, n_bndry, true)
grid:SetBoundaryIDFromLogicalVolume(s_vol, s_bndry, true)

diff_options = {
  boundary_conditions = {
    {
      boundary = e_bndry,
      type = "dirichlet",
      coeffs = { 0.0 },
    },
    {
      boundary = n_bndry,
      type = "dirichlet",
      coeffs = { 0.0 },
    },
    {
      boundary = s_bndry,
      type = "dirichlet",
      coeffs = { 0.0 },
    },
    {
      boundary = w_bndry,
      type = "dirichlet",
      coeffs = { 0.0 },
    },
  },
}

d_coef_fn = LuaScalarSpatialMaterialFunction.Create({ function_name = "D_coef" })
Q_ext_fn = LuaScalarSpatialMaterialFunction.Create({ function_name = "Q_ext" })
Sigma_a_fn = LuaScalarSpatialMaterialFunction.Create({ function_name = "Sigma_a" })

-- DFEM solver
phys1 = diffusion.DFEMDiffusionSolver.Create({
  name = "DFEMDiffusionSolver",
  mesh = grid,
  residual_tolerance = 1e-8,
})
phys1:SetOptions(diff_options)
phys1:SetDCoefFunction(d_coef_fn)
phys1:SetQExtFunction(Q_ext_fn)
phys1:SetSigmaAFunction(Sigma_a_fn)

phys1:Initialize()
phys1:Execute()

-- Get field functions
fflist = phys1:GetFieldFunctions()

-- Export VTU
if master_export == nil then
  fieldfunc.ExportToVTK(fflist[1], "CFEMDiff2D_analytic_coef", "flux")
end

-- Volume integrations
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })

ffvol = fieldfunc.FieldFunctionInterpolationVolume.Create()
ffvol:SetOperationType(OP_MAX)
ffvol:SetLogicalVolume(vol0)
ffvol:AddFieldFunction(fflist[1])

ffvol:Initialize(ffvol)
ffvol:Execute(ffvol)
maxval = ffvol:GetValue()

log.Log(LOG_0, string.format("Max-value=%.6f", maxval))
