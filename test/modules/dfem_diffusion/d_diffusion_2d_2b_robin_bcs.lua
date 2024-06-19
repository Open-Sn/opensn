--############################################### Setup mesh
nodes = {}
N = 40
L = 2
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
mesh.SetMaterialIDFromLogicalVolume(vol0, 0)

vol1 =
  logvol.RPPLogicalVolume.Create({ xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5, infz = true })
mesh.SetMaterialIDFromLogicalVolume(vol1, 1)

D = { 1.0, 5.0 }
Q = { 0.0, 1.0 }
XSa = { 1.0, 1.0 }
function D_coef(i, pt)
  return D[i + 1] -- + x
end
function Q_ext(i, pt)
  return Q[i + 1] -- x*x
end
function Sigma_a(i, pt)
  return XSa[i + 1]
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

mesh.SetBoundaryIDFromLogicalVolume(e_vol, e_bndry)
mesh.SetBoundaryIDFromLogicalVolume(w_vol, w_bndry)
mesh.SetBoundaryIDFromLogicalVolume(n_vol, n_bndry)
mesh.SetBoundaryIDFromLogicalVolume(s_vol, s_bndry)

diff_options = {
  boundary_conditions = {
    {
      boundary = e_bndry,
      type = "robin",
      coeffs = { 0.25, 0.5, 0.0 }
    },
    {
      boundary = n_bndry,
      type = "reflecting"
    },
    {
      boundary = s_bndry,
      type = "reflecting"
    },
    {
      boundary = w_bndry,
      type = "robin",
      coeffs = { 0.25, 0.5, 0.1 }
    }
  }
}

-- DFEM solver
phys1 = diffusion.DFEMSolver.Create({
  name = "DFEMSolver",
  residual_tolerance = 1e-8
})
diffusion.SetOptions(phys1, diff_options)

solver.Initialize(phys1)
solver.Execute(phys1)

--############################################### Get field functions
fflist, count = solver.GetFieldFunctionList(phys1)

--############################################### Export VTU
if master_export == nil then
  fieldfunc.ExportToVTK(fflist[1], "DFEMDiff2D_RobinRefl", "flux")
end

--############################################### Volume integrations
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })

ffvol = fieldfunc.FFInterpolationCreate(VOLUME)
fieldfunc.SetProperty(ffvol, OPERATION, OP_AVG)
fieldfunc.SetProperty(ffvol, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(ffvol, ADD_FIELDFUNCTION, fflist[1])

fieldfunc.Initialize(ffvol)
fieldfunc.Execute(ffvol)
maxval = fieldfunc.GetValue(ffvol)

log.Log(LOG_0, string.format("Avg-value=%.6f", maxval))
