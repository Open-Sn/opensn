-- Setup mesh
nodes = {}
N = 40
L = 1
xmin = 0
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
meshgen1:Execute()

-- Set Material IDs
mesh.SetUniformMaterialID(0)

-- governing law: -(u_xx + u_yy) = q, on domain [0,1]x[0,1]
-- when the exact solution is chosen u(x,y) = sin(pi.x) * sin(pi.y)
-- this automatically gives:
--    boundary = zero-Dirichlet on all 4 sides
--    volumetric source term: q(,x) = 2*pi*pi * sin(pi.x) * sin(pi.y)
-- the factor 2 is the dim of the problem
function D_coef(i, pt)
  return 1.0
end
function Q_ext(i, pt)
  return 2. * math.pi * math.pi * math.sin(math.pi * pt.x) * math.sin(math.pi * pt.y)
end
function Sigma_a(i, pt)
  return 0.0
end
function MMS_phi(i, pt)
  return math.sin(math.pi * pt.x) * math.sin(math.pi * pt.y)
end

-- Set boundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = logvol.RPPLogicalVolume.Create({ xmin = 0.99999, xmax = 1000.0, infy = true, infz = true })
w_vol =
  logvol.RPPLogicalVolume.Create({ xmin = -1000.0, xmax = -0.99999, infy = true, infz = true })
n_vol = logvol.RPPLogicalVolume.Create({ ymin = 0.99999, ymax = 1000.0, infx = true, infz = true })
s_vol =
  logvol.RPPLogicalVolume.Create({ ymin = -1000.0, ymax = -0.99999, infx = true, infz = true })

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

mesh.SetBoundaryIDFromLogicalVolume(e_vol, e_bndry)
mesh.SetBoundaryIDFromLogicalVolume(w_vol, w_bndry)
mesh.SetBoundaryIDFromLogicalVolume(n_vol, n_bndry)
mesh.SetBoundaryIDFromLogicalVolume(s_vol, s_bndry)

-- Call Lua Sim Test
SimTest_IP_MMS_L2error() --simtest_IP_MMS_L2_handle becomes available here

-- Export VTU
if master_export == nil then
  fieldfunc.ExportToVTK(simtest_IP_MMS_L2_handle, "DFEMDiff2D_MMS", "flux")
end

---- Volume integrations
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
--
ffvol = fieldfunc.FFInterpolationCreate(VOLUME)
fieldfunc.SetProperty(ffvol, OPERATION, OP_MAX)
fieldfunc.SetProperty(ffvol, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(ffvol, ADD_FIELDFUNCTION, simtest_IP_MMS_L2_handle)
--
fieldfunc.Initialize(ffvol)
fieldfunc.Execute(ffvol)
maxval = fieldfunc.GetValue(ffvol)

log.Log(LOG_0, string.format("Max-value=%.6f", maxval))
