--[[# test for boolean operations on logical volumes
lv1 = sphere
lv2 = right circular cylinder (rcc)
lv3 = in rcc but not in sphere
--]]

-- set up orthogonal 3D geometry
nodes={}
N=40
L=5
xmin = -L/2
dx = L/N
for i=1,(N+1) do
  k=i-1
  nodes[i] = xmin + k*dx
end
meshgen1 = mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes,nodes,nodes},
})
mesh.MeshGenerator.Execute(meshgen1)

-- assign matID 10 to all cells
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetProperty(MATID_FROMLOGICAL,vol0,10)

-- create logical volume lv1 as an analytical sphere
lv1 = mesh.SphereLogicalVolume.Create({r = 1.3, x=1.0, y=-1.0, z=2.0})

-- create logical volume lv2 as an analytical rcc
lv2 = mesh.RCCLogicalVolume.Create({r = 1.3,
                                        x0=-0.8, y0=-0.8, z0=-1.5,
                                        vx=1.0, vy=1.0, vz=3.0})

-- create logical volume lv3 as boolean: true if cell is in lv2 and false if in lv1
lv3 = mesh.BooleanLogicalVolume.Create
({
  parts = { { op=true, lv=lv2 },
            { op=false, lv=lv1 } }
})

-- assign matID 1 to all cells in lv3 which is the part of lv2 that is not in lv1
mesh.SetProperty(MATID_FROMLOGICAL, lv3, 1)

-- assign matID 5 to all cells in lv1
mesh.SetProperty(MATID_FROMLOGICAL, lv1, 5)

-- export to vtk
mesh.ExportToVTK("lv_boolean_test1")
