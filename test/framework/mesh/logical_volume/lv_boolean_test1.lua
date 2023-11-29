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

lv1 = mesh.SphereLogicalVolume.Create({r = 1.3, x=1.0, y=-1.0, z=2.0})

lv2 = mesh.RCCLogicalVolume.Create({r = 1.3,
                                        x0=-0.8, y0=-0.8, z0=-1.5,
                                        vx=1.0, vy=1.0, vz=3.0})
lv3 = mesh.BooleanLogicalVolume.Create
({
  parts = { { op=true, lv=lv2 },
            { op=false, lv=lv1 } }
})
VolumeMesherSetProperty(MATID_FROMLOGICAL, lv3, 1)

MeshHandlerExportMeshToVTK("lv_boolean_test1")
