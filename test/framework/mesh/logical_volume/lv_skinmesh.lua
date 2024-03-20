-- test for skin (surface) mesh used as a delimiter for a logical volume
-- set up orthogonal 3D geometry
nodes={}
N=50
L=5.0
xmin = -L/2
dx = L/N
for i=1,(N+1) do
  k=i-1
  nodes[i] = xmin + k*dx
end

meshgen = mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes,nodes,nodes},
})
mesh.MeshGenerator.Execute(meshgen)

-- assign mat ID 10 to whole domain
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0, 10)

-- create a logical volume as an analytical RPP
vol1 = mesh.RPPLogicalVolume.Create
({ xmin=-0.5,xmax=0.5,ymin=0.8,ymax=1.5, zmin=-1.5,zmax=0.5,  })
-- assign mat ID 11 to lv of RPP
mesh.SetMaterialIDFromLogicalVolume(vol1, 11)

-- create a logical volume as the interior of a skin mesh
surfmesh = mesh.SurfaceMeshCreate()
skin_mesh_file = "./cube_with_normals.obj"
mesh.SurfaceMeshImportFromOBJFile(surfmesh, skin_mesh_file)
lv_skinmesh = mesh.SurfaceMeshLogicalVolume.Create({surface_mesh_handle=surfmesh})
-- assign mat ID 15 to lv of skin mesh
mesh.SetMaterialIDFromLogicalVolume(lv_skinmesh, 15)

-- export to vtk
mesh.ExportToVTK("lv_skinmesh_out")
