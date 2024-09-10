-- test for skin (surface) mesh used as a delimiter for a logical volume
-- set up orthogonal 3D geometry
nodes = {}
N = 50
L = 5.0
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { nodes, nodes, nodes },
})
meshgen:Execute()

-- assign mat ID 10 to whole domain
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
mesh.SetMaterialIDFromLogicalVolume(vol0, 10, true)

-- create a logical volume as an analytical RPP
vol1 = logvol.RPPLogicalVolume.Create({
  xmin = -0.5,
  xmax = 0.5,
  ymin = 0.8,
  ymax = 1.5,
  zmin = -1.5,
  zmax = 0.5,
})
-- assign mat ID 11 to lv of RPP
mesh.SetMaterialIDFromLogicalVolume(vol1, 11, true)
-- create a logical volume as the interior of a skin mesh
surfmesh = mesh.SurfaceMesh.Create({})
surfmesh:ImportFromOBJFile("./cube_with_normals.obj", false, Vector3(0, 0, 0))
lv_skinmesh = logvol.SurfaceMeshLogicalVolume.Create({ surface_mesh = surfmesh })
-- assign mat ID 15 to lv of skin mesh
mesh.SetMaterialIDFromLogicalVolume(lv_skinmesh, 15, true)

-- export to vtk
mesh.ExportToPVTU("lv_skinmesh_out")
