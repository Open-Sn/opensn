--[[ @doc
# 3D Orthogonal Grid

This is very similar to the previous `2D Orthogonal Grid` tutorial. We use the two partitioners again (KBA and Parmetis).

## Mesh and KBA partition
A simple orthogonal 3D mesh with KBA partitioner. We will run with 8 processes, so we will
partition the domain in 2x2x2 parts with cuts placed exactly at x=0, y=0, and z=0.


--]]
-- Setup the mesh
nodes={}
n_cells=10
length=2.
xmin = -length/2.
dx = length/n_cells
for i=1,(n_cells+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen = mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes,nodes,nodes},
  partitioner = mesh.KBAGraphPartitioner.Create
  ({
    nx = 2, ny=2, nz=2,
    xcuts = {0.}, ycuts = {0.}, zcuts = {0.}
  })
})
mesh.MeshGenerator.Execute(meshgen)

-- Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0,0)

--[[ @doc
## Export the mesh
The resulting mesh partition is shown below
![below](images/ortho_3D_KBA.png)
--]]
-- Exporting the mesh
mesh.ExportToPVTU("ortho_3D_KBA")

--[[ @doc
## Mesh (again) and Parmetis partition
A simple orthogonal 3D mesh with Parmetis partitioner.
--]]
meshgen = mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes,nodes,nodes},
  partitioner = mesh.PETScGraphPartitioner.Create({type="parmetis"})

})
mesh.MeshGenerator.Execute(meshgen)

-- Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0,0)

--[[ @doc
## Export the mesh
Note that now, both partitioners are not giving the same result. The Parmetis partition is shown below
![below](images/ortho_3D_Parmetis.png)
--]]
-- Exporting the mesh
mesh.ExportToPVTU("ortho_3D_Parmetis")
