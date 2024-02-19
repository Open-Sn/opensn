--[[ @doc
# 2D Orthogonal Grid

We have used the `OrthogonalMeshGenerator` in previous examples. Here, we use it again but also
introduce partitioning:
- KBA partitioning for regular grids that can be cut into right parallelepiped
- Parmetis partitioning, applicable with any grid type.

## Mesh
A simple orthogonal 2D mesh with KBA partitioner. We will run with 4 processes, so we will
partition the domain in 2x2 parts with cuts placed exactly at x=0 and y=0.


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
  node_sets = {nodes,nodes},
  partitioner = KBAGraphPartitioner.Create
  ({
    nx = 2, ny=2,
    xcuts = {0.}, ycuts = {0.}
  })
})
mesh.MeshGenerator.Execute(meshgen)

-- Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0,0)

--[[ @doc
## Export the mesh
The resulting mesh partition is shown below
![below](images/ortho_2D_KBA.png)
--]]
-- Exporting the mesh
mesh.ExportToVTK("ortho_2D_KBA")

--[[ @doc
## Mesh (again)
A simple orthogonal 2D mesh with Parmetis partitioner.
--]]
meshgen = mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes,nodes},
  partitioner = PETScGraphPartitioner.Create({type="parmetis"})

})
mesh.MeshGenerator.Execute(meshgen)

-- Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0,0)

--[[ @doc
## Export the mesh
On such a simple regular mesh, both partitioners are giving the same result. The Parmetis partition is shown below
![below](images/ortho_2D_Parmetis.png)
--]]
-- Exporting the mesh
mesh.ExportToVTK("ortho_2D_Parmetis")

