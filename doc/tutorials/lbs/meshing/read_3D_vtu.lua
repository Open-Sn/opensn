--[[ @doc
# Reading a 3D .vtu file

## Reading the Mesh
We start by reading a 3D vtu file. The resulting mesh and partition is shown below:

![Mesh_Partition](images/read_3d_vtu_partition.png)

--]]
-- Setup the mesh
meshgen = mesh.MeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
        filename="./GMSH_AllTets.vtu"
    }),
  },
  partitioner = PETScGraphPartitioner.Create({type="parmetis"})
})
mesh.MeshGenerator.Execute(meshgen)

mesh.ExportToVTK("Read_3D_mesh_only")

--[[ @doc
## The rest of the simulation
The following line inserts the rest of the simulation data:
+ materials and sources
+ angular quadrature
+ LBS solver options and execution (caveat: single aggregation on tetrahedral meshes, more on this later)
+ VTK post-processing

You can view its contents in [transport_simulation_part_3D_single_agg.lua](transport_simulation_part_3D_single_agg.md)
--]]
-- Rest of the simulation
dofile("transport_simulation_part_3D_single_agg.lua")

--[[ @doc
## Post-Processing via Field Functions
We extract the scalar flux (i.e., the first entry in the field function list; recall that lua indexing starts at 1)
and export it to a VTK file whose name is supplied by the user.
--]]
-- Get field functions
fflist,count = LBSGetScalarFieldFunctionList(phys)
vtk_basename = "read_3d_vtu"
ExportFieldFunctionToVTK(fflist[1],vtk_basename)

--[[ @doc
The resulting scalar flux is shown below:

![Scalar_flux](images/read_3d_vtu_scalar_flux.png)
--]]

