--[[ @doc
# Reading a 2D .obj file and extrude it

## Reading the Mesh
We start by reading a 2D obj file that we will extrude. (See a previous tutorial for reading 2D obj file). We inserts
2 layers between z=0 and z=1.1, followed by 3 layers, between z=1.1 and z=2.1.

Finally, we export the mesh to a VTU file.

The resulting mesh and material layout is shown below:

![Mesh_Material](images/extruded_mesh_material.png)

The resulting mesh and partition is shown below:

![Mesh_Partition](images/extruded_mesh_partition.png)

--]]
-- Setup the mesh
meshgen = mesh.ExtruderMeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
        filename="./tri_2mat_bc_1542.obj"
    }),
  },
  layers = {{z=1.1, n=2}, {z=2.1, n=3}},
  partitioner = PETScGraphPartitioner.Create({type="parmetis"})
})
mesh.MeshGenerator.Execute(meshgen)

mesh.ExportToVTK("Extruded_mesh_only")

--[[ @doc
## The rest of the simulation
The following line inserts the rest of the simulation data:
+ materials and sources
+ angular quadrature
+ LBS solver options and execution
+ VTK post-processing

You can view its contents in [transport_simulation_part_3D.lua](transport_simulation_part_3D.md)
--]]
-- Rest of the simulation
dofile("transport_simulation_part_3D.lua")

--[[ @doc
## Post-Processing via Field Functions
We extract the scalar flux (i.e., the first entry in the field function list; recall that lua indexing starts at 1)
and export it to a VTK file whose name is supplied by the user.
--]]
-- Get field functions
fflist,count = LBSGetScalarFieldFunctionList(phys)
vtk_basename = "extrusion"
ExportFieldFunctionToVTK(fflist[1],vtk_basename)

--[[ @doc
The resulting scalar flux is shown below:

![Scalar_flux](images/extruded_scalar_flux.png)
--]]

