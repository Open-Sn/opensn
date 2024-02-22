/** \page MeshTutorial_07 Mesh Tutorial 7: A 3D Extruded Mesh

 The idea behind extruded meshes is that the 3D complex geometrical features
 can be reasonably well projected into a 2D plane. This 2D geometry is then meshed
 (it is expected that this 2D mesh is unstructured).

 Then, the 2D meshes is extruded and appropriate material IDs/boundary IDs
 in 3D are provided to complete the 3D mesh.

Thus, in order to obtain a 3D extruded mesh, you need
 1. A 2D mesh
 2. Extrusion data

## The 2D mesh
 Currently, Chi-Tech can read in VTU and OBJ files. In the case of VTU files, you
 may provide an additional field such as 2D material IDs (i.e., material IDs that
 correspond to the 2D geometry).

 You may create your 2D mesh with our favorite mesh generator. To convert that mesh
 into an object readable by Chi-Tech, you want need to use mesh converters (in Python,
 there is meshio but there are other possibilities).

\code
 MeshHandlerCreate()

 mesh2d_file = "your_mesh.obj"
 umesh = UnpartitionedMeshFromWavefrontOBJ(mesh2d_file)
\endcode
or
\code
 MeshHandlerCreate()

 mesh2d_file = "your_mesh.vtu"
 umesh = UnpartitionedMeshFromVTU(mesh2d_file, "attribute")
\endcode

 ## The extrusion process

 \code
 SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 VolumeMesherCreate(VOLUMEMESHER_EXTRUDER,
                      ExtruderTemplateType.UNPARTITIONED_MESH,
                      umesh);

 nbr_planes_in_layer = 10
 height=0.1
 VolumeMesherSetProperty(EXTRUSION_LAYER,height,nbr_planes_in_layer,"Charlie");
 VolumeMesherSetProperty(EXTRUSION_LAYER,height,nbr_planes_in_layer,"Chuck");
 VolumeMesherSetProperty(EXTRUSION_LAYER,height,nbr_planes_in_layer,"Bob");
 VolumeMesherSetProperty(EXTRUSION_LAYER,height,nbr_planes_in_layer,"SarahConner");

 VolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)
 VolumeMesherSetKBAPartitioningPxPyPz(1,1,1)

 SurfaceMesherExecute();
 VolumeMesherExecute();

 mesh.ExportToVTK("export_mesh_without_IDs")
 \endcode

  ## Completing the 3D mesh

 To complete the 3D mesh, we need to provide material IDs and boundary IDs.
 We use the concept of a logical volume (LV). For each material zone, we will assign
 a material ID for cells (identified by their centroids) that are included in the LV.

 We have several ways of defining a Logical Volume, LV:
 1. Using pre-defined volumes, such as RPP (rectangular paralleliped)and  RCC (right circular
cylinder),
 2. Using a surface mesh read in as an```.obj``` file,
 3. Using lua functions to describe the  surface mesh of the LV.

 This should be repeated as many times as necessary to assign all 3D material IDs.

 ### LV of pre-defined type
 \code
-- Logical Volumes
my_LV = logvol.Create(RCC, 0, 0, 0.1, 0, 0, 0.2, 0.4)
mesh.SetMaterialIDFromLogicalVolume(Air, 1)
mesh.ExportToVTK("export_mesh_with_IDs")
 \endcode

  ### LV defined as read-in surfaces
 \code
surf_LV = SurfaceMeshCreate()
SurfaceMeshImportFromOBJFile(surf_LV, "LV_file.obj", false)
my_LV = logvol.Create(SURFACE, surf_LV)

mesh.SetMaterialIDFromLogicalVolume(Air, 1)
mesh.ExportToVTK("export_mesh_with_IDs")
 \endcode

  ### LV using a Lua function
 \code
mesh.SetMaterialIDFromFunction("my_LV_func.lua")
mesh.ExportToVTK("export_mesh_with_IDs")
 \endcode

 where the Lua function was previously defined, e.g.,
  \code
aaa
 \endcode
*/
