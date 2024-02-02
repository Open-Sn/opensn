#pragma once

#include "framework/lua.h"

/** Creates a new volume mesher.
 *
 * \param Type int Volume Remesher type.
 * \param OtherArgs varying Additional arguments depending on `Type`.
 *
 * Remesher types:\n
 *  VOLUMEMESHER_EXTRUDER = Creates an extruded mesh from a 2D template mesh.
 *  Requires two additional arguments. `TemplateType` and `handle`. See below.\n
 *  VOLUMEMESHER_UNPARTITIONED = Create the mesh from the latest UnpartitionedMesh.
 *  Requires a single additional argument, `handle`, which is a handle to
 *  a valid unpartitioned mesh.\n
 *
 * ##_
 *
 * ###Extruder parameters
 *
 * When the mesher type is specified to be VOLUMEMESHER_EXTRUDER then two
 * additional arguments are required. `TemplateType` and `handle`.\n
 *
 * - `TemplateType` can for now only be
 *  `ExtruderTemplateType.UNPARTITIONED_MESH`.\n
 * - `handle` is a handle to the template mesh. When `TemplateType` is
 *   set to `ExtruderTemplateType.UNPARTITIONED_MESH` then the handle must point
 *   to a valid unpartitioned mesh.
 *
 * \ingroup LuaVolumeMesher
 * \author Jan
 */
int VolumeMesherCreate(lua_State* L);

/** Executes the volume meshing pipeline.
 *
 * \ingroup LuaVolumeMesher
 * \author Jan
 */
int VolumeMesherExecute(lua_State* L);

/** Sets a volume mesher property.
 *
 * \param PropertyIndex int Index of the property to change. See below
 * \param PropertyValue varying Value of the property.
 *
 * ##_
 *
 * ###PropertyIndex:
 *  FORCE_POLYGONS = <B>PropertyValue:[bool]</B> Forces the 2D Meshing to use
 *                   polygon cells even if the underlying surface mesh is
 *                   triangles. Expects a boolean value.\n
 *  MESH_GLOBAL = <B>PropertyValue:[bool]</B> Generate/Read the full mesh at each
 *                location. Expects a boolean value [Default=true].\n
 *  PARTITION_TYPE = <B>PartitionType</B>. See below.\n
 *  EXTRUSION_LAYER = <B>PropertyValue:[double,(int),(char)]</B> Adds a layer to
 * the extruder volume mesher if it exists. Expects 1 required parameter, the layer
 * height, followed by 2 optional parameters: number of subdivisions (defaults to
 * 1), and layer id (char)(defaults to nothing). Only supported if partition-type
 * is
 *                    ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```. \n
 *  CUTS_X = Adds a cut at the given x-value. Only supported if partition-type is
 *                    ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 *  CUTS_Y = Adds a cut at the given y-value. Only supported if partition-type is
 *                    ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 *  CUTS_Z = Adds a cut at the given z-value. Only supported if partition-type is
 *                    ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 *  PARTITION_X   = <B>PropertyValue:[int]</B> Number of partitions in X.
 *                     Only supported if partition-type is
 *                    ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 *  PARTITION_Y   = <B>PropertyValue:[int]</B> Number of partitions in Y.
 *                     Only supported if partition-type is
 *                    ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 *  PARTITION_Z   = <B>PropertyValue:[int]</B> Number of partitions in Z.
 *                     Only supported if partition-type is
 *                    ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 *  MATID_FROMLOGICAL = <B>LogicalVolumeHandle:[int],Mat_id:[int],
 *                      Sense:[bool](Optional, default:true)</B> Sets the material
 *                      id of cells that meet the sense requirement for the given
 *                      logical volume.\n
 *  BNDRYID_FROMLOGICAL = <B>LogicalVolumeHandle:[int],Bndry_name:[string],
 *                      Sense:[bool](Optional, default:true)</B> Sets the cell
 *                      boundary id to the specified value for cells
 *                      that meet the sense requirement for the given
 *                      logical volume.\n
 *  MATID_FROM_LUA_FUNCTION = <B>LuaFunctionName:[string]</B>. For each cell, will
 *                            call a lua function that can change the material id.
 *                            The lua function must have 4 parameters,
 *                            the cell's centroid x,y,z values (doubles) and the
 *                            current cell-material id (int). The function must
 *                            return a material id.
 *  BNDRYID_FROM_LUA_FUNCTION = <B>LuaFunctionName:[string]</B>. For each boundary
 *                            face, will call a lua function that can change the
 *                            boundary id.
 *                            The lua function must have 7 parameters,
 *                            the face's centroid x,y,z values (doubles), the
 *                            face's normal x,y,z values (double), and the
 *                            current face-boundary id (int). The function must
 *                            return a boundary id.
 * ## _
 *
 * ### PartitionType
 * Can be any of the following:
 *  - KBA_STYLE_XYZ
 *  - PARMETIS
 *
 * \ingroup LuaVolumeMesher
 * \author Jan
 */
int VolumeMesherSetProperty(lua_State* L);

/**Sets the Px, Py and Pz partititioning parameters for a
 * KBA-type partitioning. This also fixes the process count required to
 * a total of Px*Py*Pz.
 *
 * \param Px int Number partitions in x.
 * \param Py int Number partitions in y.
 * \param Pz int Number partitions in z.
 * \ingroup LuaVolumeMesher
 */
int VolumeMesherSetKBAPartitioningPxPyPz(lua_State* L);

/**Sets the x-cuts for KBA type partitioning with a lua array.
 * \ingroup LuaVolumeMesher
 */
int VolumeMesherSetKBACutsX(lua_State* L);

/**Sets the y-cuts for KBA type partitioning with a lua array.
 * \ingroup LuaVolumeMesher
 */
int VolumeMesherSetKBACutsY(lua_State* L);

/**Sets the z-cuts for KBA type partitioning with a lua array.
 * \ingroup LuaVolumeMesher
 */
int VolumeMesherSetKBACutsZ(lua_State* L);

/** Sets all cell-material id's to the supplied value.
 * \param material_id int The id.
 * \ingroup LuaVolumeMesher
 */
int VolumeMesherSetMatIDToAll(lua_State* L);

/** Sets boundary numbers on boundaries orthogonal to the cardinal directions
 * as xmax=0, xmin=1, ymax=2, ymin=3, zmax=4, zmin=5.
 *
 * \ingroup LuaVolumeMesher
 */
int VolumeMesherSetupOrthogonalBoundaries(lua_State* L);
