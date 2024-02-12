#pragma once

#include "framework/lua.h"

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

/** Sets all cell-material id's to the supplied value.
 * \param material_id int The id.
 * \ingroup LuaVolumeMesher
 */
int VolumeMesherSetMatIDToAll(lua_State* L);
