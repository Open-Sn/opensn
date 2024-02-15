#pragma once

#include "framework/lua.h"

/**
 * Sets mesh property.
 *
 * \param PropertyIndex int Index of the property to change. See below
 * \param PropertyValue varying Value of the property.
 *
 * ##_
 *
 * ###PropertyIndex:
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
 * \ingroup LuaVolumeMesher
 * \author Jan
 */
int MeshSetProperty(lua_State* L);
