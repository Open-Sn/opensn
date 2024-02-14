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
