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
 * ## _
 *
 * \ingroup LuaVolumeMesher
 * \author Jan
 */
int MeshSetProperty(lua_State* L);
