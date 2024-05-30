// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/lua.h"

namespace opensnlua
{

/**
 * Evaluates whether a point is within the logical volume.
 *
 * \param LVHandle int Handle to the logical volume.
 * \param Point_x double X-coordinate of the point.
 * \param Point_y double Y-coordinate of the point.
 * \param Point_z double Z-coordinate of the point.
 *
 * \return Sense true if inside the logical volume and false if outside.
 * \ingroup LuaLogicVolumes
 */
int LogVolPointSense(lua_State* L);

} // namespace opensnlua
