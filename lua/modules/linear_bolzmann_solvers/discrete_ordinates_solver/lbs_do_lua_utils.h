// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua::lbs
{

/**
 * Computes balance tables and prints it to the console.
 *
 * \param SolverIndex int Handle to the solver for which the list is to be
 * obtained.
 *
 * \ingroup LBSLuaFunctions
 * \author Jan
 */
int LBSComputeBalance(lua_State* L);

/**
 * Computes the group-wise leakage on all boundaries.
 *
 * \param SolverIndex int Handle to the solver.
 * \param BoundaryNames array List of boundary names. These must be the
 *      standard boundary names used in OpenSn:
 *      xmax, xmin, ymax, ymin, zmax, zmin
 *
 * \return A table mapping boundary names to group-wise leakage arrays.
 *
 * \ingroup LBSLuaFunctions
 * \author Zachary K. Hardy
 */
int ComputeLeakage(lua_State* L);

} // namespace opensnlua::lbs
