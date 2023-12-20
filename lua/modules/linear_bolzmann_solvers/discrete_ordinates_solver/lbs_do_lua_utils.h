#pragma once

#include "framework/lua.h"

namespace opensnlua::lbs
{
/**Computes balance tables and prints it to the console.
 *
 * \param SolverIndex int Handle to the solver for which the list is to be
 * obtained.
 *
 * \ingroup LBSLuaFunctions
 * \author Jan
 */
int chiLBSComputeBalance(lua_State* L);

/**Computes the leakage for the specified groupset and boundary id.
 *
 * \param SolverIndex int Handle to the solver.
 * \param GroupSetHandle int Handle to the groupset.
 * \param BoundaryID int Id of the boundary for which leakage is to be computed.
 *
 * \return The leakage on a per group basis.
 *
 * \ingroup LBSLuaFunctions
 * \author Jan
 */
int chiLBSComputeLeakage(lua_State* L);
} // namespace opensnlua::lbs
