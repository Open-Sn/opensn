#pragma once

#include "framework/lua.h"
#include "framework/parameters/input_parameters.h"

namespace lbs
{
class LBSSolver;
}

namespace opensnlua::lbs
{

// void SetBoundaryOptions(LBSSolver& lbs_solver,
//                         const objects::InputParameters& params);

opensn::InputParameters GetSyntax_SetOptions();
opensn::ParameterBlock SetOptions(const opensn::InputParameters& params);

int LBSSetOptions(lua_State* L);

/**Sets the internal phi vector to the value in the associated
 * field function.
 * \param handle int Handle to the lbs-based object.
 * \param specs Table Various options in a table. Detailed below.
 *
 * ##_
 *
 * ### Example usage
 * Example:
 * \code
 * LBSSetPhiFromFieldFunction(phys1,
 * {
 *   which_phi = "old",  --Optional
 *   m_ids = {0,1,2,3},  --Optional Empty means all of them
 *   g_ids = {}          --Optional Empty means all of them
 * })
 * \endcode
 *
 * ### Table keys
 * `which_phi`\n
 * <I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
 * Optional. Can be "old" or "new". Denotes which phi version to copy to.
 * Default: `"old"`.\n\n
 *
 * `m_ids`\n
 * <I>type=</I><span style="color: blue;"><TT>ARRAY</TT></span>
 * Optional. Array of moment IDs. If this is empty all the moments will be copied.
 * Default: `{}`.\n\n
 *
 * `g_ids`\n
 * <I>type=</I><span style="color: blue;"><TT>ARRAY</TT></span>
 * Optional. Array of group IDs. If this is empty all the groups will be copied.
 * Default: `{}`.\n\n
 *
 */
int LBSSetPhiFromFieldFunction(lua_State* L);

/**
 * Adds a point source to an LBS solver.
 *
 * \param SolverIndex int Handle to the solver.
 * \param PointSourceIndex int Handle to the point source
 * \ingroup LBSLuaFunctions
 */
int AddPointSource(lua_State* L);

/**
 * Clears the point sources within an LBS solver.
 * This is primarily useful for adjoint problems.
 *
 * \param int Handle to LBS solver.
 * \ingroup LBSLuaFunctions
 */
int ClearPointSources(lua_State* L);

/**
 * Adds a distributed source to an LBS solver.
 *
 * \param SolverIndex int Handle to the solver.
 * \param DistributedSourceHandle int Handle to the distributed source
 * \ingroup LBSLuaFunctions
 */
int AddDistributedSource(lua_State* L);

/**
 * Clears the distributed sources within an LBS solver.
 *
 * \param int Handle to LBS solver.
 * \ingroup LBSLuaFunctions
 */
int ClearDistributedSources(lua_State* L);

void RegisterLuaEntities(lua_State* L);
} // namespace opensnlua::lbs
