#pragma once

#include "framework/chi_lua.h"
#include "framework/parameters/input_parameters.h"

namespace lbs
{
class LBSSolver;
}

namespace lbs::common_lua_utils
{

// void SetBoundaryOptions(LBSSolver& lbs_solver,
//                         const chi_objects::InputParameters& params);

chi::InputParameters GetSyntax_SetOptions();
chi::ParameterBlock SetOptions(const chi::InputParameters& params);

int chiLBSSetOptions(lua_State* L);

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
 * chiLBSSetPhiFromFieldFunction(phys1,
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
int chiLBSSetPhiFromFieldFunction(lua_State* L);
void RegisterLuaEntities(lua_State* L);
} // namespace lbs::common_lua_utils
