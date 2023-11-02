#include "lbs_lua_utils.h"

#include "chi_lua.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "console/chi_console.h"
#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs::common_lua_utils
{

RegisterLuaFunctionAsIs(chiLBSSetPhiFromFieldFunction);

/**Sets the internal phi vector to the value in the associated
field function.
\param handle int Handle to the lbs-based object.
\param specs Table Various options in a table. Detailed below.

##_

### Example usage
Example:
\code
chiLBSSetPhiFromFieldFunction(phys1,
{
  which_phi = "old",  --Optional
  m_ids = {0,1,2,3},  --Optional Empty means all of them
  g_ids = {}          --Optional Empty means all of them
})
\endcode

### Table keys
`which_phi`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
Optional. Can be "old" or "new". Denotes which phi version to copy to.
Default: `"old"`.\n\n

`m_ids`\n
<I>type=</I><span style="color: blue;"><TT>ARRAY</TT></span>
Optional. Array of moment IDs. If this is empty all the moments will be copied.
Default: `{}`.\n\n

`g_ids`\n
<I>type=</I><span style="color: blue;"><TT>ARRAY</TT></span>
Optional. Array of group IDs. If this is empty all the groups will be copied.
Default: `{}`.\n\n

*/
int chiLBSSetPhiFromFieldFunction(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, /*expected=*/2, /*given=*/num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 2);

  const size_t handle = lua_tointeger(L, 1);

  auto& lbs_solver =
    Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, handle, fname);

  auto specs = chi_lua::TableParserAsParameterBlock::ParseTable(L, 2);

  lbs::PhiSTLOption phi_option = PhiSTLOption::PHI_OLD;
  std::vector<size_t> moment_indices;
  std::vector<size_t> group_indices;

  specs.SetErrorOriginScope(fname);
  for (const auto& spec : specs.Parameters())
  {
    if (spec.Name() == "which_phi")
    {
      const auto phi_str = spec.GetValue<std::string>();
      if (phi_str == "old") phi_option = PhiSTLOption::PHI_OLD;
      else if (phi_str == "new")
        phi_option = PhiSTLOption::PHI_NEW;
      else
        ChiInvalidArgument(std::string("Parameter \"which_phi\" can only be"
                                       " \"old\" or \"new\". ") +
                           "\"" + phi_str + "\" is not allowed.");
    }
    else if (spec.Name() == "m_ids")
    {
      moment_indices = spec.GetVectorValue<size_t>();
    }
    else if (spec.Name() == "g_ids")
    {
      group_indices = spec.GetVectorValue<size_t>();
    }
    else
      ChiInvalidArgument(std::string("Unsupported option ") + spec.Name());

  } // for each specification

  // ============================================ Now call the function
  lbs_solver.SetPhiFromFieldFunctions(
    phi_option, moment_indices, group_indices);

  return 0;
}

} // namespace lbs::common_lua_utils