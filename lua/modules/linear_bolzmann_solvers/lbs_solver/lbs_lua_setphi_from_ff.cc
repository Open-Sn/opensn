#include "lbs_lua_utils.h"

#include "framework/chi_lua.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "framework/console/chi_console.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace lbs::common_lua_utils
{

RegisterLuaFunctionAsIs(chiLBSSetPhiFromFieldFunction);

int
chiLBSSetPhiFromFieldFunction(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, /*expected=*/2, /*given=*/num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 2);

  const size_t handle = lua_tointeger(L, 1);

  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, handle, fname);

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
    else if (spec.Name() == "m_ids") { moment_indices = spec.GetVectorValue<size_t>(); }
    else if (spec.Name() == "g_ids") { group_indices = spec.GetVectorValue<size_t>(); }
    else
      ChiInvalidArgument(std::string("Unsupported option ") + spec.Name());

  } // for each specification

  // Now call the function
  lbs_solver.SetPhiFromFieldFunctions(phi_option, moment_indices, group_indices);

  return 0;
}

} // namespace lbs::common_lua_utils
