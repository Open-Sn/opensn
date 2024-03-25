#include "lbs_lua_utils.h"

#include "framework/lua.h"

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"

#include "framework/console/console.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionNamespace(LBSSetPhiFromFieldFunction, lbs, SetPhiFromFieldFunction);

int
LBSSetPhiFromFieldFunction(lua_State* L)
{
  const std::string fname = "lbs.SetPhiFromFieldFunction";
  LuaCheckArgs<size_t, ParameterBlock>(L, fname);

  const auto handle = LuaArg<size_t>(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, handle, fname);

  auto specs = LuaArg<ParameterBlock>(L, 2);

  opensn::lbs::PhiSTLOption phi_option = opensn::lbs::PhiSTLOption::PHI_OLD;
  std::vector<size_t> moment_indices;
  std::vector<size_t> group_indices;

  specs.SetErrorOriginScope(fname);
  for (const auto& spec : specs.Parameters())
  {
    if (spec.Name() == "which_phi")
    {
      const auto phi_str = spec.GetValue<std::string>();
      if (phi_str == "old")
        phi_option = opensn::lbs::PhiSTLOption::PHI_OLD;
      else if (phi_str == "new")
        phi_option = opensn::lbs::PhiSTLOption::PHI_NEW;
      else
        OpenSnInvalidArgument(std::string("Parameter \"which_phi\" can only be"
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
      OpenSnInvalidArgument(std::string("Unsupported option ") + spec.Name());

  } // for each specification

  // Now call the function
  lbs_solver.SetPhiFromFieldFunctions(phi_option, moment_indices, group_indices);

  return LuaReturn(L);
}

} // namespace opensnlua::lbs
