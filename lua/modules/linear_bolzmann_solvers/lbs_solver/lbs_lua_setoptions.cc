#include "lbs_lua_utils.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "framework/console/chi_console.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace lbs::common_lua_utils
{

// ##################################################################
RegisterWrapperFunction(/*namespace_in_lua=*/lbs,
                        /*name_in_lua=*/SetOptions,
                        /*syntax_function=*/GetSyntax_SetOptions,
                        /*actual_function=*/SetOptions);

chi::InputParameters
GetSyntax_SetOptions()
{
  chi::InputParameters params;

  // clang-format off
  params.SetGeneralDescription("Set options from a large list of parameters");
  params.SetDocGroup("LBSLuaFunctions");

  params.AddRequiredParameter<size_t>(
    "arg0", "Handle to a <TT>lbs::LBSSolver</TT> object.");
  params.AddRequiredParameterBlock(
    "arg1", "Block of parameters for <TT>lbs::OptionsBlock</TT>");
  params.LinkParameterToBlock("arg1", "lbs::OptionsBlock");

  // clang-format on

  return params;
}

chi::ParameterBlock
SetOptions(const chi::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const size_t handle = params.GetParamValue<size_t>("arg0");
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, handle, fname);

  auto options_params = LBSSolver::OptionsBlock();
  options_params.AssignParameters(params.GetParam("arg1"));

  lbs_solver.SetOptions(options_params);

  return chi::ParameterBlock();
}

} // namespace lbs::common_lua_utils
