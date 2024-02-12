#include "lbs_lua_utils.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

#include "framework/console/console.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensnlua::lbs
{

RegisterWrapperFunctionNamespace(lbs, SetOptions, GetSyntax_SetOptions, SetOptions);

opensn::InputParameters
GetSyntax_SetOptions()
{
  opensn::InputParameters params;

  params.SetGeneralDescription("Set options from a large list of parameters");
  params.SetDocGroup("LBSLuaFunctions");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>lbs::LBSSolver</TT> object.");
  params.AddRequiredParameterBlock("arg1", "Block of parameters for <TT>lbs::OptionsBlock</TT>");
  params.LinkParameterToBlock("arg1", "lbs::OptionsBlock");

  return params;
}

opensn::ParameterBlock
SetOptions(const opensn::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const size_t handle = params.GetParamValue<size_t>("arg0");
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, handle, fname);

  auto options_params = opensn::lbs::LBSSolver::OptionsBlock();
  options_params.AssignParameters(params.GetParam("arg1"));

  lbs_solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

} // namespace opensnlua::lbs
