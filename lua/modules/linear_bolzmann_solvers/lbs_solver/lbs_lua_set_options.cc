// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/lbs_solver/lbs_lua_utils.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensnlua
{

namespace
{

opensn::InputParameters
GetSyntax_SetOptions()
{
  opensn::InputParameters params;

  params.SetGeneralDescription("Set options from a large list of parameters");
  params.SetDocGroup("LBSLuaFunctions");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>LBSSolver</TT> object.");
  params.AddRequiredParameterBlock("arg1", "Block of parameters for <TT>OptionsBlock</TT>");
  params.LinkParameterToBlock("arg1", "OptionsBlock");

  return params;
}

opensn::ParameterBlock
SetOptions(const opensn::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const size_t handle = params.GetParamValue<size_t>("arg0");
  auto& lbs_solver = opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, handle, fname);

  auto options_params = opensn::LBSSolver::OptionsBlock();
  options_params.AssignParameters(params.GetParam("arg1"));

  lbs_solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

} // namespace

RegisterWrapperFunctionInNamespace(lbs, SetOptions, GetSyntax_SetOptions, SetOptions);

} // namespace opensnlua
