// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/cfem_diffusion/cfem_diffusion_solver.h"
#include "framework/console/console.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/framework/math/functions/lua_scalar_spatial_material_function.h"

using namespace opensn;

namespace opensnlua
{

namespace
{

std::shared_ptr<LuaScalarSpatialMaterialFunction>
CreateFunction(const std::string& function_name)
{
  ParameterBlock blk;
  blk.AddParameter("lua_function_name", function_name);
  InputParameters params = LuaScalarSpatialMaterialFunction::GetInputParameters();
  params.AssignParameters(blk);
  return std::make_shared<LuaScalarSpatialMaterialFunction>(params);
}

opensn::InputParameters
GetSyntax_SetOptions()
{
  opensn::InputParameters params;
  params.SetGeneralDescription("Set options from a large list of parameters");
  params.AddRequiredParameter<size_t>("arg0", "Handle to a `diffusion::CFEMSolver` object.");
  params.AddRequiredParameterBlock("arg1", "Block of parameters for `diffusion::OptionsBlock`");
  params.LinkParameterToBlock("arg1", "diffusion::OptionsBlock");
  return params;
}

opensn::ParameterBlock
SetOptions(const opensn::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  auto d_coef_function = CreateFunction("D_coef");
  opensn::function_stack.push_back(d_coef_function);

  auto q_ext_function = CreateFunction("Q_ext");
  opensn::function_stack.push_back(q_ext_function);

  auto sigma_a_function = CreateFunction("Sigma_a");
  opensn::function_stack.push_back(sigma_a_function);

  const size_t handle = params.GetParamValue<size_t>("arg0");
  auto& solver =
    opensn::GetStackItem<opensn::diffusion::CFEMSolver>(opensn::object_stack, handle, fname);
  solver.SetDCoefFunction(d_coef_function);
  solver.SetQExtFunction(q_ext_function);
  solver.SetSigmaAFunction(sigma_a_function);

  auto options_params = opensn::diffusion::CFEMSolver::OptionsBlock();
  options_params.AssignParameters(params.GetParam("arg1"));

  solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

} // namespace

RegisterWrapperFunctionInNamespace(diffusion, SetOptions, GetSyntax_SetOptions, SetOptions);

} // namespace opensnlua
