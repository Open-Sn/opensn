// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/console/console.h"
#include "lua/framework/math/functions/lua_scalar_spatial_material_function.h"
#include "framework/parameters/input_parameters.h"
#include "framework/runtime.h"
#include "framework/physics/solver.h"
#include "modules/diffusion/cfem_diffusion_solver.h"
#include "modules/diffusion/dfem_diffusion_solver.h"
#include "modules/diffusion/fv_diffusion_solver.h"
#include "modules/diffusion/mg_diffusion_solver.h"

namespace opensnlua
{

namespace
{

std::shared_ptr<LuaScalarSpatialMaterialFunction>
CreateFunction(const std::string& function_name)
{
  opensn::ParameterBlock blk;
  blk.AddParameter("lua_function_name", function_name);
  opensn::InputParameters params = LuaScalarSpatialMaterialFunction::GetInputParameters();
  params.AssignParameters(blk);
  return std::make_shared<LuaScalarSpatialMaterialFunction>(params);
}

} // namespace

// CFEM diffusion solver

opensn::ParameterBlock
CFEMSetOptions(const opensn::InputParameters& params)
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

  const size_t handle = params.ParamValue<size_t>("arg0");
  auto& solver =
    opensn::GetStackItem<opensn::CFEMDiffusionSolver>(opensn::object_stack, handle, fname);
  solver.SetDCoefFunction(d_coef_function);
  solver.SetQExtFunction(q_ext_function);
  solver.SetSigmaAFunction(sigma_a_function);

  auto options_params = opensn::CFEMDiffusionSolver::OptionsBlock();
  options_params.AssignParameters(params.Param("arg1"));

  solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

// CFEM diffusion solver

opensn::ParameterBlock
DFEMSetOptions(const opensn::InputParameters& params)
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

  const size_t handle = params.ParamValue<size_t>("arg0");
  auto& solver =
    opensn::GetStackItem<opensn::DFEMDiffusionSolver>(opensn::object_stack, handle, fname);
  solver.SetDCoefFunction(d_coef_function);
  solver.SetQExtFunction(q_ext_function);
  solver.SetSigmaAFunction(sigma_a_function);

  auto options_params = opensn::DFEMDiffusionSolver::OptionsBlock();
  options_params.AssignParameters(params.Param("arg1"));

  solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

// FV diffusion solver

opensn::ParameterBlock
FVSetOptions(const opensn::InputParameters& params)
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

  const size_t handle = params.ParamValue<size_t>("arg0");
  auto& solver =
    opensn::GetStackItem<opensn::FVDiffusionSolver>(opensn::object_stack, handle, fname);
  solver.SetDCoefFunction(d_coef_function);
  solver.SetQExtFunction(q_ext_function);
  solver.SetSigmaAFunction(sigma_a_function);

  auto options_params = opensn::FVDiffusionSolver::OptionsBlock();
  options_params.AssignParameters(params.Param("arg1"));

  solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

// CFEM MG diffusion solver

opensn::ParameterBlock
CFEMMGSetOptions(const opensn::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const size_t handle = params.ParamValue<size_t>("arg0");
  auto& solver =
    opensn::GetStackItem<opensn::MGDiffusionSolver>(opensn::object_stack, handle, fname);

  auto options_params = opensn::MGDiffusionSolver::OptionsBlock();
  options_params.AssignParameters(params.Param("arg1"));

  solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

//

namespace
{

opensn::InputParameters
GetSyntax_SetOptions()
{
  opensn::InputParameters params;
  params.SetGeneralDescription("Set options from a large list of parameters");
  params.AddRequiredParameter<size_t>("arg0", "Handle to a `CFEMSolver` object.");
  params.AddRequiredParameterBlock("arg1", "Block of parameters for `OptionsBlock`");
  params.LinkParameterToBlock("arg1", "OptionsBlock");
  return params;
}

opensn::ParameterBlock
SetOptions(const opensn::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const size_t handle = params.ParamValue<size_t>("arg0");
  auto& solver = opensn::GetStackItem<opensn::Solver>(opensn::object_stack, handle, fname);

  // FIXME: dispatch to the right solver until there is a common diffusion solver class,
  // FIXME: then this will be taken care of by polymorphism
  if (dynamic_cast<opensn::CFEMDiffusionSolver*>(&solver))
    return CFEMSetOptions(params);
  else if (dynamic_cast<opensn::DFEMDiffusionSolver*>(&solver))
    return DFEMSetOptions(params);
  else if (dynamic_cast<opensn::FVDiffusionSolver*>(&solver))
    return FVSetOptions(params);
  else if (dynamic_cast<opensn::MGDiffusionSolver*>(&solver))
    return CFEMMGSetOptions(params);
  else
    throw std::runtime_error("Unknown solver type");
}

} // namespace

RegisterWrapperFunctionInNamespace(diffusion, SetOptions, GetSyntax_SetOptions, SetOptions);

} // namespace opensnlua
