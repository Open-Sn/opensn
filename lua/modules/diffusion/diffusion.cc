// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/console/console.h"
#include "framework/runtime.h"
#include "framework/physics/solver_base/solver.h"
#include "modules/cfem_diffusion/cfem_diffusion_solver.h"
#include "lua/modules/cfem_diffusion/cfem_diff_solver.h"
#include "modules/dfem_diffusion/dfem_diffusion_solver.h"
#include "lua/modules/dfem_diffusion/dfem_diff_solver.h"
#include "modules/fv_diffusion/fv_diffusion_solver.h"
#include "lua/modules/fv_diffusion/fv_diff_solver.h"
#include "modules/mg_diffusion/mg_diffusion_solver.h"
#include "lua/modules/mg_diffusion/mg_diff_solver.h"

namespace opensnlua
{

namespace
{

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

  const size_t handle = params.GetParamValue<size_t>("arg0");
  auto& solver = opensn::GetStackItem<opensn::Solver>(opensn::object_stack, handle, fname);

  // FIXME: dispatch to the right solver until there is a common diffusion solver class,
  // FIXME: then this will be taken care of by polymorphism
  if (dynamic_cast<opensn::diffusion::CFEMSolver*>(&solver))
    return CFEMSetOptions(params);
  else if (dynamic_cast<opensn::diffusion::DFEMSolver*>(&solver))
    return DFEMSetOptions(params);
  else if (dynamic_cast<opensn::diffusion::FVSolver*>(&solver))
    return FVSetOptions(params);
  else if (dynamic_cast<opensn::diffusion::MGSolver*>(&solver))
    return CFEMMGSetOptions(params);
  else
    throw std::runtime_error("Unknown solver type");
}

} // namespace

RegisterWrapperFunctionInNamespace(diffusion, SetOptions, GetSyntax_SetOptions, SetOptions);

} // namespace opensnlua
