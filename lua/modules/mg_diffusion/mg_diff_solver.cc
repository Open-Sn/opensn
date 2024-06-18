// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/mg_diffusion/mg_diff_solver.h"
#include "lua/framework/console/console.h"
#include "modules/mg_diffusion/mg_diffusion_solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua
{

opensn::ParameterBlock
CFEMMGSetOptions(const opensn::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const size_t handle = params.GetParamValue<size_t>("arg0");
  auto& solver =
    opensn::GetStackItem<opensn::diffusion::MGSolver>(opensn::object_stack, handle, fname);

  auto options_params = opensn::diffusion::MGSolver::OptionsBlock();
  options_params.AssignParameters(params.GetParam("arg1"));

  solver.SetOptions(options_params);

  return opensn::ParameterBlock();
}

} // namespace opensnlua
