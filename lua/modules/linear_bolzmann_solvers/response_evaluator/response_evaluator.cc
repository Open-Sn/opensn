// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/response_evaluator/response_evaluator.h"
#include "lua/framework/console/console.h"
#include "modules/linear_boltzmann_solvers/response_evaluator/response_evaluator.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(ClearResponseSources, lbs, ClearResponseSources);

int
ClearResponseSources(lua_State* L)
{
  const std::string fname = "lbs.ClearResponseSources";
  LuaCheckArgs<size_t>(L, fname);

  // Get the response evaluator
  const auto handle = LuaArg<size_t>(L, 1);
  auto& response_evaluator =
    GetStackItem<opensn::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  // Clear the sources
  response_evaluator.ClearForwardSources();
  return LuaReturn(L);
}

RegisterWrapperFunctionInNamespace(lbs,
                                   AddResponseBuffers,
                                   GetResponseBufferSyntax,
                                   AddResponseBuffers);

InputParameters
GetResponseBufferSyntax()
{
  InputParameters params;
  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>ResponseEvaluator</TT> object.");
  params.AddRequiredParameterArray("arg1", "Array of blocks for <TT>BufferOptionsBlock</TT>.");
  params.LinkParameterToBlock("arg1", "response::BufferOptionsBlock");
  return params;
}

ParameterBlock
AddResponseBuffers(const InputParameters& params)
{
  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const auto handle = params.GetParamValue<size_t>("arg0");
  auto& response_evaluator =
    GetStackItem<opensn::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  const auto buffer_params = params.GetParam("arg1");
  for (size_t p = 0; p < buffer_params.NumParameters(); ++p)
  {
    auto spec = opensn::ResponseEvaluator::BufferOptionsBlock();
    spec.AssignParameters(buffer_params.GetParam(p));
    response_evaluator.SetBufferOptions(spec);
  }

  return ParameterBlock();
}

RegisterWrapperFunctionInNamespace(lbs,
                                   AddResponseSources,
                                   GetResponseSourceSyntax,
                                   AddResponseSources);

InputParameters
GetResponseSourceSyntax()
{
  InputParameters params;
  params.SetGeneralDescription("Add sources to the response evaluator.");
  params.SetDocGroup("LBSResponseEvaluator");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>ResponseEvaluator</TT> object.");
  params.AddRequiredParameterBlock("arg1",
                                   "A block with the syntax of <TT>SourceOptionsBlock</TT>.");
  params.LinkParameterToBlock("arg1", "response::SourceOptionsBlock");

  return params;
}

ParameterBlock
AddResponseSources(const InputParameters& params)
{
  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const auto handle = params.GetParamValue<size_t>("arg0");
  auto& response_evaluator =
    GetStackItem<opensn::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  auto spec = opensn::ResponseEvaluator::SourceOptionsBlock();
  spec.AssignParameters(params.GetParam("arg1"));
  response_evaluator.SetSourceOptions(spec);

  return ParameterBlock();
}

RegisterLuaFunctionInNamespace(EvaluateResponse, lbs, EvaluateResponse);

int
EvaluateResponse(lua_State* L)
{
  const std::string fname = "lbs.EvaluateResponse";
  LuaCheckArgs<size_t, std::string>(L, fname);

  // Get the response evaluator
  const auto handle = LuaArg<size_t>(L, 1);
  auto& response_evaluator =
    GetStackItem<opensn::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  // Get the buffer name
  const auto buffer = LuaArg<std::string>(L, 2);

  // Compute the response
  double val = response_evaluator.EvaluateResponse(buffer);
  return LuaReturn(L, val);
}

} // namespace opensnlua
