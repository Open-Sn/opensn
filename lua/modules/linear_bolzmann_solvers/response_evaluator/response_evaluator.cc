#include "response_evaluator.h"
#include "modules/linear_boltzmann_solvers/response_evaluator/response_evaluator.h"
#include "lua/framework/console/console.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionNamespace(ClearResponseSources, lbs, ClearResponseSources);

int
ClearResponseSources(lua_State* L)
{
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);

  // Get the response evaluator
  const int handle = lua_tonumber(L, 1);
  auto& response_evaluator =
    GetStackItem<opensn::lbs::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  // Clear the sources
  response_evaluator.ClearForwardSources();
  return 0;
}

RegisterWrapperFunctionNamespace(lbs,
                                 AddResponseBuffers,
                                 GetResponseBufferSyntax,
                                 AddResponseBuffers);

InputParameters
GetResponseBufferSyntax()
{
  InputParameters params;
  params.AddRequiredParameter<size_t>("arg0",
                                      "Handle to a <TT>lbs::ResponseEvaluator</TT> object.");
  params.AddRequiredParameterArray("arg1", "Array of blocks for <TT>BufferOptionsBlock</TT>.");
  params.LinkParameterToBlock("arg1", "lbs::response::BufferOptionsBlock");
  return params;
}

ParameterBlock
AddResponseBuffers(const InputParameters& params)
{
  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const auto handle = params.GetParamValue<size_t>("arg0");
  auto& response_evaluator =
    GetStackItem<opensn::lbs::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  const auto buffer_params = params.GetParam("arg1");
  for (size_t p = 0; p < buffer_params.NumParameters(); ++p)
  {
    auto spec = opensn::lbs::ResponseEvaluator::BufferOptionsBlock();
    spec.AssignParameters(buffer_params.GetParam(p));
    response_evaluator.SetBufferOptions(spec);
  }

  return ParameterBlock();
}

RegisterWrapperFunctionNamespace(lbs,
                                 AddResponseSources,
                                 GetResponseSourceSyntax,
                                 AddResponseSources);

InputParameters
GetResponseSourceSyntax()
{
  InputParameters params;
  params.SetGeneralDescription("Add sources to the response evaluator.");
  params.SetDocGroup("LBSResponseEvaluator");

  params.AddRequiredParameter<size_t>("arg0",
                                      "Handle to a <TT>lbs::ResponseEvaluator</TT> object.");
  params.AddRequiredParameterBlock("arg1",
                                   "A block with the syntax of <TT>SourceOptionsBlock</TT>.");
  params.LinkParameterToBlock("arg1", "lbs::response::SourceOptionsBlock");

  return params;
}

ParameterBlock
AddResponseSources(const InputParameters& params)
{
  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const auto handle = params.GetParamValue<size_t>("arg0");
  auto& response_evaluator =
    GetStackItem<opensn::lbs::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  auto spec = opensn::lbs::ResponseEvaluator::SourceOptionsBlock();
  spec.AssignParameters(params.GetParam("arg1"));
  response_evaluator.SetSourceOptions(spec);

  return ParameterBlock();
}

RegisterLuaFunctionNamespace(EvaluateResponse, lbs, EvaluateResponse);

int
EvaluateResponse(lua_State* L)
{
  const auto num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__, 2, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);
  LuaCheckStringValue(__FUNCTION__, L, 2);

  // Get the response evaluator
  const auto handle = lua_tointeger(L, 1);
  auto& response_evaluator =
    GetStackItem<opensn::lbs::ResponseEvaluator>(object_stack, handle, __FUNCTION__);

  // Get the buffer name
  const auto buffer = lua_tostring(L, 2);

  // Compute the response
  double val = response_evaluator.EvaluateResponse(buffer);
  lua_pushnumber(L, static_cast<lua_Number>(val));
  return 1;
}

} // namespace opensnlua::lbs
