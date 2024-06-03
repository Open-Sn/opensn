// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/math/functions/lua_dim_a_to_dim_b.h"
#include "framework/console/console.h"
#include "framework/logging/log_exceptions.h"
#include "framework/lua.h"

#include "framework/object_factory.h"

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(math::functions, LuaDimAToDimB);

InputParameters
LuaDimAToDimB::GetInputParameters()
{
  InputParameters params = FunctionDimAToDimB::GetInputParameters();

  // Inherits input_dimension and output_dimension

  params.SetGeneralDescription("Lua based parsed function");
  params.SetDocGroup("DocMathFunctions");

  params.AddRequiredParameter<std::string>("lua_function_name", "Name of the lua function");

  return params;
}

LuaDimAToDimB::LuaDimAToDimB(const InputParameters& params)
  : FunctionDimAToDimB(params),
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{
}

std::vector<double>
LuaDimAToDimB::Evaluate(const std::vector<double>& vals) const
{
  const std::string fname = __PRETTY_FUNCTION__;

  const size_t num_vals = vals.size();
  OpenSnInvalidArgumentIf(num_vals != InputDimension(),
                          std::string("Number of inputs do not match. ") +
                            "Attempted to evaluate with " + std::to_string(num_vals) +
                            " parameters but requires " + std::to_string(InputDimension()));

  lua_State* L = console.GetConsoleState();
  auto result = LuaCall<std::vector<double>>(L, lua_function_name_, vals);
  OpenSnLogicalErrorIf(result.size() != OutputDimension(),
                       "Number of outputs after the function was called does not match the "
                       "function specifications. A table is expected with " +
                         std::to_string(OutputDimension()) + " entries.");
  return result;
}

} // namespace opensnlua
