#include "function_lua_dimA_to_dimB.h"

#include "opensn/framework/console/chi_console.h"
#include "opensn/framework/logging/chi_log_exceptions.h"
#include "opensn/framework/chi_lua.h"

#include "opensn/framework/ChiObjectFactory.h"

namespace chi_math::functions
{

RegisterChiObject(chi_math::functions, LuaDimAToDimB);

chi::InputParameters
LuaDimAToDimB::GetInputParameters()
{
  chi::InputParameters params = FunctionDimAToDimB::GetInputParameters();

  // Inherits input_dimension and output_dimension

  // clang-format off
  params.SetGeneralDescription("Lua based parsed function");
  params.SetDocGroup("DocMathFunctions");
  // clang-format on

  params.AddRequiredParameter<std::string>("lua_function_name", "Name of the lua function");

  return params;
}

LuaDimAToDimB::LuaDimAToDimB(const chi::InputParameters& params)
  : FunctionDimAToDimB(params),
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{
}

std::vector<double>
LuaDimAToDimB::Evaluate(const std::vector<double>& vals) const
{
  const std::string fname = __PRETTY_FUNCTION__;
  lua_State* L = Chi::console.GetConsoleState();
  lua_getglobal(L, lua_function_name_.c_str());

  ChiLogicalErrorIf(not lua_isfunction(L, -1),
                    std::string("Attempted to access lua-function, ") + lua_function_name_ +
                      ", but it seems the function could "
                      "not be retrieved.");

  const size_t num_vals = vals.size();

  ChiInvalidArgumentIf(num_vals != InputDimension(),
                       std::string("Number of inputs do not match. ") +
                         "Attempted to evaluate with " + std::to_string(num_vals) +
                         " parameters but requires " + std::to_string(InputDimension()));

  lua_newtable(L);
  lua_Integer k = 0;
  for (double val : vals)
  {
    lua_pushinteger(L, ++k);
    lua_pushnumber(L, val);
    lua_settable(L, -3);
  }

  std::vector<double> result;
  // 1 arguments, 1 result (table), 0=original error object
  if (lua_pcall(L, 1, 1, 0) == 0)
  {
    LuaCheckTableValue(fname, L, -1);
    size_t table_length = lua_rawlen(L, -1);
    result.reserve(table_length);
    for (size_t i = 0; i < table_length; ++i)
    {
      lua_pushinteger(L, static_cast<lua_Integer>(i) + 1);
      lua_gettable(L, -2);
      result.push_back(lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
  }
  else
    throw std::logic_error(fname + " attempted to call lua-function, " + lua_function_name_ +
                           ", but the call failed. " + lua_tostring(L, -1));

  ChiLogicalErrorIf(result.size() != OutputDimension(),
                    std::string("Number of outputs after the function was ") +
                      "called does not "
                      "match the function specifications. A table is expected with " +
                      std::to_string(OutputDimension()) + " entries.");

  return result;
}

} // namespace chi_math::functions
