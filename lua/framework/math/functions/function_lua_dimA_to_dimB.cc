#include "function_lua_dimA_to_dimB.h"

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
  lua_State* L = console.GetConsoleState();
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

} // namespace opensnlua
