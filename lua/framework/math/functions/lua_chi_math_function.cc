#include "opensn/framework/chi_lua.h"

#include "opensn/framework/console/chi_console.h"

#include "opensn/framework/math/Functions/function_dimA_to_dimB.h"

#include "opensn/framework/logging/chi_log_exceptions.h"

namespace chi_math::functions::lua_utils
{

int chiFunctionDimAToDimBEvaluate(lua_State* L);

RegisterLuaFunctionAsIs(chiFunctionDimAToDimBEvaluate);

/**Evaluates a function of base type `FunctionXYZDimAToDimB`.
\param handle int. Handle to the function to evaluate.
\param params Varying. Table or individual arguments.

\return Varying Either a single number or a table of output values.*/
int
chiFunctionDimAToDimBEvaluate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  // ================================== Getting function object
  LuaCheckNilValue(fname, L, 1);
  const size_t handle = lua_tointeger(L, 1);

  const auto& function =
    Chi::GetStackItem<chi_math::FunctionDimAToDimB>(Chi::object_stack, handle, fname);

  // ================================== Getting params
  std::vector<double> params;
  if (lua_istable(L, 2))
  {
    auto table_block = chi_lua::TableParserAsParameterBlock::ParseTable(L, 2);
    ChiInvalidArgumentIf(table_block.Type() != chi::ParameterBlockType::ARRAY,
                         fname + ": Only an array type is allowed. Table can "
                                 "not have string keys.");
    params = table_block.GetVectorValue<double>();
  }
  else
  {
    for (int p = 2; p <= num_args; ++p)
    {
      LuaCheckNumberValue(fname, L, p);
      params.push_back(lua_tonumber(L, p));
    }
  }

  // ================================== Calling function
  const std::vector<double> values = function.Evaluate(params);

  //=================================== Parse outputs
  if (values.size() == 1)
  {
    lua_pushnumber(L, values.front());
    return 1;
  }
  // else

  lua_newtable(L);
  for (size_t k = 0; k < values.size(); ++k)
  {
    lua_pushinteger(L, static_cast<lua_Integer>(k) + 1);
    lua_pushnumber(L, values[k]);
    lua_settable(L, -3);
  }
  return 1;
}

} // namespace chi_math::functions::lua_utils
