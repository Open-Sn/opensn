#include "framework/lua.h"

#include "framework/console/console.h"

#include "framework/math/functions/function_dimA_to_dimB.h"

#include "framework/logging/log_exceptions.h"

using namespace opensn;

namespace opensnlua
{

/**Evaluates a function of base type `FunctionXYZDimAToDimB`.
 * \param handle int. Handle to the function to evaluate.
 * \param params Varying. Table or individual arguments.
 *
 * \return Varying Either a single number or a table of output values.
 */
int FunctionDimAToDimBEvaluate(lua_State* L);

RegisterLuaFunctionAsIs(FunctionDimAToDimBEvaluate);

int
FunctionDimAToDimBEvaluate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  LuaCheckArgs<size_t, std::vector<double>>(L, "FunctionDimAToDimBEvaluate");

  // Getting function object
  const auto handle = LuaArg<size_t>(L, 1);
  const auto& function =
    opensn::GetStackItem<FunctionDimAToDimB>(opensn::object_stack, handle, fname);

  // Getting params
  std::vector<double> params;
  if (lua_istable(L, 2))
  {
    auto table_block = TableParserAsParameterBlock::ParseTable(L, 2);
    OpenSnInvalidArgumentIf(table_block.Type() != ParameterBlockType::ARRAY,
                            fname + ": Only an array type is allowed. Table can "
                                    "not have string keys.");
    params = table_block.GetVectorValue<double>();
  }
  else
  {
    const int num_args = lua_gettop(L);
    for (int p = 2; p <= num_args; ++p)
      params.push_back(LuaArg<double>(L, p));
  }

  // Calling function
  const std::vector<double> values = function.Evaluate(params);

  // Parse outputs
  if (values.size() == 1)
  {
    LuaPush(L, values.front());
    return 1;
  }
  // else

  LuaPush(L, values);
  return 1;
}

} // namespace opensnlua
