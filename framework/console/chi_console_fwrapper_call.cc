#include "chi_console.h"

#include "chi_log_exceptions.h"
#include "chi_lua.h"

namespace chi
{

RegisterLuaFunction(Console::LuaWrapperCall, chi_console, LuaWrapperCall);

/** This function expects at least 1 parameter. The first parameter must
 * be a string indicating the registered name of the function wrapper to
 * call. All other parameters will be forwarded to the function wrapper.*/
int Console::LuaWrapperCall(lua_State* L)
{
  const int num_args = lua_gettop(L);
  // We do not check for the required parameters here because we want
  // to make this function call as fast as possible. Besides, via the
  // static registration we should never run into an issue here.

  auto& console = Console::GetInstance();

  const auto& registry = console.function_wrapper_registry_;

  const std::string fname = lua_tostring(L, 1);

  ChiLogicalErrorIf(registry.count(fname) == 0,
                    std::string("Wrapper with name \"") + fname +
                      "\" not in console registry.");

  const auto& reg_entry = registry.at(fname);

  auto input_params = reg_entry.get_in_params_func();

  ParameterBlock main_arguments_block;
  for (int p = 2; p <= num_args; ++p)
  {
    const std::string arg_name = "arg" + std::to_string(p-2);

    if (lua_isboolean(L, p))
      main_arguments_block.AddParameter(arg_name, lua_toboolean(L, p));
    else if (lua_isinteger(L, p))
      main_arguments_block.AddParameter(arg_name, lua_tointeger(L, p));
    else if (lua_isnumber(L, p))
      main_arguments_block.AddParameter(arg_name, lua_tonumber(L, p));
    else if (lua_isstring(L, p))
      main_arguments_block.AddParameter(arg_name, lua_tostring(L, p));
    else if (lua_istable(L, p))
    {
      auto block = chi_lua::TableParserAsParameterBlock::ParseTable(L, p);
      block.SetBlockName(arg_name);
      std::string scope = fname + ":";
      scope.append(arg_name + " ");
      block.SetErrorOriginScope(scope);
      main_arguments_block.AddParameter(block);
    }
    else
      ChiInvalidArgument("In call to \"" + fname +
                         "\": Unsupported argument "
                         "type \"" +
                         lua_typename(L, lua_type(L, p)) + "\" encountered.");
  }
  // Set input parameters here
  input_params.SetErrorOriginScope(fname + "()");
  input_params.AssignParameters(main_arguments_block);

  auto output_params = reg_entry.call_func(input_params);

  output_params.SetErrorOriginScope(fname + ":output:");
  chi_lua::PushParameterBlock(L, output_params);

  const int num_sub_params = static_cast<int>(output_params.NumParameters());

  return output_params.IsScalar() ? 1 : num_sub_params;
}

} // namespace chi_objects