#include "framework/console/console.h"
#ifdef OPENSN_WITH_LUA
#include "lua/modules/modules_lua.h"
#include "framework/lua.h"
#endif
#include "config.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/logging/log_exceptions.h"
#include "framework/mpi/mpi.h"
#include "framework/utils/utils.h"
#include <iostream>

namespace chi::lua_utils
{
int chiMakeObject(lua_State* L);
}

namespace chi
{

RegisterLuaFunction(Console::LuaWrapperCall, chi_console, LuaWrapperCall);

Console&
Console::GetInstance() noexcept
{
  static Console singleton;
  return singleton;
}

Console::Console() noexcept
#ifdef OPENSN_WITH_LUA
  : console_state_(luaL_newstate())
#endif
{
}

void
Console::FlushConsole()
{
#ifdef OPENSN_WITH_LUA
  try
  {
    for (auto& command : command_buffer_)
    {
      bool error = luaL_dostring(console_state_, command.c_str());
      if (error)
      {
        Chi::log.LogAll() << lua_tostring(console_state_, -1);
        lua_pop(console_state_, 1);
      }
    }
  }
  catch (const std::exception& e)
  {
    Chi::log.LogAllError() << e.what();
    Chi::Exit(EXIT_FAILURE);
  }
#endif
}

#ifdef OPENSN_WITH_LUA
int
Console::LuaWrapperCall(lua_State* L)
{
  const int num_args = lua_gettop(L);
  // We do not check for the required parameters here because we want
  // to make this function call as fast as possible. Besides, via the
  // static registration we should never run into an issue here.

  auto& console = Console::GetInstance();

  const auto& registry = console.function_wrapper_registry_;

  const std::string fname = lua_tostring(L, 1);

  ChiLogicalErrorIf(registry.count(fname) == 0,
                    std::string("Wrapper with name \"") + fname + "\" not in console registry.");

  const auto& reg_entry = registry.at(fname);

  auto input_params = reg_entry.get_in_params_func();

  ParameterBlock main_arguments_block;
  for (int p = 2; p <= num_args; ++p)
  {
    const std::string arg_name = "arg" + std::to_string(p - 2);

    if (lua_isboolean(L, p)) main_arguments_block.AddParameter(arg_name, lua_toboolean(L, p));
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
#endif

void
Console::RunConsoleLoop(char*) const
{
  Chi::log.Log() << "Console loop started. "
                 << "Type \"exit\" to quit (or Ctl-C).";

  /** Wrapper to an MPI_Bcast call for a single integer
   * broadcast from location 0. */
  auto BroadcastSingleInteger = [](int* int_being_bcast)
  { MPI_Bcast(int_being_bcast, 1, MPI_INT, 0, Chi::mpi.comm); };

  /** Wrapper to an MPI_Bcast call for an array of characters
   * broadcast from location 0. */
  auto HomeBroadcastStringAsRaw = [](std::string string_to_bcast, int length)
  {
    char* raw_string_to_bcast = string_to_bcast.data();
    MPI_Bcast(raw_string_to_bcast, length, MPI_CHAR, 0, Chi::mpi.comm);
  };

  /** Wrapper to an MPI_Bcast call for an array of characters
   * broadcast from location 0. This call is for non-home locations. */
  auto NonHomeBroadcastStringAsRaw = [](std::string& string_to_bcast, int length)
  {
    std::vector<char> raw_chars(length + 1, '\0');
    MPI_Bcast(raw_chars.data(), length, MPI_CHAR, 0, Chi::mpi.comm);

    string_to_bcast = std::string(raw_chars.data());
  };

#ifdef OPENSN_WITH_LUA
  /** Executes a string within the lua-console. */
  auto LuaDoString = [this](const std::string& the_string)
  {
    bool error = luaL_dostring(console_state_, the_string.c_str());
    if (error)
    {
      Chi::log.LogAll() << lua_tostring(console_state_, -1);
      lua_pop(console_state_, 1);
    }
  };
#endif

  auto ConsoleInputNumChars = [](const std::string& input)
  {
    int L = static_cast<int>(input.size());
    if (input == std::string("exit")) L = -1;

    return L;
  };

  const bool HOME = Chi::mpi.location_id == 0;

  while (not Chi::run_time::termination_posted_)
  {
    std::string console_input;

    if (HOME) std::cin >> console_input; // Home will be waiting here

    int console_input_len = ConsoleInputNumChars(console_input);

    BroadcastSingleInteger(&console_input_len); // Non-Home locs wait here

    if (console_input_len < 0) break;
    else if (HOME)
      HomeBroadcastStringAsRaw(console_input, console_input_len);
    else
      NonHomeBroadcastStringAsRaw(console_input, console_input_len);

#ifdef OPENSN_WITH_LUA
    try
    {
      LuaDoString(console_input);
    }
    catch (const Chi::RecoverableException& e)
    {
      Chi::log.LogAllError() << e.what();
    }
    catch (const std::exception& e)
    {
      Chi::log.LogAllError() << e.what();
      Chi::Exit(EXIT_FAILURE);
    }
#endif
  } // while not termination posted

  Chi::run_time::termination_posted_ = true;

  Chi::log.Log() << "Console loop stopped successfully.";
}

int
chi::Console::ExecuteFile(const std::string& fileName, int argc, char** argv) const
{
#ifdef OPENSN_WITH_LUA
  lua_State* L = this->console_state_;
  if (not fileName.empty())
  {
    if (argc > 0)
    {
      lua_newtable(L);
      for (int i = 1; i <= argc; i++)
      {
        lua_pushnumber(L, i);
        lua_pushstring(L, argv[i - 1]);
        lua_settable(L, -3);
      }
      lua_setglobal(L, "chiArgs");
    }
    int error = luaL_dofile(this->console_state_, fileName.c_str());

    if (error > 0)
    {
      Chi::log.LogAllError() << "LuaError: " << lua_tostring(this->console_state_, -1);
      return EXIT_FAILURE;
    }
  }
#endif
  return EXIT_SUCCESS;
}

void
chi::Console::PostMPIInfo(int location_id, int number_of_processes) const
{
#ifdef OPENSN_WITH_LUA
  lua_State* L = this->console_state_;

  lua_pushinteger(L, location_id);
  lua_setglobal(L, "chi_location_id");

  lua_pushinteger(L, number_of_processes);
  lua_setglobal(L, "chi_number_of_processes");
#endif
}

void
chi::Console::AddFunctionToRegistry(const std::string& name_in_lua, lua_CFunction function_ptr)
{
#ifdef OPENSN_WITH_LUA
  auto& console = GetInstance();

  // Check if the function name is already there
  if (console.lua_function_registry_.count(name_in_lua) > 0)
  {
    const auto& current_entry = console.lua_function_registry_.at(name_in_lua);

    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           ": Attempted "
                           "to register lua function \"" +
                           name_in_lua +
                           "\" but the function "
                           "is already taken by " +
                           current_entry.function_raw_name);
  }

  console.lua_function_registry_.insert(
    std::make_pair(name_in_lua, LuaFunctionRegistryEntry{function_ptr, name_in_lua}));
#endif
}

#ifdef OPENSN_WITH_LUA

char
chi::Console::AddFunctionToRegistryGlobalNamespace(const std::string& raw_name_in_lua,
                                                   lua_CFunction function_ptr)
{
  // Filter out namespace from the raw name
  const std::string name_in_lua = chi::StringUpToFirstReverse(raw_name_in_lua, "::");

  AddFunctionToRegistry(name_in_lua, function_ptr);

  return 0;
}
#endif

#ifdef OPENSN_WITH_LUA

char
chi::Console::AddFunctionToRegistryInNamespaceWithName(lua_CFunction function_ptr,
                                                       const std::string& namespace_name,
                                                       const std::string& function_name)
{
  const std::string name_in_lua = namespace_name + "::" + function_name;

  AddFunctionToRegistry(name_in_lua, function_ptr);

  return 0;
}
#endif

#ifdef OPENSN_WITH_LUA

char
chi::Console::AddLuaConstantToRegistry(const std::string& namespace_name,
                                       const std::string& constant_name,
                                       const chi_data_types::Varying& value)
{
  const std::string name_in_lua = namespace_name + "::" + constant_name;

  // Check if the constant name is already there
  auto& console = Console::GetInstance();
  if (console.lua_constants_registry_.count(name_in_lua) > 0)
  {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           ": Attempted "
                           "to register lua const  \"" +
                           name_in_lua +
                           "\" but the value "
                           "is already taken.");
  }

  console.lua_constants_registry_.insert(std::make_pair(name_in_lua, value));

  return 0;
}
#endif

chi::InputParameters
chi::Console::DefaultGetInParamsFunc()
{
  return InputParameters();
}

#ifdef OPENSN_WITH_LUA

char
chi::Console::AddWrapperToRegistryInNamespaceWithName(const std::string& namespace_name,
                                                      const std::string& name_in_lua,
                                                      WrapperGetInParamsFunc syntax_function,
                                                      WrapperCallFunc actual_function,
                                                      bool ignore_null_call_func)
{
  const std::string name =
    (namespace_name.empty()) ? name_in_lua : namespace_name + "::" + name_in_lua;

  auto& console = GetInstance();
  auto& registry = console.function_wrapper_registry_;

  ChiLogicalErrorIf(registry.count(name) > 0,
                    std::string("Attempted to register lua-function wrapper \"") + name +
                      "\" but a wrapper with the same name already exists");

  if (not syntax_function) syntax_function = DefaultGetInParamsFunc;

  if (not ignore_null_call_func)
    ChiLogicalErrorIf(not actual_function, "Problem with get_in_params_func");

  LuaFuncWrapperRegEntry reg_entry;
  reg_entry.get_in_params_func = syntax_function;
  reg_entry.call_func = actual_function;

  registry.insert(std::make_pair(name, reg_entry));

  return 0;
}
#endif

#ifdef OPENSN_WITH_LUA

void
chi::Console::SetLuaFuncNamespaceTableStructure(const std::string& full_lua_name,
                                                lua_CFunction function_ptr)
{
  auto L = GetInstance().console_state_;
  const auto lua_name_split = chi::StringSplit(full_lua_name, "::");

  if (lua_name_split.size() == 1)
  {
    lua_pushcfunction(L, function_ptr);
    lua_setglobal(L, lua_name_split.back().c_str());
    return;
  }

  const std::vector<std::string> table_names(lua_name_split.begin(), lua_name_split.end() - 1);

  FleshOutLuaTableStructure(table_names);

  lua_pushstring(L, lua_name_split.back().c_str());
  lua_pushcfunction(L, function_ptr);
  lua_settable(L, -3);

  lua_pop(L, lua_gettop(L));
}
#endif

#ifdef OPENSN_WITH_LUA

void
chi::Console::SetLuaFuncWrapperNamespaceTableStructure(const std::string& full_lua_name)
{
  auto L = GetInstance().console_state_;

  /**Lambda for making a chunk*/
  auto MakeChunk = [&L, &full_lua_name]()
  {
    std::string chunk_code = "local params = ...; ";
    chunk_code += "return chi_console.LuaWrapperCall(\"" + full_lua_name + "\", ...)";

    luaL_loadstring(L, chunk_code.c_str());
  };

  const auto table_names = chi::StringSplit(full_lua_name, "::");
  std::vector<std::string> namespace_names;
  for (const auto& table_name : table_names)
    if (table_name != table_names.back()) namespace_names.push_back(table_name);

  const auto& function_name = table_names.back();

  if (not namespace_names.empty())
  {
    FleshOutLuaTableStructure(namespace_names);
    lua_pushstring(L, function_name.c_str());
    MakeChunk();
    lua_settable(L, -3);
  }
  else
  {
    MakeChunk();
    lua_setglobal(L, function_name.c_str());
  }

  lua_pop(L, lua_gettop(L));
}
#endif

#ifdef OPENSN_WITH_LUA

void
chi::Console::SetObjectNamespaceTableStructure(const std::string& full_lua_name)
{
  auto L = GetInstance().console_state_;

  /**Lambda for registering object type and creation function.*/
  auto RegisterObjectItems = [&L](const std::string& full_name)
  {
    lua_pushstring(L, "type");
    lua_pushstring(L, full_name.c_str());
    lua_settable(L, -3);

    lua_pushstring(L, "Create");
    std::string chunk_code = "local params = ...; ";
    chunk_code += "return chiMakeObjectType(\"" + full_name + "\", ...)";

    luaL_loadstring(L, chunk_code.c_str());
    lua_settable(L, -3);
  };

  const auto table_names = chi::StringSplit(full_lua_name, "::");

  FleshOutLuaTableStructure(table_names);

  RegisterObjectItems(full_lua_name);

  lua_pop(L, lua_gettop(L));
}
#endif

#ifdef OPENSN_WITH_LUA

void
chi::Console::FleshOutLuaTableStructure(const std::vector<std::string>& table_names)
{
  auto L = GetInstance().console_state_;

  for (const auto& table_name : table_names)
  {
    // The first entry needs to be in lua's global scope,
    // so it looks a little different
    if (table_name == table_names.front())
    {
      lua_getglobal(L, table_name.c_str());
      if (not lua_istable(L, -1))
      {
        lua_pop(L, 1);
        lua_newtable(L);
        lua_setglobal(L, table_name.c_str());
        lua_getglobal(L, table_name.c_str());
      }
    }
    else
    {
      lua_getfield(L, -1, table_name.c_str());
      if (not lua_istable(L, -1))
      {
        lua_pop(L, 1);
        lua_pushstring(L, table_name.c_str());
        lua_newtable(L);
        lua_settable(L, -3);
        lua_getfield(L, -1, table_name.c_str());
      }
    }
  } // for table_key in table_keys
}
#endif

#ifdef OPENSN_WITH_LUA

void
chi::Console::SetLuaConstant(const std::string& constant_name, const chi_data_types::Varying& value)
{
  auto& console = GetInstance();
  auto L = console.console_state_;
  const auto path_names = chi::StringSplit(constant_name, "::");

  auto PushVaryingValue = [&L](const chi_data_types::Varying& var_value)
  {
    if (var_value.Type() == chi_data_types::VaryingDataType::BOOL)
      lua_pushboolean(L, var_value.BoolValue());
    else if (var_value.Type() == chi_data_types::VaryingDataType::STRING)
      lua_pushstring(L, var_value.StringValue().c_str());
    else if (var_value.Type() == chi_data_types::VaryingDataType::INTEGER)
      lua_pushinteger(L, static_cast<lua_Integer>(var_value.IntegerValue()));
    else if (var_value.Type() == chi_data_types::VaryingDataType::FLOAT)
      lua_pushnumber(L, var_value.FloatValue());
    else
      ChiInvalidArgument("Unsupported value type. Only bool, string, int and "
                         "double is supported");
  };

  if (path_names.size() == 1)
  {
    PushVaryingValue(value);
    lua_setglobal(L, path_names.front().c_str());
  }
  else
  {
    std::vector<std::string> namespace_names;
    for (const auto& table_name : path_names)
      if (table_name != path_names.back()) { namespace_names.push_back(table_name); }

    FleshOutLuaTableStructure(namespace_names);
    lua_pushstring(L, path_names.back().c_str());
    PushVaryingValue(value);
    lua_settable(L, -3);
  }

  lua_pop(L, lua_gettop(L));
}
#endif

#ifdef OPENSN_WITH_LUA

void
chi::Console::DumpRegister() const
{
  Chi::log.Log() << "\n\n";
  for (const auto& [key, entry] : function_wrapper_registry_)
  {
    if (Chi::log.GetVerbosity() == 0)
    {
      Chi::log.Log() << key;
      continue;
    }

    Chi::log.Log() << "LUA_FUNCWRAPPER_BEGIN " << key;

    if (not entry.call_func) Chi::log.Log() << "SYNTAX_BLOCK";

    const auto in_params = entry.get_in_params_func();
    in_params.DumpParameters();

    Chi::log.Log() << "LUA_FUNCWRAPPER_END\n\n";
  }
  Chi::log.Log() << "\n\n";
}
#endif

#ifdef OPENSN_WITH_LUA
void
Console::UpdateConsoleBindings(const chi::RegistryStatuses& old_statuses)
{
  auto ListHasValue = [](const std::vector<std::string>& list, const std::string& value)
  { return std::find(list.begin(), list.end(), value) != list.end(); };

  const auto& object_factory = ChiObjectFactory::GetInstance();
  for (const auto& [key, _] : object_factory.Registry())
    if (not ListHasValue(old_statuses.objfactory_keys_, key)) SetObjectNamespaceTableStructure(key);

  for (const auto& [key, entry] : lua_function_registry_)
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      SetLuaFuncNamespaceTableStructure(key, entry.function_ptr);

  for (const auto& [key, entry] : function_wrapper_registry_)
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      if (entry.call_func) SetLuaFuncWrapperNamespaceTableStructure(key);
}
#endif

} // namespace chi
