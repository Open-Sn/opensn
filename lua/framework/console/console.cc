// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/console/console.h"
#include "lua/modules/modules.h"
#include "lua/framework/lua.h"
#include "lua/framework/interfaces/plugin.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/logging/log_exceptions.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include "config.h"
#include <iostream>

using namespace opensn;

namespace opensnlua
{

Console& console = Console::GetInstance();

RegisterLuaFunctionInNamespace(Console::LuaWrapperCall, console, LuaWrapperCall);

Console&
Console::GetInstance() noexcept
{
  static Console singleton;
  return singleton;
}

Console::Console() noexcept : console_state_(luaL_newstate())
{
}

void
Console::FlushConsole()
{
  try
  {
    for (auto& command : command_buffer_)
    {
      bool error = luaL_dostring(console_state_, command.c_str());
      if (error)
      {
        opensn::log.LogAll() << lua_tostring(console_state_, -1);
        lua_pop(console_state_, 1);
      }
    }
  }
  catch (const std::exception& e)
  {
    opensn::log.LogAllError() << e.what();
    opensn::Exit(EXIT_FAILURE);
  }
}

int
Console::LuaWrapperCall(lua_State* L)
{
  const int num_args = lua_gettop(L);
  // We do not check for the required parameters here because we want
  // to make this function call as fast as possible. Besides, via the
  // static registration we should never run into an issue here.

  auto& console = Console::GetInstance();

  const auto& registry = console.function_wrapper_registry_;

  const auto fname = LuaArg<std::string>(L, 1);

  OpenSnLogicalErrorIf(registry.count(fname) == 0,
                       std::string("Wrapper with name \"") + fname + "\" not in console registry.");

  const auto& reg_entry = registry.at(fname);

  auto input_params = reg_entry.get_in_params_func();

  ParameterBlock main_arguments_block;
  for (int p = 2; p <= num_args; ++p)
  {
    const std::string arg_name = "arg" + std::to_string(p - 2);

    if (lua_isboolean(L, p))
      main_arguments_block.AddParameter(arg_name, LuaArg<bool>(L, p));
    else if (lua_isinteger(L, p))
      main_arguments_block.AddParameter(arg_name, LuaArg<int>(L, p));
    else if (lua_isnumber(L, p))
      main_arguments_block.AddParameter(arg_name, LuaArg<double>(L, p));
    else if (lua_isstring(L, p))
      main_arguments_block.AddParameter(arg_name, LuaArg<std::string>(L, p));
    else if (lua_istable(L, p))
    {
      auto block = LuaArg<ParameterBlock>(L, p);
      block.SetBlockName(arg_name);
      std::string scope = fname + ":";
      scope.append(arg_name + " ");
      block.SetErrorOriginScope(scope);
      main_arguments_block.AddParameter(block);
    }
    else
      OpenSnInvalidArgument("In call to \"" + fname +
                            "\": Unsupported argument "
                            "type \"" +
                            lua_typename(L, lua_type(L, p)) + "\" encountered.");
  }
  // Set input parameters here
  input_params.SetErrorOriginScope(fname + "()");
  input_params.AssignParameters(main_arguments_block);

  auto output_params = reg_entry.call_func(input_params);

  output_params.SetErrorOriginScope(fname + ":output:");
  PushParameterBlock(L, output_params);

  const int num_sub_params = static_cast<int>(output_params.NumParameters());

  return output_params.IsScalar() ? 1 : num_sub_params;
}

void
Console::RunConsoleLoop(char*) const
{
  opensn::log.Log() << "Console loop started. "
                    << "Type \"exit\" to quit (or Ctl-C).";

  /** Executes a string within the lua-console. */
  auto LuaDoString = [this](const std::string& the_string)
  {
    bool error = luaL_dostring(console_state_, the_string.c_str());
    if (error)
    {
      opensn::log.LogAll() << lua_tostring(console_state_, -1);
      lua_pop(console_state_, 1);
    }
  };

  const bool HOME = opensn::mpi_comm.rank() == 0;

  while (true)
  {
    std::string console_input;

    if (HOME)
      std::cin >> console_input; // Home will be waiting here

    mpi_comm.broadcast(console_input, 0);
    if (console_input == "exit")
      break;

    try
    {
      LuaDoString(console_input);
    }
    catch (const opensn::RecoverableException& e)
    {
      opensn::log.LogAllError() << e.what();
    }
    catch (const std::exception& e)
    {
      opensn::log.LogAllError() << e.what();
      Exit(EXIT_FAILURE);
    }
  } // while not termination posted

  opensn::log.Log() << "Console loop stopped successfully.";
}

int
Console::ExecuteFile(const std::string& fileName, int argc, char** argv) const
{
  lua_State* L = this->console_state_;
  if (not fileName.empty())
  {
    if (argc > 0)
    {
      std::vector<std::string> args;
      args.resize(argc);
      for (int i = 0; i < argc; ++i)
        args[i] = std::string(argv[i]);
      LuaSetGlobal(L, "Args", args);
    }
    int error = luaL_dofile(this->console_state_, fileName.c_str());

    if (error > 0)
    {
      opensn::log.LogAllError() << "LuaError: " << lua_tostring(this->console_state_, -1);
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

void
Console::PostMPIInfo(int location_id, int number_of_processes) const
{
  lua_State* L = this->console_state_;

  LuaSetGlobal(L, "location_id", location_id);
  LuaSetGlobal(L, "number_of_processes", number_of_processes);
}

void
Console::AddFunctionToRegistry(const std::string& name_in_lua, lua_CFunction function_ptr)
{
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
}

char
Console::AddFunctionToRegistryGlobalNamespace(const std::string& raw_name_in_lua,
                                              lua_CFunction function_ptr)
{
  // Filter out namespace from the raw name
  const std::string name_in_lua = StringUpToFirstReverse(raw_name_in_lua, "::");

  AddFunctionToRegistry(name_in_lua, function_ptr);

  return 0;
}

char
Console::AddFunctionToRegistryInNamespaceWithName(lua_CFunction function_ptr,
                                                  const std::string& namespace_name,
                                                  const std::string& function_name)
{
  const std::string name_in_lua = namespace_name + "::" + function_name;

  AddFunctionToRegistry(name_in_lua, function_ptr);

  return 0;
}

char
Console::AddLuaConstantToRegistry(const std::string& namespace_name,
                                  const std::string& constant_name,
                                  const Varying& value)
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

InputParameters
Console::DefaultGetInParamsFunc()
{
  return InputParameters();
}

char
Console::AddWrapperToRegistryInNamespaceWithName(const std::string& namespace_name,
                                                 const std::string& name_in_lua,
                                                 WrapperGetInParamsFunc syntax_function,
                                                 WrapperCallFunc actual_function,
                                                 bool ignore_null_call_func)
{
  return AddWrapperToRegistryInNamespaceWithName(
    namespace_name + "::" + name_in_lua, syntax_function, actual_function, ignore_null_call_func);
}

char
Console::AddWrapperToRegistryInNamespaceWithName(const std::string& name_in_lua,
                                                 WrapperGetInParamsFunc syntax_function,
                                                 WrapperCallFunc actual_function,
                                                 bool ignore_null_call_func)
{
  auto& console = GetInstance();
  auto& registry = console.function_wrapper_registry_;

  OpenSnLogicalErrorIf(registry.count(name_in_lua) > 0,
                       std::string("Attempted to register lua-function wrapper \"") + name_in_lua +
                         "\" but a wrapper with the same name already exists");

  if (not syntax_function)
    syntax_function = DefaultGetInParamsFunc;

  if (not ignore_null_call_func)
    OpenSnLogicalErrorIf(not actual_function, "Problem with get_in_params_func");

  LuaFuncWrapperRegEntry reg_entry;
  reg_entry.get_in_params_func = syntax_function;
  reg_entry.call_func = actual_function;

  registry.insert(std::make_pair(name_in_lua, reg_entry));

  return 0;
}

void
Console::SetLuaFuncNamespaceTableStructure(const std::string& full_lua_name,
                                           lua_CFunction function_ptr)
{
  auto L = GetInstance().console_state_;
  const auto lua_name_split = StringSplit(full_lua_name, "::");

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

void
Console::SetLuaFuncWrapperNamespaceTableStructure(const std::string& full_lua_name)
{
  auto L = GetInstance().console_state_;

  /**Lambda for making a chunk*/
  auto MakeChunk = [&L, &full_lua_name]()
  {
    std::string chunk_code = "local params = ...; ";
    chunk_code += "return console.LuaWrapperCall(\"" + full_lua_name + "\", ...)";

    luaL_loadstring(L, chunk_code.c_str());
  };

  const auto table_names = StringSplit(full_lua_name, "::");
  std::vector<std::string> namespace_names;
  for (const auto& table_name : table_names)
    if (table_name != table_names.back())
      namespace_names.push_back(table_name);

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

void
Console::SetObjectNamespaceTableStructure(const std::string& full_lua_name)
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
    chunk_code += "return MakeObjectType(\"" + full_name + "\", ...)";

    luaL_loadstring(L, chunk_code.c_str());
    lua_settable(L, -3);
  };

  const auto table_names = StringSplit(full_lua_name, "::");

  FleshOutLuaTableStructure(table_names);

  RegisterObjectItems(full_lua_name);

  lua_pop(L, lua_gettop(L));
}

void
Console::FleshOutLuaTableStructure(const std::vector<std::string>& table_names)
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

void
Console::SetLuaConstant(const std::string& constant_name, const Varying& value)
{
  auto& console = GetInstance();
  auto L = console.console_state_;
  const auto path_names = StringSplit(constant_name, "::");

  auto PushVaryingValue = [&L](const Varying& var_value)
  {
    if (var_value.Type() == VaryingDataType::BOOL)
      lua_pushboolean(L, var_value.BoolValue());
    else if (var_value.Type() == VaryingDataType::STRING)
      lua_pushstring(L, var_value.StringValue().c_str());
    else if (var_value.Type() == VaryingDataType::INTEGER)
      lua_pushinteger(L, static_cast<lua_Integer>(var_value.IntegerValue()));
    else if (var_value.Type() == VaryingDataType::FLOAT)
      lua_pushnumber(L, var_value.FloatValue());
    else
      OpenSnInvalidArgument("Unsupported value type. Only bool, string, int and "
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
      if (table_name != path_names.back())
      {
        namespace_names.push_back(table_name);
      }

    FleshOutLuaTableStructure(namespace_names);
    lua_pushstring(L, path_names.back().c_str());
    PushVaryingValue(value);
    lua_settable(L, -3);
  }

  lua_pop(L, lua_gettop(L));
}

void
Console::DumpRegister() const
{
  opensn::log.Log() << "\n\n";
  for (const auto& [key, entry] : function_wrapper_registry_)
  {
    if (opensn::log.GetVerbosity() == 0)
    {
      opensn::log.Log() << key;
      continue;
    }

    opensn::log.Log() << "LUA_FUNCWRAPPER_BEGIN " << key;

    if (not entry.call_func)
      opensn::log.Log() << "SYNTAX_BLOCK";

    const auto in_params = entry.get_in_params_func();
    in_params.DumpParameters();

    opensn::log.Log() << "LUA_FUNCWRAPPER_END\n\n";
  }
  opensn::log.Log() << "\n\n";
}

void
Console::UpdateConsoleBindings(const RegistryStatuses& old_statuses)
{
  auto ListHasValue = [](const std::vector<std::string>& list, const std::string& value)
  { return std::find(list.begin(), list.end(), value) != list.end(); };

  const auto& object_factory = ObjectFactory::GetInstance();
  for (const auto& [key, _] : object_factory.Registry())
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      SetObjectNamespaceTableStructure(key);

  for (const auto& [key, entry] : lua_function_registry_)
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      SetLuaFuncNamespaceTableStructure(key, entry.function_ptr);

  for (const auto& [key, entry] : function_wrapper_registry_)
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      if (entry.call_func)
        SetLuaFuncWrapperNamespaceTableStructure(key);
}

} // namespace opensnlua
