#include "console/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_utils.h"

#include "ChiObjectFactory.h"

namespace chi::lua_utils
{
int chiMakeObject(lua_State* L);
}

// #############################################################################
// Execute file
/** Executes the given file in the Lua engine.
\author Jan*/
int
chi::Console::ExecuteFile(const std::string& fileName, int argc, char** argv) const
{
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
  return EXIT_SUCCESS;
}

// ###################################################################
/**Pushes location id and number of processes to lua state.*/
void
chi::Console::PostMPIInfo(int location_id, int number_of_processes) const
{
  lua_State* L = this->console_state_;

  lua_pushinteger(L, location_id);
  lua_setglobal(L, "chi_location_id");

  lua_pushinteger(L, number_of_processes);
  lua_setglobal(L, "chi_number_of_processes");
}

// ###################################################################
/**Basic addition to registry. Used by the other public methods
 * to registry a text-key to a lua function.*/
void
chi::Console::AddFunctionToRegistry(const std::string& name_in_lua, lua_CFunction function_ptr)
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

// ###################################################################
/**Adds a lua_CFunction to the registry. The registry of functions gets
 * parsed into the lua console when `chi::Initialize` is called. This
 * particular function will strip the namespace from the the parameter
 * `raw_name_in_lua` and cause the function to be registered in the
 * global namespace of the lua console.*/
char
chi::Console::AddFunctionToRegistryGlobalNamespace(const std::string& raw_name_in_lua,
                                                   lua_CFunction function_ptr)
{
  // Filter out namespace from the raw name
  const std::string name_in_lua = chi::StringUpToFirstReverse(raw_name_in_lua, "::");

  AddFunctionToRegistry(name_in_lua, function_ptr);

  return 0;
}

// ###################################################################
/**Adds a lua_CFunction to the registry. The registry of functions gets
 * parsed into the lua console when `chi::Initialize` is called. The full
 * path of the function will be derived from `namespace_name` + "::" +
 * `function_name`.*/
char
chi::Console::AddFunctionToRegistryInNamespaceWithName(lua_CFunction function_ptr,
                                                       const std::string& namespace_name,
                                                       const std::string& function_name)
{
  const std::string name_in_lua = namespace_name + "::" + function_name;

  AddFunctionToRegistry(name_in_lua, function_ptr);

  return 0;
}

// ###################################################################
/**\brief Adds a constant to the lua state. Prepending the constant
 * within a namespace is optional.*/
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

// ###################################################################
chi::InputParameters
chi::Console::DefaultGetInParamsFunc()
{
  return InputParameters();
}

// ###################################################################
/**Wrapper functions operate with input and output parameters, essentially
 * hiding the lua interface.*/
char
chi::Console::AddWrapperToRegistryInNamespaceWithName(const std::string& namespace_name,
                                                      const std::string& name_in_lua,
                                                      WrapperGetInParamsFunc syntax_function,
                                                      WrapperCallFunc actual_function,
                                                      bool ignore_null_call_func /*=false*/)
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

// ###################################################################
/**Sets/Forms a lua function in the state using a namespace structure.*/
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

// ###################################################################
/**Sets/Forms a table structure that mimics the namespace structure of
 * a string. For example the string "sing::sob::nook::Tigger" will be
 * assigned a table structure
 * `sing.sob.nook.Tigger = "sing::sob::nook::Tigger"`. Then finally assigns
 * lua call to this table.*/
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

// ###################################################################
/**Sets/Forms a table structure that mimics the namespace structure of
 * a string. For example the string "sing::sob::nook::Tigger" will be
 * assigned a table structure
 * `sing.sob.nook.Tigger = "sing::sob::nook::Tigger"`.*/
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

// ##################################################################
/**Fleshes out a path in a table tree. For example, given
 * "fee::foo::fah::koo, this routine will make sure that
 * fee.foo.fah.koo is defined as a table tree structure. The routine will
 * create a table structure where one is needed and leave existing ones alone.
 *
 * At the end of the routine the last table in the structure will be on top
 * of the stack.*/
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

// ##################################################################
/**Sets a lua constant in the lua state.*/
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

// ##################################################################
/**Makes a formatted output, readible by the documentation scripts,
 * of all the lua wrapper functions.*/
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

// ##################################################################
/**Given an old status, will update the bindings for only newly registered
 * items.*/
void
chi::Console::UpdateConsoleBindings(const chi::RegistryStatuses& old_statuses)
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