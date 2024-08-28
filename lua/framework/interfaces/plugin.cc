// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/interfaces/plugin.h"
#include "lua/framework/console/console.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"
#include <dlfcn.h>

using namespace opensn;

namespace opensnlua
{

RegistryStatuses
GetStatusOfRegistries()
{
  RegistryStatuses stats;

  const auto& object_factory = ObjectFactory::GetInstance();
  for (const auto& [key, _] : object_factory.Registry())
    stats.objfactory_keys_.push_back(key);

  for (const auto& [key, _] : console.GetLuaFunctionRegistry())
    stats.console_lua_func_keys_.push_back(key);

  for (const auto& [key, _] : console.GetFunctionWrapperRegistry())
    stats.console_lua_wrapper_keys_.push_back(key);

  return stats;
}

OpenSnRegisterObject(Plugin);

InputParameters
Plugin::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("Object to handle the loading of shared libraries as plug-ins");
  params.SetDocGroup("DocInterfaces");

  params.AddRequiredParameter<std::string>("plugin_path",
                                           "Path to the shared library containing the plug-in.");
  params.AddOptionalParameter("entry_function", "", "Entry function to call.");

  return params;
}

Plugin::Plugin(const InputParameters& params)
  : Object(params), plugin_path_(params.GetParamValue<std::string>("plugin_path"))
{
  opensn::log.Log0Verbose1() << "Loading plugin \"" << plugin_path_ << "\"";
  RegistryStatuses registry_statuses = GetStatusOfRegistries();

  AssertReadableFile(plugin_path_);
  library_handle_ = dlopen(plugin_path_.c_str(), RTLD_LAZY);

  OpenSnLogicalErrorIf(not library_handle_, "Failure loading \"" + plugin_path_ + "\"");

  const auto& user_params = params.ParametersAtAssignment();
  if (user_params.Has("entry_function"))
  {
    const auto& entry_function = user_params.GetParamValue<std::string>("entry_function");

    using func_ptr = void (*)();
    auto func = (func_ptr)dlsym(library_handle_, entry_function.c_str());

    OpenSnLogicalErrorIf(not func, "Failed to call entry function \"" + entry_function + "\"");

    // Calling the function
    func();
  }

  console.UpdateConsoleBindings(registry_statuses);
}

Plugin::~Plugin()
{
  dlclose(library_handle_);
}

} // namespace opensnlua
