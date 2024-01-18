#pragma once

#include "framework/object.h"

namespace opensnlua
{

/**Stores all the keys currently in the registries.*/
struct RegistryStatuses
{
  std::vector<std::string> objfactory_keys_;
  std::vector<std::string> console_lua_func_keys_;
  std::vector<std::string> console_lua_wrapper_keys_;
};

/**Builds a `RegistryStatuses` structure*/
RegistryStatuses GetStatusOfRegistries();

class Plugin : public opensn::Object
{
public:
  static opensn::InputParameters GetInputParameters();
  explicit Plugin(const opensn::InputParameters& params);

  ~Plugin();

protected:
  const std::string plugin_path_;
  void* library_handle_ = nullptr;
};

} // namespace opensnlua
