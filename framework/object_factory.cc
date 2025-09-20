// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

ObjectFactory&
ObjectFactory::GetInstance() noexcept
{
  static ObjectFactory singleton;
  return singleton;
}

const std::map<std::string, ObjectFactory::ObjectRegistryEntry>&
ObjectFactory::GetRegistry() const
{
  return object_registry_;
}

bool
ObjectFactory::RegistryHasKey(const std::string& key) const
{
  return object_registry_.count(key) > 0;
}

InputParameters
ObjectFactory::GetRegisteredObjectParameters(const std::string& type) const
{
  auto iter = object_registry_.find(type);
  OpenSnInvalidArgumentIf(iter == object_registry_.end(),
                          "Object type \"" + type + "\" is not registered in ObjectFactory.");

  const auto& reg_entry = iter->second;

  return reg_entry.get_in_params_func();
}

void
ObjectFactory::DumpRegister() const
{
  log.Log() << "\n\n";
  for (const auto& [key, entry] : object_registry_)
  {
    if (log.GetVerbosity() == 0)
    {
      log.Log() << key;
      continue;
    }

    log.Log() << "OBJECT_BEGIN " << key;

    if (entry.get_in_params_func == nullptr)
      log.Log() << "NOT_CONSTRUCTIBLE";
    else
    {
      const auto in_params = entry.get_in_params_func();
      in_params.DumpParameters();
    }

    log.Log() << "OBJECT_END\n\n";
  }
  log.Log() << "\n\n";
}

void
ObjectFactory::AssertRegistryKeyAvailable(const std::string& key,
                                          const std::string& calling_function) const
{
  if (RegistryHasKey(key))
    OpenSnLogicalError(calling_function + ": Attempted to register Object \"" + key +
                       "\" but an object with the same name is already registered.");
}

} // namespace opensn
