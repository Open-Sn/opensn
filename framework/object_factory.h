// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include "framework/logging/log_exceptions.h"
#include "framework/utils/utils.h"
#include <memory>

/**
 * Macro for registering an object within the ObjectFactory singleton.
 * \param namespace_name Name of the namespace within which the object is.
 * \param object_name Name of the object in the registry.
 * Example:
 * \code
 * OpenSnRegisterObjectInNamespace(kaka, Zorba);
 * \endcode
 * \note Remember to include the header "framework/object_factory.h".*/
#define OpenSnRegisterObjectInNamespace(namespace_name, object_name)                               \
  static char OpenSnJoinWords(unique_var_name_object_##object_name##_, __COUNTER__) =              \
    opensn::ObjectFactory::AddObjectToRegistry<object_name>(#namespace_name, #object_name)

#define OpenSnRegisterObject(object_name)                                                          \
  static char OpenSnJoinWords(unique_var_name_object_##object_name##_, __COUNTER__) =              \
    opensn::ObjectFactory::AddObjectToRegistry<object_name>(#object_name)

/**
 * Macro for registering an object alias within the ObjectFactory
 *
 * \param namespace_name Namespace name
 * \param alias Name of the object
 * \param object_name C++ class name to register.
 *
 * \note This will register a C++ class `object_name` such that it will show up as
 * `namespace_name.alias`
 */
#define OpenSnRegisterObjectAliasInNamespace(namespace_name, alias, object_name)                   \
  static char OpenSnJoinWords(unique_var_name_object_##object_name##_, __COUNTER__) =              \
    opensn::ObjectFactory::AddObjectToRegistry<object_name>(#namespace_name, #object_name)

/**
 * Macro for registering an object (parameters only) within the
 * ObjectFactory singleton.
 * \param namespace_name Name of the namespace within which the object is.
 * \param object_name Name of the object in the registry.
 * Example:
 * \code
 * OpenSnRegisterObjectParametersOnlyInNamespace(kaka, Zorba);
 * \endcode
 *
 * \note Remember to include the header "framework/object_factory.h"*/
#define OpenSnRegisterObjectParametersOnlyInNamespace(namespace_name, object_name)                 \
  static char OpenSnJoinWords(unique_var_name_object_##object_name##_, __COUNTER__) =              \
    opensn::ObjectFactory::AddObjectToRegistryParamsOnly<object_name>(#namespace_name,             \
                                                                      #object_name)

#define OpenSnRegisterObjectParametersOnly(object_name)                                            \
  static char OpenSnJoinWords(unique_var_name_object_##object_name##_, __COUNTER__) =              \
    opensn::ObjectFactory::AddObjectToRegistryParamsOnly<object_name>(#object_name)

namespace opensn
{

/// Singleton object for handling the registration and making of `Object`s.
class ObjectFactory
{
public:
  using ObjectGetInParamsFunc = InputParameters (*)();

  /// Structure storing the entities necessary for creating an object
  struct ObjectRegistryEntry
  {
    ObjectGetInParamsFunc get_in_params_func = nullptr;
  };

  // Deleted copy, move constructors and copy assignment operator
  ObjectFactory(const ObjectFactory&) = delete;
  ObjectFactory(const ObjectFactory&&) = delete;
  ObjectFactory& operator=(const ObjectFactory&) = delete;

  /// Access to the singleton
  static ObjectFactory& GetInstance() noexcept;

  /// Returns a constant reference to the object registry.
  const std::map<std::string, ObjectRegistryEntry>& GetRegistry() const;

  /// Checks if the object registry has a specific text key.
  bool RegistryHasKey(const std::string& key) const;

  template <typename T>
  static char AddObjectToRegistry(const std::string& namespace_name, const std::string& object_name)
  {
    return AddObjectToRegistry<T>(namespace_name + "::" + object_name);
  }

  template <typename T>
  static char AddObjectToRegistry(const std::string& object_name)
  {
    auto& object_maker = GetInstance();

    const std::string name = object_name;
    object_maker.AssertRegistryKeyAvailable(name, __PRETTY_FUNCTION__);

    ObjectRegistryEntry reg_entry;
    reg_entry.get_in_params_func = &CallGetInputParamsFunction<T>;
    object_maker.object_registry_.insert(std::make_pair(name, reg_entry));

    return 0;
  }

  template <typename T>
  static char AddObjectToRegistryParamsOnly(const std::string& namespace_name,
                                            const std::string& object_name)
  {
    return AddObjectToRegistryParamsOnly<T>(namespace_name + "::" + object_name);
  }

  template <typename T>
  static char AddObjectToRegistryParamsOnly(const std::string& object_name)
  {
    auto& object_maker = GetInstance();

    const std::string name = object_name;
    object_maker.AssertRegistryKeyAvailable(name, __PRETTY_FUNCTION__);

    ObjectRegistryEntry reg_entry;
    reg_entry.get_in_params_func = &CallGetInputParamsFunction<T>;
    object_maker.object_registry_.insert(std::make_pair(name, reg_entry));

    return 0;
  }

  static char AddSyntaxBlockToRegistry(const std::string& namespace_name,
                                       const std::string& block_name,
                                       ObjectGetInParamsFunc syntax_function)
  {
    return AddSyntaxBlockToRegistry(namespace_name + "::" + block_name, syntax_function);
  }

  static char AddSyntaxBlockToRegistry(const std::string& block_name,
                                       ObjectGetInParamsFunc syntax_function)
  {
    auto& object_maker = GetInstance();

    const std::string name = block_name;
    object_maker.AssertRegistryKeyAvailable(name, __PRETTY_FUNCTION__);

    ObjectRegistryEntry reg_entry;
    reg_entry.get_in_params_func = syntax_function;
    object_maker.object_registry_.insert(std::make_pair(name, reg_entry));

    return 0;
  }

  template <class TYPE>
  std::shared_ptr<TYPE> Create(const std::string& type, const ParameterBlock& params) const
  {
    if (object_registry_.count(type) == 0)
      throw std::logic_error("No registered type \"" + type + "\" found.");

    auto object_entry = object_registry_.at(type);
    if (not object_entry.get_in_params_func)
      throw std::runtime_error(
        "Object is not constructable since it has no registered constructor");

    auto input_params = object_entry.get_in_params_func();
    input_params.SetObjectType(type);
    input_params.SetErrorOriginScope(type);
    input_params.AssignParameters(params);
    auto obj = std::make_shared<TYPE>(input_params);
    return obj;
  }

  /// Returns the input parameters of a registered object.
  InputParameters GetRegisteredObjectParameters(const std::string& type) const;

  /// Dumps the object registry to stdout.
  void DumpRegister() const;

private:
  std::map<std::string, ObjectRegistryEntry> object_registry_;

  /// Private constructor because this is a singleton.
  ObjectFactory() = default;

  /// Utility redirection to call an object's static `GetInputParameters` function.
  template <typename T>
  static InputParameters CallGetInputParamsFunction()
  {
    return T::GetInputParameters();
  }

  /// Utility redirection to call an object's constructor with a specified list of input parameters.
  template <typename T, typename base_T>
  static std::shared_ptr<base_T> CallObjectConstructor(const InputParameters& params)
  {
    static_assert(std::is_base_of_v<base_T, T>, "Is not a base");
    return std::make_shared<T>(params);
  }

  /// Checks that the registry key is available and throws a `std::logical_error` if it is not.
  void AssertRegistryKeyAvailable(const std::string& key,
                                  const std::string& calling_function) const;
};

} // namespace opensn
