// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/materials/material.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "framework/materials/material.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/scalar_value.h"
#include "framework/materials/isotropic_multigroup_source.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <iostream>

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(MatAddMaterial, mat, AddMaterial);
RegisterLuaFunctionInNamespace(MatSetProperty, mat, SetProperty);
RegisterLuaFunctionInNamespace(MatGetProperty, mat, GetProperty);

RegisterLuaConstant(SCALAR_VALUE, Varying(1));
RegisterLuaConstant(TRANSPORT_XSECTIONS, Varying(10));
RegisterLuaConstant(ISOTROPIC_MG_SOURCE, Varying(11));

namespace
{

void
ScalarPropertyPushTable(lua_State* L, std::shared_ptr<MaterialProperty> property)
{
  lua_newtable(L);
  LuaPushTableKey(L, "is_empty", false);
  LuaPushTableKey(L, "value", property->GetScalarValue());
}

void
IsotropicMGSourcePropertyPushTable(lua_State* L,
                                   std::shared_ptr<IsotropicMultiGroupSource> property)
{
  lua_newtable(L);
  LuaPushTableKey(L, "is_empty", false);
  LuaPushTableKey(L, "G", property->source_value_g.size());
  LuaPushTableKey(L, "source_value_g", property->source_value_g);
}

void
MaterialPropertyPushLuaTable(lua_State* L)
{
  lua_newtable(L);
  LuaPushTableKey(L, "is_empty", true);
}

void
PropertyPushLuaTable(lua_State* L, std::shared_ptr<MaterialProperty> property)
{
  if (property->Type() == PropertyType::SCALAR_VALUE)
    ScalarPropertyPushTable(L, property);
  else if (property->Type() == PropertyType::ISOTROPIC_MG_SOURCE)
    IsotropicMGSourcePropertyPushTable(
      L, std::dynamic_pointer_cast<IsotropicMultiGroupSource>(property));
  else
    MaterialPropertyPushLuaTable(L);
}

} // namespace

int
MatAddMaterial(lua_State* L)
{
  auto new_material = std::make_shared<Material>();
  new_material->name = LuaArgOptional<std::string>(L, 1, "");

  opensn::material_stack.push_back(new_material);

  const size_t index = opensn::material_stack.size() - 1;

  opensn::log.Log0Verbose1() << "New material added at index " << index << " with name \""
                             << new_material->name << "\"";

  return LuaReturn(L, index);
}

int
MatSetProperty(lua_State* L)
{
  const std::string fname = "mat.SetProperty";
  const int num_args = lua_gettop(L);

  if (num_args < 3)
    LuaPostArgAmountError(fname, L, 3, num_args);

  const auto material_handle = LuaArg<int>(L, 1);
  const auto property_type = LuaArg<int>(L, 2);
  const auto op_type = LuaArg<int>(L, 3);

  // Get a pointer to the material
  auto material = opensn::StackItemPtr(opensn::material_stack, material_handle, fname);

  // Find the index of the material property on the material
  int property_index = -1;
  for (int p = 0; p < material->properties.size(); ++p)
  {
    const auto& property = material->properties.at(p);
    if (static_cast<int>(property->Type()) == property_type)
    {
      property_index = p;
      break;
    }
  }

  // Process property
  if (property_type == static_cast<int>(PropertyType::SCALAR_VALUE))
  {
    // Create the property, if no location was found
    if (property_index == -1)
    {
      auto property = std::make_shared<ScalarValue>();
      material->properties.push_back(property);
      property_index = static_cast<int>(material->properties.size()) - 1;
    }
    OpenSnLogicalErrorIf(property_index < 0, "Error creating or finding ScalarValue property.");

    // Get the property
    auto property = std::static_pointer_cast<ScalarValue>(material->properties.at(property_index));

    // Process operation
    if (op_type == static_cast<int>(OperationType::SINGLE_VALUE))
    {
      const auto value = LuaArg<double>(L, 4);
      property->Set(value);
      opensn::log.Log0Verbose1() << "Scalar value for material at index " << material_handle
                                 << " set to " << value;
    }
    else
      OpenSnLogicalError("Unsupported operation for ScalarValue property.");
  } // if scalar value

  else if (property_type == static_cast<int>(PropertyType::TRANSPORT_XSECTIONS))
  {
    // Create the property, if no location was found
    if (property_index == -1)
    {
      auto property = std::make_shared<MultiGroupXS>();
      material->properties.push_back(property);
      property_index = static_cast<int>(material->properties.size()) - 1;
    }
    OpenSnLogicalErrorIf(property_index < 0, "Error creating or finding MultiGroupXS property.");

    // Get the property
    auto property = std::static_pointer_cast<MultiGroupXS>(material->properties.at(property_index));

    // Process operation
    if (op_type == static_cast<int>(OperationType::SIMPLE_ONE_GROUP))
    {
      if (num_args != 5)
        LuaPostArgAmountError(fname, L, 5, num_args);

      const auto sigma_t = LuaArg<double>(L, 4);
      const auto c = LuaArg<double>(L, 5);
      property->Initialize(sigma_t, c);
      opensn::log.Log0Verbose1() << "Simple one group cross sections with sigma_t=" << sigma_t
                                 << ", c=" << c << " set on the material at index "
                                 << material_handle << ".";
    }

    else if (op_type == static_cast<int>(OperationType::OPENSN_XSFILE))
    {
      if (num_args != 4)
        LuaPostArgAmountError(fname, L, 4, num_args);

      const auto file_name = LuaArg<std::string>(L, 4);
      property->Initialize(file_name);
      opensn::log.Log0Verbose1() << "Cross sections from OpenSn XS file " << file_name
                                 << " set on the material at index " << material_handle << ".";
    }

    else if (op_type == static_cast<int>(OperationType::OPENMC_XSLIB))
    {
      if (num_args < 5)
        LuaPostArgAmountError(fname, L, 5, num_args);

      const auto file_name = LuaArg<std::string>(L, 4);
      const auto temperature = LuaArg<double>(L, 5);
      const auto dataset_name = LuaArgOptional<std::string>(L, 6, "set1");
      property->Initialize(file_name, dataset_name, temperature);
      opensn::log.Log0Verbose1() << "Cross sections from OpenMC library " << file_name
                                 << ", dataset " << dataset_name << ", and temperature "
                                 << temperature << " set on the material at index "
                                 << material_handle << ".";
    }

    else if (op_type == static_cast<int>(OperationType::EXISTING))
    {
      if (num_args != 4)
        LuaPostArgAmountError(fname, L, 4, num_args);

      const auto xs_handle = LuaArg<int>(L, 4);
      material->properties.at(property_index) = std::dynamic_pointer_cast<MultiGroupXS>(
        StackItemPtr(multigroup_xs_stack, xs_handle, fname));

      opensn::log.Log0Verbose1() << "Cross sections at index " << xs_handle
                                 << " set on material at index " << material_handle << ".";
    }

    else
      OpenSnLogicalError("Unsupported operation for MultiGroupXS property.");
  } // if transport xsections

  else if (property_type == static_cast<int>(PropertyType::ISOTROPIC_MG_SOURCE))
  {
    // Create the property, if no location was found
    if (property_index == -1)
    {
      auto property = std::make_shared<IsotropicMultiGroupSource>();
      material->properties.push_back(property);
      property_index = static_cast<int>(material->properties.size()) - 1;
    }
    OpenSnLogicalErrorIf(property_index < 0,
                         "Error creating or finding IsotropicMultiGrpSource property.");

    // Get the property
    auto property =
      std::static_pointer_cast<IsotropicMultiGroupSource>(material->properties.at(property_index));

    // Process operation
    if (op_type == static_cast<int>(OperationType::SINGLE_VALUE))
    {
      if (num_args != 4)
        LuaPostArgAmountError(fname, L, 4, num_args);

      const auto value = LuaArg<double>(L, 4);
      property->source_value_g.resize(1, value);
      opensn::log.Log0Verbose1() << "Isotropic multi-group source value for material at index "
                                 << material_handle << " set to " << value << ".";
    }
    else if (op_type == static_cast<int>(OperationType::FROM_ARRAY))
    {
      if (num_args != 4)
        LuaPostArgAmountError("MatSetProperty", L, 4, num_args);

      property->source_value_g = LuaArg<std::vector<double>>(L, 4);
      opensn::log.Log0Verbose1() << "Isotropic Multi-group Source populated with "
                                 << property->source_value_g.size() << " values.";
    }
    else
      OpenSnInvalidArgument("Unsupported operation for IsotropicMultiGrpSource");
  } // if isotropic multi-group source

  else
    OpenSnInvalidArgument("Unsupported material property type specified in call to " + fname + ".");

  return LuaReturn(L);
}

int
MatGetProperty(lua_State* L)
{
  const std::string fname = "mat.GetProperty";
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("MatGetProperty", L, 2, num_args);

  auto material_index = LuaArg<int>(L, 1);
  int property_index = -1;
  std::string property_index_name;
  if (lua_isnumber(L, 2))
    property_index = LuaArg<int>(L, 2);
  else
    property_index_name = LuaArg<std::string>(L, 2);

  // Get reference to material
  auto cur_material = opensn::StackItemPtr(opensn::material_stack, material_index, fname);

  // If user supplied name then find property index
  if (not lua_isnumber(L, 2))
  {
    for (auto& property : cur_material->properties)
      if (property->property_name == property_index_name)
        property_index = static_cast<int>(property->Type());
  }

  // Process property
  bool property_populated = false;
  for (auto& property : cur_material->properties)
  {
    if (static_cast<int>(property->Type()) == property_index)
    {
      PropertyPushLuaTable(L, property);
      property_populated = true;
    }
  }

  if (not property_populated)
  {
    opensn::log.LogAllError() << "Invalid material property specified in call to MatGetProperty."
                              << property_index << std::endl;
    opensn::Exit(EXIT_FAILURE);
  }

  return 1;
}

} // namespace opensnlua
