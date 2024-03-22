#include "lua_material.h"
#include "framework/lua.h"
#include "framework/physics/physics_material/physics_material.h"
#include "framework/physics/physics_material/material_property_scalar_value.h"
#include "framework/physics/physics_material/multi_group_xs/single_state_mgxs.h"
#include "framework/physics/physics_material/material_property_isotropic_mg_src.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include <iostream>

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(MatAddMaterial, mat, AddMaterial);
RegisterLuaFunctionNamespace(MatAddProperty, mat, AddProperty);
RegisterLuaFunctionNamespace(MatSetProperty, mat, SetProperty);
RegisterLuaFunctionNamespace(MatGetProperty, mat, GetProperty);

RegisterLuaConstantAsIs(SCALAR_VALUE, Varying(1));
RegisterLuaConstantAsIs(TRANSPORT_XSECTIONS, Varying(10));
RegisterLuaConstantAsIs(ISOTROPIC_MG_SOURCE, Varying(11));

namespace
{

void
ScalarPropertyPushTable(lua_State* L, std::shared_ptr<PhysicsMaterialProperty> property)
{
  lua_newtable(L);
  LuaPushTableKey(L, "is_empty", false);
  LuaPushTableKey(L, "value", property->GetScalarValue());
}

void
IsotropicMGSourcePropertyPushTable(lua_State* L, std::shared_ptr<IsotropicMultiGrpSource> property)
{
  lua_newtable(L);
  LuaPushTableKey(L, "is_empty", false);
  LuaPushTableKey(L, "G", property->source_value_g_.size());
  LuaPushTableKey(L, "source_value_g", property->source_value_g_);
}

void
MaterialPropertyPushLuaTable(lua_State* L)
{
  lua_newtable(L);
  LuaPushTableKey(L, "is_empty", true);
}

void
PropertyPushLuaTable(lua_State* L, std::shared_ptr<PhysicsMaterialProperty> property)
{
  if (property->Type() == PropertyType::SCALAR_VALUE)
    ScalarPropertyPushTable(L, property);
  else if (property->Type() == PropertyType::ISOTROPIC_MG_SOURCE)
    IsotropicMGSourcePropertyPushTable(
      L, std::dynamic_pointer_cast<IsotropicMultiGrpSource>(property));
  else
    MaterialPropertyPushLuaTable(L);
}

} // namespace

int
MatAddMaterial(lua_State* L)
{
  auto new_material = std::make_shared<Material>();
  new_material->name_ = LuaArgOptional<std::string>(L, 1, "");

  opensn::material_stack.push_back(new_material);

  const size_t index = opensn::material_stack.size() - 1;
  LuaPush(L, static_cast<lua_Integer>(index));

  opensn::log.Log0Verbose1() << "New material added at index " << index << " with name \""
                             << new_material->name_ << "\"";

  return 1;
}

int
MatAddProperty(lua_State* L)
{
  const std::string fname = "mat.AddProperty";
  LuaCheckArgs<int, int>(L, fname);

  auto material_index = LuaArg<int>(L, 1);
  auto property_index = LuaArg<int>(L, 2);

  // Get reference to material
  auto cur_material = opensn::GetStackItemPtr(opensn::material_stack, material_index, fname);

  auto provided_name = LuaArgOptional<std::string>(
    L, 3, std::string("Property ") + std::to_string(cur_material->properties_.size()));

  // Process property
  using MatProperty = PropertyType;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALAR_VALUE
  if (property_index == static_cast<int>(MatProperty::SCALAR_VALUE))
  {
    // Duplicates are allowed

    auto prop = std::make_shared<ScalarValue>();

    prop->property_name = provided_name;

    cur_material->properties_.push_back(prop);
    opensn::log.Log0Verbose1() << "Scalar Value Property added to material at index "
                               << material_index;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSPORT_XSECTIONS
  else if (property_index == static_cast<int>(MatProperty::TRANSPORT_XSECTIONS))
  {
    // Check for duplicate
    for (int p = 0; p < cur_material->properties_.size(); p++)
    {
      if (cur_material->properties_[p]->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        opensn::log.Log0Error() << "Material " << material_index << " \"" << cur_material->name_
                                << "\""
                                << " already has property TRANSPORT_XSECTIONS" << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }

    auto prop = std::make_shared<SingleStateMGXS>();

    prop->property_name = provided_name;

    cur_material->properties_.push_back(prop);
    opensn::log.Log0Verbose1() << "Transport cross-sections added to material at index "
                               << material_index;

    opensn::multigroup_xs_stack.push_back(prop);

    const size_t index = opensn::multigroup_xs_stack.size() - 1;

    LuaPush(L, index);
    return 1;
  }
  else if (property_index == static_cast<int>(MatProperty::ISOTROPIC_MG_SOURCE))
  {
    // Check for duplicate
    for (int p = 0; p < cur_material->properties_.size(); p++)
    {
      if (cur_material->properties_[p]->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        opensn::log.Log0Error() << "Material " << material_index << " \"" << cur_material->name_
                                << "\""
                                << " already has property ISOTROPIC_MG_SOURCE " << property_index
                                << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }

    auto prop = std::make_shared<IsotropicMultiGrpSource>();

    prop->property_name = provided_name;

    cur_material->properties_.push_back(prop);
    opensn::log.Log0Verbose1() << "Isotropic Multigroup Source added to material at index "
                               << material_index;
  }
  else
  {
    opensn::log.Log0Error() << "Unsupported property type in call to mat.AddProperty.";
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}

int
MatSetProperty(lua_State* L)
{
  const std::string fname = "mat.SetProperty";
  const int num_args = lua_gettop(L);

  if (num_args < 3)
  {
    opensn::log.Log0Error() << "Invalid number of arguments when calling mat.SetProperty";
    opensn::Exit(EXIT_FAILURE);
  }

  auto material_index = LuaArg<int>(L, 1);
  int property_index = -1;
  std::string property_index_name;
  if (lua_isnumber(L, 2))
    property_index = LuaArg<int>(L, 2);
  else
    property_index_name = LuaArg<std::string>(L, 2);

  auto operation_index = LuaArg<int>(L, 3);

  // Get reference to material
  auto cur_material = opensn::GetStackItemPtr(opensn::material_stack, material_index, fname);

  // If user supplied name then find property index
  if (not lua_isnumber(L, 2))
  {
    for (auto& property : cur_material->properties_)
      if (property->property_name == property_index_name)
        property_index = static_cast<int>(property->Type());
  }

  // Process property
  using MatProperty = PropertyType;
  using OpType = OperationType;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALAR_VALUE
  if (property_index == static_cast<int>(MatProperty::SCALAR_VALUE))
  {
    int location_of_prop = -1;
    // Check if the material has this property
    if (lua_isnumber(L, 2))
    {
      for (int p = 0; p < cur_material->properties_.size(); p++)
        if (cur_material->properties_[p]->Type() == MatProperty::SCALAR_VALUE)
          location_of_prop = p;
    }
    else
    {
      for (int p = 0; p < cur_material->properties_.size(); p++)
        if (cur_material->properties_[p]->property_name == property_index_name)
          location_of_prop = p;
    }

    // If the property is valid
    if (location_of_prop >= 0)
    {
      auto prop =
        std::static_pointer_cast<ScalarValue>(cur_material->properties_[location_of_prop]);

      // Process operation
      if (operation_index == static_cast<int>(OpType::SINGLE_VALUE))
      {
        auto value = LuaArg<double>(L, 4);
        prop->value_ = value;
        opensn::log.Log0Verbose1()
          << "Scalar value for material at index " << material_index << " set to " << value;
      }
      else
      {
        opensn::log.Log0Error() << "ERROR: Unsupported operation for SCALAR_VALUE." << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }
    else
    {
      opensn::log.Log0Error() << "ERROR: Material has no property SCALAR_VALUE." << std::endl;
      opensn::Exit(EXIT_FAILURE);
    }
  } // if scalar value
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSPORT_XSECTIONS
  else if (property_index == static_cast<int>(MatProperty::TRANSPORT_XSECTIONS))
  {
    int location_of_prop = -1;
    // Check if the material has this property
    if (lua_isnumber(L, 2))
    {
      for (int p = 0; p < cur_material->properties_.size(); p++)
      {
        if (cur_material->properties_[p]->Type() == MatProperty::TRANSPORT_XSECTIONS)
        {
          location_of_prop = p;
        }
      }
    }
    else
    {
      for (int p = 0; p < cur_material->properties_.size(); p++)
      {
        if (cur_material->properties_[p]->property_name == property_index_name)
        {
          location_of_prop = p;
        }
      }
    }

    // If the property is valid
    if (location_of_prop >= 0)
    {
      auto prop =
        std::static_pointer_cast<SingleStateMGXS>(cur_material->properties_[location_of_prop]);

      // Process operation
      if (operation_index == static_cast<int>(OpType::SIMPLEXS0))
      {
        if (num_args != 5)
          LuaPostArgAmountError("MatSetProperty", 5, num_args);

        auto G = LuaArg<int>(L, 4);
        auto sigma_t = LuaArg<double>(L, 5);

        prop->MakeSimple0(G, sigma_t);
      }
      else if (operation_index == static_cast<int>(OpType::SIMPLEXS1))
      {
        if (num_args != 6)
          LuaPostArgAmountError("MatSetProperty", 6, num_args);

        auto G = LuaArg<int>(L, 4);
        auto sigma_t = LuaArg<double>(L, 5);
        auto c = LuaArg<double>(L, 6);

        prop->MakeSimple1(G, sigma_t, c);
      }
      else if (operation_index == static_cast<int>(OpType::OPENSN_XSFILE))
      {
        if (num_args != 4)
          LuaPostArgAmountError("MatSetProperty", 4, num_args);

        const auto file_name = LuaArg<std::string>(L, 4);

        prop->MakeFromOpenSnXSFile(file_name);
      }
      else if (operation_index == static_cast<int>(OpType::EXISTING))
      {
        if (num_args != 4)
          LuaPostArgAmountError("MatSetProperty", 4, num_args);

        auto handle = LuaArg<int>(L, 4);

        std::shared_ptr<SingleStateMGXS> xs;
        try
        {
          xs = std::dynamic_pointer_cast<SingleStateMGXS>(
            opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle, fname));
        }
        catch (const std::out_of_range& o)
        {
          opensn::log.LogAllError()
            << "ERROR: Invalid cross-section handle in call to MatSetProperty." << std::endl;
          opensn::Exit(EXIT_FAILURE);
        }
        //        auto old_prop = prop;
        prop = xs;

        cur_material->properties_[location_of_prop] = prop;

        //        delete old_prop; //Still debating if this should be deleted
      }
      else
      {
        opensn::log.LogAllError() << "Unsupported operation for TRANSPORT_XSECTIONS." << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }
    else
    {
      opensn::log.LogAllError() << "Material has no property TRANSPORT_XSECTIONS." << std::endl;
      opensn::Exit(EXIT_FAILURE);
    }
  } // if thermal conductivity
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ISOTROPIC_MG_SOURCE
  else if (property_index == static_cast<int>(MatProperty::ISOTROPIC_MG_SOURCE))
  {
    int location_of_prop = -1;
    // Check if the material has this property
    if (lua_isnumber(L, 2))
    {
      for (int p = 0; p < cur_material->properties_.size(); p++)
      {
        if (cur_material->properties_[p]->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
        {
          location_of_prop = p;
        }
      }
    }
    else
    {
      for (int p = 0; p < cur_material->properties_.size(); p++)
      {
        if (cur_material->properties_[p]->property_name == property_index_name)
        {
          location_of_prop = p;
        }
      }
    }

    // If the property is valid
    if (location_of_prop >= 0)
    {
      auto prop = std::static_pointer_cast<IsotropicMultiGrpSource>(
        cur_material->properties_[location_of_prop]);

      if (operation_index == static_cast<int>(OpType::SINGLE_VALUE))
      {
        if (num_args != 4)
          LuaPostArgAmountError("MatSetProperty", 4, num_args);

        auto value = LuaArg<double>(L, 4);

        prop->source_value_g_.resize(1, value);
        opensn::log.Log0Verbose1() << "Isotropic Multigroup Source value for material at index "
                                   << material_index << " set to " << value;
      }
      else if (operation_index == static_cast<int>(OpType::FROM_ARRAY))
      {
        if (num_args != 4)
          LuaPostArgAmountError("MatSetProperty", 4, num_args);
        prop->source_value_g_ = LuaArg<std::vector<double>>(L, 4);
        opensn::log.Log0Verbose1() << "Isotropic Multigroup Source populated with "
                                   << prop->source_value_g_.size() << " values";
      }
      else
      {
        opensn::log.LogAllError() << "Unsupported operation for ISOTROPIC_MG_SOURCE." << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }
    else
    {
      opensn::log.LogAllError() << "Material \"" << cur_material->name_
                                << "\" has no property ISOTROPIC_MG_SOURCE." << std::endl;
      opensn::Exit(EXIT_FAILURE);
    }
  }
  else
  {
    opensn::log.LogAllError()
      << "Unsupported material property specified in call to MatSetProperty." << property_index
      << std::endl;
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}

int
MatGetProperty(lua_State* L)
{
  const std::string fname = "mat.GetProperty";
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("MatGetProperty", 2, num_args);

  auto material_index = LuaArg<int>(L, 1);
  int property_index = -1;
  std::string property_index_name;
  if (lua_isnumber(L, 2))
    property_index = LuaArg<int>(L, 2);
  else
    property_index_name = LuaArg<std::string>(L, 2);

  // Get reference to material
  auto cur_material = opensn::GetStackItemPtr(opensn::material_stack, material_index, fname);

  // If user supplied name then find property index
  if (not lua_isnumber(L, 2))
  {
    for (auto& property : cur_material->properties_)
      if (property->property_name == property_index_name)
        property_index = static_cast<int>(property->Type());
  }

  // Process property
  bool property_populated = false;
  for (auto& property : cur_material->properties_)
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
