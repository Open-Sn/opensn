#include "framework/lua.h"
#include <iostream>

#include "framework/physics/physics_material/physics_material.h"
#include "framework/physics/physics_material/material_property_scalar_value.h"
#include "framework/physics/physics_material/multi_group_xs/single_state_mgxs.h"
#include "framework/physics/physics_material/material_property_isotropic_mg_src.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "physics_lua_utils.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(PhysicsMaterialAddProperty);
RegisterLuaFunctionAsIs(PhysicsMaterialSetProperty);
RegisterLuaFunctionAsIs(chiPhysicsMaterialGetProperty);

RegisterLuaConstantAsIs(SCALAR_VALUE, Varying(1));
RegisterLuaConstantAsIs(TRANSPORT_XSECTIONS, Varying(10));
RegisterLuaConstantAsIs(ISOTROPIC_MG_SOURCE, Varying(11));

namespace
{

void
ScalarPropertyPushTable(lua_State* L, std::shared_ptr<MaterialProperty> property)
{
  lua_newtable(L);
  lua_pushstring(L, "is_empty");
  lua_pushboolean(L, false);
  lua_settable(L, -3);

  lua_pushstring(L, "value");
  lua_pushnumber(L, property->GetScalarValue());
  lua_settable(L, -3);
}

void
IsotropicMGSourcePropertyPushTable(lua_State* L, std::shared_ptr<IsotropicMultiGrpSource> property)
{
  lua_newtable(L);
  lua_pushstring(L, "is_empty");
  lua_pushboolean(L, false);
  lua_settable(L, -3);

  lua_pushstring(L, "G");
  lua_pushnumber(L, property->source_value_g_.size());
  lua_settable(L, -3);

  lua_pushstring(L, "source_value_g");
  lua_newtable(L);
  int g = 0;
  for (auto val : property->source_value_g_)
  {
    ++g;
    lua_pushnumber(L, g);
    lua_pushnumber(L, val);
    lua_settable(L, -3);
  }
  lua_settable(L, -3);
}

void
MaterialPropertyPushLuaTable(lua_State* L)
{
  lua_newtable(L);
  lua_pushstring(L, "is_empty");
  lua_pushboolean(L, true);
  lua_settable(L, -3);
}

void
PropertyPushLuaTable(lua_State* L, std::shared_ptr<MaterialProperty> property)
{
  if (property->Type() == PropertyType::SCALAR_VALUE) ScalarPropertyPushTable(L, property);
  else if (property->Type() == PropertyType::ISOTROPIC_MG_SOURCE)
    IsotropicMGSourcePropertyPushTable(
      L, std::dynamic_pointer_cast<IsotropicMultiGrpSource>(property));
  else
    MaterialPropertyPushLuaTable(L);
}

} // namespace

int
PhysicsMaterialAddProperty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int numArgs = lua_gettop(L);

  if (!((numArgs >= 2) && (numArgs <= 3)))
  {
    opensn::log.Log0Error() << "Incorrect amount of arguments "
                               "in PhysicsMaterialAddProperty";
    opensn::Exit(EXIT_FAILURE);
  }

  int material_index = lua_tonumber(L, 1);
  int property_index = lua_tonumber(L, 2);

  const char* provided_name = "";
  if (numArgs == 3) { provided_name = lua_tostring(L, 3); }

  // Get reference to material
  auto cur_material = opensn::GetStackItemPtr(opensn::material_stack, material_index, fname);

  // Process property
  using MatProperty = PropertyType;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALAR_VALUE
  if (property_index == static_cast<int>(MatProperty::SCALAR_VALUE))
  {
    // Duplicates are allowed

    auto prop = std::make_shared<ScalarValue>();

    prop->property_name =
      std::string("Property ") + std::to_string(cur_material->properties_.size());

    if (numArgs == 3) prop->property_name = std::string(provided_name);

    cur_material->properties_.push_back(prop);
    opensn::log.Log0Verbose1() << "Scalar Value Property added to material"
                                  " at index "
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
                                << " already has property "
                                   "TRANSPORT_XSECTIONS"
                                << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }

    auto prop = std::make_shared<SingleStateMGXS>();

    prop->property_name =
      std::string("Property ") + std::to_string(cur_material->properties_.size());

    if (numArgs == 3) prop->property_name = std::string(provided_name);

    cur_material->properties_.push_back(prop);
    opensn::log.Log0Verbose1() << "Transport cross-sections added to material"
                                  " at index "
                               << material_index;

    opensn::multigroup_xs_stack.push_back(prop);

    const size_t index = opensn::multigroup_xs_stack.size() - 1;

    lua_pushnumber(L, static_cast<lua_Number>(index));
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
                                << " already has property "
                                   "ISOTROPIC_MG_SOURCE "
                                << property_index << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }

    auto prop = std::make_shared<IsotropicMultiGrpSource>();

    prop->property_name =
      std::string("Property ") + std::to_string(cur_material->properties_.size());

    if (numArgs == 3) prop->property_name = std::string(provided_name);

    cur_material->properties_.push_back(prop);
    opensn::log.Log0Verbose1() << "Isotropic Multigroup Source added to material"
                                  " at index "
                               << material_index;
  }
  else
  {
    opensn::log.Log0Error() << "Unsupported property type in call to PhysicsMaterialAddProperty.";
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}

int
PhysicsMaterialSetProperty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int numArgs = lua_gettop(L);

  if (numArgs < 3)
  {
    opensn::log.Log0Error() << "Incorrect amount of arguments "
                               "in PhysicsMaterialSetProperty";
    opensn::Exit(EXIT_FAILURE);
  }

  int material_index = lua_tonumber(L, 1);
  int property_index = -1;
  std::string property_index_name;
  if (lua_isnumber(L, 2)) property_index = lua_tonumber(L, 2);
  else
  {
    const char* temp_name = lua_tostring(L, 2);
    property_index_name = std::string(temp_name);
  }

  int operation_index = lua_tonumber(L, 3);

  // Get reference to material
  auto cur_material = opensn::GetStackItemPtr(opensn::material_stack, material_index, fname);

  // If user supplied name then find property index
  if (!lua_isnumber(L, 2))
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
        if (cur_material->properties_[p]->Type() == MatProperty::SCALAR_VALUE) location_of_prop = p;
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
        double value = lua_tonumber(L, 4);
        prop->value_ = value;
        opensn::log.Log0Verbose1() << "Scalar value for material"
                                      " at index "
                                   << material_index << " set to " << value;
      }
      else
      {
        opensn::log.Log0Error() << "ERROR: Unsupported operation for "
                                   "SCALAR_VALUE."
                                << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }
    else
    {
      opensn::log.Log0Error() << "ERROR: Material has no property "
                                 "SCALAR_VALUE."
                              << std::endl;
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
        if (numArgs != 5) LuaPostArgAmountError("PhysicsMaterialSetProperty", 5, numArgs);

        int G = lua_tonumber(L, 4);
        double sigma_t = lua_tonumber(L, 5);

        prop->MakeSimple0(G, sigma_t);
      }
      else if (operation_index == static_cast<int>(OpType::SIMPLEXS1))
      {
        if (numArgs != 6) LuaPostArgAmountError("PhysicsMaterialSetProperty", 6, numArgs);

        int G = lua_tonumber(L, 4);
        double sigma_t = lua_tonumber(L, 5);
        double c = lua_tonumber(L, 6);

        prop->MakeSimple1(G, sigma_t, c);
      }
      else if (operation_index == static_cast<int>(OpType::CHI_XSFILE))
      {
        if (numArgs != 4) LuaPostArgAmountError("PhysicsMaterialSetProperty", 4, numArgs);

        const char* file_name_c = lua_tostring(L, 4);

        prop->MakeFromChiXSFile(std::string(file_name_c));
      }
      else if (operation_index == static_cast<int>(OpType::EXISTING))
      {
        if (numArgs != 4) LuaPostArgAmountError("PhysicsMaterialSetProperty", 4, numArgs);

        LuaCheckNilValue("PhysicsMaterialSetProperty", L, 4);
        int handle = lua_tonumber(L, 4);

        std::shared_ptr<SingleStateMGXS> xs;
        try
        {
          xs = std::dynamic_pointer_cast<SingleStateMGXS>(
            opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle, fname));
        }
        catch (const std::out_of_range& o)
        {
          opensn::log.LogAllError() << "ERROR: Invalid cross-section handle"
                                    << " in call to PhysicsMaterialSetProperty." << std::endl;
          opensn::Exit(EXIT_FAILURE);
        }
        //        auto old_prop = prop;
        prop = xs;

        cur_material->properties_[location_of_prop] = prop;

        //        delete old_prop; //Still debating if this should be deleted
      }
      else
      {
        opensn::log.LogAllError() << "Unsupported operation for "
                                     "TRANSPORT_XSECTIONS."
                                  << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }
    else
    {
      opensn::log.LogAllError() << "Material has no property "
                                   "TRANSPORT_XSECTIONS."
                                << std::endl;
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
        if (numArgs != 4) LuaPostArgAmountError("PhysicsMaterialSetProperty", 4, numArgs);

        double value = lua_tonumber(L, 4);

        prop->source_value_g_.resize(1, value);
        opensn::log.Log0Verbose1() << "Isotropic Multigroup Source value "
                                      "for material"
                                      " at index "
                                   << material_index << " set to " << value;
      }
      else if (operation_index == static_cast<int>(OpType::FROM_ARRAY))
      {
        if (numArgs != 4) LuaPostArgAmountError("PhysicsMaterialSetProperty", 4, numArgs);

        if (!lua_istable(L, 4))
        {
          opensn::log.LogAllError() << "In call to PhysicsMaterialSetProperty: "
                                    << "Material \"" << cur_material->name_ << "\", when setting "
                                    << "ISOTROPIC_MG_SOURCE using operation FROM_ARRAY, the fourth "
                                       "argument was detected not to be a lua table.";
          opensn::Exit(EXIT_FAILURE);
        }

        const size_t table_len = lua_rawlen(L, 4);

        std::vector<double> values(table_len, 0.0);
        for (int g = 0; g < table_len; g++)
        {
          lua_pushnumber(L, g + 1);
          lua_gettable(L, 4);
          values[g] = lua_tonumber(L, -1);
          lua_pop(L, 1);
        }

        prop->source_value_g_.resize(table_len, 0.0);
        std::copy(values.begin(), values.end(), prop->source_value_g_.begin());
        opensn::log.Log0Verbose1() << "Isotropic Multigroup Source populated "
                                   << " with " << table_len << " values";
      }
      else
      {
        opensn::log.LogAllError() << "Unsupported operation for "
                                     "ISOTROPIC_MG_SOURCE."
                                  << std::endl;
        opensn::Exit(EXIT_FAILURE);
      }
    }
    else
    {
      opensn::log.LogAllError() << "Material \"" << cur_material->name_
                                << "\" has no property "
                                   "ISOTROPIC_MG_SOURCE."
                                << std::endl;
      opensn::Exit(EXIT_FAILURE);
    }
  }
  else
  {
    opensn::log.LogAllError() << "Unsupported material property specified in "
                                 "call to PhysicsMaterialSetProperty."
                              << property_index << std::endl;
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}

int
chiPhysicsMaterialGetProperty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError("chiPhysicsMaterialGetProperty", 2, num_args);

  int material_index = lua_tonumber(L, 1);
  int property_index = -1;
  std::string property_index_name;
  if (lua_isnumber(L, 2)) property_index = lua_tonumber(L, 2);
  else
  {
    const char* temp_name = lua_tostring(L, 2);
    property_index_name = std::string(temp_name);
  }

  // Get reference to material
  auto cur_material = opensn::GetStackItemPtr(opensn::material_stack, material_index, fname);

  // If user supplied name then find property index
  if (!lua_isnumber(L, 2))
  {
    for (auto& property : cur_material->properties_)
      if (property->property_name == property_index_name)
        property_index = static_cast<int>(property->Type());
  }

  // Process property
  bool property_polulated = false;
  for (auto& property : cur_material->properties_)
  {
    if (static_cast<int>(property->Type()) == property_index)
    {
      PropertyPushLuaTable(L, property);
      property_polulated = true;
    }
  }

  if (not property_polulated)
  {
    opensn::log.LogAllError() << "Invalid material property specified in "
                                 "call to chiPhysicsMaterialGetProperty."
                              << property_index << std::endl;
    opensn::Exit(EXIT_FAILURE);
  }

  return 1;
}
