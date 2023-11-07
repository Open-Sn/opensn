#include "framework/lua.h"
#include <iostream>

#include "framework/physics/physics_namespace.h"
#include "framework/physics/physics_material/multi_group_xs/single_state_mgxs.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "multigroup_xs_lua_utils.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiPhysicsTransportXSCreate);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSSet);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSMakeCombined);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSSetCombined);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSGet);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSExportToChiTechFormat);

RegisterLuaConstantAsIs(SINGLE_VALUE, chi_data_types::Varying(0));
RegisterLuaConstantAsIs(FROM_ARRAY, chi_data_types::Varying(1));
RegisterLuaConstantAsIs(SIMPLEXS0, chi_data_types::Varying(20));
RegisterLuaConstantAsIs(SIMPLEXS1, chi_data_types::Varying(21));
RegisterLuaConstantAsIs(EXISTING, chi_data_types::Varying(22));
RegisterLuaConstantAsIs(CHI_XSFILE, chi_data_types::Varying(23));

int
chiPhysicsTransportXSCreate(lua_State* L)
{
  auto xs = std::make_shared<chi_physics::SingleStateMGXS>();

  Chi::multigroup_xs_stack.push_back(xs);

  const size_t index = Chi::multigroup_xs_stack.size() - 1;

  lua_pushinteger(L, static_cast<lua_Integer>(index));
  return 1;
}

int
chiPhysicsTransportXSSet(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args < 3)
  {
    LuaPostArgAmountError("chiPhysicsTransportXSSet", 3, num_args);
    Chi::Exit(EXIT_FAILURE);
  }

  LuaCheckNilValue("chiPhysicsTransportXSSet", L, 1);
  LuaCheckNilValue("chiPhysicsTransportXSSet", L, 2);

  int handle = lua_tonumber(L, 1);
  int operation_index = lua_tonumber(L, 2);

  std::shared_ptr<chi_physics::SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<chi_physics::SingleStateMGXS>(
      Chi::GetStackItemPtr(Chi::multigroup_xs_stack, handle));
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "ERROR: Invalid cross section handle"
                           << " in call to chiPhysicsTransportXSSet." << std::endl;
    Chi::Exit(EXIT_FAILURE);
  }

  // Process operation
  using OpType = chi_physics::OperationType;
  if (operation_index == static_cast<int>(OpType::SIMPLEXS0))
  {
    if (num_args != 4) LuaPostArgAmountError("chiPhysicsTransportXSSet", 4, num_args);

    int G = lua_tonumber(L, 3);
    double sigma_t = lua_tonumber(L, 4);

    xs->MakeSimple0(G, sigma_t);
  }
  else if (operation_index == static_cast<int>(OpType::SIMPLEXS1))
  {
    if (num_args != 5) LuaPostArgAmountError("chiPhysicsTransportXSSet", 5, num_args);

    int G = lua_tonumber(L, 3);
    double sigma_t = lua_tonumber(L, 4);
    double c = lua_tonumber(L, 5);

    xs->MakeSimple1(G, sigma_t, c);
  }
  else if (operation_index == static_cast<int>(OpType::CHI_XSFILE))
  {
    if (num_args != 3) LuaPostArgAmountError("chiPhysicsTransportXSSet", 3, num_args);

    const char* file_name_c = lua_tostring(L, 3);

    xs->MakeFromChiXSFile(std::string(file_name_c));
  }
  else
  {
    Chi::log.LogAllError() << "Unsupported operation in "
                           << "chiPhysicsTransportXSSet. " << operation_index << std::endl;
    Chi::Exit(EXIT_FAILURE);
  }
  return 0;
}

int
chiPhysicsTransportXSGet(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args < 1)
  {
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);
    Chi::Exit(EXIT_FAILURE);
  }

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int handle = lua_tonumber(L, 1);

  std::shared_ptr<chi_physics::SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<chi_physics::SingleStateMGXS>(
      Chi::GetStackItemPtr(Chi::multigroup_xs_stack, handle));
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "ERROR: Invalid cross section handle"
                           << " in call to " << __FUNCTION__ << "." << std::endl;
    Chi::Exit(EXIT_FAILURE);
  }

  xs->PushLuaTable(L);

  return 1;
}

int
chiPhysicsTransportXSMakeCombined(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError("chiPhysicsMakeCombinedTransportXS", 1, num_args);

  if (!lua_istable(L, 1))
  {
    Chi::log.LogAllError() << "In call to chiPhysicsMakeCombinedTransportXS: "
                           << "Argument must be a lua table.";
    Chi::Exit(EXIT_FAILURE);
  }

  size_t table_len = lua_rawlen(L, 1);

  std::vector<std::pair<int, double>> combinations;
  combinations.reserve(table_len);

  // Process table
  for (int v = 0; v < table_len; ++v)
  {
    lua_pushnumber(L, v + 1);
    lua_gettable(L, 1);

    if (!lua_istable(L, -1))
    {
      Chi::log.LogAllError() << "In call to chiPhysicsMakeCombinedTransportXS: "
                             << "The elements of the supplied table must themselves also"
                                "be lua tables of the xs handle and its scalar multiplier.";
      Chi::Exit(EXIT_FAILURE);
    }

    lua_pushinteger(L, 1);
    lua_gettable(L, -2);
    LuaCheckNilValue("chiPhysicsMakeCombinedTransportXS:A1:E1", L, -1);

    int handle = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_pushinteger(L, 2);
    lua_gettable(L, -2);
    LuaCheckNilValue("chiPhysicsMakeCombinedTransportXS:A1:E2", L, -1);

    double scalar = lua_tonumber(L, -1);
    lua_pop(L, 1);

    combinations.emplace_back(handle, scalar);
    lua_pop(L, 1); // pop off table
  }

  // Print out table
  Chi::log.Log() << "Generating XS with following combination:";
  for (auto& elem : combinations)
    Chi::log.Log() << "  Element handle: " << elem.first << " scalar value: " << elem.second;

  // Make the new cross section
  auto new_xs = std::make_shared<chi_physics::SingleStateMGXS>();

  new_xs->MakeCombined(combinations);

  Chi::multigroup_xs_stack.push_back(new_xs);
  lua_pushinteger(L, static_cast<lua_Integer>(Chi::multigroup_xs_stack.size()) - 1);

  return 1;
}

int
chiPhysicsTransportXSSetCombined(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args < 2)
  {
    LuaPostArgAmountError(__FUNCTION__, 2, num_args);
    Chi::Exit(EXIT_FAILURE);
  }

  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);
  LuaCheckTableValue(__FUNCTION__, L, 2);

  // Process handle
  int xs_handle = lua_tonumber(L, 1);

  std::shared_ptr<chi_physics::SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<chi_physics::SingleStateMGXS>(
      Chi::GetStackItemPtr(Chi::multigroup_xs_stack, xs_handle));
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "ERROR: Invalid cross section handle"
                           << " in call to " << __FUNCTION__ << "." << std::endl;
    Chi::Exit(EXIT_FAILURE);
  }

  // Process table
  size_t table_len = lua_rawlen(L, 2);

  std::vector<std::pair<int, double>> combinations;
  combinations.reserve(table_len);

  for (int v = 0; v < table_len; ++v)
  {
    lua_pushnumber(L, v + 1);
    lua_gettable(L, 1);

    if (!lua_istable(L, -1))
    {
      Chi::log.LogAllError() << "In call to " << __FUNCTION__ << ": "
                             << "The elements of the supplied table must themselves also"
                                "be lua tables of the xs handle and its scalar multiplier.";
      Chi::Exit(EXIT_FAILURE);
    }

    lua_pushinteger(L, 1);
    lua_gettable(L, -2);
    LuaCheckNilValue((std::string(__FUNCTION__) + ":A1:E1").c_str(), L, -1);

    int handle = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_pushinteger(L, 2);
    lua_gettable(L, -2);
    LuaCheckNilValue((std::string(__FUNCTION__) + ":A1:E2").c_str(), L, -1);

    double scalar = lua_tonumber(L, -1);
    lua_pop(L, 1);

    combinations.emplace_back(handle, scalar);
    lua_pop(L, 1); // pop off table
  }

  // Print out table
  Chi::log.Log() << "Setting XS with following combination:";
  for (auto& elem : combinations)
    Chi::log.Log() << "  Element handle: " << elem.first << " scalar value: " << elem.second;

  xs->MakeCombined(combinations);

  return 0;
}

int
chiPhysicsTransportXSExportToChiTechFormat(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2) LuaPostArgAmountError(__FUNCTION__, 2, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);

  // Process handle
  int handle = lua_tonumber(L, 1);

  std::shared_ptr<chi_physics::MultiGroupXS> xs;
  try
  {
    xs = Chi::GetStackItemPtr(Chi::multigroup_xs_stack, handle);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "ERROR: Invalid cross section handle"
                           << " in call to " << __FUNCTION__ << "." << std::endl;
    Chi::Exit(EXIT_FAILURE);
  }

  std::string file_name = lua_tostring(L, 2);

  xs->ExportToChiXSFile(file_name);

  return 0;
}
