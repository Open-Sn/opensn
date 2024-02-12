#include "framework/lua.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "lua/framework/mesh/lua_mesh_set_property.h"
#include "lua/framework/mesh/lua_mesh_set_ids.h"
#include "framework/console/console.h"
#include <iostream>

using namespace opensn;

RegisterLuaFunctionNamespace(MeshSetProperty, mesh, SetProperty);

enum MeshProperty
{
  MATID_FROMLOGICAL = 11,
  BNDRYID_FROMLOGICAL = 12,
  MATID_FROM_LUA_FUNCTION = 13,
  BNDRYID_FROM_LUA_FUNCTION = 14
};

RegisterLuaConstantAsIs(MATID_FROMLOGICAL, Varying(11));
RegisterLuaConstantAsIs(BNDRYID_FROMLOGICAL, Varying(12));
RegisterLuaConstantAsIs(MATID_FROM_LUA_FUNCTION, Varying(13));
RegisterLuaConstantAsIs(BNDRYID_FROM_LUA_FUNCTION, Varying(14));

int
MeshSetProperty(lua_State* L)
{
  const std::string fname = "mesh.SetProperty";

  // Get property index
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  int property_index = lua_tonumber(L, 1);

  if (property_index == MeshProperty::MATID_FROMLOGICAL)
  {
    if (not((num_args == 3) or (num_args == 4)))
    {
      opensn::log.LogAllError()
        << "Invalid number of arguments when calling mesh.SetProperty(MATID_FROMLOGICAL...";
      opensn::Exit(EXIT_FAILURE);
    }
    int volume_hndl = lua_tonumber(L, 2);
    int mat_id = lua_tonumber(L, 3);
    int sense = true;
    if (num_args == 4)
      sense = lua_toboolean(L, 4);

    const auto& log_vol =
      opensn::GetStackItem<LogicalVolume>(opensn::object_stack, volume_hndl, fname);

    opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                               << " Setting material id from logical volume.";
    std::shared_ptr<MeshContinuum> vol_cont = GetCurrentMesh();
    vol_cont->SetMaterialIDFromLogical(log_vol, sense, mat_id);
  }

  else if (property_index == MeshProperty::BNDRYID_FROMLOGICAL)
  {
    if (not((num_args == 3) or (num_args == 4)))
    {
      opensn::log.LogAllError() << "Invalid amount of arguments used for "
                                   "VolumeMesherSetProperty("
                                   "BNDRYID_FROMLOGICAL...";
      opensn::Exit(EXIT_FAILURE);
    }
    LuaCheckNilValue(fname, L, 2);
    LuaCheckStringValue(fname, L, 3);
    int volume_hndl = lua_tonumber(L, 2);
    std::string bndry_name = lua_tostring(L, 3);
    int sense = true;
    if (num_args == 4)
      sense = lua_toboolean(L, 4);

    ChiLogicalErrorIf(bndry_name.empty(), "argument 3 must not be an empty string.");

    const auto& log_vol =
      opensn::GetStackItem<LogicalVolume>(opensn::object_stack, volume_hndl, fname);

    opensn::log.Log() << program_timer.GetTimeString()
                      << " Setting boundary id from logical volume.";
    std::shared_ptr<MeshContinuum> vol_cont = GetCurrentMesh();
    vol_cont->SetBoundaryIDFromLogical(log_vol, sense, bndry_name);
  }
  else if (property_index == MeshProperty::MATID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    SetMatIDFromLuaFunction(lua_fname);
  }
  else if (property_index == MeshProperty::BNDRYID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    SetBoundaryIDFromLuaFunction(lua_fname);
  }
  else
  {
    opensn::log.LogAllError() << "Invalid property specified " << property_index
                              << " in call to VolumeMesherSetProperty().";
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}
