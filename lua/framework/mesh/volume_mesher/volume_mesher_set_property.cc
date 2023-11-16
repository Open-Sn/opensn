#include "framework/lua.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/surface_mesher/surface_mesher.h"
#include "framework/mesh/volume_mesher/extruder/volmesher_extruder.h"

#include "framework/runtime.h"
#include "framework/mesh/logical_volume/logical_volume.h"

#include <iostream>
#include "volume_mesher_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiVolumeMesherSetProperty);

RegisterLuaConstantAsIs(FORCE_POLYGONS, Varying(1));
RegisterLuaConstantAsIs(MESH_GLOBAL, Varying(2));
RegisterLuaConstantAsIs(PARTITION_Z, Varying(3));
RegisterLuaConstantAsIs(VOLUMEPARTITION_Y, Varying(4));
RegisterLuaConstantAsIs(VOLUMEPARTITION_X, Varying(5));
RegisterLuaConstantAsIs(CUTS_Z, Varying(6));
RegisterLuaConstantAsIs(CUTS_Y, Varying(7));
RegisterLuaConstantAsIs(CUTS_X, Varying(8));
RegisterLuaConstantAsIs(PARTITION_TYPE, Varying(9));
RegisterLuaConstantAsIs(KBA_STYLE_XYZ, Varying(2));
RegisterLuaConstantAsIs(PARMETIS, Varying(3));
RegisterLuaConstantAsIs(EXTRUSION_LAYER, Varying(10));
RegisterLuaConstantAsIs(MATID_FROMLOGICAL, Varying(11));
RegisterLuaConstantAsIs(BNDRYID_FROMLOGICAL, Varying(12));
RegisterLuaConstantAsIs(MATID_FROM_LUA_FUNCTION, Varying(13));
RegisterLuaConstantAsIs(BNDRYID_FROM_LUA_FUNCTION, Varying(14));

RegisterLuaFunctionAsIs(chiVolumeMesherSetKBAPartitioningPxPyPz);
RegisterLuaFunctionAsIs(chiVolumeMesherSetKBACutsX);
RegisterLuaFunctionAsIs(chiVolumeMesherSetKBACutsY);
RegisterLuaFunctionAsIs(chiVolumeMesherSetKBACutsZ);

int
chiVolumeMesherSetProperty(lua_State* L)
{
  const std::string fname = "chiVolumeMesherSetProperty";
  // Get current mesh handler
  auto& cur_hndlr = GetCurrentHandler();
  auto& volume_mesher = cur_hndlr.GetVolumeMesher();

  // Get property index
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  int property_index = lua_tonumber(L, 1);

  typedef VolumeMesherProperty VMP;

  // Selects property
  if (property_index == VMP::FORCE_POLYGONS)
  {
    bool force_condition = lua_toboolean(L, 2);
    volume_mesher.options.force_polygons = force_condition;
  }

  else if (property_index == VMP::MESH_GLOBAL)
  {
    bool mesh_global = lua_toboolean(L, 2);
    volume_mesher.options.mesh_global = mesh_global;
  }

  else if (property_index == VMP::PARTITION_Z)
  {
    int pz = lua_tonumber(L, 2);
    volume_mesher.options.partition_z = pz;
    opensn::log.LogAllVerbose1() << "Partition z set to " << pz;
  }
  else if (property_index == VMP::PARTITION_Y)
  {
    int p = lua_tonumber(L, 2);
    volume_mesher.options.partition_y = p;
    opensn::log.LogAllVerbose1() << "Partition y set to " << p;
  }
  else if (property_index == VMP::PARTITION_X)
  {
    int p = lua_tonumber(L, 2);
    volume_mesher.options.partition_x = p;
    opensn::log.LogAllVerbose1() << "Partition x set to " << p;
  }
  else if (property_index == VMP::CUTS_Z)
  {
    double p = lua_tonumber(L, 2);
    volume_mesher.options.zcuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_Y)
  {
    double p = lua_tonumber(L, 2);
    cur_hndlr.GetSurfaceMesher().AddYCut(p);
    volume_mesher.options.ycuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_X)
  {
    double p = lua_tonumber(L, 2);
    cur_hndlr.GetSurfaceMesher().AddXCut(p);
    volume_mesher.options.xcuts.push_back(p);
  }
  else if (property_index == VMP::PARTITION_TYPE)
  {
    int p = lua_tonumber(L, 2);
    if (p >= VolumeMesher::PartitionType::KBA_STYLE_XYZ and
        p <= VolumeMesher::PartitionType::PARMETIS)
      volume_mesher.options.partition_type = (VolumeMesher::PartitionType)p;
    else
    {
      opensn::log.LogAllError() << "Unsupported partition type used in call to " << fname << ".";
      opensn::Chi::Exit(EXIT_FAILURE);
    }
  }

  else if (property_index == VMP::EXTRUSION_LAYER)
  {
    if (volume_mesher.Type() == VolumeMesherType::EXTRUDER)
    {
      auto& mesher = dynamic_cast<VolumeMesherExtruder&>(volume_mesher);

      double layer_height = lua_tonumber(L, 2);
      int subdivisions = 1;

      if (num_args >= 3) { subdivisions = lua_tonumber(L, 3); }
      VolumeMesherExtruder::MeshLayer new_layer;
      new_layer.height = layer_height;
      new_layer.sub_divisions = subdivisions;

      if (num_args == 4) { new_layer.name = std::string(lua_tostring(L, 4)); }
      mesher.AddLayer(new_layer);
    }
    else
    {
      fprintf(stdout,
              "VolumeMesherExtruder is not the current volume mesher"
              " therefore the z-layer property is ignored.\n");
    }
  }

  else if (property_index == VMP::MATID_FROMLOGICAL)
  {
    if (!((num_args == 3) || (num_args == 4)))
    {
      opensn::log.LogAllError() << "Invalid amount of arguments used for "
                                   "chiVolumeMesherSetProperty("
                                   "MATID_FROMLOGICAL...";
      opensn::Chi::Exit(EXIT_FAILURE);
    }
    int volume_hndl = lua_tonumber(L, 2);
    int mat_id = lua_tonumber(L, 3);
    int sense = true;
    if (num_args == 4) sense = lua_toboolean(L, 4);

    const auto& log_vol =
      opensn::Chi::GetStackItem<LogicalVolume>(opensn::Chi::object_stack, volume_hndl, fname);

    VolumeMesher::SetMatIDFromLogical(log_vol, sense, mat_id);
  }

  else if (property_index == VMP::BNDRYID_FROMLOGICAL)
  {
    if (!((num_args == 3) || (num_args == 4)))
    {
      opensn::log.LogAllError() << "Invalid amount of arguments used for "
                                   "chiVolumeMesherSetProperty("
                                   "BNDRYID_FROMLOGICAL...";
      opensn::Chi::Exit(EXIT_FAILURE);
    }
    LuaCheckNilValue(fname, L, 2);
    LuaCheckStringValue(fname, L, 3);
    int volume_hndl = lua_tonumber(L, 2);
    std::string bndry_name = lua_tostring(L, 3);
    int sense = true;
    if (num_args == 4) sense = lua_toboolean(L, 4);

    if (bndry_name.empty())
      throw std::invalid_argument(fname + ": argument 3 must not be"
                                          "an empty string.");

    const auto& log_vol =
      opensn::Chi::GetStackItem<LogicalVolume>(opensn::Chi::object_stack, volume_hndl, fname);

    VolumeMesher::SetBndryIDFromLogical(log_vol, sense, bndry_name);
  }
  else if (property_index == VMP::MATID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    VolumeMesher::SetMatIDFromLuaFunction(lua_fname);
  }
  else if (property_index == VMP::BNDRYID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    VolumeMesher::SetBndryIDFromLuaFunction(lua_fname);
  }
  else
  {
    opensn::log.LogAllError() << "Invalid property specified " << property_index
                              << " in call to chiVolumeMesherSetProperty().";
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  return 0;
}

int
chiVolumeMesherSetKBAPartitioningPxPyPz(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(__FUNCTION__, 3, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);
  LuaCheckNilValue(__FUNCTION__, L, 3);

  // Get current mesh handler
  auto& cur_hndlr = GetCurrentHandler();
  auto& vol_mesher = cur_hndlr.GetVolumeMesher();

  int px = lua_tonumber(L, 1);
  int py = lua_tonumber(L, 2);
  int pz = lua_tonumber(L, 3);

  vol_mesher.options.partition_x = px;
  vol_mesher.options.partition_y = py;
  vol_mesher.options.partition_z = pz;

  return 0;
}

int
chiVolumeMesherSetKBACutsX(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  std::vector<double> cuts;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 1, cuts);

  auto& mesh_handler = GetCurrentHandler();
  mesh_handler.GetVolumeMesher().options.xcuts = cuts;

  return 0;
}

int
chiVolumeMesherSetKBACutsY(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  std::vector<double> cuts;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 1, cuts);

  auto& mesh_handler = GetCurrentHandler();
  mesh_handler.GetVolumeMesher().options.ycuts = cuts;

  return 0;
}

int
chiVolumeMesherSetKBACutsZ(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  std::vector<double> cuts;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 1, cuts);

  auto& mesh_handler = GetCurrentHandler();
  mesh_handler.GetVolumeMesher().options.zcuts = cuts;

  return 0;
}
