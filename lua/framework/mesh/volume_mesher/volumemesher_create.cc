#include "framework/lua.h"
#include "framework/mesh/volume_mesher/extruder/volmesher_extruder.h"
#include "framework/mesh/volume_mesher/predefined_unpartitioned/volmesher_predefunpart.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"

#include <iostream>

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "volumemesher_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiVolumeMesherCreate);
RegisterLuaConstantAsIs(VOLUMEMESHER_EXTRUDER, chi_data_types::Varying(4));
RegisterLuaConstantAsIs(VOLUMEMESHER_UNPARTITIONED, chi_data_types::Varying(6));

RegisterLuaConstant(ExtruderTemplateType, SURFACE_MESH, chi_data_types::Varying(1));
RegisterLuaConstant(ExtruderTemplateType, UNPARTITIONED_MESH, chi_data_types::Varying(2));

int
chiVolumeMesherCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  // Arguments check
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  // Mesher type
  const auto mesher_type = static_cast<chi_mesh::VolumeMesherType>(lua_tointeger(L, 1));

  std::shared_ptr<chi_mesh::VolumeMesher> new_mesher = nullptr;

  if (mesher_type == chi_mesh::VolumeMesherType::EXTRUDER)
  {
    if (num_args != 3)
    {
      Chi::log.LogAllError() << fname +
                                  ": "
                                  "When specifying VOLUMEMESHER_EXTRUDER, the template type and "
                                  "handle must also be supplied.";
      Chi::Exit(EXIT_FAILURE);
    }

    LuaCheckNilValue(fname, L, 2);
    LuaCheckNilValue(fname, L, 3);

    int template_type = lua_tonumber(L, 2);
    int template_handle = lua_tonumber(L, 3);

    const auto UNPART_MESH_TEMPLATE =
      chi_mesh::VolumeMesherExtruder::TemplateType::UNPARTITIONED_MESH;

    if (template_type == (int)UNPART_MESH_TEMPLATE)
    {
      auto p_umesh = Chi::GetStackItemPtr(Chi::unpartitionedmesh_stack, template_handle, fname);

      new_mesher = std::make_shared<chi_mesh::VolumeMesherExtruder>(p_umesh);
    }
    else
    {
      Chi::log.LogAllError() << "In call to " << __FUNCTION__
                             << ". Invalid template type specified.";
      Chi::Exit(EXIT_FAILURE);
    }
  }
  else if (mesher_type == chi_mesh::VolumeMesherType::UNPARTITIONED)
  {
    if (num_args != 2)
    {
      Chi::log.LogAllError() << fname + ": "
                                        "When specifying VOLUMEMESHER_UNPARTITIONED, the "
                                        "handle must also be supplied.";
      Chi::Exit(EXIT_FAILURE);
    }

    LuaCheckNilValue(fname, L, 2);
    const int template_handle = lua_tonumber(L, 2);

    auto p_umesh = Chi::GetStackItemPtr(Chi::unpartitionedmesh_stack, template_handle, fname);

    new_mesher = std::make_shared<chi_mesh::VolumeMesherPredefinedUnpartitioned>(p_umesh);
  }
  else
  {
    Chi::log.Log0Error() << "Invalid Volume mesher type in function "
                            "chiVolumeMesherCreate. Allowed options are"
                            "VOLUMEMESHER_EXTRUDER or "
                            "VOLUMEMESHER_UNPARTITIONED";
    Chi::Exit(EXIT_FAILURE);
  }

  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  cur_hndlr.SetVolumeMesher(new_mesher);

  Chi::log.LogAllVerbose2() << "chiVolumeMesherCreate: Volume mesher created." << std::endl;

  return 0;
}
