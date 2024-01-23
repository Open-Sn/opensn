#include "framework/lua.h"
#include "framework/mesh/volume_mesher/predefined_unpartitioned/volmesher_predefunpart.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"

#include <iostream>

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "volume_mesher_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiVolumeMesherCreate);
RegisterLuaConstantAsIs(VOLUMEMESHER_EXTRUDER, Varying(4));
RegisterLuaConstantAsIs(VOLUMEMESHER_UNPARTITIONED, Varying(6));

RegisterLuaConstant(ExtruderTemplateType, SURFACE_MESH, Varying(1));
RegisterLuaConstant(ExtruderTemplateType, UNPARTITIONED_MESH, Varying(2));

int
chiVolumeMesherCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  // Arguments check
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  // Mesher type
  const auto mesher_type = static_cast<opensn::VolumeMesherType>(lua_tointeger(L, 1));

  std::shared_ptr<opensn::VolumeMesher> new_mesher = nullptr;

  if (mesher_type == opensn::VolumeMesherType::UNPARTITIONED)
  {
    if (num_args != 2)
    {
      opensn::log.LogAllError() << fname + ": "
                                           "When specifying VOLUMEMESHER_UNPARTITIONED, the "
                                           "handle must also be supplied.";
      opensn::Exit(EXIT_FAILURE);
    }

    LuaCheckNilValue(fname, L, 2);
    const int template_handle = lua_tonumber(L, 2);

    auto p_umesh =
      opensn::Chi::GetStackItemPtr(opensn::Chi::unpartitionedmesh_stack, template_handle, fname);

    new_mesher = std::make_shared<opensn::VolumeMesherPredefinedUnpartitioned>(p_umesh);
  }
  else
  {
    opensn::log.Log0Error() << "Invalid Volume mesher type in function "
                               "chiVolumeMesherCreate. Allowed options are"
                               "VOLUMEMESHER_EXTRUDER or "
                               "VOLUMEMESHER_UNPARTITIONED";
    opensn::Exit(EXIT_FAILURE);
  }

  auto& cur_hndlr = opensn::GetCurrentHandler();
  cur_hndlr.SetVolumeMesher(new_mesher);

  opensn::log.LogAllVerbose2() << "chiVolumeMesherCreate: Volume mesher created." << std::endl;

  return 0;
}
