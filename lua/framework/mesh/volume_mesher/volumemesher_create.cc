#include "opensn/framework/chi_lua.h"
#include "opensn/framework/mesh/VolumeMesher/Extruder/volmesher_extruder.h"
#include "opensn/framework/mesh/VolumeMesher/PredefinedUnpartitioned/volmesher_predefunpart.h"

#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"
#include "opensn/framework/mesh/UnpartitionedMesh/unpartitioned_mesh.h"

#include <iostream>

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"
#include "volumemesher_lua.h"
#include "opensn/framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiVolumeMesherCreate);
RegisterLuaConstantAsIs(VOLUMEMESHER_EXTRUDER, chi_data_types::Varying(4));
RegisterLuaConstantAsIs(VOLUMEMESHER_UNPARTITIONED, chi_data_types::Varying(6));

RegisterLuaConstant(ExtruderTemplateType, SURFACE_MESH, chi_data_types::Varying(1));
RegisterLuaConstant(ExtruderTemplateType, UNPARTITIONED_MESH, chi_data_types::Varying(2));

// #############################################################################
/** Creates a new volume mesher.
 *
\param Type int Volume Remesher type.
\param OtherArgs varying Additional arguments depending on `Type`.

Remesher types:\n
 VOLUMEMESHER_EXTRUDER = Creates an extruded mesh from a 2D template mesh.
 Requires two additional arguments. `TemplateType` and `handle`. See below.\n
 VOLUMEMESHER_UNPARTITIONED = Create the mesh from the latest UnpartitionedMesh.
 Requires a single additional argument, `handle`, which is a handle to
 a valid unpartitioned mesh.\n

##_

###Extruder parameters

When the mesher type is specified to be VOLUMEMESHER_EXTRUDER then two
additional arguments are required. `TemplateType` and `handle`.\n

- `TemplateType` can for now only be
 `ExtruderTemplateType.UNPARTITIONED_MESH`.\n
- `handle` is a handle to the template mesh. When `TemplateType` is
  set to `ExtruderTemplateType.UNPARTITIONED_MESH` then the handle must point
  to a valid unpartitioned mesh.

\ingroup LuaVolumeMesher
\author Jan*/
int
chiVolumeMesherCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  //============================================= Arguments check
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Mesher type
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
