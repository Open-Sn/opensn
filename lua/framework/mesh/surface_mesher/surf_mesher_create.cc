#include "framework/lua.h"
#include <iostream>
#include "framework/mesh/surface_mesher/predefined/surfmesher_predefined.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "surf_mesher_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiSurfaceMesherCreate);
RegisterLuaConstantAsIs(SURFACEMESHER_PREDEFINED, Varying(1));

int
chiSurfaceMesherCreate(lua_State* L)
{
  auto& cur_hndlr = opensn::GetCurrentHandler();

  // Get argument
  LuaCheckNilValue("chiSurfaceMesherCreate", L, 1);
  int type = lua_tonumber(L, 1);

  // Create the surface mesher
  std::shared_ptr<opensn::SurfaceMesher> new_mesher = nullptr;
  if (type == (int)opensn::SurfaceMesherType::Predefined)
  {
    new_mesher = std::make_shared<opensn::SurfaceMesherPredefined>();
  }
  else
  {
    std::cerr << "ERROR: Illegal surface mesher specified"
                 "in chiSurfaceMesherCreate"
              << std::endl;
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  cur_hndlr.SetSurfaceMesher(new_mesher);

  opensn::Chi::log.LogAllVerbose2()
    << "chiSurfaceMesherCreate: Surface remesher created." << std::endl;

  return 0;
}
