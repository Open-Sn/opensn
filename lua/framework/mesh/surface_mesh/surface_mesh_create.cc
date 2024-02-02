#include "framework/lua.h"

#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/runtime.h"

#include "framework/logging/log.h"
#include "lua_surface_mesh.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(SurfaceMeshCreate);

int
SurfaceMeshCreate(lua_State* L)
{
  auto new_mesh = new SurfaceMesh;

  opensn::surface_mesh_stack.emplace_back(new_mesh);

  size_t index = opensn::surface_mesh_stack.size() - 1;
  lua_pushnumber(L, static_cast<lua_Number>(index));

  opensn::log.LogAllVerbose2() << "SurfaceMeshCreate: "
                                  "Empty SurfaceMesh object, "
                               << index << ", created" << std::endl;

  return 1;
}
