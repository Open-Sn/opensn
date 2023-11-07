#include "framework/lua.h"

#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/runtime.h"

#include "framework/logging/log.h"
#include "lua_surface_mesh.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiSurfaceMeshCreate);

int
chiSurfaceMeshCreate(lua_State* L)
{
  auto new_mesh = new chi_mesh::SurfaceMesh;

  Chi::surface_mesh_stack.emplace_back(new_mesh);

  size_t index = Chi::surface_mesh_stack.size() - 1;
  lua_pushnumber(L, static_cast<lua_Number>(index));

  Chi::log.LogAllVerbose2() << "chiSurfaceMeshCreate: "
                               "Empty SurfaceMesh object, "
                            << index << ", created" << std::endl;

  return 1;
}
