#include "framework/lua.h"

#include "framework/mesh/surface_mesh/surface_mesh.h"

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

  opensn::object_stack.emplace_back(new_mesh);

  size_t index = opensn::object_stack.size() - 1;
  lua_pushinteger(L, static_cast<lua_Integer>(index));

  opensn::log.LogAllVerbose2() << "SurfaceMeshCreate: "
                                  "Empty SurfaceMesh object, "
                               << index << ", created" << std::endl;

  return 1;
}
