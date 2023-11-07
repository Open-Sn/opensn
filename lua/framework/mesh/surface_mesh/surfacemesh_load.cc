#include "framework/lua.h"

#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua_surface_mesh.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiSurfaceMeshImportFromOBJFile);
RegisterLuaFunctionAsIs(chiSurfaceMeshImportFromTriangleFiles);

int
chiSurfaceMeshImportFromOBJFile(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const std::string file_name = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args >= 3) as_poly = lua_toboolean(L, 3);

  auto& surface_mesh =
    Chi::GetStackItem<chi_mesh::SurfaceMesh>(Chi::surface_mesh_stack, handle, __FUNCTION__);

  Chi::log.Log0Verbose2() << fname << ": Loading Wavefront .obj file: " << std::endl;

  // Transform if necessary
  chi_mesh::Vector3 Tvec(0.0, 0.0, 0.0);
  if (num_args == 4)
  {
    std::vector<double> T;
    LuaPopulateVectorFrom1DArray(fname, L, 4, T);
    if (T.size() != 3) throw std::invalid_argument(fname + ": Argument 4. Table length not 3.");
    Tvec = chi_mesh::Vector3(T[0], T[1], T[2]);
    Chi::log.Log0Verbose2() << "Transform vector: " << Tvec.PrintStr();
  }

  surface_mesh.ImportFromOBJFile(file_name, as_poly, Tvec);

  return 1;
}

int
chiSurfaceMeshImportFromTriangleFiles(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  // Get arguments
  int num_args = lua_gettop(L);
  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args == 3) { as_poly = lua_toboolean(L, 3); }

  auto& surface_mesh =
    Chi::GetStackItem<chi_mesh::SurfaceMesh>(Chi::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ImportFromTriangleFiles(temp, as_poly);

  return 1;
}

int
chiSurfaceMeshImportFromMshFiles(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  // Get arguments
  int num_args = lua_gettop(L);
  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args == 3) { as_poly = lua_toboolean(L, 3); }

  auto& surface_mesh =
    Chi::GetStackItem<chi_mesh::SurfaceMesh>(Chi::surface_mesh_stack, handle, __FUNCTION__);

  std::stringstream outtext;
  outtext << "chiSurfaceMeshImportFromMshFiles: "
             "Loading a gmsh ascii file: ";
  outtext << temp << std::endl;
  Chi::log.LogAllVerbose2() << outtext.str();
  surface_mesh.ImportFromMshFiles(temp, as_poly);

  return 1;
}
