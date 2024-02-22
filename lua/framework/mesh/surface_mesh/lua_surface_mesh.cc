#include "lua_surface_mesh.h"
#include "framework/lua.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionNamespace(MeshSurfaceMeshCreate, mesh, SurfaceMeshCreate);
RegisterLuaFunctionNamespace(MeshComputeLoadBalancing, mesh, ComputeLoadBalancing);
RegisterLuaFunctionNamespace(MeshSurfaceMeshImportFromOBJFile, mesh, SurfaceMeshImportFromOBJFile);
RegisterLuaFunctionNamespace(MeshSurfaceMeshImportFromTriangleFiles,
                             mesh,
                             SurfaceMeshImportFromTriangleFiles);

int
MeshSurfaceMeshCreate(lua_State* L)
{
  auto new_mesh = new SurfaceMesh;

  opensn::object_stack.emplace_back(new_mesh);

  size_t index = opensn::object_stack.size() - 1;
  lua_pushinteger(L, static_cast<lua_Integer>(index));

  opensn::log.LogAllVerbose2() << "mesh.SurfaceMeshCreate: Empty SurfaceMesh object, " << index
                               << ", created" << std::endl;

  return 1;
}

int
MeshComputeLoadBalancing(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("mesh.ComputeLoadBalancing", 3, num_args);

  // Get reference surface mesh
  int surf_handle = lua_tonumber(L, 1);

  auto& cur_surf =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, surf_handle, __FUNCTION__);

  // Extract x-cuts
  if (not lua_istable(L, 2))
  {
    opensn::log.LogAllError() << "In call to mesh.ComputeLoadBalancing: expected table for "
                                 "argument 2. Incompatible value supplied.";
    opensn::Exit(EXIT_FAILURE);
  }

  int x_table_len = lua_rawlen(L, 2);

  std::vector<double> x_cuts(x_table_len, 0.0);
  for (int g = 0; g < x_table_len; g++)
  {
    lua_pushnumber(L, g + 1);
    lua_gettable(L, 2);
    x_cuts[g] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  // Extract y-cuts
  if (not lua_istable(L, 3))
  {
    opensn::log.LogAllError() << "In call to mesh.ComputeLoadBalancing: expected table for "
                                 "argument 3. Incompatible value supplied.";
    opensn::Exit(EXIT_FAILURE);
  }

  int y_table_len = lua_rawlen(L, 3);

  std::vector<double> y_cuts(y_table_len, 0.0);
  for (int g = 0; g < y_table_len; g++)
  {
    lua_pushnumber(L, g + 1);
    lua_gettable(L, 3);
    y_cuts[g] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  // Call compute balance
  std::stable_sort(x_cuts.begin(), x_cuts.end());
  std::stable_sort(y_cuts.begin(), y_cuts.end());
  cur_surf.ComputeLoadBalancing(x_cuts, y_cuts);

  return 0;
}

/** Exports all open edges of a surface mesh to file. This is used mostly
 * for graphical error checking.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param FileName char Filename to which the edges are to be exported.
 *
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int
SurfaceMeshExtractOpenEdgesToObj(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("SurfaceMeshExtractOpenEdgesToObj", 2, num_args);

  int surf_handle = lua_tonumber(L, 1);
  const char* file_name = lua_tostring(L, 2);

  auto& surface_mesh =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, surf_handle, __FUNCTION__);

  surface_mesh.ExtractOpenEdgesToObj(file_name);
  return 0;
}

int
MeshSurfaceMeshImportFromOBJFile(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, 2, num_args);

  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const std::string file_name = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args >= 3)
    as_poly = lua_toboolean(L, 3);

  auto& surface_mesh =
    opensn::GetStackItem<SurfaceMesh>(opensn::object_stack, handle, __FUNCTION__);

  opensn::log.Log0Verbose2() << fname << ": Loading Wavefront .obj file: " << std::endl;

  // Transform if necessary
  Vector3 Tvec(0.0, 0.0, 0.0);
  if (num_args == 4)
  {
    std::vector<double> T;
    LuaPopulateVectorFrom1DArray(fname, L, 4, T);
    if (T.size() != 3)
      throw std::invalid_argument(fname + ": Argument 4. Table length not 3.");
    Tvec = Vector3(T[0], T[1], T[2]);
    opensn::log.Log0Verbose2() << "Transform vector: " << Tvec.PrintStr();
  }

  surface_mesh.ImportFromOBJFile(file_name, as_poly, Tvec);

  return 1;
}

int
MeshSurfaceMeshImportFromTriangleFiles(lua_State* L)
{
  // Get arguments
  int num_args = lua_gettop(L);
  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args == 3)
  {
    as_poly = lua_toboolean(L, 3);
  }

  auto& surface_mesh =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ImportFromTriangleFiles(temp, as_poly);

  return 1;
}

int
SurfaceMeshImportFromMshFiles(lua_State* L)
{
  // Get arguments
  int num_args = lua_gettop(L);
  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args == 3)
  {
    as_poly = lua_toboolean(L, 3);
  }

  auto& surface_mesh =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, handle, __FUNCTION__);

  std::stringstream outtext;
  outtext << "SurfaceMeshImportFromMshFiles: "
             "Loading a gmsh ascii file: ";
  outtext << temp << std::endl;
  opensn::log.LogAllVerbose2() << outtext.str();
  surface_mesh.ImportFromMshFiles(temp, as_poly);

  return 1;
}
