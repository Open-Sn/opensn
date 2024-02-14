#include "lua/framework/mesh/lua_mesh_set_ids.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/console/console.h"
#include "framework/mesh/logical_volume/logical_volume.h"

RegisterLuaFunctionNamespace(MeshSetUniformMaterialID, mesh, SetUniformMaterialID);
RegisterLuaFunctionNamespace(MeshSetMaterialIDFromLogicalVolume,
                             mesh,
                             SetMaterialIDFromLogicalVolume);

using namespace opensn;

int
MeshSetUniformMaterialID(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int mat_id = lua_tonumber(L, 1);

  auto vol_cont = GetCurrentMesh();
  vol_cont->SetUniformMaterialID(mat_id);
  mpi_comm.barrier();
  opensn::log.Log() << program_timer.GetTimeString() << " Done setting material id " << mat_id
                    << " to all cells";

  return 0;
}

int
MeshSetMaterialIDFromLogicalVolume(lua_State* L)
{
  const std::string fname = "mesh.SetMaterialIDFromLogicalVolume";

  int num_args = lua_gettop(L);
  if ((num_args == 2) or (num_args == 3))
  {
    int volume_handle = lua_tonumber(L, 1);
    int mat_id = lua_tonumber(L, 2);
    int sense = true;
    if (num_args == 3)
      sense = lua_toboolean(L, 3);

    const auto& lv =
      opensn::GetStackItem<LogicalVolume>(opensn::object_stack, volume_handle, fname);

    opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                               << " Setting material id from logical volume.";
    std::shared_ptr<MeshContinuum> mesh = GetCurrentMesh();
    mesh->SetMaterialIDFromLogical(lv, sense, mat_id);
  }
  else
  {
    opensn::log.LogAllError()
      << "Invalid number of arguments when calling 'mesh.SetMaterialIDFromLogicalVolume'";
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}

void
SetMatIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "SetMatIDFromLuaFunction";

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting material id from lua function.";

  // Define console call
  auto L = opensnlua::console.GetConsoleState();
  auto CallLuaXYZFunction = [&L, &lua_fname, &fname](const Cell& cell)
  {
    // Load lua function
    lua_getglobal(L, lua_fname.c_str());

    // Error check lua function
    if (not lua_isfunction(L, -1))
      ChiLogicalError("Attempted to access lua-function, " + lua_fname +
                      ", but it seems the function could not be retrieved.");

    const auto& xyz = cell.centroid_;

    // Push arguments
    lua_pushnumber(L, xyz.x);
    lua_pushnumber(L, xyz.y);
    lua_pushnumber(L, xyz.z);
    lua_pushinteger(L, cell.material_id_);

    // Call lua function
    // 4 arguments, 1 result (double), 0=original error object
    int lua_return;
    if (lua_pcall(L, 4, 1, 0) == 0)
    {
      LuaCheckNumberValue(fname, L, -1);
      lua_return = lua_tointeger(L, -1);
    }
    else
      ChiLogicalError("Attempted to call lua-function, " + lua_fname + ", but the call failed.");

    lua_pop(L, 1); // pop the int, or error code

    return lua_return;
  };

  // Get back mesh
  MeshContinuum& grid = *GetCurrentMesh();

  int local_num_cells_modified = 0;
  for (auto& cell : grid.local_cells)
  {
    int new_matid = CallLuaXYZFunction(cell);

    if (cell.material_id_ != new_matid)
    {
      cell.material_id_ = new_matid;
      ++local_num_cells_modified;
    }
  } // for local cell

  const auto& ghost_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid.cells[ghost_id];
    int new_matid = CallLuaXYZFunction(cell);

    if (cell.material_id_ != new_matid)
    {
      cell.material_id_ = new_matid;
      ++local_num_cells_modified;
    }
  } // for ghost cell id

  int globl_num_cells_modified;
  mpi_comm.all_reduce(local_num_cells_modified, globl_num_cells_modified, mpi::op::sum<int>());

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Done setting material id from lua function. "
                             << "Number of cells modified = " << globl_num_cells_modified << ".";
}

void
SetBoundaryIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "SetBoundaryIDFromLuaFunction";

  ChiLogicalErrorIf(opensn::mpi_comm.size() != 1, "Can for now only be used in serial.");

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting boundary id from lua function.";

  // Define console call
  auto L = opensnlua::console.GetConsoleState();
  auto CallLuaXYZFunction = [&L, &lua_fname, &fname](const CellFace& face)
  {
    // Load lua function
    lua_getglobal(L, lua_fname.c_str());

    // Error check lua function
    ChiLogicalErrorIf(not lua_isfunction(L, -1),
                      "Attempted to access lua-function, " + lua_fname +
                        ", but it seems the function could not be retrieved.");

    const auto& xyz = face.centroid_;
    const auto& n = face.normal_;

    // Push arguments
    lua_pushnumber(L, xyz.x);
    lua_pushnumber(L, xyz.y);
    lua_pushnumber(L, xyz.z);
    lua_pushnumber(L, n.x);
    lua_pushnumber(L, n.y);
    lua_pushnumber(L, n.z);
    lua_pushinteger(L, static_cast<lua_Integer>(face.neighbor_id_));

    // Call lua function
    // 7 arguments, 1 result (string), 0=original error object
    std::string lua_return_bname;
    if (lua_pcall(L, 7, 1, 0) == 0)
    {
      LuaCheckNumberValue(fname, L, -1);
      LuaCheckStringValue(fname, L, -2);
      lua_return_bname = lua_tostring(L, -1);
    }
    else
      ChiLogicalError("Attempted to call lua-function, " + lua_fname + ", but the call failed.");

    lua_pop(L, 1); // pop the string, or error code

    return lua_return_bname;
  };

  MeshContinuum& grid = *GetCurrentMesh();

  // Check if name already has id
  auto& grid_bndry_id_map = grid.GetBoundaryIDMap();

  int local_num_faces_modified = 0;
  for (auto& cell : grid.local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        const std::string bndry_name = CallLuaXYZFunction(face);
        const uint64_t bndry_id = grid.MakeBoundaryID(bndry_name);

        if (face.neighbor_id_ != bndry_id)
        {
          face.neighbor_id_ = bndry_id;
          ++local_num_faces_modified;

          if (grid_bndry_id_map.count(bndry_id) == 0)
            grid_bndry_id_map[bndry_id] = bndry_name;
        }
      } // for bndry face

  const auto& ghost_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid.cells[ghost_id];
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        const std::string bndry_name = CallLuaXYZFunction(face);
        const uint64_t bndry_id = grid.MakeBoundaryID(bndry_name);

        if (face.neighbor_id_ != bndry_id)
        {
          face.neighbor_id_ = bndry_id;
          ++local_num_faces_modified;

          if (grid_bndry_id_map.count(bndry_id) == 0)
            grid_bndry_id_map[bndry_id] = bndry_name;
        }
      } // for bndry face
  }     // for ghost cell id

  int globl_num_faces_modified;
  mpi_comm.all_reduce(local_num_faces_modified, globl_num_faces_modified, mpi::op::sum<int>());

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Done setting boundary id from lua function. "
                             << "Number of cells modified = " << globl_num_faces_modified << ".";
}
