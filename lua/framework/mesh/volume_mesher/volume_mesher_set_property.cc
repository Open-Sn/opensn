#include "framework/lua.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "framework/mesh/logical_volume/logical_volume.h"

#include <iostream>
#include "volume_mesher_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(VolumeMesherSetProperty);

RegisterLuaConstantAsIs(FORCE_POLYGONS, Varying(1));
RegisterLuaConstantAsIs(MESH_GLOBAL, Varying(2));
RegisterLuaConstantAsIs(PARTITION_Z, Varying(3));
RegisterLuaConstantAsIs(VOLUMEPARTITION_Y, Varying(4));
RegisterLuaConstantAsIs(VOLUMEPARTITION_X, Varying(5));
RegisterLuaConstantAsIs(CUTS_Z, Varying(6));
RegisterLuaConstantAsIs(CUTS_Y, Varying(7));
RegisterLuaConstantAsIs(CUTS_X, Varying(8));
RegisterLuaConstantAsIs(PARTITION_TYPE, Varying(9));
RegisterLuaConstantAsIs(KBA_STYLE_XYZ, Varying(2));
RegisterLuaConstantAsIs(PARMETIS, Varying(3));
RegisterLuaConstantAsIs(EXTRUSION_LAYER, Varying(10));
RegisterLuaConstantAsIs(MATID_FROMLOGICAL, Varying(11));
RegisterLuaConstantAsIs(BNDRYID_FROMLOGICAL, Varying(12));
RegisterLuaConstantAsIs(MATID_FROM_LUA_FUNCTION, Varying(13));
RegisterLuaConstantAsIs(BNDRYID_FROM_LUA_FUNCTION, Varying(14));

RegisterLuaFunctionAsIs(VolumeMesherSetKBAPartitioningPxPyPz);
RegisterLuaFunctionAsIs(VolumeMesherSetKBACutsZ);

using namespace opensn;

namespace
{

/**
 * Sets material id's using a lua function. The lua function is called with for each cell with 4
 * arguments, the cell's centroid x,y,z values and the cell's current material id.
 *
 * The lua function's prototype should be:
 * \code
 * function LuaFuncName(x,y,z,id)
 *   --stuff
 * end
 * \endcode
 */
void
SetMatIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "VolumeMesher::SetMatIDFromLuaFunction";

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

  // Get current mesh handler
  auto& handler = GetCurrentHandler();

  // Get back mesh
  MeshContinuum& grid = *handler.GetGrid();

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

/**Sets boundary id's using a lua function. The lua function is called for each boundary face
 * with 7 arguments, the face's centroid x,y,z values, the face's normal x,y,z values and the
 * face's current boundary id. The function must return a new_bndry_name (string).
 *
 * The lua function's prototype should be:
 * \code
 * function LuaFuncName(x,y,z,nx,ny,nz,id)
 * --stuff
 * end
 * \endcode
 */
void
SetBndryIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "VolumeMesher::SetBndryIDFromLuaFunction";

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

  // Get current mesh handler
  auto& handler = GetCurrentHandler();

  // Get back mesh
  MeshContinuum& grid = *handler.GetGrid();

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

} // namespace

int
VolumeMesherSetProperty(lua_State* L)
{
  const std::string fname = "VolumeMesherSetProperty";
  // Get current mesh handler
  auto& cur_hndlr = GetCurrentHandler();
  auto& volume_mesher = cur_hndlr.GetVolumeMesher();

  // Get property index
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  int property_index = lua_tonumber(L, 1);

  typedef VolumeMesherProperty VMP;

  // Selects property
  if (property_index == VMP::FORCE_POLYGONS)
  {
    bool force_condition = lua_toboolean(L, 2);
    volume_mesher.options.force_polygons = force_condition;
  }

  else if (property_index == VMP::MESH_GLOBAL)
  {
    bool mesh_global = lua_toboolean(L, 2);
    volume_mesher.options.mesh_global = mesh_global;
  }

  else if (property_index == VMP::PARTITION_Z)
  {
    int pz = lua_tonumber(L, 2);
    volume_mesher.options.partition_z = pz;
    opensn::log.LogAllVerbose1() << "Partition z set to " << pz;
  }
  else if (property_index == VMP::PARTITION_Y)
  {
    int p = lua_tonumber(L, 2);
    volume_mesher.options.partition_y = p;
    opensn::log.LogAllVerbose1() << "Partition y set to " << p;
  }
  else if (property_index == VMP::PARTITION_X)
  {
    int p = lua_tonumber(L, 2);
    volume_mesher.options.partition_x = p;
    opensn::log.LogAllVerbose1() << "Partition x set to " << p;
  }
  else if (property_index == VMP::CUTS_Z)
  {
    double p = lua_tonumber(L, 2);
    volume_mesher.options.zcuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_Y)
  {
    double p = lua_tonumber(L, 2);
    volume_mesher.options.ycuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_X)
  {
    double p = lua_tonumber(L, 2);
    volume_mesher.options.xcuts.push_back(p);
  }
  else if (property_index == VMP::PARTITION_TYPE)
  {
    int p = lua_tonumber(L, 2);
    if (p >= VolumeMesher::PartitionType::KBA_STYLE_XYZ and
        p <= VolumeMesher::PartitionType::PARMETIS)
      volume_mesher.options.partition_type = (VolumeMesher::PartitionType)p;
    else
    {
      opensn::log.LogAllError() << "Unsupported partition type used in call to " << fname << ".";
      opensn::Exit(EXIT_FAILURE);
    }
  }

  else if (property_index == VMP::MATID_FROMLOGICAL)
  {
    if (not((num_args == 3) or (num_args == 4)))
    {
      opensn::log.LogAllError() << "Invalid amount of arguments used for "
                                   "VolumeMesherSetProperty("
                                   "MATID_FROMLOGICAL...";
      opensn::Exit(EXIT_FAILURE);
    }
    int volume_hndl = lua_tonumber(L, 2);
    int mat_id = lua_tonumber(L, 3);
    int sense = true;
    if (num_args == 4)
      sense = lua_toboolean(L, 4);

    const auto& log_vol =
      opensn::GetStackItem<LogicalVolume>(opensn::object_stack, volume_hndl, fname);

    VolumeMesher::SetMatIDFromLogical(log_vol, sense, mat_id);
  }

  else if (property_index == VMP::BNDRYID_FROMLOGICAL)
  {
    if (not((num_args == 3) or (num_args == 4)))
    {
      opensn::log.LogAllError() << "Invalid amount of arguments used for "
                                   "VolumeMesherSetProperty("
                                   "BNDRYID_FROMLOGICAL...";
      opensn::Exit(EXIT_FAILURE);
    }
    LuaCheckNilValue(fname, L, 2);
    LuaCheckStringValue(fname, L, 3);
    int volume_hndl = lua_tonumber(L, 2);
    std::string bndry_name = lua_tostring(L, 3);
    int sense = true;
    if (num_args == 4)
      sense = lua_toboolean(L, 4);

    ChiLogicalErrorIf(bndry_name.empty(), "argument 3 must not be an empty string.");

    const auto& log_vol =
      opensn::GetStackItem<LogicalVolume>(opensn::object_stack, volume_hndl, fname);

    VolumeMesher::SetBndryIDFromLogical(log_vol, sense, bndry_name);
  }
  else if (property_index == VMP::MATID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    SetMatIDFromLuaFunction(lua_fname);
  }
  else if (property_index == VMP::BNDRYID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    SetBndryIDFromLuaFunction(lua_fname);
  }
  else
  {
    opensn::log.LogAllError() << "Invalid property specified " << property_index
                              << " in call to VolumeMesherSetProperty().";
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}

int
VolumeMesherSetKBAPartitioningPxPyPz(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError(__FUNCTION__, 3, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);
  LuaCheckNilValue(__FUNCTION__, L, 3);

  // Get current mesh handler
  auto& cur_hndlr = GetCurrentHandler();
  auto& vol_mesher = cur_hndlr.GetVolumeMesher();

  int px = lua_tonumber(L, 1);
  int py = lua_tonumber(L, 2);
  int pz = lua_tonumber(L, 3);

  vol_mesher.options.partition_x = px;
  vol_mesher.options.partition_y = py;
  vol_mesher.options.partition_z = pz;

  return 0;
}

int
VolumeMesherSetKBACutsZ(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  std::vector<double> cuts;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 1, cuts);

  auto& mesh_handler = GetCurrentHandler();
  mesh_handler.GetVolumeMesher().options.zcuts = cuts;

  return 0;
}
