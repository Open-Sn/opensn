#include "chi_volumemesher.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/LogicalVolume/LogicalVolume.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "console/chi_console.h"

#include "utils/chi_timer.h"
#include "chi_lua.h"


//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
  SetMatIDFromLogical(const chi_mesh::LogicalVolume& log_vol,
                      bool sense, int mat_id)
{
  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Setting material id from logical volume.";
  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuumPtr vol_cont = handler.GetGrid();

  int num_cells_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    if (log_vol.Inside(cell.centroid_) && sense){
      cell.material_id_ = mat_id;
      ++num_cells_modified;
    }
  }

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = vol_cont->cells[ghost_id];
    if (log_vol.Inside(cell.centroid_) && sense)
      cell.material_id_ = mat_id;
  }

  int global_num_cells_modified;
  MPI_Allreduce(&num_cells_modified,        //sendbuf
                &global_num_cells_modified, //recvbuf
                1, MPI_INT,                 //count + datatype
                MPI_SUM,                    //operation
                Chi::mpi.comm);            //comm

  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Done setting material id from logical volume. "
    << "Number of cells modified = " << global_num_cells_modified << ".";
}

//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
  SetBndryIDFromLogical(const chi_mesh::LogicalVolume& log_vol,
                        bool sense,
                        const std::string& bndry_name)
{
  Chi::log.Log()
    << Chi::program_timer.GetTimeString()
    << " Setting boundary id from logical volume.";
  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuumPtr vol_cont = handler.GetGrid();

  //============================================= Check if name already has id
  auto& grid_bndry_id_map = vol_cont->GetBoundaryIDMap();
  uint64_t bndry_id = vol_cont->MakeBoundaryID(bndry_name);

  //============================================= Loop over cells
  int num_faces_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_) continue;
      if (log_vol.Inside(face.centroid_) && sense){
        face.neighbor_id_ = bndry_id;
        ++num_faces_modified;
      }
    }
  }

  int global_num_faces_modified;
  MPI_Allreduce(&num_faces_modified,        //sendbuf
                &global_num_faces_modified, //recvbuf
                1, MPI_INT,                 //count + datatype
                MPI_SUM,                    //operation
                Chi::mpi.comm);            //comm


  if (global_num_faces_modified > 0 and
      grid_bndry_id_map.count(bndry_id) == 0)
      grid_bndry_id_map[bndry_id] = bndry_name;

  Chi::log.Log()
    << Chi::program_timer.GetTimeString()
    << " Done setting boundary id from logical volume. "
    << "Number of faces modified = " << global_num_faces_modified << ".";
}

//###################################################################
/**Sets material id's for all cells to the specified material id.*/
void chi_mesh::VolumeMesher::SetMatIDToAll(int mat_id)
{
  Chi::log.Log()
    << Chi::program_timer.GetTimeString()
    << " Setting material id " << mat_id << " to all cells.";

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  auto vol_cont = handler.GetGrid();

  for (auto& cell : vol_cont->local_cells)
    cell.material_id_ = mat_id;

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    vol_cont->cells[ghost_id].material_id_ = mat_id;

  Chi::mpi.Barrier();
  Chi::log.Log()
    << Chi::program_timer.GetTimeString()
    << " Done setting material id " << mat_id << " to all cells";
}


//###################################################################
/**Sets material id's using a lua function. The lua function is called
with for each cell with 4 arguments, the cell's centroid x,y,z values
and the cell's current material id.

The lua function's prototype should be:
\code
function LuaFuncName(x,y,z,id)
 --stuff
end
\endcode
*/
void chi_mesh::VolumeMesher::
  SetMatIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "chi_mesh::VolumeMesher::SetMatIDFromLuaFunction";

  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Setting material id from lua function.";

  //============================================= Define console call
  auto L = Chi::console.GetConsoleState();
  auto CallLuaXYZFunction = [&L,&lua_fname,&fname](const chi_mesh::Cell& cell)
  {
    //============= Load lua function
    lua_getglobal(L, lua_fname.c_str());

    //============= Error check lua function
    if (not lua_isfunction(L, -1))
      throw std::logic_error(fname + " attempted to access lua-function, " +
                               lua_fname + ", but it seems the function"
                                             " could not be retrieved.");

    const auto& xyz = cell.centroid_;

    //============= Push arguments
    lua_pushnumber(L, xyz.x);
    lua_pushnumber(L, xyz.y);
    lua_pushnumber(L, xyz.z);
    lua_pushinteger(L, cell.material_id_);

    //============= Call lua function
    //4 arguments, 1 result (double), 0=original error object
    int lua_return;
    if (lua_pcall(L,4,1,0) == 0)
    {
      LuaCheckNumberValue(fname, L, -1);
      lua_return = lua_tointeger(L,-1);
    }
    else
      throw std::logic_error(fname + " attempted to call lua-function, " +
                             lua_fname + ", but the call failed.");

    lua_pop(L,1); //pop the int, or error code

    return lua_return;
  };

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuum& grid = *handler.GetGrid();

  int local_num_cells_modified = 0;
  for (auto& cell : grid.local_cells)
  {
    int new_matid = CallLuaXYZFunction(cell);

    if (cell.material_id_ != new_matid)
    {
      cell.material_id_ = new_matid;
      ++local_num_cells_modified;
    }
  }//for local cell

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
  }//for ghost cell id

  int globl_num_cells_modified;
  MPI_Allreduce(&local_num_cells_modified, //sendbuf
                &globl_num_cells_modified, //recvbuf
                1, MPI_INT,                //count+datatype
                MPI_SUM,                   //operation
                Chi::mpi.comm);           //comm

  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Done setting material id from lua function. "
    << "Number of cells modified = " << globl_num_cells_modified << ".";
}


//###################################################################
/**Sets boundary id's using a lua function. The lua function is called
for each boundary face with 7 arguments, the face's centroid x,y,z values,
the face's normal x,y,z values and the face's current boundary id. The function
must return a new_bndry_name (string).

The lua function's prototype should be:
\code
function LuaFuncName(x,y,z,nx,ny,nz,id)
 --stuff
end
\endcode
*/
void chi_mesh::VolumeMesher::
  SetBndryIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "chi_mesh::VolumeMesher::SetBndryIDFromLuaFunction";

  if (Chi::mpi.process_count != 1)
    throw std::logic_error(fname + ": Can for now only be used in serial.");

  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Setting boundary id from lua function.";

  //============================================= Define console call
  auto L = Chi::console.GetConsoleState();
  auto CallLuaXYZFunction = [&L,&lua_fname,&fname]
    (const chi_mesh::CellFace& face)
  {
    //============= Load lua function
    lua_getglobal(L, lua_fname.c_str());

    //============= Error check lua function
    if (not lua_isfunction(L, -1))
      throw std::logic_error(fname + " attempted to access lua-function, " +
                             lua_fname + ", but it seems the function"
                                         " could not be retrieved.");

    const auto& xyz = face.centroid_;
    const auto& n   = face.normal_;

    //============= Push arguments
    lua_pushnumber(L, xyz.x);
    lua_pushnumber(L, xyz.y);
    lua_pushnumber(L, xyz.z);
    lua_pushnumber(L, n.x);
    lua_pushnumber(L, n.y);
    lua_pushnumber(L, n.z);
    lua_pushinteger(L, static_cast<lua_Integer>(face.neighbor_id_));

    //============= Call lua function
    //7 arguments, 1 result (string), 0=original error object
    std::string lua_return_bname;
    if (lua_pcall(L,7,1,0) == 0)
    {
      LuaCheckNumberValue(fname, L, -1);
      LuaCheckStringValue(fname, L, -2);
      lua_return_bname = lua_tostring(L,-1);
    }
    else
      throw std::logic_error(fname + " attempted to call lua-function, " +
                             lua_fname + ", but the call failed.");

    lua_pop(L,1); //pop the string, or error code

    return lua_return_bname;
  };

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuum& grid = *handler.GetGrid();

  //============================================= Check if name already has id
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
      }//for bndry face

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
      }//for bndry face
  }//for ghost cell id

  int globl_num_faces_modified;
  MPI_Allreduce(&local_num_faces_modified, //sendbuf
                &globl_num_faces_modified, //recvbuf
                1, MPI_INT,                //count+datatype
                MPI_SUM,                   //operation
                Chi::mpi.comm);           //comm

  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Done setting boundary id from lua function. "
    << "Number of cells modified = " << globl_num_faces_modified << ".";
}