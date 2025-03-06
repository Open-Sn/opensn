// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/runtime.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/mesh/io/mesh_io.h"
#include "lua/lib/console.h"

using namespace opensn;

namespace opensnlua
{

void
MeshSetBlockIDFromFunction(std::shared_ptr<MeshContinuum> grid, const char* lua_fname)
{
  auto L = Console::GetInstance().GetConsoleState();
  auto lua_fn = luabridge::getGlobal(L, lua_fname);

  int local_num_cells_modified = 0;
  for (auto& cell : grid->local_cells)
  {
    int new_blkid = lua_fn(cell.centroid, cell.block_id)[0];
    if (cell.block_id != new_blkid)
    {
      cell.block_id = new_blkid;
      ++local_num_cells_modified;
    }
  } // for local cell

  const auto& ghost_ids = grid->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid->cells[ghost_id];
    int new_blkid = lua_fn(cell.centroid, cell.block_id)[0];

    if (cell.block_id != new_blkid)
    {
      cell.block_id = new_blkid;
      ++local_num_cells_modified;
    }
  } // for ghost cell id

  int global_num_cells_modified;
  mpi_comm.all_reduce(local_num_cells_modified, global_num_cells_modified, mpi::op::sum<int>());

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Done setting material id from lua function. "
                             << "Number of cells modified = " << global_num_cells_modified << ".";
}

void
MeshSetBoundaryIDFromFunction(std::shared_ptr<MeshContinuum> grid, const char* lua_fname)
{
  OpenSnLogicalErrorIf(opensn::mpi_comm.size() != 1, "Can for now only be used in serial.");

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting boundary id from lua function.";

  auto L = Console::GetInstance().GetConsoleState();
  auto lua_fn = luabridge::getGlobal(L, lua_fname);

  // Check if name already has id
  auto& grid_boundary_id_map = grid->GetBoundaryIDMap();

  int local_num_faces_modified = 0;
  for (auto& cell : grid->local_cells)
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        std::string boundary_name = lua_fn(face.centroid, face.normal, face.neighbor_id)[0];

        const uint64_t boundary_id = grid->MakeBoundaryID(boundary_name);

        if (face.neighbor_id != boundary_id)
        {
          face.neighbor_id = boundary_id;
          ++local_num_faces_modified;

          if (grid_boundary_id_map.count(boundary_id) == 0)
            grid_boundary_id_map[boundary_id] = boundary_name;
        }
      }

  const auto& ghost_ids = grid->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid->cells[ghost_id];
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        std::string boundary_name = lua_fn(face.centroid, face.normal, face.neighbor_id)[0];
        const uint64_t boundary_id = grid->MakeBoundaryID(boundary_name);

        if (face.neighbor_id != boundary_id)
        {
          face.neighbor_id = boundary_id;
          ++local_num_faces_modified;

          if (grid_boundary_id_map.count(boundary_id) == 0)
            grid_boundary_id_map[boundary_id] = boundary_name;
        }
      }
  }

  int global_num_faces_modified;
  mpi_comm.all_reduce(local_num_faces_modified, global_num_faces_modified, mpi::op::sum<int>());

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Done setting boundary id from lua function. "
                             << "Number of cells modified = " << global_num_faces_modified << ".";
}

void
MeshExportToPVTU(std::shared_ptr<opensn::MeshContinuum> grid, const std::string& file_name)
{
  opensn::MeshIO::ToPVTU(grid, file_name);
}

} // namespace opensnlua
