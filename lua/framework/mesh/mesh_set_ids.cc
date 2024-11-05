// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/mesh/mesh_set_ids.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(MeshSetUniformMaterialID, mesh, SetUniformMaterialID);
RegisterLuaFunctionInNamespace(MeshSetMaterialIDFromLogicalVolume,
                               mesh,
                               SetMaterialIDFromLogicalVolume);
RegisterLuaFunctionInNamespace(MeshSetBoundaryIDFromLogicalVolume,
                               mesh,
                               SetBoundaryIDFromLogicalVolume);
RegisterLuaFunctionInNamespace(MeshSetMaterialIDFromLuaFunction, mesh, SetMaterialIDFromFunction);
RegisterLuaFunctionInNamespace(MeshSetBoundaryIDFromLuaFunction, mesh, SetBoundaryIDFromFunction);

using namespace opensn;

int
MeshSetUniformMaterialID(lua_State* L)
{
  const std::string fname = "mesh.SetUniformMaterialID";
  LuaCheckArgs<int>(L, fname);

  auto mat_id = LuaArg<int>(L, 1);

  auto vol_cont = GetCurrentMesh();
  vol_cont->SetUniformMaterialID(mat_id);
  mpi_comm.barrier();
  opensn::log.Log() << program_timer.GetTimeString() << " Done setting material id " << mat_id
                    << " to all cells";

  return LuaReturn(L);
}

int
MeshSetMaterialIDFromLogicalVolume(lua_State* L)
{
  const std::string fname = "mesh.SetMaterialIDFromLogicalVolume";
  LuaCheckArgs<size_t, int>(L, fname);

  auto volume_handle = LuaArg<size_t>(L, 1);
  auto mat_id = LuaArg<int>(L, 2);
  auto sense = LuaArgOptional<bool>(L, 3, true);

  const auto& lv = opensn::GetStackItem<LogicalVolume>(opensn::object_stack, volume_handle, fname);

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting material id from logical volume.";
  std::shared_ptr<MeshContinuum> mesh = GetCurrentMesh();
  mesh->SetMaterialIDFromLogical(lv, sense, mat_id);

  return LuaReturn(L);
}

int
MeshSetMaterialIDFromLuaFunction(lua_State* L)
{
  const std::string fname = "mesh.SetMaterialIDFromLuaFunction";
  LuaCheckArgs<std::string>(L, fname);

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting material id from lua function.";

  const auto lua_fname = LuaArg<std::string>(L, 1);

  // Get back mesh
  MeshContinuum& grid = *GetCurrentMesh();

  int local_num_cells_modified = 0;
  for (auto& cell : grid.local_cells)
  {
    int new_matid = LuaCall<int>(L, lua_fname, cell.centroid, cell.material_id);
    if (cell.material_id != new_matid)
    {
      cell.material_id = new_matid;
      ++local_num_cells_modified;
    }
  } // for local cell

  const auto& ghost_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid.cells[ghost_id];
    int new_matid = LuaCall<int>(L, lua_fname, cell.centroid, cell.material_id);

    if (cell.material_id != new_matid)
    {
      cell.material_id = new_matid;
      ++local_num_cells_modified;
    }
  } // for ghost cell id

  int global_num_cells_modified;
  mpi_comm.all_reduce(local_num_cells_modified, global_num_cells_modified, mpi::op::sum<int>());

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Done setting material id from lua function. "
                             << "Number of cells modified = " << global_num_cells_modified << ".";

  return LuaReturn(L);
}

int
MeshSetBoundaryIDFromLuaFunction(lua_State* L)
{
  const std::string fname = "mesh.SetBoundaryIDFromFunction";
  LuaCheckArgs<std::string>(L, fname);

  OpenSnLogicalErrorIf(opensn::mpi_comm.size() != 1, "Can for now only be used in serial.");

  const auto lua_fname = LuaArg<std::string>(L, 1);

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting boundary id from lua function.";

  MeshContinuum& grid = *GetCurrentMesh();

  // Check if name already has id
  auto& grid_boundary_id_map = grid.BoundaryIDMap();

  int local_num_faces_modified = 0;
  for (auto& cell : grid.local_cells)
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        auto boundary_name =
          LuaCall<std::string>(L, lua_fname, face.centroid, face.normal, face.neighbor_id);

        const uint64_t boundary_id = grid.MakeBoundaryID(boundary_name);

        if (face.neighbor_id != boundary_id)
        {
          face.neighbor_id = boundary_id;
          ++local_num_faces_modified;

          if (grid_boundary_id_map.count(boundary_id) == 0)
            grid_boundary_id_map[boundary_id] = boundary_name;
        }
      }

  const auto& ghost_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid.cells[ghost_id];
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        auto boundary_name =
          LuaCall<std::string>(L, lua_fname, face.centroid, face.normal, face.neighbor_id);
        const uint64_t boundary_id = grid.MakeBoundaryID(boundary_name);

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

  return LuaReturn(L);
}

int
MeshSetBoundaryIDFromLogicalVolume(lua_State* L)
{
  const std::string fname = "mesh.SetBoundaryIDFromLogicalVolume";
  LuaCheckArgs<size_t, std::string>(L, fname);

  auto volume_handle = LuaArg<size_t>(L, 1);
  auto boundary_name = LuaArg<std::string>(L, 2);
  auto sense = LuaArgOptional<bool>(L, 3, true);

  OpenSnLogicalErrorIf(boundary_name.empty(), "argument 2 must not be an empty string.");

  const auto& log_vol =
    opensn::GetStackItem<LogicalVolume>(opensn::object_stack, volume_handle, fname);

  opensn::log.Log() << program_timer.GetTimeString() << " Setting boundary id from logical volume.";
  std::shared_ptr<MeshContinuum> mesh = GetCurrentMesh();
  mesh->SetBoundaryIDFromLogical(log_vol, sense, boundary_name);

  return LuaReturn(L);
}

} // namespace opensnlua
