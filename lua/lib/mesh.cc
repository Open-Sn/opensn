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
MeshSetUniformMaterialID(int mat_id)
{
  auto mesh = GetCurrentMesh();
  mesh->SetUniformMaterialID(mat_id);
  mpi_comm.barrier();
  opensn::log.Log() << program_timer.GetTimeString() << " Done setting material id " << mat_id
                    << " to all cells";
}

void
MeshSetMaterialIDFromLogicalVolume(std::shared_ptr<LogicalVolume> lv, int mat_id, bool sense)
{
  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting material id from logical volume.";
  std::shared_ptr<MeshContinuum> mesh = GetCurrentMesh();
  mesh->SetMaterialIDFromLogical(*lv.get(), sense, mat_id);
}

void
MeshSetBoundaryIDFromLogicalVolume(std::shared_ptr<LogicalVolume> lv,
                                   const std::string& boundary_name,
                                   bool sense)
{
  OpenSnLogicalErrorIf(boundary_name.empty(), "argument 2 must not be an empty string.");

  opensn::log.Log() << program_timer.GetTimeString() << " Setting boundary id from logical volume.";
  std::shared_ptr<MeshContinuum> mesh = GetCurrentMesh();
  mesh->SetBoundaryIDFromLogical(*lv.get(), boundary_name, sense);
}

void
MeshSetupOrthogonalBoundaries()
{
  opensn::log.Log() << program_timer.GetTimeString() << " Setting orthogonal boundaries.";

  auto vol_cont = GetCurrentMesh();

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);
  const Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : vol_cont->local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (not face.has_neighbor)
      {
        Vector3& n = face.normal;

        std::string boundary_name;
        if (n.Dot(ihat) < -0.999)
          boundary_name = "XMIN";
        else if (n.Dot(ihat) > 0.999)
          boundary_name = "XMAX";
        else if (n.Dot(jhat) < -0.999)
          boundary_name = "YMIN";
        else if (n.Dot(jhat) > 0.999)
          boundary_name = "YMAX";
        else if (n.Dot(khat) < -0.999)
          boundary_name = "ZMIN";
        else if (n.Dot(khat) > 0.999)
          boundary_name = "ZMAX";

        uint64_t bndry_id = vol_cont->MakeBoundaryID(boundary_name);

        face.neighbor_id = bndry_id;

        vol_cont->GetBoundaryIDMap()[bndry_id] = boundary_name;
      }
    }
  }

  opensn::mpi_comm.barrier();
  opensn::log.Log() << program_timer.GetTimeString() << " Done setting orthogonal boundaries.";
}

void
MeshSetMaterialIDFromFunction(const char* lua_fname)
{
  // Get back mesh
  MeshContinuum& grid = *GetCurrentMesh();

  auto L = Console::GetInstance().GetConsoleState();
  auto lua_fn = luabridge::getGlobal(L, lua_fname);

  int local_num_cells_modified = 0;
  for (auto& cell : grid.local_cells)
  {
    int new_matid = lua_fn(cell.centroid, cell.material_id)[0];
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
    int new_matid = lua_fn(cell.centroid, cell.material_id)[0];

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
}

void
MeshSetBoundaryIDFromFunction(const char* lua_fname)
{
  OpenSnLogicalErrorIf(opensn::mpi_comm.size() != 1, "Can for now only be used in serial.");

  opensn::log.Log0Verbose1() << program_timer.GetTimeString()
                             << " Setting boundary id from lua function.";

  auto L = Console::GetInstance().GetConsoleState();
  auto lua_fn = luabridge::getGlobal(L, lua_fname);

  MeshContinuum& grid = *GetCurrentMesh();

  // Check if name already has id
  auto& grid_boundary_id_map = grid.GetBoundaryIDMap();

  int local_num_faces_modified = 0;
  for (auto& cell : grid.local_cells)
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        std::string boundary_name = lua_fn(face.centroid, face.normal, face.neighbor_id)[0];

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
        std::string boundary_name = lua_fn(face.centroid, face.normal, face.neighbor_id)[0];
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
}

void
MeshExportToPVTU(const std::string& file_name)
{
  auto grid = opensn::GetCurrentMesh();
  opensn::MeshIO::ToPVTU(grid, file_name);
}

void
MeshComputeVolumePerMaterialID()
{
  auto curr_mesh = GetCurrentMesh();
  curr_mesh->ComputeVolumePerMaterialID();
}

} // namespace opensnlua
