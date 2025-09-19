// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_volume/finite_volume_mapping.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log_exceptions.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"

namespace opensn
{

FiniteVolume::FiniteVolume(std::shared_ptr<MeshContinuum> grid)
  : SpatialDiscretization(grid, SpatialDiscretizationType::FINITE_VOLUME)
{
  CreateCellMappings();
  OrderNodes();
}

std::shared_ptr<FiniteVolume>
FiniteVolume::New(std::shared_ptr<MeshContinuum> grid)
{
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : sdm_stack)
    if (sdm->GetType() == SpatialDiscretizationType::FINITE_VOLUME and sdm->GetGrid() == grid and
        sdm->GetGrid()->GetCoordinateSystem() == grid->GetCoordinateSystem())
    {
      auto sdm_ptr = std::dynamic_pointer_cast<FiniteVolume>(sdm);

      OpenSnLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  // If no existing discretization was found then go ahead and make a
  // new one
  auto new_sdm = std::shared_ptr<FiniteVolume>(new FiniteVolume(grid));

  sdm_stack.push_back(new_sdm);

  return new_sdm;
}

void
FiniteVolume::CreateCellMappings()
{
  constexpr std::string_view fname = "SpatialDiscretization_FV::"
                                     "CreateCellMappings";

  auto MakeCellMapping = [this, fname](const Cell& cell)
  {
    using namespace std;
    using namespace opensn;
    std::unique_ptr<CellMapping> mapping;

    switch (cell.GetType())
    {
      case CellType::SLAB:
      case CellType::POLYGON:
      case CellType::POLYHEDRON:
      {
        mapping = make_unique<FiniteVolumeMapping>(
          grid_, cell, cell.centroid, std::vector<std::vector<int>>(cell.faces.size(), {-1}));
        break;
      }
      default:
        throw std::logic_error(std::string(fname) +
                               std::string(": Invalid cell type encountered."));
    }
    return mapping;
  };

  for (const auto& cell : grid_->local_cells)
    cell_mappings_.push_back(MakeCellMapping(cell));

  const auto ghost_ids = grid_->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto ghost_mapping = MakeCellMapping(grid_->cells[ghost_id]);
    nb_cell_mappings_.insert(std::make_pair(ghost_id, std::move(ghost_mapping)));
  }
}

void
FiniteVolume::OrderNodes()
{
  static std::string mapping_error = "FiniteVolume::OrderNodes: Error mapping neighbor cells";

  // Communicate node counts
  const uint64_t local_num_nodes = grid_->local_cells.size();
  mpi_comm.all_gather(local_num_nodes, locJ_block_size_);

  // Build block addresses
  locJ_block_address_.assign(opensn::mpi_comm.size(), 0);
  uint64_t global_num_nodes = 0;
  for (int j = 0; j < opensn::mpi_comm.size(); ++j)
  {
    locJ_block_address_[j] = global_num_nodes;
    global_num_nodes += locJ_block_size_[j];
  }

  local_block_address_ = locJ_block_address_[opensn::mpi_comm.rank()];

  local_base_block_size_ = local_num_nodes;
  global_base_block_size_ = global_num_nodes;

  // Sort neigbor ids
  const auto neighbor_gids = grid_->cells.GetGhostGlobalIDs();
  std::map<uint64_t, std::vector<uint64_t>> sorted_nb_gids;
  for (uint64_t gid : neighbor_gids)
  {
    const auto& cell = grid_->cells[gid];
    sorted_nb_gids[cell.partition_id].push_back(gid);
  }

  // Communicate neighbor ids requiring mapping
  const auto query_nb_gids = MapAllToAll(sorted_nb_gids, mpi_comm);

  // Map the ids
  std::map<uint64_t, std::vector<uint64_t>> mapped_query_nb_gids;
  for (const auto& pid_list_pair : query_nb_gids)
  {
    const uint64_t pid = pid_list_pair.first;
    const auto& gids = pid_list_pair.second;

    auto& local_ids = mapped_query_nb_gids[pid];
    local_ids.reserve(gids.size());
    for (uint64_t gid : gids)
    {
      if (not grid_->IsCellLocal(gid))
        throw std::logic_error(mapping_error);

      const auto& local_cell = grid_->cells[gid];
      local_ids.push_back(local_cell.local_id);
    } // for gid
  } // for pid_list_pair

  // Communicate back the mapped ids
  const auto mapped_nb_gids = MapAllToAll(mapped_query_nb_gids, mpi_comm);

  // Create the neighbor cell mapping
  neighbor_cell_local_ids_.clear();
  for (const auto& pid_list_pair : sorted_nb_gids)
  {
    try
    {
      const auto& pid = pid_list_pair.first;
      const auto& gid_list = pid_list_pair.second;
      const auto& lid_list = mapped_nb_gids.at(pid);

      if (gid_list.size() != lid_list.size())
        throw std::logic_error(mapping_error + std::string(" Size-mismatch."));

      for (size_t i = 0; i < gid_list.size(); ++i)
        neighbor_cell_local_ids_.insert(std::make_pair(gid_list[i], lid_list[i]));
    }
    catch (const std::out_of_range& oor)
    {
      throw std::logic_error(mapping_error + std::string(" OOR."));
    }
  } // for pid_list_pair

  local_base_block_size_ = grid_->local_cells.size();
  global_base_block_size_ = grid_->GetGlobalNumberOfCells();
}

void
FiniteVolume::BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                                   std::vector<int64_t>& nodal_nnz_off_diag,
                                   const UnknownManager& unknown_manager) const
{
  unsigned int num_uk = unknown_manager.unknowns.size();           // Number of unknowns
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize(); // Total number of unknowns

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  const size_t num_local_cells = grid_->local_cells.size();

  nodal_nnz_in_diag.resize(num_local_cells * N, 0.0);
  nodal_nnz_off_diag.resize(num_local_cells * N, 0.0);

  for (int uk = 0; uk < num_uk; ++uk)
  {
    const unsigned int num_comps = unknown_manager.unknowns[uk].num_components;
    for (int comp = 0; comp < num_comps; ++comp)
    {
      for (auto& cell : grid_->local_cells)
      {
        const int64_t i = MapDOFLocal(cell, 0, unknown_manager, uk, comp);

        nodal_nnz_in_diag[i] += 1;

        for (auto& face : cell.faces)
        {
          if (not face.has_neighbor)
            continue;

          if (face.IsNeighborLocal(grid_.get()))
            nodal_nnz_in_diag[i] += 1;
          else
            nodal_nnz_off_diag[i] += 1;
        }
      } // for cell
    } // for components
  } // for unknown
}

int64_t
FiniteVolume::MapDOF(const Cell& cell,
                     const unsigned int,
                     const UnknownManager& unknown_manager,
                     const unsigned int unknown_id,
                     const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type;

  const size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  const size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  const size_t num_local_cells = grid_->local_cells.size();

  if (component >= num_unknowns)
    return -1;

  int64_t address = -1;
  if (cell.partition_id == opensn::mpi_comm.rank())
  {
    if (storage == UnknownStorageType::BLOCK)
      address = static_cast<int64_t>(local_block_address_) * num_unknowns +
                num_local_cells * block_id + cell.local_id;
    else if (storage == UnknownStorageType::NODAL)
      address = static_cast<int64_t>(local_block_address_) * num_unknowns +
                cell.local_id * num_unknowns + block_id;
  }
  else
  {
    const uint64_t ghost_local_id = neighbor_cell_local_ids_.at(cell.global_id);

    if (storage == UnknownStorageType::BLOCK)
      address = static_cast<int64_t>(locJ_block_address_[cell.partition_id]) * num_unknowns +
                locJ_block_size_[cell.partition_id] * block_id + ghost_local_id;
    else if (storage == UnknownStorageType::NODAL)
      address = static_cast<int64_t>(locJ_block_address_[cell.partition_id]) * num_unknowns +
                ghost_local_id * num_unknowns + block_id;
  }

  return address;
}

int64_t
FiniteVolume::MapDOFLocal(const Cell& cell,
                          const unsigned int,
                          const UnknownManager& unknown_manager,
                          const unsigned int unknown_id,
                          const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type;

  const size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  const size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  const size_t num_local_cells = grid_->local_cells.size();

  if (component >= num_unknowns)
    return -1;

  int64_t address = -1;
  if (cell.partition_id == opensn::mpi_comm.rank())
  {
    if (storage == UnknownStorageType::BLOCK)
      address = static_cast<int64_t>(num_local_cells) * block_id + cell.local_id;
    else if (storage == UnknownStorageType::NODAL)
      address = static_cast<int64_t>(cell.local_id) * num_unknowns + block_id;
  }
  else
  {
    const size_t num_local_dofs = GetNumLocalDOFs(unknown_manager);
    const size_t num_ghost_nodes = GetNumGhostDOFs(UNITARY_UNKNOWN_MANAGER);
    const uint64_t ghost_local_id = grid_->cells.GetGhostLocalID(cell.global_id);

    if (storage == UnknownStorageType::BLOCK)
      address = static_cast<int64_t>(num_local_dofs) +
                static_cast<int64_t>(num_ghost_nodes) * block_id + ghost_local_id;
    else if (storage == UnknownStorageType::NODAL)
      address = static_cast<int64_t>(num_local_dofs) +
                num_unknowns * static_cast<int64_t>(ghost_local_id) + block_id;
  }

  return address;
}

size_t
FiniteVolume::GetNumGhostDOFs(const UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return grid_->cells.GhostCellCount() * N;
}

std::vector<int64_t>
FiniteVolume::GetGhostDOFIndices(const UnknownManager& unknown_manager) const
{
  std::vector<int64_t> dof_ids;
  dof_ids.reserve(GetNumGhostDOFs(unknown_manager));

  std::vector<uint64_t> ghost_cell_ids = grid_->cells.GetGhostGlobalIDs();

  const size_t num_uks = unknown_manager.unknowns.size();

  for (const auto cell_id : ghost_cell_ids)
  {
    const auto& cell = grid_->cells[cell_id];
    for (size_t u = 0; u < num_uks; ++u)
    {
      const auto& unkn = unknown_manager.unknowns[u];
      const size_t num_comps = unkn.num_components;
      for (size_t c = 0; c < num_comps; ++c)
      {
        const int64_t dofmap = MapDOF(cell, 0, unknown_manager, u, c);
        dof_ids.push_back(dofmap);
      } // for c
    } // for u
  }

  return dof_ids;
}

} // namespace opensn
