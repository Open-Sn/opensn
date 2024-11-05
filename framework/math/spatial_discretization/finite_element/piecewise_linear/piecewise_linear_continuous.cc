// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/mpi/mpi_utils.h"
#include <algorithm>

namespace opensn
{

PieceWiseLinearContinuous::PieceWiseLinearContinuous(const MeshContinuum& grid,
                                                     QuadratureOrder q_order,
                                                     CoordinateSystemType cs_type)
  : PieceWiseLinearBase(
      grid, q_order, SpatialDiscretizationType::PIECEWISE_LINEAR_CONTINUOUS, cs_type)
{
  CreateCellMappings();

  OrderNodes();
}

std::shared_ptr<PieceWiseLinearContinuous>
PieceWiseLinearContinuous::New(const MeshContinuum& grid,
                               QuadratureOrder q_order,
                               CoordinateSystemType cs_type)

{
  const auto PWLC = SpatialDiscretizationType::PIECEWISE_LINEAR_CONTINUOUS;
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : sdm_stack)
    if (sdm->Type() == PWLC and std::addressof(sdm->Grid()) == std::addressof(grid) and
        sdm->CoordinateSystem() == cs_type)
    {
      auto fe_ptr = std::dynamic_pointer_cast<FiniteElementBase>(sdm);

      OpenSnLogicalErrorIf(not fe_ptr, "Casting failure to FE");

      if (fe_ptr->GetQuadratureOrder() != q_order)
        break;

      auto sdm_ptr = std::dynamic_pointer_cast<PieceWiseLinearContinuous>(fe_ptr);

      OpenSnLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  auto new_sdm = std::shared_ptr<PieceWiseLinearContinuous>(
    new PieceWiseLinearContinuous(grid, q_order, cs_type));

  sdm_stack.push_back(new_sdm);

  return new_sdm;
}

void
PieceWiseLinearContinuous::OrderNodes()
{
  const std::string fname = __FUNCTION__;
  // Build set of local scope nodes
  // ls_node_id = local scope node id
  std::set<uint64_t> ls_node_ids_set;
  for (const auto& cell : ref_grid_.local_cells)
    for (uint64_t node_id : cell.vertex_ids)
      ls_node_ids_set.insert(node_id);

  // Build node partition subscriptions
  // psub = partition subscription
  // Multiple partitions can subscribe to a given
  // node. We build this list here.
  // We start by adding the current location id
  // as the first subscription
  std::map<uint64_t, std::set<uint64_t>> ls_node_ids_psubs;
  for (const uint64_t node_id : ls_node_ids_set)
    ls_node_ids_psubs[node_id] = {static_cast<uint64_t>(opensn::mpi_comm.rank())};

  // Now we add the partitions associated with the
  // ghost cells.
  const auto ghost_cell_ids = ref_grid_.cells.GetGhostGlobalIDs();
  for (const uint64_t ghost_id : ghost_cell_ids)
  {
    const auto& ghost_cell = ref_grid_.cells[ghost_id];
    for (const uint64_t vid : ghost_cell.vertex_ids)
      ls_node_ids_psubs[vid].insert(ghost_cell.partition_id);
  } // for ghost_id

  // Build lists of local- and non-local nodes
  // The lowest partition-# owns a node.
  std::vector<uint64_t> local_node_ids;
  std::map<uint64_t, std::vector<uint64_t>> nonlocal_node_ids_map;
  for (const uint64_t node_id : ls_node_ids_set)
  {
    uint64_t smallest_partition_id = opensn::mpi_comm.rank();
    for (const uint64_t pid : ls_node_ids_psubs[node_id]) // pid = partition id
      smallest_partition_id = std::min(smallest_partition_id, pid);

    if (smallest_partition_id == opensn::mpi_comm.rank())
      local_node_ids.push_back(node_id);
    else
      nonlocal_node_ids_map[smallest_partition_id].push_back(node_id);
  }

  // Communicate node counts
  const uint64_t local_num_nodes = local_node_ids.size();
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
  globl_base_block_size_ = global_num_nodes;

  // Build node mapping for local
  //                                              nodes
  node_mapping_.clear();
  for (uint64_t i = 0; i < local_num_nodes; ++i)
    node_mapping_[local_node_ids[i]] = static_cast<int64_t>(local_block_address_ + i);

  // Communicate nodes in need
  //                                              of mapping
  std::map<uint64_t, std::vector<uint64_t>> query_node_ids = MapAllToAll(nonlocal_node_ids_map);

  // Map the query nodes
  std::map<uint64_t, std::vector<int64_t>> mapped_node_ids;
  for (const auto& key_value : query_node_ids)
  {
    const uint64_t& pid = key_value.first;
    const auto& node_list = key_value.second;

    for (const uint64_t node_id : node_list)
      if (node_mapping_.count(node_id) == 0)
        throw std::logic_error("Error mapping query node.");
      else
      {
        const int64_t mapping = node_mapping_.at(node_id);
        mapped_node_ids[pid].push_back(mapping);
      }
  } // for query location and nodes

  // Communicate back the mappings
  std::map<uint64_t, std::vector<int64_t>> nonlocal_node_ids_map_mapped =
    MapAllToAll(mapped_node_ids);

  // Processing the mapping for non-local nodes
  ghost_node_mapping_.clear();
  try
  {
    for (const auto& pid_node_ids : nonlocal_node_ids_map)
    {
      const uint64_t& pid = pid_node_ids.first;
      const auto& node_list = pid_node_ids.second;
      const auto& mappings = nonlocal_node_ids_map_mapped.at(pid);

      if (mappings.size() != node_list.size())
        throw std::logic_error("mappings.size() != node_list.size()");

      const size_t num_nodes = node_list.size();
      for (size_t i = 0; i < num_nodes; ++i)
      {
        node_mapping_[node_list[i]] = mappings[i];
        ghost_node_mapping_[node_list[i]] = mappings[i];
      }
    } // for pid and non-local id
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(fname + ": Processing non-local mapping failed.");
  }
  catch (const std::logic_error& lerr)
  {
    throw std::logic_error(fname + ": Processing non-local mapping failed." + lerr.what());
  }
}

void
PieceWiseLinearContinuous::BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                                                std::vector<int64_t>& nodal_nnz_off_diag,
                                                const UnknownManager& unknown_manager) const
{
  //**************************************** DEFINE UTILITIES

  // Dof-handler
  /**Utility mappings*/
  struct DOFHandler
  {
    const int64_t local_block_start = 0;
    const int64_t local_block_end = 0;
    const std::vector<uint64_t>& locI_block_addr;

    DOFHandler(int64_t block_start,
               int64_t block_end,
               const std::vector<uint64_t>& locJ_block_address)
      : local_block_start(block_start),
        local_block_end(block_end),
        locI_block_addr(locJ_block_address)
    {
    }

    bool IsMapLocal(int64_t ir) const { return (ir >= local_block_start and ir < local_block_end); }

    int64_t MapIRLocal(int64_t ir) const { return ir - local_block_start; }

    int64_t GetLocFromIR(int64_t ir)
    {
      int64_t locI = std::upper_bound(locI_block_addr.begin(), locI_block_addr.end(), ir) -
                     locI_block_addr.begin() - 1;
      return locI;
    }

  } dof_handler(static_cast<int64_t>(local_block_address_),
                static_cast<int64_t>(local_block_address_ + local_base_block_size_),
                locJ_block_address_);

  // Writes a message on ir error
  auto IR_MAP_ERROR = []()
  {
    log.LogAllError() << "PWL-MapCFEMDOF: ir Mapping error node ";
    Exit(EXIT_FAILURE);
  };

  // Writes a message on jr error
  auto JR_MAP_ERROR = []()
  {
    log.LogAllError() << "PWL-MapCFEMDOF: jr Mapping error node ";
    Exit(EXIT_FAILURE);
  };

  // Checks whether an integer is already in a vector
  auto IS_VALUE_IN_VECTOR = [](const std::vector<int64_t>& vec, int64_t val)
  {
    bool already_there = false;
    for (auto check_val : vec)
      if (check_val == val)
      {
        already_there = true;
        break;
      }
    return already_there;
  };

  //**************************************** END OF UTILITIES

  // Build local sparsity pattern
  log.Log0Verbose1() << "Building local sparsity pattern.";
  std::vector<std::vector<int64_t>> nodal_connections(local_base_block_size_);

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_base_block_size_, 0);
  nodal_nnz_off_diag.resize(local_base_block_size_, 0);

  for (auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);
    for (unsigned int i = 0; i < cell_mapping.NumNodes(); ++i)
    {
      const int64_t ir = MapDOF(cell, i);
      if (ir < 0)
        IR_MAP_ERROR();

      if (dof_handler.IsMapLocal(ir))
      {
        const int64_t il = dof_handler.MapIRLocal(ir);
        std::vector<int64_t>& node_links = nodal_connections[il];

        for (unsigned int j = 0; j < cell_mapping.NumNodes(); ++j)
        {
          const int64_t jr = MapDOF(cell, j);
          if (jr < 0)
            JR_MAP_ERROR();

          if (IS_VALUE_IN_VECTOR(node_links, jr))
            continue;

          node_links.push_back(jr);
          if (dof_handler.IsMapLocal(jr))
            nodal_nnz_in_diag[il] += 1;
          else
            nodal_nnz_off_diag[il] += 1;
        } // for j
      }   // if i local
    }     // for i
  }       // for cell

  // Build non-local sparsity pattern
  log.Log0Verbose1() << "Building non-local sparsity pattern.";

  // In this process we build a list
  // of ir-nodes that are not local. Each ir-node needs to
  // be furnished with the jr-nodes it links to.

  using ROWJLINKS = std::pair<int64_t, std::vector<int64_t>>;
  std::vector<ROWJLINKS> ir_links;

  for (auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);

    for (unsigned int i = 0; i < cell_mapping.NumNodes(); ++i)
    {
      const int64_t ir = MapDOF(cell, i);
      if (ir < 0)
        IR_MAP_ERROR();

      if (not dof_handler.IsMapLocal(ir))
      {
        ROWJLINKS new_ir_link;
        ROWJLINKS* cur_ir_link = &new_ir_link;

        // Check if ir already there
        bool ir_already_there = false;
        for (auto& ir_link : ir_links)
          if (ir == ir_link.first)
          {
            ir_already_there = true;
            cur_ir_link = &ir_link;
            break;
          }

        // Now add links
        auto& node_links = cur_ir_link->second;
        for (unsigned int j = 0; j < cell_mapping.NumNodes(); ++j)
        {
          const int64_t jr = MapDOF(cell, j);
          if (jr < 0)
            JR_MAP_ERROR();

          if (IS_VALUE_IN_VECTOR(node_links, jr))
            continue;
          else
            node_links.push_back(jr);
        } // for j

        // If its not add it
        if (not ir_already_there)
        {
          cur_ir_link->first = ir;
          ir_links.push_back(*cur_ir_link);
        }

      } // if i not local
    }   // for i
  }     // for cell

  // Build communication structure
  log.Log0Verbose1() << "Building communication structure.";

  // Step 1
  // We now serialize the non-local data
  std::vector<std::vector<int64_t>> locI_serialized(opensn::mpi_comm.size());

  for (const auto& ir_linkage : ir_links)
  {
    const int64_t locI = dof_handler.GetLocFromIR(ir_linkage.first);

    // row cols amount
    locI_serialized[locI].push_back(static_cast<int64_t>(ir_linkage.second.size()));
    // row num
    locI_serialized[locI].push_back(static_cast<int64_t>(ir_linkage.first));
    for (int64_t jr : ir_linkage.second)
      locI_serialized[locI].push_back(jr); // col num
  }

  std::vector<int64_t> recvbuf;
  mpi_comm.all_to_all(locI_serialized, recvbuf);

  // Deserialze data
  log.Log0Verbose1() << "Deserialize data.";

  std::vector<ROWJLINKS> foreign_ir_links;

  for (size_t k = 0; k < recvbuf.size();)
  {
    const int64_t num_values = recvbuf[k++];
    const int64_t ir = recvbuf[k++];

    ROWJLINKS new_links;
    new_links.first = ir;
    new_links.second.reserve(num_values);
    for (int i = 0; i < num_values; ++i)
      new_links.second.push_back(recvbuf[k++]);

    foreign_ir_links.push_back(std::move(new_links));
  }

  // Adding to sparsity pattern
  for (const auto& ir_linkage : foreign_ir_links)
  {
    const int64_t ir = ir_linkage.first;

    if (not dof_handler.IsMapLocal(ir))
      IR_MAP_ERROR();

    int64_t il = dof_handler.MapIRLocal(ir);
    std::vector<int64_t>& node_links = nodal_connections[il];

    for (int64_t jr : ir_linkage.second)
    {
      if (IS_VALUE_IN_VECTOR(node_links, jr))
        continue;

      node_links.push_back(jr);
      if (dof_handler.IsMapLocal(jr))
        nodal_nnz_in_diag[il] += 1;
      else
        nodal_nnz_off_diag[il] += 1;
    }
  }

  opensn::mpi_comm.barrier();

  // Spacing according to unknown manager
  auto backup_nnz_in_diag = nodal_nnz_in_diag;
  auto backup_nnz_off_diag = nodal_nnz_off_diag;

  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_base_block_size_ * N, 0);
  nodal_nnz_off_diag.resize(local_base_block_size_ * N, 0);

  if (unknown_manager.dof_storage_type == UnknownStorageType::NODAL)
  {
    int ir = -1;
    for (int i = 0; i < local_base_block_size_; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        ++ir;
        nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
        nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
      } // for j
    }   // for i
  }
  else if (unknown_manager.dof_storage_type == UnknownStorageType::BLOCK)
  {
    int ir = -1;
    for (int j = 0; j < N; ++j)
    {
      for (int i = 0; i < local_base_block_size_; ++i)
      {
        ++ir;
        nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
        nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
      } // for i
    }   // for j
  }
}

int64_t
PieceWiseLinearContinuous::MapDOF(const Cell& cell,
                                  const unsigned int node,
                                  const UnknownManager& unknown_manager,
                                  const unsigned int unknown_id,
                                  const unsigned int component) const
{
  const uint64_t vertex_id = cell.vertex_ids[node];

  OpenSnLogicalErrorIf(node_mapping_.count(vertex_id) == 0,
                       std::string("Bad trouble mapping vertex ") + std::to_string(vertex_id));
  const int64_t global_id = node_mapping_.at(vertex_id);

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  auto storage = unknown_manager.dof_storage_type;

  int64_t address = -1;
  if (storage == UnknownStorageType::BLOCK)
  {
    for (int locJ = 0; locJ < opensn::mpi_comm.size(); ++locJ)
    {
      const int64_t local_id = global_id - static_cast<int64_t>(locJ_block_address_[locJ]);

      if (local_id < 0 or local_id >= locJ_block_size_[locJ])
        continue;

      address = static_cast<int64_t>(locJ_block_address_[locJ] * num_unknowns) +
                static_cast<int64_t>(locJ_block_size_[locJ] * block_id) + local_id;
      break;
    }
  }
  else if (storage == UnknownStorageType::NODAL)
    address = global_id * static_cast<int64_t>(num_unknowns) + static_cast<int64_t>(block_id);

  return address;
}

int64_t
PieceWiseLinearContinuous::MapDOFLocal(const Cell& cell,
                                       const unsigned int node,
                                       const UnknownManager& unknown_manager,
                                       const unsigned int unknown_id,
                                       const unsigned int component) const
{
  const uint64_t vertex_id = cell.vertex_ids[node];

  OpenSnLogicalErrorIf(node_mapping_.count(vertex_id) == 0, "Bad trouble");
  const int64_t node_global_id = node_mapping_.at(vertex_id);

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  auto storage = unknown_manager.dof_storage_type;

  const int64_t local_id = node_global_id - static_cast<int64_t>(local_block_address_);
  const bool is_local = not(local_id < 0 or local_id >= local_base_block_size_);

  int64_t address = -1;
  if (is_local)
  {
    if (storage == UnknownStorageType::BLOCK)
    {
      address = static_cast<int64_t>(local_base_block_size_ * block_id) + local_id;
    }
    else if (storage == UnknownStorageType::NODAL)
      address = local_id * static_cast<int64_t>(num_unknowns) + static_cast<int64_t>(block_id);
  } // if is_local
  else
  {
    const size_t num_local_dofs = NumLocalDOFs(unknown_manager);
    int64_t ghost_local_node_id = -1;
    int64_t counter = 0;
    for (const auto& vid_gnid : ghost_node_mapping_)
    {
      if (node_global_id == vid_gnid.second)
      {
        ghost_local_node_id = counter;
        break;
      }
      ++counter;
    }
    if (storage == UnknownStorageType::BLOCK)
    {
      address = static_cast<int64_t>(ghost_node_mapping_.size() * block_id) + ghost_local_node_id;
    }
    else if (storage == UnknownStorageType::NODAL)
      address =
        ghost_local_node_id * static_cast<int64_t>(num_unknowns) + static_cast<int64_t>(block_id);

    address += static_cast<int64_t>(num_local_dofs);
  }

  return address;
}

size_t
PieceWiseLinearContinuous::NumGhostDOFs(const UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return ghost_node_mapping_.size() * N;
}

std::vector<int64_t>
PieceWiseLinearContinuous::GhostDOFIndices(const UnknownManager& unknown_manager) const
{
  std::vector<int64_t> dof_ids;
  dof_ids.reserve(NumGhostDOFs(unknown_manager));

  const size_t num_unknown_comps = unknown_manager.GetTotalUnknownStructureSize();
  const auto storage = unknown_manager.dof_storage_type;
  const size_t num_unknowns = unknown_manager.unknowns.size();

  for (const auto& vid_gnid : ghost_node_mapping_)
  {
    const int64_t global_id = vid_gnid.second;

    for (size_t u = 0; u < num_unknowns; ++u)
    {
      const auto& unkn = unknown_manager.unknowns[u];
      const size_t num_comps = unkn.num_components;
      for (size_t c = 0; c < num_comps; ++c)
      {
        size_t block_id = unknown_manager.MapUnknown(u, c);
        int64_t address = -1;
        if (storage == UnknownStorageType::BLOCK)
        {
          for (int locJ = 0; locJ < opensn::mpi_comm.size(); ++locJ)
          {
            const int64_t local_id = global_id - static_cast<int64_t>(locJ_block_address_[locJ]);

            if (local_id < 0 or local_id >= locJ_block_size_[locJ])
              continue;

            address = static_cast<int64_t>(locJ_block_address_[locJ] * num_unknown_comps) +
                      static_cast<int64_t>(locJ_block_size_[locJ] * block_id) + local_id;
            break;
          }
        }
        else if (storage == UnknownStorageType::NODAL)
          address =
            global_id * static_cast<int64_t>(num_unknown_comps) + static_cast<int64_t>(block_id);

        dof_ids.push_back(address);
      } // for c
    }   // for u
  }

  return dof_ids;
}

} // namespace opensn
