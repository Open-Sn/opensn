// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_continuum/grid_face_histogram.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

using LockBox = std::vector<std::pair<int, short>>;

AAH_FLUDSCommonData::AAH_FLUDSCommonData(
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SPDS& spds,
  const GridFaceHistogram& grid_face_histogram)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
  this->InitializeAlphaElements(spds, grid_face_histogram);
  this->InitializeBetaElements(spds);
}

void
AAH_FLUDSCommonData::InitializeAlphaElements(const SPDS& spds,
                                             const GridFaceHistogram& grid_face_histogram)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::InitializeAlphaElements");

  const auto grid = spds.GetGrid();
  const std::vector<int>& spls = spds.GetLocalSubgrid();

  // Initialize face categorization
  num_face_categories_ = grid_face_histogram.GetNumberOfFaceHistogramBins();
  // TODO: Check if we can move this down
  local_psi_stride_.resize(num_face_categories_, 0);
  local_psi_max_elements_.resize(num_face_categories_, 0);
  local_psi_n_block_stride_.resize(num_face_categories_, 0);
  local_psi_Gn_block_strideG_.resize(num_face_categories_, 0);

  // Initialize dependent locations
  size_t num_of_deplocs = spds.GetLocationSuccessors().size();
  deplocI_face_dof_count_.resize(num_of_deplocs, 0);
  deplocI_cell_views_.resize(num_of_deplocs);

  // PERFORM SLOT DYNAMICS
  // Loop over cells in sweep order

  // Given a local cell index, gives the so index
  std::vector<int> local_so_cell_mapping;
  local_so_cell_mapping.resize(grid->local_cells.size(), 0);

  largest_face_ = 0;                                     // Will contain the max dofs per face
  std::vector<LockBox> lock_boxes(num_face_categories_); // cell,face index pairs
  LockBox delayed_lock_box;
  std::set<int> location_boundary_dependency_set;

  // csoi = cell sweep order index
  so_cell_inco_face_face_category_.reserve(spls.size());
  so_cell_outb_face_slot_indices_.reserve(spls.size());
  so_cell_outb_face_face_category_.reserve(spls.size());
  for (auto csoi = 0; csoi < spls.size(); ++csoi)
  {
    auto cell_local_id = spls[csoi];
    const auto& cell = grid->local_cells[cell_local_id];

    local_so_cell_mapping[cell.local_id] = csoi; // Set mapping

    SlotDynamics(cell, spds, grid_face_histogram, lock_boxes, delayed_lock_box);

  } // for csoi

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Done with Slot Dynamics.";
  opensn::mpi_comm.barrier();

  // PERFORM INCIDENT MAPPING
  // Loop over cells in sweep order
  so_cell_inco_face_dof_indices_.reserve(spls.size());
  for (auto csoi = 0; csoi < spls.size(); ++csoi)
  {
    auto cell_local_id = spls[csoi];
    const auto& cell = grid->local_cells[cell_local_id];

    LocalIncidentMapping(cell, spds, local_so_cell_mapping);

  } // for csoi

  for (auto fc = 0; fc < num_face_categories_; ++fc)
  {
    local_psi_stride_[fc] = grid_face_histogram.GetFaceHistogramBinDOFSize(fc);
    local_psi_max_elements_[fc] = lock_boxes[fc].size();
    local_psi_n_block_stride_[fc] = local_psi_stride_[fc] * lock_boxes[fc].size();
    local_psi_Gn_block_strideG_[fc] = local_psi_n_block_stride_[fc] * /*G=*/1;
  }
  delayed_local_psi_stride_ = largest_face_;
  delayed_local_psi_max_elements_ = delayed_lock_box.size();
  delayed_local_psi_Gn_block_stride_ = largest_face_ * delayed_lock_box.size();
  delayed_local_psi_Gn_block_strideG_ = delayed_local_psi_Gn_block_stride_ * /*G=*/1;

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Done with Local Incidence mapping.";
  opensn::mpi_comm.barrier();

  // Clean up
  so_cell_outb_face_slot_indices_.shrink_to_fit();

  local_so_cell_mapping.clear();
  local_so_cell_mapping.shrink_to_fit();

  so_cell_inco_face_dof_indices_.shrink_to_fit();

  nonlocal_outb_face_deplocI_slot_.shrink_to_fit();
}

void
AAH_FLUDSCommonData::SlotDynamics(const Cell& cell,
                                  const SPDS& spds,
                                  const GridFaceHistogram& grid_face_histogram,
                                  std::vector<std::vector<std::pair<int, short>>>& lock_boxes,
                                  std::vector<std::pair<int, short>>& delayed_lock_box)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::SlotDynamics");

  auto grid = spds.GetGrid();

  // Loop over faces but process only incident faces
  std::vector<short> inco_face_face_category;
  inco_face_face_category.reserve(cell.faces.size());
  for (auto f = 0; f < cell.faces.size(); ++f)
  {
    const CellFace& face = cell.faces[f];
    const auto& orientation = spds.GetCellFaceOrientations()[cell.local_id][f];

    // Incident face
    if (orientation == FaceOrientation::INCOMING)
    {

      // Local cell dependence
      if (face.IsNeighborLocal(grid.get()))
      {
        size_t num_face_dofs = face.vertex_ids.size();
        size_t face_categ = grid_face_histogram.MapFaceHistogramBins(num_face_dofs);

        inco_face_face_category.push_back(static_cast<short>(face_categ));

        LockBox& lock_box = lock_boxes[face_categ];

        // Check if part of cyclic dependency
        bool is_cyclic = false;
        for (auto cyclic_dependency : spds.GetLocalSweepFAS())
        {
          auto edge_v1_id = cyclic_dependency.first;
          auto edge_v2_id = cyclic_dependency.second;
          auto cell_id = cell.local_id;
          auto neighbor_cell_id = face.GetNeighborLocalID(grid.get());

          if ((edge_v1_id == cell_id) and (edge_v2_id == neighbor_cell_id))
          {
            is_cyclic = true;
            inco_face_face_category.back() *= -1;
            inco_face_face_category.back() -= 1;
          }

          if ((edge_v1_id == neighbor_cell_id) and (edge_v2_id == cell_id))
          {
            is_cyclic = true;
            inco_face_face_category.back() *= -1;
            inco_face_face_category.back() -= 1;
          }
        }
        if (is_cyclic)
          continue;

        // Find associated face for dof mapping and lock box
        auto adj_face_idx = (short)face.GetNeighborAdjacentFaceIndex(grid.get());

        // Now find the cell (index,face) pair in the lock box and empty slot
        bool found = false;
        for (auto& lock_box_slot : lock_box)
        {
          if ((static_cast<uint64_t>(lock_box_slot.first) == face.neighbor_id) and
              (lock_box_slot.second == adj_face_idx))
          {
            lock_box_slot.first = -1;
            lock_box_slot.second = -1;
            found = true;
            break;
          }
        }
        if (not found)
        {
          std::ostringstream oss;
          oss << "AAH_FLUDSCommonData: Lock-box location not found in call to "
                 "InitializeAlphaElements.\n"
              << "Local cell: " << std::to_string(cell.local_id) << ", "
              << "Face: " << std::to_string(f) << ", "
              << "Looking for cell: " << std::to_string(face.GetNeighborLocalID(grid.get())) << ", "
              << "Adjacent face: " << std::to_string(adj_face_idx) << ", "
              << "Category: " << std::to_string(face_categ) << ", "
              << "Omega: " << spds.GetOmega().PrintStr() << ", "
              << "Lock-box size: " << std::to_string(lock_box.size());
          throw std::runtime_error(oss.str());
        }
      } // if local
    } // if incident

  } // for f

  so_cell_inco_face_face_category_.push_back(inco_face_face_category);

  // Loop over faces but process only outgoing faces
  std::vector<int> outb_face_slot_indices;
  std::vector<short> outb_face_face_category;
  outb_face_slot_indices.reserve(cell.faces.size());
  outb_face_face_category.reserve(cell.faces.size());
  for (auto f = 0; f < cell.faces.size(); ++f)
  {
    const CellFace& face = cell.faces[f];
    auto cell_g_index = cell.global_id;
    const auto& orientation = spds.GetCellFaceOrientations()[cell.local_id][f];

    // Outgoing face
    if (orientation == FaceOrientation::OUTGOING)
    {
      size_t num_face_dofs = face.vertex_ids.size();
      size_t face_categ = grid_face_histogram.MapFaceHistogramBins(num_face_dofs);

      outb_face_face_category.push_back(static_cast<short>(face_categ));

      LockBox* temp_lock_box = &lock_boxes[face_categ];

      // Check if part of cyclic dependency
      if (face.IsNeighborLocal(grid.get()))
      {
        for (auto cyclic_dependency : spds.GetLocalSweepFAS())
        {
          auto a = cyclic_dependency.first;
          auto b = cyclic_dependency.second;
          auto c = cell.local_id;
          auto d = face.GetNeighborLocalID(grid.get());

          if ((a == c) and (b == d))
          {
            temp_lock_box = &delayed_lock_box;
            outb_face_face_category.back() *= -1;
            outb_face_face_category.back() -= 1;
          }

          if ((a == d) and (b == c))
          {
            temp_lock_box = &delayed_lock_box;
            outb_face_face_category.back() *= -1;
            outb_face_face_category.back() -= 1;
          }
        }
      }

      LockBox& lock_box = *temp_lock_box;

      // Check if this face is the max size
      if (num_face_dofs > largest_face_)
        largest_face_ = static_cast<int>(num_face_dofs);

      // Find a open slot
      bool slot_found = false;
      for (auto k = 0; k < lock_box.size(); ++k)
      {
        if (lock_box[k].first < 0)
        {
          outb_face_slot_indices.push_back(k);
          lock_box[k].first = cell_g_index;
          lock_box[k].second = static_cast<short>(f);
          slot_found = true;
          break;
        }
      }

      // If an open slot was not found push a new one
      if (not slot_found)
      {
        outb_face_slot_indices.push_back(lock_box.size());
        lock_box.push_back(std::pair<int, short>(cell_g_index, f));
      }

      // Non-local outgoing
      if (face.has_neighbor and (not face.IsNeighborLocal(grid.get())))
      {
        auto locJ = face.GetNeighborPartitionID(grid.get());
        auto deplocI = spds.MapLocJToDeplocI(locJ);
        auto face_slot = deplocI_face_dof_count_[deplocI];

        deplocI_face_dof_count_[deplocI] += face.vertex_ids.size();

        nonlocal_outb_face_deplocI_slot_.emplace_back(deplocI, face_slot);

        // The following function is defined below
        AddFaceViewToDepLocI(deplocI, cell_g_index, face_slot, face);

      } // non-local neighbor
    } // if outgoing

  } // for f

  so_cell_outb_face_slot_indices_.push_back(outb_face_slot_indices);
  so_cell_outb_face_face_category_.push_back(outb_face_face_category);
}

void
AAH_FLUDSCommonData::AddFaceViewToDepLocI(int deplocI,
                                          int cell_g_index,
                                          int face_slot,
                                          const CellFace& face)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::AddFaceViewToDepLocI");

  // Check if cell is already there
  bool cell_already_there = false;
  for (auto& cell_view : deplocI_cell_views_[deplocI])
  {
    if (cell_view.first == cell_g_index)
    {
      cell_already_there = true;
      cell_view.second.emplace_back(face_slot, face.vertex_ids);
      break;
    }
  }

  // If the cell is not there yet
  if (not cell_already_there)
  {
    CompactCellView new_cell_view;
    new_cell_view.first = cell_g_index;
    new_cell_view.second.emplace_back(face_slot, face.vertex_ids);

    deplocI_cell_views_[deplocI].push_back(new_cell_view);
  }
}

void
AAH_FLUDSCommonData::LocalIncidentMapping(const Cell& cell,
                                          const SPDS& spds,
                                          std::vector<int>& local_so_cell_mapping)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::LocalIncidentMapping");

  const auto grid = spds.GetGrid();
  auto& cell_nodal_mapping = grid_nodal_mappings_[cell.local_id];
  std::vector<std::pair<int, std::vector<short>>> inco_face_dof_mapping;

  // Loop over faces but process only incident faces
  for (auto f = 0; f < cell.faces.size(); ++f)
  {
    const CellFace& face = cell.faces[f];
    const auto& orienation = spds.GetCellFaceOrientations()[cell.local_id][f];

    // Incident face
    if (orienation == FaceOrientation::INCOMING)
    {
      if (face.IsNeighborLocal(grid.get()))
      {
        // Find associated face for dof mapping
        auto adj_face_idx = cell_nodal_mapping[f].associated_face_;

        std::pair<int, std::vector<short>> dof_mapping;
        dof_mapping.second = cell_nodal_mapping[f].face_node_mapping_;

        // Find associated face counter for slot lookup
        const auto& adj_cell = grid->cells[face.neighbor_id];
        const auto adj_so_index = local_so_cell_mapping[adj_cell.local_id];
        const auto& face_oris = spds.GetCellFaceOrientations()[adj_cell.local_id];
        int adj_f_counter = -1;

        int out_f = -1;
        for (auto af = 0; af < adj_cell.faces.size(); ++af)
        {
          if (face_oris[af] == FaceOrientation::OUTGOING)
          {
            ++out_f;
          }

          if (af == adj_face_idx)
          {
            adj_f_counter = out_f;
            break;
          }
        }

        dof_mapping.first = /*local_psi_stride*G**/
          so_cell_outb_face_slot_indices_[adj_so_index][adj_f_counter];

        dof_mapping.second.shrink_to_fit();
        inco_face_dof_mapping.push_back(dof_mapping);
      } // if local
    } // if incident
  } // for incindent f

  std::vector<INCOMING_FACE_INFO> inco_face_info_array(inco_face_dof_mapping.size());
  for (auto i = 0; i < inco_face_dof_mapping.size(); ++i)
    inco_face_info_array[i].Setup(inco_face_dof_mapping[i]);

  so_cell_inco_face_dof_indices_.push_back(inco_face_info_array);
}

void
AAH_FLUDSCommonData::InitializeBetaElements(const SPDS& spds, int tag_index /*=0*/)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::InitializeBetaElements");

  const auto grid = spds.GetGrid();
  const std::vector<int>& spls = spds.GetLocalSubgrid();

  // The first two major steps here are: Send delayed successor information
  // and Receive delayed predecessor information. The send portion is done
  // first because the delayed information does not follow the
  // Task Dependency Graph and hence when a location receives its delayed
  // information, the information might not have been sent yet.

  // Send delayed successor information
  const auto& location_successors = spds.GetLocationSuccessors();
  const auto& delayed_location_successors = spds.GetDelayedLocationSuccessors();
  std::vector<mpi::Request> send_requests;
  send_requests.resize(location_successors.size());
  std::vector<std::vector<int>> multi_face_indices(location_successors.size(), std::vector<int>());
  for (auto deplocI = 0; deplocI < location_successors.size(); ++deplocI)
  {
    auto locJ = location_successors[deplocI];

    std::vector<int>::const_iterator delayed_successor =
      std::find(delayed_location_successors.begin(), delayed_location_successors.end(), locJ);
    if ((delayed_successor == delayed_location_successors.end()))
      continue;

    std::vector<CompactCellView> cell_views = deplocI_cell_views_[deplocI];

    SerializeCellInfo(cell_views, multi_face_indices[deplocI], deplocI_face_dof_count_[deplocI]);

    send_requests[deplocI] = mpi_comm.isend(locJ, 101 + tag_index, multi_face_indices[deplocI]);

    // TODO: Watch eager limits on sent data

    deplocI_cell_views_[deplocI].clear();
    deplocI_cell_views_[deplocI].shrink_to_fit();
  }

  // Receive delayed predecessor information
  const auto& delayed_location_dependencies = spds.GetDelayedLocationDependencies();
  delayed_prelocI_cell_views_.resize(delayed_location_dependencies.size(),
                                     std::vector<CompactCellView>());
  delayed_prelocI_face_dof_count_.resize(delayed_location_dependencies.size(), 0);
  for (auto prelocI = 0; prelocI < delayed_location_dependencies.size(); ++prelocI)
  {
    auto locJ = delayed_location_dependencies[prelocI];

    std::vector<int> face_indices;
    mpi_comm.recv(locJ, 101 + tag_index, face_indices);

    DeSerializeCellInfo(delayed_prelocI_cell_views_[prelocI],
                        &face_indices,
                        delayed_prelocI_face_dof_count_[prelocI]);
  }

  // The next two operations is to receive predecessor information followed
  // by the sending of successor information. The receives are blocking but
  // will cause a dead lock because the successor/predecessor combination
  // follows the TDG.

  // Receive predecessor information
  const auto& location_dependencies = spds.GetLocationDependencies();
  prelocI_cell_views_.resize(location_dependencies.size(), std::vector<CompactCellView>());
  prelocI_face_dof_count_.resize(location_dependencies.size(), 0);
  for (auto prelocI = 0; prelocI < location_dependencies.size(); ++prelocI)
  {
    auto locJ = location_dependencies[prelocI];

    std::vector<int> face_indices;
    mpi_comm.recv(locJ, 101 + tag_index, face_indices);

    DeSerializeCellInfo(
      prelocI_cell_views_[prelocI], &face_indices, prelocI_face_dof_count_[prelocI]);
  }

  // Send successor information
  for (auto deplocI = 0; deplocI < location_successors.size(); ++deplocI)
  {
    auto locJ = location_successors[deplocI];

    auto delayed_successor =
      std::find(delayed_location_successors.begin(), delayed_location_successors.end(), locJ);
    if ((delayed_successor != delayed_location_successors.end()))
      continue;

    std::vector<CompactCellView> cell_views = deplocI_cell_views_[deplocI];

    SerializeCellInfo(cell_views, multi_face_indices[deplocI], deplocI_face_dof_count_[deplocI]);

    send_requests[deplocI] = mpi_comm.isend(locJ, 101 + tag_index, multi_face_indices[deplocI]);

    // TODO: Watch eager limits on sent data

    deplocI_cell_views_[deplocI].clear();
    deplocI_cell_views_[deplocI].shrink_to_fit();
  }

  // Verify sends completed
  for (auto deplocI = 0; deplocI < location_successors.size(); ++deplocI)
    mpi::wait(send_requests[deplocI]);
  multi_face_indices.clear();
  multi_face_indices.shrink_to_fit();

  // In the next process we loop over cells in the sweep order and perform
  // the non-local face mappings. This is dependent on having the compact
  // cellviews on the partition interfaces.

  // Loop over cells in sorder
  for (auto csoi = 0; csoi < spls.size(); ++csoi)
  {
    auto cell_local_index = spls[csoi];
    const auto& cell = grid->local_cells[cell_local_index];

    NonLocalIncidentMapping(cell, spds);
  } // for csoi

  deplocI_cell_views_.clear();
  deplocI_cell_views_.shrink_to_fit();

  prelocI_cell_views_.clear();
  prelocI_cell_views_.shrink_to_fit();

  delayed_prelocI_cell_views_.clear();
  delayed_prelocI_cell_views_.shrink_to_fit();

  // Clear unneccesary data
  auto empty_vector = std::vector<std::vector<CompactCellView>>(0);
  deplocI_cell_views_.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  prelocI_cell_views_.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  delayed_prelocI_cell_views_.swap(empty_vector);
}

void
AAH_FLUDSCommonData::SerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                       std::vector<int>& face_indices,
                                       int num_face_dofs)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::SerializeCellInfo");

  const size_t num_cells = cell_views.size();

  // First entry is number of face dofs
  face_indices.push_back(num_face_dofs);

  // Second entry is amount of cells
  face_indices.push_back(static_cast<int>(num_cells));

  // Third entry is negative global cell index
  // Each time a negative entry occurs it denotes a cell face but
  // the actual number is -cell_g_index-1. The offset is necessary
  // for evaluating the negative. The offset is restored during the
  // deserialization process.
  // It is followed by a positive number which is the store location
  // of the face
  for (size_t c = 0; c < num_cells; ++c)
  {
    int glob_index = -cell_views[c].first - 1;

    std::vector<CompactFaceView>& cell_face_views = cell_views[c].second;

    size_t num_faces = cell_face_views.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      face_indices.push_back(glob_index);
      face_indices.push_back(cell_face_views[f].first);
      std::vector<uint64_t>& face_vertices = cell_face_views[f].second;

      size_t num_verts = face_vertices.size();
      for (auto fi = 0; fi < num_verts; ++fi)
      {
        face_indices.push_back(static_cast<int>(face_vertices[fi]));
      }
    }
  }
}

void
AAH_FLUDSCommonData::DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                         std::vector<int>* face_indices,
                                         int& num_face_dofs)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::DeSerializeCellInfo");

  num_face_dofs = (*face_indices)[0];
  int num_cells = (*face_indices)[1];

  cell_views.resize(num_cells);

  int k = 2;
  int last_cell = -1;
  int c = -1; // cell counter
  int f = -1;
  while (k < face_indices->size())
  {
    int entry = (*face_indices)[k];
    // Cell/Face indicator
    if (entry < 0)
    {
      if (-entry != last_cell)
      {
        cell_views.emplace_back();
        c++;
        cell_views[c].first = -entry - 1;

        cell_views[c].second.emplace_back();
        f = 0;

        last_cell = -entry;

        cell_views[c].second[f].first = (*face_indices)[k + 1];
        k++;
      }
      else
      {
        cell_views[c].second.emplace_back();
        f++;

        cell_views[c].second[f].first = (*face_indices)[k + 1];
        k++;
      }
    }
    // Face vertex
    else
    {
      cell_views[c].second[f].second.push_back(entry);
    }
    k++;
  } // while k
}

void
AAH_FLUDSCommonData::NonLocalIncidentMapping(const Cell& cell, const SPDS& spds)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDSCommonData::NonLocalIncidentMapping");

  const auto grid = spds.GetGrid();

  // Loop over faces but process only incident faces
  for (auto f = 0; f < cell.faces.size(); ++f)
  {
    const CellFace& face = cell.faces[f];
    const auto& orientation = spds.GetCellFaceOrientations()[cell.local_id][f];

    // Incident face
    if (orientation == FaceOrientation::INCOMING)
    {
      if ((face.has_neighbor) and (!face.IsNeighborLocal(grid.get())))
      {
        // Find prelocI
        auto locJ = face.GetNeighborPartitionID(grid.get());
        auto prelocI = spds.MapLocJToPrelocI(locJ);

        if (prelocI >= 0)
        {
          // Find the cell in prelocI cell views
          int adj_cell = -1;
          for (auto c = 0; c < prelocI_cell_views_[prelocI].size(); ++c)
          {
            if (prelocI_cell_views_[prelocI][c].first == face.neighbor_id)
            {
              adj_cell = c;
              break;
            }
          }
          if (adj_cell < 0)
          {
            std::ostringstream oss;
            oss << "AAH_FLUDSCommonData: Required predecessor cell not found in call to "
                   "InitializeBetaElements ("
                << "locJ = " << locJ << ", prelocI = " << prelocI << ", cell = " << face.neighbor_id
                << ")";
            throw std::runtime_error(oss.str());
          }

          // Find associated face
          std::set<int> cfvids(face.vertex_ids.begin(), face.vertex_ids.end());
          CompactCellView* adj_cell_view = &prelocI_cell_views_[prelocI][adj_cell];
          int adj_face_idx = -1, af = -1;
          for (auto& adj_face : adj_cell_view->second)
          {
            ++af;
            bool face_matches = true;

            std::set<int> afvids(adj_face.second.begin(), adj_face.second.end());

            if (cfvids != afvids)
              face_matches = false;

            if (face_matches)
            {
              adj_face_idx = af;
              break;
            }
          }
          if (adj_face_idx < 0)
            throw std::runtime_error(
              "AAH_FLUDSCommonData: Associated face not found in call to InitializeBetaElements");

          // Map dofs
          std::pair<int, std::vector<int>> dof_mapping;
          dof_mapping.first = adj_cell_view->second[adj_face_idx].first;
          std::vector<uint64_t>* adj_face_verts = &adj_cell_view->second[adj_face_idx].second;
          for (auto fv = 0; fv < face.vertex_ids.size(); ++fv)
          {
            bool match_found = false;
            for (auto afv = 0; afv < adj_face_verts->size(); ++afv)
            {
              if (face.vertex_ids[fv] == adj_face_verts->operator[](afv))
              {
                match_found = true;
                dof_mapping.second.push_back(afv);
                break;
              }
            }

            if (not match_found)
              throw std::runtime_error("AAH_FLUDSCommonData: Associated vertex not found in call "
                                       "to InitializeBetaElements");
          }

          // Push back final face info
          std::pair<int, std::pair<int, std::vector<int>>> inc_face_prelocI_info;
          inc_face_prelocI_info.first = prelocI;
          inc_face_prelocI_info.second = std::pair<int, std::vector<int>>(dof_mapping);

          std::pair<int, std::pair<int, std::vector<int>>> empty_delayed_info;
          empty_delayed_info.first = prelocI;

          nonlocal_inc_face_prelocI_slot_dof_.push_back(inc_face_prelocI_info);
          delayed_nonlocal_inc_face_prelocI_slot_dof_.push_back(empty_delayed_info);
        } // If not delayed predecessor

        else
        {
          int delayed_preLocI = abs(prelocI) - 1;
          // Find the cell in prelocI cell views
          int ass_cell = -1;
          for (auto c = 0; c < delayed_prelocI_cell_views_[delayed_preLocI].size(); ++c)
          {
            if (delayed_prelocI_cell_views_[delayed_preLocI][c].first == face.neighbor_id)
            {
              ass_cell = c;
              break;
            }
          }
          if (ass_cell < 0)
          {
            std::ostringstream oss;
            oss << "AAH_FLUDSCommonData: Required predecessor cell not located in call to "
                   "InitializeBetaElements ("
                << "locJ = " << locJ << ", delayed prelocI = " << delayed_preLocI
                << ", cell = " << face.neighbor_id << ")";
            throw std::runtime_error(oss.str());
          }

          // Find associated face
          CompactCellView* adj_cell_view = &delayed_prelocI_cell_views_[delayed_preLocI][ass_cell];
          int adj_face_idx = -1;
          for (auto af = 0; af < adj_cell_view->second.size(); ++af)
          {
            bool face_matches = true;
            for (auto afv = 0; afv < adj_cell_view->second[af].second.size(); ++afv)
            {
              bool match_found = false;
              for (auto fv = 0; fv < face.vertex_ids.size(); ++fv)
              {
                if (adj_cell_view->second[af].second[afv] == face.vertex_ids[fv])
                {
                  match_found = true;
                  break;
                }
              }

              if (not match_found)
              {
                face_matches = false;
                break;
              }
            }

            if (face_matches)
            {
              adj_face_idx = af;
              break;
            }
          }
          if (adj_face_idx < 0)
            throw std::runtime_error(
              "AAH_FLUDSCommonData: Associated face not found in call to InitializeBetaElements");

          // Map dofs
          std::pair<int, std::vector<int>> dof_mapping;
          dof_mapping.first = adj_cell_view->second[adj_face_idx].first;
          std::vector<uint64_t>* adj_face_verts = &adj_cell_view->second[adj_face_idx].second;
          for (auto fv = 0; fv < face.vertex_ids.size(); ++fv)
          {
            bool match_found = false;
            for (auto afv = 0; afv < adj_face_verts->size(); ++afv)
            {
              if (face.vertex_ids[fv] == adj_face_verts->operator[](afv))
              {
                match_found = true;
                dof_mapping.second.push_back(afv);
                break;
              }
            }

            if (not match_found)
              throw std::runtime_error("AAH_FLUDSCommonData: Associated vertex not found in call "
                                       "to InitializeBetaElements");
          }

          // Push back final face info
          std::pair<int, std::pair<int, std::vector<int>>> inc_face_prelocI_info;
          inc_face_prelocI_info.first = delayed_preLocI;
          inc_face_prelocI_info.second = std::pair<int, std::vector<int>>(dof_mapping);

          std::pair<int, std::pair<int, std::vector<int>>> delayed_info;
          delayed_info.first = -(delayed_preLocI + 1);

          delayed_nonlocal_inc_face_prelocI_slot_dof_.push_back(inc_face_prelocI_info);
          nonlocal_inc_face_prelocI_slot_dof_.push_back(delayed_info);
        } // If delayed predecessor

      } // if not local and not boundary
    } // if incident
  } // for incindent f
}

} // namespace opensn
