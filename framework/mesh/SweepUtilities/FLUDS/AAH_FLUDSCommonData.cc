#include "framework/mesh/SweepUtilities/FLUDS/AAH_FLUDSCommonData.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/mesh/SweepUtilities/SPDS/SPDS.h"
#include "framework/mesh/MeshContinuum/chi_grid_face_histogram.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include <algorithm>

typedef std::vector<std::pair<int, short>> LockBox;

namespace chi_mesh::sweep_management
{

AAH_FLUDSCommonData::AAH_FLUDSCommonData(
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SPDS& spds,
  const chi_mesh::GridFaceHistogram& grid_face_histogram)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
  this->InitializeAlphaElements(spds, grid_face_histogram);
  this->InitializeBetaElements(spds);
}

void
AAH_FLUDSCommonData::InitializeAlphaElements(const SPDS& spds,
                                             const GridFaceHistogram& grid_face_histogram)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();
  const chi_mesh::sweep_management::SPLS& spls = spds.GetSPLS();

  //================================================== Initialize face
  //                                                   categorization
  num_face_categories = grid_face_histogram.NumberOfFaceHistogramBins();
  // TODO: Check if we can move this down
  local_psi_stride.resize(num_face_categories, 0);
  local_psi_max_elements.resize(num_face_categories, 0);
  local_psi_n_block_stride.resize(num_face_categories, 0);
  local_psi_Gn_block_strideG.resize(num_face_categories, 0);

  //================================================== Initialize dependent
  //                                                   locations
  size_t num_of_deplocs = spds.GetLocationSuccessors().size();
  deplocI_face_dof_count.resize(num_of_deplocs, 0);
  deplocI_cell_views.resize(num_of_deplocs);

  //                      PERFORM SLOT DYNAMICS
  //================================================== Loop over cells in
  //                                                   sweep order

  // Given a local cell index, gives the so index
  std::vector<int> local_so_cell_mapping;
  local_so_cell_mapping.resize(grid.local_cells.size(), 0);

  largest_face = 0;                                     // Will contain the max dofs per face
  std::vector<LockBox> lock_boxes(num_face_categories); // cell,face index pairs
  LockBox delayed_lock_box;
  std::set<int> location_boundary_dependency_set;

  // csoi = cell sweep order index
  so_cell_inco_face_face_category.reserve(spls.item_id.size());
  so_cell_outb_face_slot_indices.reserve(spls.item_id.size());
  so_cell_outb_face_face_category.reserve(spls.item_id.size());
  for (int csoi = 0; csoi < spls.item_id.size(); csoi++)
  {
    int cell_local_id = spls.item_id[csoi];
    const auto& cell = grid.local_cells[cell_local_id];

    local_so_cell_mapping[cell.local_id_] = csoi; // Set mapping

    SlotDynamics(cell,
                 spds,
                 grid_face_histogram,
                 lock_boxes,
                 delayed_lock_box,
                 location_boundary_dependency_set);

  } // for csoi

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Done with Slot Dynamics.";
  Chi::mpi.Barrier();

  //================================================== Populate boundary
  //                                                   dependencies
  for (auto bndry : location_boundary_dependency_set)
    boundary_dependencies.push_back(bndry);

  //                      PERFORM INCIDENT MAPPING
  //================================================== Loop over cells in
  //                                                   sweep order
  so_cell_inco_face_dof_indices.reserve(spls.item_id.size());
  for (int csoi = 0; csoi < spls.item_id.size(); csoi++)
  {
    int cell_local_id = spls.item_id[csoi];
    const auto& cell = grid.local_cells[cell_local_id];

    LocalIncidentMapping(cell, spds, local_so_cell_mapping);

  } // for csoi

  for (size_t fc = 0; fc < num_face_categories; ++fc)
  {
    local_psi_stride[fc] = grid_face_histogram.GetFaceHistogramBinDOFSize(fc);
    local_psi_max_elements[fc] = lock_boxes[fc].size();
    local_psi_n_block_stride[fc] = local_psi_stride[fc] * lock_boxes[fc].size();
    local_psi_Gn_block_strideG[fc] = local_psi_n_block_stride[fc] * /*G=*/1;
  }
  delayed_local_psi_stride = largest_face;
  delayed_local_psi_max_elements = delayed_lock_box.size();
  delayed_local_psi_Gn_block_stride = largest_face * delayed_lock_box.size();
  delayed_local_psi_Gn_block_strideG = delayed_local_psi_Gn_block_stride * /*G=*/1;

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Done with Local Incidence mapping.";
  Chi::mpi.Barrier();

  //================================================== Clean up
  so_cell_outb_face_slot_indices.shrink_to_fit();

  local_so_cell_mapping.clear();
  local_so_cell_mapping.shrink_to_fit();

  so_cell_inco_face_dof_indices.shrink_to_fit();

  nonlocal_outb_face_deplocI_slot.shrink_to_fit();
}

void
AAH_FLUDSCommonData::SlotDynamics(const chi_mesh::Cell& cell,
                                  const SPDS& spds,
                                  const GridFaceHistogram& grid_face_histogram,
                                  std::vector<std::vector<std::pair<int, short>>>& lock_boxes,
                                  std::vector<std::pair<int, short>>& delayed_lock_box,
                                  std::set<int>& location_boundary_dependency_set)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();

  chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  std::vector<short> inco_face_face_category;
  inco_face_face_category.reserve(cell.faces_.size());
  for (int f = 0; f < cell.faces_.size(); f++)
  {
    const CellFace& face = cell.faces_[f];
    const auto& orientation = spds.CellFaceOrientations()[cell.local_id_][f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (orientation == FaceOrientation::INCOMING)
    {

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ LOCAL CELL DEPENDENCE
      if (face.IsNeighborLocal(grid))
      {
        size_t num_face_dofs = face.vertex_ids_.size();
        size_t face_categ = grid_face_histogram.MapFaceHistogramBins(num_face_dofs);

        inco_face_face_category.push_back(static_cast<short>(face_categ));

        LockBox& lock_box = lock_boxes[face_categ];

        //========================================== Check if part of cyclic
        //                                           dependency
        bool is_cyclic = false;
        for (auto cyclic_dependency : spds.GetLocalCyclicDependencies())
        {
          int a = cyclic_dependency.first;
          int b = cyclic_dependency.second;
          int c = cell.local_id_;
          int d = face.GetNeighborLocalID(grid);

          if ((a == c) && (b == d))
          {
            is_cyclic = true;
            inco_face_face_category.back() *= -1;
            inco_face_face_category.back() -= 1;
          }

          if ((a == d) && (b == c))
          {
            is_cyclic = true;
            inco_face_face_category.back() *= -1;
            inco_face_face_category.back() -= 1;
          }
        }
        if (is_cyclic) continue;

        //======================================== Find associated face for
        //                                         dof mapping and lock box
        auto ass_face = (short)face.GetNeighborAssociatedFace(grid);

        // Now find the cell (index,face) pair in the lock box and empty slot
        bool found = false;
        for (auto& lock_box_slot : lock_box)
        {
          if ((lock_box_slot.first == face.neighbor_id_) && (lock_box_slot.second == ass_face))
          {
            lock_box_slot.first = -1;
            lock_box_slot.second = -1;
            found = true;
            break;
          }
        }
        if (!found)
        {
          Chi::log.LogAllError() << "Lock-box location not found in call to "
                                 << "InitializeAlphaElements. Local Cell " << cell.local_id_
                                 << " face " << f << " looking for cell "
                                 << face.GetNeighborLocalID(grid) << " face " << ass_face
                                 << " cat: " << face_categ << " omg=" << spds.Omega().PrintS()
                                 << " lbsize=" << lock_box.size();
          Chi::Exit(EXIT_FAILURE);
        }

      } // if local
      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ BOUNDARY DEPENDENCE
      else if (not face.has_neighbor_)
      {
        const chi_mesh::Vector3& face_norm = face.normal_;

        if (face_norm.Dot(ihat) > 0.999) location_boundary_dependency_set.insert(0);
        else if (face_norm.Dot(ihat) < -0.999)
          location_boundary_dependency_set.insert(1);
        else if (face_norm.Dot(jhat) > 0.999)
          location_boundary_dependency_set.insert(2);
        else if (face_norm.Dot(jhat) < -0.999)
          location_boundary_dependency_set.insert(3);
        else if (face_norm.Dot(khat) > 0.999)
          location_boundary_dependency_set.insert(4);
        else if (face_norm.Dot(khat) < -0.999)
          location_boundary_dependency_set.insert(5);
      }
    } // if incident

  } // for f

  auto raw_inco_face_face_category = new short[inco_face_face_category.size()];
  std::copy(
    inco_face_face_category.begin(), inco_face_face_category.end(), raw_inco_face_face_category);

  so_cell_inco_face_face_category.push_back(raw_inco_face_face_category);

  //=================================================== Loop over faces
  //                OUTGOING                            but process
  //                                                    only outgoing faces
  std::vector<int> outb_face_slot_indices;
  std::vector<short> outb_face_face_category;
  outb_face_slot_indices.reserve(cell.faces_.size());
  outb_face_face_category.reserve(cell.faces_.size());
  for (int f = 0; f < cell.faces_.size(); f++)
  {
    const CellFace& face = cell.faces_[f];
    int cell_g_index = cell.global_id_;
    const auto& orientation = spds.CellFaceOrientations()[cell.local_id_][f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outgoing face
    if (orientation == FaceOrientation::OUTGOING)
    {
      size_t num_face_dofs = face.vertex_ids_.size();
      size_t face_categ = grid_face_histogram.MapFaceHistogramBins(num_face_dofs);

      outb_face_face_category.push_back(static_cast<short>(face_categ));

      LockBox* temp_lock_box = &lock_boxes[face_categ];

      //========================================== Check if part of cyclic
      //                                           dependency
      if (face.IsNeighborLocal(grid))
      {
        for (auto cyclic_dependency : spds.GetLocalCyclicDependencies())
        {
          int a = cyclic_dependency.first;
          int b = cyclic_dependency.second;
          int c = cell.local_id_;
          int d = face.GetNeighborLocalID(grid);

          if ((a == c) && (b == d))
          {
            temp_lock_box = &delayed_lock_box;
            outb_face_face_category.back() *= -1;
            outb_face_face_category.back() -= 1;
          }

          if ((a == d) && (b == c))
          {
            temp_lock_box = &delayed_lock_box;
            outb_face_face_category.back() *= -1;
            outb_face_face_category.back() -= 1;
          }
        }
      }

      LockBox& lock_box = *temp_lock_box;

      //========================================== Check if this face is
      //                                           the max size
      if (num_face_dofs > largest_face) largest_face = static_cast<int>(num_face_dofs);

      //========================================== Find a open slot
      bool slot_found = false;
      for (int k = 0; k < lock_box.size(); k++)
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

      //========================================= If an open slot was not found
      //                                          push a new one
      if (!slot_found)
      {
        outb_face_slot_indices.push_back(lock_box.size());
        lock_box.push_back(std::pair<int, short>(cell_g_index, f));
      }

      //========================================== Non-local outgoing
      if (face.has_neighbor_ and (not face.IsNeighborLocal(grid)))
      {
        int locJ = face.GetNeighborPartitionID(grid);
        int deplocI = spds.MapLocJToDeplocI(locJ);
        int face_slot = deplocI_face_dof_count[deplocI];

        deplocI_face_dof_count[deplocI] += face.vertex_ids_.size();

        nonlocal_outb_face_deplocI_slot.emplace_back(deplocI, face_slot);

        // The following function is defined below
        AddFaceViewToDepLocI(deplocI, cell_g_index, face_slot, face);

      } // non-local neighbor
    }   // if outgoing

  } // for f

  auto raw_outb_face_slot_indices = new int[outb_face_slot_indices.size()];
  std::copy(
    outb_face_slot_indices.begin(), outb_face_slot_indices.end(), raw_outb_face_slot_indices);

  so_cell_outb_face_slot_indices.push_back(raw_outb_face_slot_indices);

  auto raw_outb_face_face_category = new short[outb_face_face_category.size()];
  std::copy(
    outb_face_face_category.begin(), outb_face_face_category.end(), raw_outb_face_face_category);

  so_cell_outb_face_face_category.push_back(raw_outb_face_face_category);
}

void
AAH_FLUDSCommonData::AddFaceViewToDepLocI(int deplocI,
                                          int cell_g_index,
                                          int face_slot,
                                          const chi_mesh::CellFace& face)
{
  //======================================== Check if cell is already there
  bool cell_already_there = false;
  for (auto& cell_view : deplocI_cell_views[deplocI])
  {
    if (cell_view.first == cell_g_index)
    {
      cell_already_there = true;
      cell_view.second.emplace_back(face_slot, face.vertex_ids_);
      break;
    }
  }

  //======================================== If the cell is not there yet
  if (!cell_already_there)
  {
    CompactCellView new_cell_view;
    new_cell_view.first = cell_g_index;
    new_cell_view.second.emplace_back(face_slot, face.vertex_ids_);

    deplocI_cell_views[deplocI].push_back(new_cell_view);
  }
}

void
AAH_FLUDSCommonData::LocalIncidentMapping(const chi_mesh::Cell& cell,
                                          const SPDS& spds,
                                          std::vector<int>& local_so_cell_mapping)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();
  auto& cell_nodal_mapping = grid_nodal_mappings_[cell.local_id_];
  std::vector<std::pair<int, std::vector<short>>> inco_face_dof_mapping;

  short incoming_face_count = -1;

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  for (short f = 0; f < cell.faces_.size(); f++)
  {
    const CellFace& face = cell.faces_[f];
    const auto& orienation = spds.CellFaceOrientations()[cell.local_id_][f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (orienation == FaceOrientation::INCOMING)
    {
      if (face.IsNeighborLocal(grid))
      {
        incoming_face_count++;
        //======================================== Find associated face for
        //                                         dof mapping
        int ass_face = cell_nodal_mapping[f].associated_face_;

        std::pair<int, std::vector<short>> dof_mapping;
        dof_mapping.second = cell_nodal_mapping[f].face_node_mapping_;

        //======================================== Find associated face
        //                                         counter for slot lookup
        const auto& adj_cell = grid.cells[face.neighbor_id_];
        const int adj_so_index = local_so_cell_mapping[adj_cell.local_id_];
        const auto& face_oris = spds.CellFaceOrientations()[adj_cell.local_id_];
        int ass_f_counter = -1;

        int out_f = -1;
        for (size_t af = 0; af < adj_cell.faces_.size(); ++af)
        {
          if (face_oris[af] == FaceOrientation::OUTGOING) { ++out_f; }

          if (af == ass_face)
          {
            ass_f_counter = out_f;
            break;
          }
        }

        dof_mapping.first = /*local_psi_stride*G**/
          so_cell_outb_face_slot_indices[adj_so_index][ass_f_counter];

        dof_mapping.second.shrink_to_fit();
        inco_face_dof_mapping.push_back(dof_mapping);
      } // if local
    }   // if incident
  }     // for incindent f

  auto inco_face_info_array = new INCOMING_FACE_INFO[inco_face_dof_mapping.size()];
  for (int i = 0; i < inco_face_dof_mapping.size(); ++i)
    inco_face_info_array[i].Setup(inco_face_dof_mapping[i]);

  so_cell_inco_face_dof_indices.push_back(inco_face_info_array);
}

void
AAH_FLUDSCommonData::InitializeBetaElements(const SPDS& spds, int tag_index /*=0*/)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();
  const chi_mesh::sweep_management::SPLS& spls = spds.GetSPLS();

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // The first two major steps here are: Send delayed successor information
  // and Receive delayed predecessor information. The send portion is done
  // first because the delayed information does not follow the
  // Task Dependency Graph and hence when a location receives its delayed
  // information, the information might not have been sent yet.

  //=============================================== Send delayed successor
  // information
  const auto& location_successors = spds.GetLocationSuccessors();
  const auto& delayed_location_successors = spds.GetDelayedLocationSuccessors();
  std::vector<MPI_Request> send_requests;
  send_requests.resize(location_successors.size(), MPI_Request());
  std::vector<std::vector<int>> multi_face_indices(location_successors.size(), std::vector<int>());
  for (int deplocI = 0; deplocI < location_successors.size(); deplocI++)
  {
    int locJ = location_successors[deplocI];

    std::vector<int>::const_iterator delayed_successor =
      std::find(delayed_location_successors.begin(), delayed_location_successors.end(), locJ);
    if ((delayed_successor == delayed_location_successors.end())) continue;

    std::vector<CompactCellView> cell_views = deplocI_cell_views[deplocI];

    SerializeCellInfo(cell_views, multi_face_indices[deplocI], deplocI_face_dof_count[deplocI]);

    MPI_Isend(multi_face_indices[deplocI].data(),
              static_cast<int>(multi_face_indices[deplocI].size()),
              MPI_INT,
              locJ,
              101 + tag_index,
              Chi::mpi.comm,
              &send_requests[deplocI]);

    // TODO: Watch eager limits on sent data

    deplocI_cell_views[deplocI].clear();
    deplocI_cell_views[deplocI].shrink_to_fit();
  }

  //=============================================== Receive delayed predecessor
  //                                                information
  const auto& delayed_location_dependencies = spds.GetDelayedLocationDependencies();
  delayed_prelocI_cell_views.resize(delayed_location_dependencies.size(),
                                    std::vector<CompactCellView>());
  delayed_prelocI_face_dof_count.resize(delayed_location_dependencies.size(), 0);
  for (int prelocI = 0; prelocI < delayed_location_dependencies.size(); prelocI++)
  {
    int locJ = delayed_location_dependencies[prelocI];

    MPI_Status probe_status;
    MPI_Probe(locJ, 101 + tag_index, Chi::mpi.comm, &probe_status);

    int amount_to_receive = 0;
    MPI_Get_count(&probe_status, MPI_INT, &amount_to_receive);

    std::vector<int> face_indices;
    face_indices.resize(amount_to_receive, 0);

    MPI_Recv(face_indices.data(),
             amount_to_receive,
             MPI_INT,
             locJ,
             101 + tag_index,
             Chi::mpi.comm,
             MPI_STATUS_IGNORE);

    DeSerializeCellInfo(
      delayed_prelocI_cell_views[prelocI], &face_indices, delayed_prelocI_face_dof_count[prelocI]);
  }

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // The next two operations is to receive predecessor information followed
  // by the sending of successor information. The receives are blocking but
  // will cause a dead lock because the successor/predecessor combination
  // follows the TDG.

  //=============================================== Receive predecessor
  //                                                information
  const auto& location_dependencies = spds.GetLocationDependencies();
  prelocI_cell_views.resize(location_dependencies.size(), std::vector<CompactCellView>());
  prelocI_face_dof_count.resize(location_dependencies.size(), 0);
  for (int prelocI = 0; prelocI < location_dependencies.size(); prelocI++)
  {
    int locJ = location_dependencies[prelocI];

    MPI_Status probe_status;
    MPI_Probe(locJ, 101 + tag_index, Chi::mpi.comm, &probe_status);

    int amount_to_receive = 0;
    MPI_Get_count(&probe_status, MPI_INT, &amount_to_receive);

    std::vector<int> face_indices;
    face_indices.resize(amount_to_receive, 0);

    MPI_Recv(face_indices.data(),
             amount_to_receive,
             MPI_INT,
             locJ,
             101 + tag_index,
             Chi::mpi.comm,
             MPI_STATUS_IGNORE);

    DeSerializeCellInfo(
      prelocI_cell_views[prelocI], &face_indices, prelocI_face_dof_count[prelocI]);
  }

  //=============================================== Send successor information
  for (int deplocI = 0; deplocI < location_successors.size(); deplocI++)
  {
    int locJ = location_successors[deplocI];

    auto delayed_successor =
      std::find(delayed_location_successors.begin(), delayed_location_successors.end(), locJ);
    if ((delayed_successor != delayed_location_successors.end())) continue;

    std::vector<CompactCellView> cell_views = deplocI_cell_views[deplocI];

    SerializeCellInfo(cell_views, multi_face_indices[deplocI], deplocI_face_dof_count[deplocI]);

    MPI_Isend(multi_face_indices[deplocI].data(),
              static_cast<int>(multi_face_indices[deplocI].size()),
              MPI_INT,
              locJ,
              101 + tag_index,
              Chi::mpi.comm,
              &send_requests[deplocI]);

    // TODO: Watch eager limits on sent data

    deplocI_cell_views[deplocI].clear();
    deplocI_cell_views[deplocI].shrink_to_fit();
  }

  //================================================== Verify sends completed
  for (int deplocI = 0; deplocI < location_successors.size(); deplocI++)
    MPI_Wait(&send_requests[deplocI], MPI_STATUS_IGNORE);
  multi_face_indices.clear();
  multi_face_indices.shrink_to_fit();

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // In the next process we loop over cells in the sweep order and perform
  // the non-local face mappings. This is dependent on having the compact
  // cellviews on the partition interfaces.

  //================================================== Loop over cells in sorder
  for (int csoi = 0; csoi < spls.item_id.size(); csoi++)
  {
    int cell_local_index = spls.item_id[csoi];
    const auto& cell = grid.local_cells[cell_local_index];

    NonLocalIncidentMapping(cell, spds);
  } // for csoi

  deplocI_cell_views.clear();
  deplocI_cell_views.shrink_to_fit();

  prelocI_cell_views.clear();
  prelocI_cell_views.shrink_to_fit();

  delayed_prelocI_cell_views.clear();
  delayed_prelocI_cell_views.shrink_to_fit();

  //================================================== Clear unneccesary data
  auto empty_vector = std::vector<std::vector<CompactCellView>>(0);
  deplocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  prelocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  delayed_prelocI_cell_views.swap(empty_vector);
}

void
AAH_FLUDSCommonData::SerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                       std::vector<int>& face_indices,
                                       int num_face_dofs)
{
  const size_t num_cells = cell_views.size();

  //======================== First entry is number of face dofs
  face_indices.push_back(num_face_dofs);

  //======================== Second entry is amount of cells
  face_indices.push_back(static_cast<int>(num_cells));

  //======================== Third entry is negative global cell index
  // Each time a negative entry occurs it denotes a cell face but
  // the actual number is -cell_g_index-1. The offset is necessary
  // for evaluating the negative. The offset is restored during the
  // deserialization process.
  // It is followed by a positive number which is the store location
  // of the face
  for (size_t c = 0; c < num_cells; c++)
  {
    int glob_index = -cell_views[c].first - 1;

    std::vector<CompactFaceView>& cell_face_views = cell_views[c].second;

    size_t num_faces = cell_face_views.size();
    for (size_t f = 0; f < num_faces; f++)
    {
      face_indices.push_back(glob_index);
      face_indices.push_back(cell_face_views[f].first);
      std::vector<uint64_t>& face_vertices = cell_face_views[f].second;

      size_t num_verts = face_vertices.size();
      for (int fi = 0; fi < num_verts; fi++)
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
  num_face_dofs = (*face_indices)[0];
  int num_cells = (*face_indices)[1];

  cell_views.resize(num_cells);

  int k = 2;
  int last_cell = -1;
  int c = -1; // cell counter
  int f = -1;
  int v = -1;
  while (k < face_indices->size())
  {
    int entry = (*face_indices)[k];
    //================================= Cell/Face indicator
    if (entry < 0)
    {
      if (-entry != last_cell)
      {
        cell_views.emplace_back();
        c++;
        cell_views[c].first = -entry - 1;

        cell_views[c].second.emplace_back();
        f = 0;

        v = 0;
        last_cell = -entry;

        cell_views[c].second[f].first = (*face_indices)[k + 1];
        k++;
      }
      else
      {
        cell_views[c].second.emplace_back();
        f++;
        v = 0;

        cell_views[c].second[f].first = (*face_indices)[k + 1];
        k++;
      }
    }
    //================================= Face vertex
    else
    {
      cell_views[c].second[f].second.push_back(entry);
      v++;
    }
    k++;
  } // while k
}

void
AAH_FLUDSCommonData::NonLocalIncidentMapping(const chi_mesh::Cell& cell, const SPDS& spds)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  for (short f = 0; f < cell.faces_.size(); f++)
  {
    const CellFace& face = cell.faces_[f];
    const auto& orientation = spds.CellFaceOrientations()[cell.local_id_][f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (orientation == FaceOrientation::INCOMING)
    {
      if ((face.has_neighbor_) and (!face.IsNeighborLocal(grid)))
      {
        //============================== Find prelocI
        int locJ = face.GetNeighborPartitionID(grid);
        int prelocI = spds.MapLocJToPrelocI(locJ);

        // ###########################################################
        if (prelocI >= 0)
        {
          //============================== Find the cell in prelocI cell views
          int ass_cell = -1;
          for (int c = 0; c < prelocI_cell_views[prelocI].size(); c++)
          {
            if (prelocI_cell_views[prelocI][c].first == face.neighbor_id_)
            {
              ass_cell = c;
              break;
            }
          }
          if (ass_cell < 0)
          {
            Chi::log.LogAll() << "Required predecessor cell not located in call to"
                              << " InitializeBetaElements. locJ=" << locJ << " prelocI=" << prelocI
                              << " cell=" << face.neighbor_id_;
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Find associated face
          std::set<int> cfvids(face.vertex_ids_.begin(), face.vertex_ids_.end());
          CompactCellView* adj_cell_view = &prelocI_cell_views[prelocI][ass_cell];
          int ass_face = -1, af = -1;
          for (auto& adj_face : adj_cell_view->second)
          {
            ++af;
            bool face_matches = true;

            std::set<int> afvids(adj_face.second.begin(), adj_face.second.end());

            if (cfvids != afvids) face_matches = false;

            if (face_matches)
            {
              ass_face = af;
              break;
            }
          }
          if (ass_face < 0)
          {
            Chi::log.LogAll() << "Associated face not found in call to InitializeBetaElements";
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Map dofs
          std::pair<int, std::vector<int>> dof_mapping;
          dof_mapping.first = adj_cell_view->second[ass_face].first;
          std::vector<uint64_t>* ass_face_verts = &adj_cell_view->second[ass_face].second;
          for (int fv = 0; fv < face.vertex_ids_.size(); fv++)
          {
            bool match_found = false;
            for (int afv = 0; afv < ass_face_verts->size(); afv++)
            {
              if (face.vertex_ids_[fv] == ass_face_verts->operator[](afv))
              {
                match_found = true;
                dof_mapping.second.push_back(afv);
                break;
              }
            }

            if (!match_found)
            {
              Chi::log.LogAll() << "Associated vertex not found in call to "
                                   "InitializeBetaElements";
              Chi::Exit(EXIT_FAILURE);
            }
          }

          //============================== Push back final face info
          std::pair<int, std::pair<int, std::vector<int>>> inc_face_prelocI_info;
          inc_face_prelocI_info.first = prelocI;
          inc_face_prelocI_info.second = std::pair<int, std::vector<int>>(dof_mapping);

          std::pair<int, std::pair<int, std::vector<int>>> empty_delayed_info;
          empty_delayed_info.first = prelocI;

          nonlocal_inc_face_prelocI_slot_dof.push_back(inc_face_prelocI_info);
          delayed_nonlocal_inc_face_prelocI_slot_dof.push_back(empty_delayed_info);
        } // If not delayed predecessor
        // ###########################################################
        else
        {
          int delayed_preLocI = abs(prelocI) - 1;
          //============================== Find the cell in prelocI cell views
          int ass_cell = -1;
          for (int c = 0; c < delayed_prelocI_cell_views[delayed_preLocI].size(); c++)
          {
            if (delayed_prelocI_cell_views[delayed_preLocI][c].first == face.neighbor_id_)
            {
              ass_cell = c;
              break;
            }
          }
          if (ass_cell < 0)
          {
            Chi::log.LogAll() << "Required predecessor cell not located in call to"
                              << " InitializeBetaElements. locJ=" << locJ
                              << " delayed prelocI=" << delayed_preLocI
                              << " cell=" << face.neighbor_id_;
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Find associated face
          CompactCellView* adj_cell_view = &delayed_prelocI_cell_views[delayed_preLocI][ass_cell];
          int ass_face = -1;
          for (int af = 0; af < adj_cell_view->second.size(); af++)
          {
            bool face_matches = true;
            for (int afv = 0; afv < adj_cell_view->second[af].second.size(); afv++)
            {
              bool match_found = false;
              for (int fv = 0; fv < face.vertex_ids_.size(); fv++)
              {
                if (adj_cell_view->second[af].second[afv] == face.vertex_ids_[fv])
                {
                  match_found = true;
                  break;
                }
              }

              if (!match_found)
              {
                face_matches = false;
                break;
              }
            }

            if (face_matches)
            {
              ass_face = af;
              break;
            }
          }
          if (ass_face < 0)
          {
            Chi::log.LogAll() << "Associated face not found in call to InitializeBetaElements";
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Map dofs
          std::pair<int, std::vector<int>> dof_mapping;
          dof_mapping.first = adj_cell_view->second[ass_face].first;
          std::vector<uint64_t>* ass_face_verts = &adj_cell_view->second[ass_face].second;
          for (int fv = 0; fv < face.vertex_ids_.size(); fv++)
          {
            bool match_found = false;
            for (int afv = 0; afv < ass_face_verts->size(); afv++)
            {
              if (face.vertex_ids_[fv] == ass_face_verts->operator[](afv))
              {
                match_found = true;
                dof_mapping.second.push_back(afv);
                break;
              }
            }

            if (!match_found)
            {
              Chi::log.LogAll() << "Associated vertex not found in call to "
                                   "InitializeBetaElements";
              Chi::Exit(EXIT_FAILURE);
            }
          }

          //============================== Push back final face info
          std::pair<int, std::pair<int, std::vector<int>>> inc_face_prelocI_info;
          inc_face_prelocI_info.first = delayed_preLocI;
          inc_face_prelocI_info.second = std::pair<int, std::vector<int>>(dof_mapping);

          std::pair<int, std::pair<int, std::vector<int>>> delayed_info;
          delayed_info.first = -(delayed_preLocI + 1);

          delayed_nonlocal_inc_face_prelocI_slot_dof.push_back(inc_face_prelocI_info);
          nonlocal_inc_face_prelocI_slot_dof.push_back(delayed_info);
        } // If delayed predecessor

      } // if not local and not boundary
    }   // if incident
  }     // for incindent f
}

} // namespace chi_mesh::sweep_management
