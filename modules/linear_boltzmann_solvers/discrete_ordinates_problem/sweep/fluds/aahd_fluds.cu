// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"

namespace opensn
{

template <typename T>
inline static char*
write_to(char* buffer, const T& src)
{
  T* dest = reinterpret_cast<T*>(buffer);
  *(dest++) = src;
  return reinterpret_cast<char*>(dest);
}

static void
record_incoming_flux(int f,
                     std::size_t num_face_nodes,
                     FaceOrientation face_orientation,
                     FloodFace& fluds_face,
                     const Cell& cell,
                     const CellLBSView& cell_transport_view,
                     int cell_local_idx,
                     AngleSet& angle_set,
                     int spls_idx,
                     int gs_gi,
                     bool is_surface_source_active,
                     int& in_face_counter,
                     int& preloc_face_counter,
                     AAH_FLUDS& fluds)
{
  std::size_t nangles = angle_set.GetAngleIndices().size();
  fluds_face.flux = NDArray<FloodPtr, 2>({num_face_nodes, nangles}, FloodPtr());
  // increase counters
  bool is_local_face = cell_transport_view.IsFaceLocal(f);
  bool is_boundary_face = not cell.faces[f].has_neighbor;
  if (is_local_face)
  {
    ++in_face_counter;
  }
  else if (not is_boundary_face)
  {
    ++preloc_face_counter;
  }
  // record pointer to FLUDS of a given angle-set index and a given face node
  for (std::size_t as_ss_idx = 0; as_ss_idx < nangles; ++as_ss_idx)
  {
    const std::uint32_t& direction_num = angle_set.GetAngleIndices()[as_ss_idx];
    for (int fi = 0; fi < num_face_nodes; ++fi)
    {
      const double* psi = nullptr;
      if (is_local_face)
      {
        psi = fluds.UpwindPsi(spls_idx, in_face_counter, fi, 0, as_ss_idx);
      }
      else if (not is_boundary_face)
      {
        psi = fluds.NLUpwindPsi(preloc_face_counter, fi, 0, as_ss_idx);
      }
      else
      {
        psi = angle_set.PsiBoundary(cell.faces[f].neighbor_id,
                                    direction_num,
                                    cell_local_idx,
                                    f,
                                    fi,
                                    gs_gi,
                                    is_surface_source_active);
      }
      fluds_face.flux(fi, as_ss_idx).incoming = psi;
    }
  }
}

static void
record_outgoing_flux(int f,
                     std::size_t num_face_nodes,
                     FaceOrientation face_orientation,
                     FloodFace& fluds_face,
                     const Cell& cell,
                     const CellLBSView& cell_transport_view,
                     int cell_local_idx,
                     AngleSet& angle_set,
                     int spls_idx,
                     int gs_gi,
                     bool is_surface_source_active,
                     int& out_face_counter,
                     int& deploc_face_counter,
                     AAH_FLUDS& fluds)
{
  std::size_t nangles = angle_set.GetAngleIndices().size();
  fluds_face.flux = NDArray<FloodPtr, 2>({num_face_nodes, nangles}, FloodPtr());
  // increase counters
  out_face_counter++;
  bool is_local_face = cell_transport_view.IsFaceLocal(f);
  bool is_boundary_face = not cell.faces[f].has_neighbor;
  bool is_reflecting_boundary_face =
    (is_boundary_face and angle_set.GetBoundaries()[cell.faces[f].neighbor_id]->IsReflecting());
  if (not is_boundary_face and not is_local_face)
  {
    ++deploc_face_counter;
  }
  // record pointer to FLUDS of a given angle-set index and a given face node
  for (std::size_t as_ss_idx = 0; as_ss_idx < nangles; ++as_ss_idx)
  {
    const std::uint32_t& direction_num = angle_set.GetAngleIndices()[as_ss_idx];
    for (int fi = 0; fi < num_face_nodes; ++fi)
    {
      double* psi = nullptr;
      if (is_local_face)
      {
        psi = fluds.OutgoingPsi(spls_idx, out_face_counter, fi, as_ss_idx);
      }
      else if (not is_boundary_face)
      {
        psi = fluds.NLOutgoingPsi(deploc_face_counter, fi, as_ss_idx);
      }
      else if (is_reflecting_boundary_face)
      {
        psi =
          angle_set.PsiReflected(cell.faces[f].neighbor_id, direction_num, cell_local_idx, f, fi);
      }
      else
      {
        psi = reinterpret_cast<double*>(UINTPTR_MAX);
      }
      fluds_face.flux(fi, as_ss_idx).outgoing = psi;
    }
  }
}

AAHD_FLUDS::AAHD_FLUDS(LBSProblem& lbs_problem,
                       const LBSGroupset& group_set,
                       AngleSet& angle_set,
                       bool is_surface_source_active)
{
  // get SPDS
  const SPDS& spds = angle_set.GetSPDS();
  const std::vector<std::vector<int>>& levelized_spls = spds.GetLevelizedLocalSubgrid();
  // get FLUDS
  AAH_FLUDS& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  // get groupset and angleset details
  groupset_size_ = group_set.groups.size();
  angleset_size_ = angle_set.GetAngleIndices().size();
  int gs_gi = group_set.groups.front().id;
  // get mesh and cell transportview
  MeshContinuum& mesh = *(lbs_problem.GetGrid());
  const SpatialDiscretization& discretization = lbs_problem.GetSpatialDiscretization();
  // allocate memory
  fluds_map_.reserve(levelized_spls.size());
  // loop over each level
  int preloc_face_counter = -1, deploc_face_counter = -1;
  for (int spls_idx = 0; const std::vector<int>& level : levelized_spls)
  {
    // compute the start position of cell data
    std::uint64_t cell_data_size = 0, face_data_size = 0, flux_data_size = 0;
    for (const int& cell_local_idx : level)
    {
      const Cell& cell = mesh.local_cells[cell_local_idx];
      const CellMapping& cell_mapping = discretization.GetCellMapping(cell);
      const std::vector<std::vector<int>>& face_node_mappings = cell_mapping.GetFaceNodeMappings();
      cell_data_size += sizeof(std::uint64_t);
      face_data_size += cell.faces.size() * sizeof(std::uint64_t);
      const std::vector<FaceOrientation>& face_orientations =
        spds.GetCellFaceOrientations()[cell_local_idx];
      for (int f = 0; f < cell.faces.size(); ++f)
      {
        std::size_t num_face_nodes = face_node_mappings[f].size();
        if (face_orientations[f] != FaceOrientation::PARALLEL)
        {
          flux_data_size += num_face_nodes * angleset_size_ * groupset_size_ * sizeof(double);
        }
      }
    }
    size_ = std::max(size_, cell_data_size + face_data_size + flux_data_size);
    // record offset for each group, angles and face node
    FloodLevel fluds_level;
    fluds_level.precomputed_offset.resize(cell_data_size + face_data_size);
    std::uint64_t face_offset = cell_data_size;
    std::int64_t flux_offset = cell_data_size + face_data_size;
    char* cell_data = fluds_level.precomputed_offset.data();
    char* face_data = cell_data + cell_data_size;
    for (const int& cell_local_idx : level)
    {
      // record cell data
      cell_data = write_to<std::uint64_t>(cell_data, face_offset);
      // update face_offset
      const Cell& cell = mesh.local_cells[cell_local_idx];
      const CellMapping& cell_mapping = discretization.GetCellMapping(cell);
      const std::vector<std::vector<int>>& face_node_mappings = cell_mapping.GetFaceNodeMappings();
      const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_idx];
      face_offset += cell.faces.size() * sizeof(std::uint64_t);
      // record face data
      for (int f = 0; f < cell.faces.size(); ++f)
      {
        std::size_t num_face_nodes = face_node_mappings[f].size();
        if (face_orientations[f] == FaceOrientation::INCOMING)
        {
          face_data = write_to<std::int64_t>(face_data, -flux_offset);
          flux_offset += num_face_nodes * angleset_size_ * groupset_size_ * sizeof(double);
        }
        else if (face_orientations[f] == FaceOrientation::OUTGOING)
        {
          face_data = write_to<std::int64_t>(face_data, flux_offset);
          flux_offset += num_face_nodes * angleset_size_ * groupset_size_ * sizeof(double);
        }
        else
        {
          face_data = write_to<std::int64_t>(face_data, 0);
        }
      }
    }
    // record pointer to the flux in the FLUDS of each face
    fluds_level.cells.reserve(level.size());
    for (const int& cell_local_idx : level)
    {
      FloodCell fluds_cell;
      const Cell& cell = mesh.local_cells[cell_local_idx];
      const CellMapping& cell_mapping = discretization.GetCellMapping(cell);
      const std::vector<std::vector<int>>& face_node_mappings = cell_mapping.GetFaceNodeMappings();
      const CellLBSView& cell_transport_view = lbs_problem.GetCellTransportViews()[cell_local_idx];
      const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_idx];
      fluds_cell.reserve(cell.faces.size());
      int in_face_counter = -1, out_face_counter = -1;
      for (int f = 0; f < cell.faces.size(); ++f)
      {
        FloodFace fluds_face;
        if (face_orientations[f] == FaceOrientation::INCOMING)
        {
          record_incoming_flux(f,
                               face_node_mappings[f].size(),
                               face_orientations[f],
                               fluds_face,
                               cell,
                               cell_transport_view,
                               cell_local_idx,
                               angle_set,
                               spls_idx,
                               gs_gi,
                               is_surface_source_active,
                               in_face_counter,
                               preloc_face_counter,
                               fluds);
        }
        else if (face_orientations[f] == FaceOrientation::OUTGOING)
        {
          record_outgoing_flux(f,
                               face_node_mappings[f].size(),
                               face_orientations[f],
                               fluds_face,
                               cell,
                               cell_transport_view,
                               cell_local_idx,
                               angle_set,
                               spls_idx,
                               gs_gi,
                               is_surface_source_active,
                               out_face_counter,
                               deploc_face_counter,
                               fluds);
        }
        fluds_cell.push_back(std::move(fluds_face));
      }
      fluds_level.cells.push_back(std::move(fluds_cell));
      ++spls_idx;
    }
    // record FLUDS pointers for each level
    fluds_map_.push_back(std::move(fluds_level));
  }
  // allocate memory for host vector
  host_buffer_.resize(size_);
  device_buffer_ = crb::DeviceMemory<char>(size_);
}

void
AAHD_FLUDS::CopyToDevice(int level)
{
  // write the offset data
  char* host_data = host_buffer_.data();
  char* precomputed_offset = host_data;
  FloodLevel& level_data = fluds_map_[level];
  precomputed_offset = std::copy(
    level_data.precomputed_offset.begin(), level_data.precomputed_offset.end(), precomputed_offset);
  host_data = precomputed_offset;
  // copy flux of each cells in the level to the host buffer
  double* flux_data = reinterpret_cast<double*>(host_data);
  for (FloodCell& cell_data : level_data.cells)
  {
    for (FloodFace& face_data : cell_data)
    {
      std::size_t num_face_nodes = face_data.flux.dimension()[0];
      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        for (std::size_t as_ss_idx = 0; as_ss_idx < angleset_size_; ++as_ss_idx)
        {
          const double* incoming_flux_ptr = face_data.flux(fi, as_ss_idx).incoming;
          double* outgoing_flux_ptr = face_data.flux(fi, as_ss_idx).outgoing;
          if (incoming_flux_ptr != nullptr)
          {
            flux_data = std::copy(incoming_flux_ptr, incoming_flux_ptr + groupset_size_, flux_data);
          }
          else if (outgoing_flux_ptr != nullptr)
          {
            std::fill(flux_data, flux_data + groupset_size_, 0.0);
            flux_data += groupset_size_;
          }
        }
      }
    }
  }
  host_data = reinterpret_cast<char*>(flux_data);
  // copy from host buffer to device buffer
  crb::copy(device_buffer_, host_buffer_, size_);
}

void
AAHD_FLUDS::CopyFromDevice(int level)
{
  // copy data from device buffer back to the host buffer
  crb::copy(host_buffer_, device_buffer_, size_);
  // copy flux of each cells in the level to the host buffer
  FloodLevel& level_data = fluds_map_[level];
  const double* flux_data =
    reinterpret_cast<const double*>(host_buffer_.data() + level_data.precomputed_offset.size());
  for (FloodCell& cell_data : level_data.cells)
  {
    for (std::uint32_t face_idx = 0; FloodFace & face_data : cell_data)
    {
      std::size_t num_face_nodes = face_data.flux.dimension()[0];
      for (std::size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        for (std::size_t as_ss_idx = 0; as_ss_idx < angleset_size_; ++as_ss_idx)
        {
          const double* incoming_flux_ptr = face_data.flux(fi, as_ss_idx).incoming;
          double* outgoing_flux_ptr = face_data.flux(fi, as_ss_idx).outgoing;
          if (incoming_flux_ptr != nullptr ||
              outgoing_flux_ptr == reinterpret_cast<double*>(UINTPTR_MAX))
          {
            flux_data += groupset_size_;
          }
          else if (outgoing_flux_ptr != nullptr)
          {
            std::copy(flux_data, flux_data + groupset_size_, outgoing_flux_ptr);
            flux_data += groupset_size_;
          }
        }
      }
      face_idx++;
    }
  }
}

} // namespace opensn
