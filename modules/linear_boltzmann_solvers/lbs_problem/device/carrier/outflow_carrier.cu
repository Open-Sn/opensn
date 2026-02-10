// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/outflow_carrier.h"

namespace opensn
{

static inline std::uint64_t
merge(const std::uint32_t& high, const std::uint32_t& low)
{
  return (static_cast<std::uint64_t>(high) << 32) | low;
}

static inline std::pair<std::uint32_t, std::uint32_t>
unpack(std::uint64_t merged)
{
  return {static_cast<std::uint32_t>(merged >> 32),
          static_cast<std::uint32_t>(merged & 0xFFFFFFFF)};
}

OutflowCarrier::OutflowCarrier(LBSProblem& lbs_problem)
{
  // get information
  MeshContinuum& mesh = *(lbs_problem.GetGrid());
  num_groups = lbs_problem.GetNumGroups();
  // compute the size and offset for the outflux of boundary faces
  std::uint64_t alloc_size = 0;
  for (const Cell& cell : mesh.local_cells)
  {
    // loop over each face
    std::size_t cell_num_faces = cell.faces.size();
    for (int f = 0; f < cell_num_faces; ++f)
    {
      // skip if face is not boundary face
      const CellFace& face = cell.faces[f];
      bool is_boundary_face = not face.has_neighbor;
      if (not is_boundary_face)
      {
        continue;
      }
      // record cell index and face index to the map
      std::uint64_t cell_face_idx = merge(cell.local_id, f);
      outflow_map[cell_face_idx] = alloc_size;
      alloc_size += num_groups;
    }
  }
  // allocate memory on CPU and GPU
  std::uint64_t size = sizeof(double) * alloc_size;
  host_memory_.reserve(size);
  host_memory_.resize(size);
  device_memory_ = crb::DeviceMemory<char>(size);
  Reset();
}

void
OutflowCarrier::AccumulateBack(std::vector<CellLBSView>& cell_transport_views)
{
  // copy data back from the GPU
  crb::copy(host_memory_, device_memory_, host_memory_.size());
  // add outflow to each boundary face
  for (auto [cell_face_idx, offset] : outflow_map)
  {
    // get corresponding cell view
    auto [cell_idx, face_idx] = unpack(cell_face_idx);
    CellLBSView& cell_view = cell_transport_views[cell_idx];
    double* src_outflow = reinterpret_cast<double*>(host_memory_.data()) + offset;
    // add outflow for each group
    for (std::uint32_t g = 0; g < num_groups; ++g)
    {
      cell_view.AddOutflow(face_idx, g, src_outflow[g]);
    }
  }
}

void
OutflowCarrier::Reset(void)
{
  std::fill(host_memory_.begin(), host_memory_.end(), '\0');
  crb::copy(device_memory_, host_memory_, host_memory_.size());
}

std::uint64_t
OutflowCarrier::GetOffset(const std::uint32_t& cell_local_idx, const std::uint32_t& face_idx)
{
  std::uint64_t cell_face_idx = merge(cell_local_idx, face_idx);
  return outflow_map.at(cell_face_idx);
}

} // namespace opensn
