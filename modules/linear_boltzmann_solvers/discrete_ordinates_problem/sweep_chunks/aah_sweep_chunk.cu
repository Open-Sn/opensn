// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "caliper/cali.h"
#include "caribou/caribou.h"

namespace crb = caribou;

namespace opensn
{

namespace gpu_kernel
{

/// Kernel performing the sweep
__global__ void
UnsaturatedKernel(Arguments args, const std::uint32_t* level, double* saved_psi)
{
  // reference to indexes
  const unsigned int& cell_idx = blockIdx.x;
  const unsigned int& angle_idx = threadIdx.y;
  const unsigned int& group_idx = threadIdx.x;
  unsigned int angle_group_idx = angle_idx * args.groupset_size + group_idx;
  // get cell view
  std::uint32_t cell_local_idx = level[cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  // skip when cell has no nodes
  if (cell.num_nodes == 0)
    return;
  // get direction view and number of moments
  std::uint32_t num_moments;
  std::uint32_t direction_num = args.directions[angle_idx];
  DirectionView direction;
  {
    QuadratureView quadrature(args.quad_data);
    num_moments = quadrature.num_moments;
    quadrature.GetDirectionView(direction, direction_num);
  }
  // get pointer to the corresponding FLUDS index
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  // launch the kernel
  sweep_spec_map[cell.num_nodes - 1](
    args, cell, direction, cell_edge_data, angle_group_idx, group_idx, num_moments, saved_psi);
}

/// Kernel performing the sweep
__global__ void
SaturatedKernel(Arguments args, const std::uint32_t* level, double* saved_psi)
{
  // reference to indexes
  const unsigned int& cell_idx = blockIdx.x;
  // get cell view
  std::uint32_t cell_local_idx = level[cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  // skip when cell has no nodes
  if (cell.num_nodes == 0)
    return;
  // get pointer to the corresponding FLUDS index
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  // loop on each group and angle
  for (unsigned int angle_group_idx = threadIdx.x; angle_group_idx < args.flud_data.stride_size;
       angle_group_idx += blockDim.x)
  {
    // get group and angle index
    std::uint32_t group_idx = angle_group_idx % args.groupset_size;
    std::uint32_t angle_idx = angle_group_idx / args.groupset_size;
    // get direction view and number of moments
    std::uint32_t num_moments;
    std::uint32_t direction_num = args.directions[angle_idx];
    DirectionView direction;
    {
      QuadratureView quadrature(args.quad_data);
      num_moments = quadrature.num_moments;
      quadrature.GetDirectionView(direction, direction_num);
    }
    // launch the kernel
    sweep_spec_map[cell.num_nodes - 1](
      args, cell, direction, cell_edge_data, angle_group_idx, group_idx, num_moments, saved_psi);
  }
}

const unsigned int threshold = 256;

} // namespace gpu_kernel

void
AAHSweepChunk::GPUSweep(AngleSet& angle_set)
{
  // prepare arguments
  AAHD_FLUDS& fluds = dynamic_cast<AAHD_FLUDS&>(angle_set.GetFLUDS());
  gpu_kernel::Arguments args(problem_, groupset_, angle_set, fluds, surface_source_active_);
  // allocate memory for saved angular flux
  crb::DeviceMemory<double> saved_psi;
  if (save_angular_flux_)
  {
    auto* mesh_carrier_ptr = reinterpret_cast<MeshCarrier*>(problem_.GetCarrier(2));
    saved_psi =
      crb::DeviceMemory<double>(mesh_carrier_ptr->num_nodes_total * fluds.GetStrideSize());
  }
  // retrieve SPDS levels
  const AAH_SPDS& spds = dynamic_cast<const AAH_SPDS&>(angle_set.GetSPDS());
  const std::vector<std::vector<std::uint32_t>>& levelized_spls = spds.GetLevelizedLocalSubgrid();
  // loop over each level based on saturation status
  unsigned int block_size_x = groupset_.groups.size(), block_size_y = angle_set.GetNumAngles();
  if (block_size_x * block_size_y <= gpu_kernel::threshold)
  {
    for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
    {
      // perform the sweep on device
      std::size_t level_size = levelized_spls[level].size();
      const std::uint32_t* level_data = spds.GetDeviceLevelVector(level);
      gpu_kernel::UnsaturatedKernel<<<level_size, {block_size_x, block_size_y}>>>(
        args, level_data, saved_psi.get());
    }
  }
  else
  {
    for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
    {
      // perform the sweep on device
      std::size_t level_size = levelized_spls[level].size();
      const std::uint32_t* level_data = spds.GetDeviceLevelVector(level);
      gpu_kernel::SaturatedKernel<<<level_size, 256>>>(args, level_data, saved_psi.get());
    }
  }
  // clean up memory after sweep and copy data back to the host
  fluds.CleanUpAfterSweep(*problem_.GetGrid(), angle_set);
  auto* phi = reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(1));
  phi->CopyFromDevice();
  auto* outflow = reinterpret_cast<OutflowCarrier*>(problem_.GetCarrier(1));
  outflow->AccumulateBack(cell_transport_views_);
  outflow->Reset();
  // save angular flux to unknowns manager
  if (save_angular_flux_)
  {
    // copy saved psi back to host
    crb::HostVector<double> psi_host(saved_psi.size());
    crb::copy(psi_host, saved_psi, psi_host.size());
    // loop for each cell in the mesh
    auto* mesh_carrier = reinterpret_cast<MeshCarrier*>(problem_.GetCarrier(2));
    for (const Cell& cell : grid_->local_cells)
    {
      // get pointer to the cell's angular fluxes
      double* dst_psi =
        &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
      double* src_psi =
        psi_host.data() + mesh_carrier->saved_psi_offset[cell.local_id] * fluds.GetStrideSize();
      // get number of cell nodes
      std::uint32_t cell_num_nodes = discretization_.GetCellMapping(cell).GetNumNodes();
      // loop for each cell node
      for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
      {
        // loop for each angle
        for (std::uint32_t as_ss_idx = 0; as_ss_idx < angle_set.GetNumAngles(); ++as_ss_idx)
        {
          auto direction_num = angle_set.GetAngleIndices()[as_ss_idx];
          // compute dst and src corresponding to the direction
          double* dst = dst_psi + direction_num * groupset_group_stride_;
          double* src = src_psi + as_ss_idx * groupset_.groups.size();
          // copy the flux for each group
          std::memcpy(dst, src, groupset_group_stride_ * sizeof(double));
        }
        // move src and dts to next node
        dst_psi += groupset_angle_group_stride_;
        src_psi += fluds.GetStrideSize();
      }
    }
  }
}

} // namespace opensn
