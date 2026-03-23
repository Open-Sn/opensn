// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/buffer.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include <utility>
#include <type_traits>

namespace opensn::gpu_kernel
{

template <std::uint32_t... I, class F>
__device__ inline void
ForDOFs1ToN(std::integer_sequence<std::uint32_t, I...>, F&& f)
{
  (f(std::integral_constant<std::uint32_t, I + 1>{}), ...);
}

template <class F>
__device__ inline void
ForDOFs1ToMax(F&& f)
{
  ForDOFs1ToN(std::make_integer_sequence<std::uint32_t, LBSProblem::max_dofs_gpu>{},
              std::forward<F>(f));
}

template <SweepType t, class... Args>
__device__ inline void
SweepDispatch(std::uint32_t n, Args&&... args)
{
  bool done = false;
  ForDOFs1ToMax(
    [&](auto dof_c)
    {
      constexpr std::uint32_t dof = decltype(dof_c)::value;
      if (!done && n == dof)
      {
        gpu_kernel::Sweep<dof, t>(std::forward<Args>(args)...);
        done = true;
      }
    });
}

template <SweepType t>
__global__ void
SweepKernel(Arguments<t> args,
            const std::uint32_t* cell_local_ids,
            unsigned int num_ready_cells,
            double* saved_psi)
{
  unsigned int cell_idx = threadIdx.y + blockDim.y * blockIdx.y;
  unsigned int angle_group_idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (cell_idx >= num_ready_cells || angle_group_idx >= args.flud_data.stride_size)
    return;
  unsigned int angle_idx = angle_group_idx / args.groupset_size;
  unsigned int group_idx = angle_group_idx - angle_idx * args.groupset_size;
  const std::uint32_t cell_local_idx = cell_local_ids[cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  if (cell.num_nodes == 0)
    return;
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  std::uint32_t num_moments;
  std::uint32_t direction_num = args.directions[angle_idx];
  DirectionView direction;
  {
    QuadratureView quadrature(args.quad_data);
    num_moments = quadrature.num_moments;
    quadrature.GetDirectionView(direction, direction_num);
  }
  opensn::gpu_kernel::SweepDispatch<t>(cell.num_nodes,
                                       args,
                                       cell,
                                       direction,
                                       cell_edge_data,
                                       angle_group_idx,
                                       group_idx,
                                       num_moments,
                                       saved_psi);
}

} // namespace opensn::gpu_kernel