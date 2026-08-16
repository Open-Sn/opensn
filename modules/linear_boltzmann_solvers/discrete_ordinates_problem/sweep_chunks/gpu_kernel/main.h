// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/buffer.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "caribou/main.hpp"
#include <utility>
#include <type_traits>

namespace crb = caribou;

namespace opensn::gpu_kernel
{

template <std::uint32_t... I, class F>
__CRB_DEVICE_FUNC__ void
ForDOFs1ToN(std::integer_sequence<std::uint32_t, I...> sequence, F& f)
{
  static_cast<void>(sequence);
  (f(std::integral_constant<std::uint32_t, I + 1>{}), ...);
}

template <class F>
__CRB_DEVICE_FUNC__ void
ForDOFs1ToMax(F& f)
{
  ForDOFs1ToN(std::make_integer_sequence<std::uint32_t, LBSProblem::max_dofs_gpu>{}, f);
}

template <SweepKind k, class... Args>
__CRB_DEVICE_FUNC__ void
SweepDispatch(std::uint32_t n, Args&... args)
{
  bool done = false;
  auto dispatch = [&](auto dof_c)
  {
    constexpr std::uint32_t dof = decltype(dof_c)::value;
    if (!done && n == dof)
    {
      gpu_kernel::Sweep<dof, k>(args...);
      done = true;
    }
  };
  ForDOFs1ToMax(dispatch);
}

template <SweepKind k>
__CRB_GLOBAL_FUNC__ void
SweepKernel(Arguments<k> args,
            const std::uint32_t* cells_to_sweep,
            unsigned int num_cells,
            double* saved_psi,
            bool cells_on_fast_axis)
{
  // Whichever of {num_cells, stride_size} is larger drives the warp-contiguous fast axis,
  // so warps stay full whether the level or the groupset is the small one. The host side
  // picks `cells_on_fast_axis` and sizes the block/grid to match.
#if defined(__NVCC__) || defined(__HIPCC__)
  unsigned int fast_idx = threadIdx.x + blockDim.x * blockIdx.x;
  unsigned int slow_idx = threadIdx.y + blockDim.y * blockIdx.y;
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  auto work_index = ::sycl::ext::oneapi::this_work_item::get_nd_item<3>();
  unsigned int fast_idx = work_index.get_global_id(2);
  unsigned int slow_idx = work_index.get_global_id(1);
#endif
  unsigned int cell_idx = cells_on_fast_axis ? fast_idx : slow_idx;
  unsigned int angle_group_idx = cells_on_fast_axis ? slow_idx : fast_idx;
  if (cell_idx >= num_cells || angle_group_idx >= args.flud_data.stride_size)
    return;
  unsigned int angle_idx = angle_group_idx / args.groupset_size;
  unsigned int group_idx = angle_group_idx - angle_idx * args.groupset_size;
  const std::uint32_t cell_local_idx = cells_to_sweep[cell_idx];
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
  opensn::gpu_kernel::SweepDispatch<k>(cell.num_nodes,
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
