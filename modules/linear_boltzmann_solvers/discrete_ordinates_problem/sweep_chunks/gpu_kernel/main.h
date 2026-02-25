// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/buffer.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/dof_limits.h"

namespace opensn::gpu_kernel
{

template <std::size_t From, std::size_t To, std::size_t... IntSeq>
constexpr auto
MakeIndexSequenceFromRangeImpl(std::index_sequence<IntSeq...>)
{
  return std::index_sequence<(IntSeq + From)...>{};
}

template <std::size_t From, std::size_t To>
using MakeIndexSequenceFromRange =
  decltype(MakeIndexSequenceFromRangeImpl<From, To>(std::make_index_sequence<To - From>{}));

using SweepFunc = std::add_pointer_t<void(const gpu_kernel::Arguments&,
                                          CellView&,
                                          DirectionView&,
                                          const std::uint64_t*,
                                          const unsigned int&,
                                          const unsigned int&,
                                          const std::uint32_t&,
                                          double*,
                                          double*)>;
template <std::size_t... IntSequence>
__device__ constexpr std::array<SweepFunc, sizeof...(IntSequence)>
MakeSweepKernelSpecificationMap(std::index_sequence<IntSequence...>)
{
  return std::array<SweepFunc, sizeof...(IntSequence)>{&gpu_kernel::Sweep<IntSequence>...};
}
__device__ std::array<SweepFunc, max_dof_gpu> sweep_spec_map =
  MakeSweepKernelSpecificationMap(MakeIndexSequenceFromRange<1, max_dof_gpu + 1>{});

} // namespace opensn::gpu_kernel
