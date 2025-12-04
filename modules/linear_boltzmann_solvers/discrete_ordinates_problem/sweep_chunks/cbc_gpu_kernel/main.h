// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_gpu_kernel/buffer.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_gpu_kernel/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn::cbc_gpu_kernel
{

/// Helper to generate an index sequence from From to To-1
template <std::size_t From, std::size_t To, std::size_t... IntSeq>
constexpr auto
MakeIndexSequenceFromRangeImpl(std::index_sequence<IntSeq...>)
{
  return std::index_sequence<(IntSeq + From)...>{};
}

/// Generate an index sequence from From to To-1
template <std::size_t From, std::size_t To>
using MakeIndexSequenceFromRange =
  decltype(MakeIndexSequenceFromRangeImpl<From, To>(std::make_index_sequence<To - From>{}));

/// Function pointer type for CBC sweep kernels
using SweepFunc = std::add_pointer_t<void(const Arguments&,
                                          CellView&,
                                          DirectionView&,
                                          const std::uint64_t*,
                                          const std::uint64_t&,
                                          const std::uint32_t&,
                                          const std::uint32_t&)>;

/// Create compile-time array of CBC sweep kernel function pointers for cell spatial DOF counts from
/// 1 to N
template <std::size_t... IntSequence>
__device__ constexpr std::array<SweepFunc, sizeof...(IntSequence)>
MakeCBCSweepSpecMap(std::index_sequence<IntSequence...>)
{
  return std::array{&Sweep<IntSequence>...};
}

/// Device map of CBC sweep kernel functions for cell spatial DOF counts from 1 to N
__device__ std::array<SweepFunc, LBSProblem::max_dofs_gpu> cbc_sweep_spec_map =
  MakeCBCSweepSpecMap(MakeIndexSequenceFromRange<1, LBSProblem::max_dofs_gpu + 1>{});

} // namespace opensn::cbc_gpu_kernel