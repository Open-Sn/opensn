// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_fused_main.h"
#include <array>
#include <cstdint>
#include <utility>

namespace opensn::gpu_kernel
{

namespace
{

inline constexpr std::array<std::uint32_t, 5> CBCD_FIXED_CELL_DOFS{2, 3, 4, 6, 8};

template <std::size_t... launch_class>
std::array<CBCD_LaunchFunction, sizeof...(launch_class)>
FixedSweepLaunches(std::uint32_t ndofs, std::index_sequence<launch_class...> /*indices*/)
{
  switch (ndofs)
  {
    case CBCD_FIXED_CELL_DOFS[0]:
      return {&LaunchCBCDFusedSweepAs<launch_class, CBCD_FIXED_CELL_DOFS[0]>...};
    case CBCD_FIXED_CELL_DOFS[1]:
      return {&LaunchCBCDFusedSweepAs<launch_class, CBCD_FIXED_CELL_DOFS[1]>...};
    case CBCD_FIXED_CELL_DOFS[2]:
      return {&LaunchCBCDFusedSweepAs<launch_class, CBCD_FIXED_CELL_DOFS[2]>...};
    case CBCD_FIXED_CELL_DOFS[3]:
      return {&LaunchCBCDFusedSweepAs<launch_class, CBCD_FIXED_CELL_DOFS[3]>...};
    case CBCD_FIXED_CELL_DOFS[4]:
      return {&LaunchCBCDFusedSweepAs<launch_class, CBCD_FIXED_CELL_DOFS[4]>...};
    default:
      return {};
  }
}

} // namespace

std::array<CBCD_LaunchFunction, CBCD_NUM_LAUNCH_CLASSES>
CBCDFixedSweepLaunches(std::uint32_t ndofs)
{
  return FixedSweepLaunches(ndofs, std::make_index_sequence<CBCD_NUM_LAUNCH_CLASSES>{});
}

} // namespace opensn::gpu_kernel
