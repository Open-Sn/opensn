// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_batch_slot.h"
#include "caribou/main.hpp"
#include <array>
#include <cstdint>

namespace crb = caribou;

namespace opensn::gpu_kernel
{

#if defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
/** Borrowed SYCL submission state, retained by its owning lane until completion. */
struct CBCD_LaunchStream
{
  crb::Stream& queue;
  const CBCD_Launch* metadata;
  std::uint32_t grid_blocks;
};
#else
using CBCD_LaunchStream = crb::Stream;
#endif

using CBCD_LaunchFunction = void (*)(const CBCD_Launch& launch,
                                     const Arguments<SweepKind::CBC>* arguments,
                                     double* const* saved_psi,
                                     CBCD_LaunchStream& stream);

/// Return fixed-DOF launch functions by size class, all null if the count is unsupported.
std::array<CBCD_LaunchFunction, CBCD_NUM_LAUNCH_CLASSES>
CBCDFixedSweepLaunches(std::uint32_t ndofs);

} // namespace opensn::gpu_kernel
