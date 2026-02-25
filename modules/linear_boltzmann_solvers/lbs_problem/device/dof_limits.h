// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include <cstdint>

namespace opensn
{

constexpr std::uint32_t max_dof_gpu_register = 4;

extern std::uint32_t max_dof_gpu_shared_mem;

std::uint32_t ComputeMaxDofSharedMem();

constexpr std::uint32_t max_dof_gpu = 10;

} // namespace opensn
