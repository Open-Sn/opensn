// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <optional>

namespace opensn
{

/// Returns this process's peak (high-water-mark) resident set size, in bytes.
std::optional<std::uint64_t> GetPeakMemoryUsageBytes();

/// Returns the total physical memory usage over all MPI ranks.
std::optional<std::uint64_t> GetTotalPeakMemoryUsageBytes();

} // namespace opensn
