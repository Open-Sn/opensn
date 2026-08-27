// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/utils/memory.h"
#include "framework/runtime.h"
#include <sys/resource.h>

namespace opensn
{

std::optional<std::uint64_t>
GetPeakMemoryUsageBytes()
{
  struct rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    return std::nullopt;
  // ru_maxrss is glibc's normal POSIX-mandated way to read this field, but some libc
  // implementations define it via an anonymous union which can trip clang-tidy
  const long max_rss = usage.ru_maxrss; // NOLINT(cppcoreguidelines-pro-type-union-access)
  // A successful call with max_rss == 0 means the platform doesn't actually fill in this
  // field rather than genuine zero peak usage
  if (max_rss <= 0)
    return std::nullopt;
  // ru_maxrss is in bytes on macOS/Darwin, but kilobytes on Linux
#ifdef __APPLE__
  return static_cast<std::uint64_t>(max_rss);
#else
  return static_cast<std::uint64_t>(max_rss) * 1024;
#endif
}

std::optional<std::uint64_t>
GetTotalPeakMemoryUsageBytes()
{
  const auto local_peak = GetPeakMemoryUsageBytes();
  const std::uint64_t local_peak_bytes = local_peak.value_or(0);
  const int local_available = local_peak.has_value() ? 1 : 0;

  std::uint64_t total_peak_bytes = 0;
  int all_available = 0;
  mpi_comm.all_reduce(local_peak_bytes, total_peak_bytes, mpi::op::sum<std::uint64_t>());
  mpi_comm.all_reduce(local_available, all_available, mpi::op::min<int>());

  if (!all_available)
    return std::nullopt;
  return total_peak_bytes;
}

} // namespace opensn
