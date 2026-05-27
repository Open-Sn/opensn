// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <string>

namespace opensn
{

/// Hardware and launch resources visible to the current process.
struct ThreadResourceInfo
{
  /// Hardware-thread count reported by the standard library.
  std::size_t hardware_threads = 1;
  /// Explicit OpenSn or OpenMP thread request, or zero when unset.
  std::size_t requested_threads = 0;
  /// Environment variable supplying `requested_threads`.
  std::string requested_source;
  /// Threads available to this process after launch constraints.
  std::size_t available_threads = 1;
  /// Allocated job nodes and the environment variable supplying the value.
  std::size_t job_nodes = 0;
  std::string job_nodes_source;
  /// Allocated job tasks and the environment variable supplying the value.
  std::size_t job_tasks = 0;
  std::string job_tasks_source;
};

/// Return the hardware and launch resources available to the current process.
ThreadResourceInfo GetThreadResourceInfo();

/// Format CPU-thread resource metadata.
std::string FormatThreadResourceInfo(const ThreadResourceInfo& info);

} // namespace opensn
