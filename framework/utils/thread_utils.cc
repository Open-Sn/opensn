// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/utils/thread_utils.h"
#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <initializer_list>
#include <sstream>
#include <string_view>
#include <thread>
#include <utility>

namespace opensn
{

namespace
{

std::size_t
ParsePositiveEnvironmentValue(const char* name)
{
  const char* value = std::getenv(name); // NOLINT(concurrency-mt-unsafe)
  if (value == nullptr or value[0] == '\0')
    return 0;

  const std::string_view text(value);
  std::size_t parsed = 0;
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), parsed);
  if (error != std::errc{} or end != text.data() + text.size() or parsed == 0)
    return 0;

  return parsed;
}

std::pair<std::size_t, std::string>
GetEnvironmentCount(const std::initializer_list<const char*> names)
{
  for (const char* name : names)
    if (const auto count = ParsePositiveEnvironmentValue(name); count > 0)
      return {count, name};

  return {0, {}};
}

} // namespace

ThreadResourceInfo
GetThreadResourceInfo()
{
  ThreadResourceInfo info;
  info.hardware_threads = std::max<std::size_t>(1, std::thread::hardware_concurrency());

  const auto [requested_threads, requested_source] = GetEnvironmentCount({"OPENSN_NUM_THREADS",
                                                                          "SLURM_CPUS_PER_TASK",
                                                                          "OMP_NUM_THREADS",
                                                                          "FLUX_CPUS_PER_TASK",
                                                                          "FLUX_TASK_CPUS"});
  info.requested_threads = requested_threads;
  info.requested_source = requested_source;

  const auto [job_nodes, job_nodes_source] =
    GetEnvironmentCount({"FLUX_JOB_NNODES", "SLURM_JOB_NUM_NODES", "SLURM_NNODES", "PBS_NNODES"});
  info.job_nodes = job_nodes;
  info.job_nodes_source = job_nodes_source;

  const auto [job_tasks, job_tasks_source] =
    GetEnvironmentCount({"FLUX_JOB_SIZE", "SLURM_NTASKS", "SLURM_NPROCS", "PBS_NP"});
  info.job_tasks = job_tasks;
  info.job_tasks_source = job_tasks_source;

  info.available_threads = info.hardware_threads;
  if (info.requested_threads > 0)
    info.available_threads = std::min(info.available_threads, info.requested_threads);
  info.available_threads = std::max<std::size_t>(1, info.available_threads);

  return info;
}

std::string
FormatThreadResourceInfo(const ThreadResourceInfo& info)
{
  std::ostringstream out;
  out << "hardware_threads=" << info.hardware_threads
      << ", requested_threads=" << info.requested_threads;
  if (not info.requested_source.empty())
    out << " from " << info.requested_source;
  out << ", available_threads=" << info.available_threads;

  if (info.job_nodes > 0)
    out << ", job_nodes=" << info.job_nodes << " from " << info.job_nodes_source;
  if (info.job_tasks > 0)
    out << ", job_tasks=" << info.job_tasks << " from " << info.job_tasks_source;

  return out.str();
}

} // namespace opensn
