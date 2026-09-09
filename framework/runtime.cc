// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/math/math.h"
#include "framework/utils/memory.h"
#include "framework/utils/timer.h"
#include "config.h"
#include "caliper/cali.h"
#include "hdf5.h"
#include <cstdlib>
#include <iostream>

namespace opensn
{
namespace
{

unsigned int
ResolveOpenSnNumThreads()
{
  // OpenSn reads this once during runtime initialization before any internal
  // worker threads are created.
  const char* env_value = std::getenv("OPENSN_NUM_THREADS"); // NOLINT(concurrency-mt-unsafe)
  if (env_value == nullptr)
    return 1;

  try
  {
    const auto parsed_value = std::stoul(env_value);
    return parsed_value > 0 ? static_cast<unsigned int>(parsed_value) : 1U;
  }
  catch (const std::exception&)
  {
    return 1;
  }
}

} // namespace

// Global variables
mpi::Communicator mpi_comm;
bool use_caliper = false;
std::string cali_config("runtime-report(calc.inclusive=true),max_column_width=80");
cali::ConfigManager cali_mgr;
Timer program_timer;
std::filesystem::path input_path;
unsigned int opensn_num_threads = 1;

namespace
{

// Check for the necessary services to support Caliper memory usage stats
std::string
SanitizeMemHighwatermarkOption(const cali::ConfigManager& mgr, std::string config)
{
  static const std::string option = "mem.highwatermark";

  const auto pos = config.find(option);
  if (pos == std::string::npos)
    return config; // not requested

  if (mgr.check(config.c_str()).empty())
    return config; // requested and supported on this build

  auto begin = pos;
  auto end = pos + option.size();
  if (begin > 0 and config[begin - 1] == ',')
    --begin;
  else if (end < config.size() and config[end] == ',')
    ++end;
  config.erase(begin, end - begin);

  log.Log() << "Note: \"mem.highwatermark\" was requested in the Caliper config, but is "
               "not available in this build.";

  return config;
}

} // namespace

int
Initialize()
{
  opensn_num_threads = ResolveOpenSnNumThreads();

  if (use_caliper)
  {
    cali_config = SanitizeMemHighwatermarkOption(cali_mgr, cali_config);
    cali_mgr.add(cali_config.c_str());
    cali_set_global_string_byname("opensn.version", GetVersionStr().c_str());
    cali_set_global_string_byname("opensn.input", input_path.c_str());
    cali_mgr.start();
  }

  CALI_MARK_PHASE_BEGIN(opensn::program.c_str());

  // Disable internal HDF error reporting
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

  return 0;
}

void
Finalize()
{
  // Flush standard streams
  std::cout.flush();
  std::cerr.flush();
  std::clog.flush();

  opensn::mpi_comm.barrier();

  // This must be set after CALI_MARK_PHASE_END to avoid spurious output in
  // the Caliper report
  CALI_MARK_PHASE_END(opensn::program.c_str());

  if (use_caliper)
  {
    if (const auto peak_bytes = GetPeakMemoryUsageBytes())
      cali_set_global_uint_byname("opensn.memory.highwatermark_bytes", *peak_bytes);
  }
}

std::string
GetVersionStr()
{
  return PROJECT_VERSION;
}

} // namespace opensn
