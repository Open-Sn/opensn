// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/runtime.h"
#include "framework/math/math.h"
#include "framework/object_factory.h"
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

int
Initialize()
{
  opensn_num_threads = ResolveOpenSnNumThreads();

  if (use_caliper)
  {
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

  CALI_MARK_PHASE_END(opensn::program.c_str());
}

std::string
GetVersionStr()
{
  return PROJECT_VERSION;
}

} // namespace opensn
