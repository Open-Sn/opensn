// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_env.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/utils/memory.h"
#include "framework/utils/timer.h"
#include "mpi4py/mpi4py.h"
#include "petscsys.h"

namespace opensn
{

PyEnv::PyEnv()
{
  // check if environment is already initialized
  if (PyEnv::p_default_env)
  {
    return;
  }
  // import mpi4py and get the communicator
  ::import_mpi4py();
  py::object mpi_module = py::module::import("mpi4py.MPI");
  py::object comm_world = mpi_module.attr("COMM_WORLD");
  mpi_comm = mpi::Communicator(*(::PyMPIComm_Get(comm_world.ptr())));
  if (mpi_comm.rank() == 0)
  {
    std::cout << program << " version " << GetVersionStr() << "\n";
    std::cout << Timer::GetLocalDateTimeString() << " Running " << program << " with "
              << mpi_comm.size() << " processes.\n\n";
  }
  // check MPI threading level
  int provided;
  MPI_Query_thread(&provided);
  if (provided != MPI_THREAD_MULTIPLE)
  {
    throw std::runtime_error("MPI must be initialized with thread_level \"multiple\". See "
                             "mpi4py.rc for more information\n");
  }
  // initialize PETSc and OpenSn
  ::PetscOptionsSetValue(NULL, "-options_left", "0");
  ::PetscOptionsInsertString(nullptr, "-no_signal_handler");
  ::PetscInitialize(nullptr, nullptr, nullptr, nullptr);
  log.SetColorEnabled(false);
  Initialize();
}

PyEnv::~PyEnv()
{
  // Finalize the run
  Finalize();
  ::PetscFinalize();
  // Flush Caliper report
  cali_mgr.flush();
  // Print memory usage
  const auto total_peak_bytes = GetTotalPeakMemoryUsageBytes();
  if (mpi_comm.rank() == 0)
  {
    std::cout << "\n";
    if (total_peak_bytes)
      std::cout << "Total peak host memory usage: "
                << static_cast<double>(*total_peak_bytes) / 1.0e9 << " GB\n";
  }
  // Print execution time
  if (mpi_comm.rank() == 0)
  {
    std::cout << "Elapsed execution time: " << program_timer.GetTimeString() << "\n";
    std::cout << Timer::GetLocalDateTimeString() << " " << program << " finished execution.\n";
  }
}

std::unique_ptr<PyEnv> PyEnv::p_default_env{};

} // namespace opensn
