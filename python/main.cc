// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_app.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include "slepc.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace mpi = mpicpp_lite;

namespace
{

std::vector<char*>
BuildSLEPcArgv(int argc, char** argv)
{
  std::vector<char*> slepc_argv;
  slepc_argv.reserve(static_cast<size_t>(argc));
  slepc_argv.push_back(argv[0]);

  for (int k = 1; k < argc; ++k)
  {
    const std::string arg = argv[k];

    // These options are consumed by the OpenSn Python app, not PETSc/SLEPc.
    if (arg == "-h" or arg == "--help" or arg == "-c" or arg == "--suppress-color")
      continue;

    if (arg == "-v" or arg == "--verbose" or arg == "-i" or arg == "--input" or arg == "-p" or
        arg == "--py")
    {
      ++k; // Skip the accompanying value.
      continue;
    }

    if (arg.starts_with("--verbose=") or arg.starts_with("--input=") or arg.starts_with("--py=") or
        arg.starts_with("--caliper="))
      continue;

    if (arg == "--caliper")
    {
      // cxxopts uses an implicit value when no explicit value is provided.
      if (k + 1 < argc and argv[k + 1][0] != '-')
        ++k;
      continue;
    }

    slepc_argv.push_back(argv[k]);
  }

  return slepc_argv;
}

} // namespace

int
main(int argc, char** argv) // NOLINT(bugprone-exception-escape)
{
  mpi::Environment env(argc, argv, mpi::ThreadSupport::MULTIPLE);

  auto slepc_argv = BuildSLEPcArgv(argc, argv);
  int slepc_argc = static_cast<int>(slepc_argv.size());
  char** slepc_argv_ptr = slepc_argv.data();
  PetscCall(SlepcInitialize(&slepc_argc, &slepc_argv_ptr, nullptr, nullptr));

  int retval = EXIT_SUCCESS;
  try
  {
    py::scoped_interpreter guard{};
    opensnpy::PyApp app(MPI_COMM_WORLD); // NOLINT(bugprone-casting-through-void)
    retval = app.Run(argc, argv);
  }
  catch (...)
  {
    std::cerr << "Unknown fatal error\n";
    retval = EXIT_FAILURE;
  }

  SlepcFinalize();

  return retval;
}
