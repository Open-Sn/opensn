// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_app.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <cstdio>
#include <cstdlib>
#include "petsc.h"

namespace mpi = mpicpp_lite;

int
main(int argc, char** argv)
{
  mpi::Environment env(argc, argv);

  PetscCall(PetscInitializeNoArguments());

  int retval = EXIT_SUCCESS;
  try
  {
    py::scoped_interpreter guard{};
    opensnpy::PyApp app(MPI_COMM_WORLD);
    retval = app.Run(argc, argv);
  }
  catch (...)
  {
    std::fprintf(stderr, "Unknown fatal error\n");
    retval = EXIT_FAILURE;
  }

  PetscFinalize();

  return retval;
}
