// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_app.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <cstdio>
#include <cstdlib>

namespace mpi = mpicpp_lite;

int
main(int argc, char** argv)
{
  try
  {
    mpi::Environment env(argc, argv);
    py::scoped_interpreter guard{};
    opensnpy::PyApp app(MPI_COMM_WORLD);
    return app.Run(argc, argv);
  }
  catch (...)
  {
    std::fprintf(stderr, "Unknown fatal error\n");
    return EXIT_FAILURE;
  }
}
