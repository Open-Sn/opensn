// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_app.h"
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;

int
main(int argc, char** argv)
{
  mpi::Environment env(argc, argv);
  py::scoped_interpreter guard{};
  opensnpy::PyApp app(MPI_COMM_WORLD);
  int error_code = app.Run(argc, argv);
  return error_code;
}
