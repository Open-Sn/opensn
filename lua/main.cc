// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/lua_app.h"
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;

int
main(int argc, char** argv)
{
  mpi::Environment env(argc, argv);
  opensnlua::LuaApp app(MPI_COMM_WORLD);
  int error_code = app.Run(argc, argv);
  return error_code;
}
