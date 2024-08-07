// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua_app.h"
#include "lua/modules/modules.h"
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;
using namespace opensn;

/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int
main(int argc, char** argv)
{
  mpi::Environment env(argc, argv);

  opensnlua::LuaApp app(MPI_COMM_WORLD);
  opensnlua::LoadRegisteredLuaItems();
  int error_code = app.Run(argc, argv);

  return error_code;
}
