#include "mpicpp-lite/mpicpp-lite.h"
#include "lua/lib/lua_app.h"

namespace mpi = mpicpp_lite;

int
main(int argc, char** argv)
{
  mpi::Environment env(argc, argv);

  opensnlua::LuaApp app(MPI_COMM_WORLD);
  int error_code = app.Run(argc, argv);

  return error_code;
}
