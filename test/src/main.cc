#include "mpicpp-lite/mpicpp-lite.h"
#include "lua/modules/modules.h"
#include "lua/framework/lua_app.h"

namespace mpi = mpicpp_lite;
using namespace opensn;

int
main(int argc, char** argv)
{
  mpi::Environment env(argc, argv);

  opensnlua::LuaApp app(MPI_COMM_WORLD);
  opensnlua::LoadRegisteredLuaItems();
  int error_code = app.Run(argc, argv);

  return error_code;
}
