#include "lua/modules/modules_lua.h"
#include "lua/framework/lua_app.h"

using namespace opensn;

int
main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  opensnlua::LuaApp app(MPI_COMM_WORLD);
  opensnlua::LoadRegisteredLuaItems();
  int error_code = app.Run(argc, argv);

  MPI_Finalize();
  return error_code;
}
