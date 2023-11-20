#include "modules/modules_lua.h"
#include "lua/framework/lua_app.h"

using namespace opensn;

/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
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
