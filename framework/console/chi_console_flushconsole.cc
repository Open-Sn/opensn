#include "console/chi_console.h"
#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/** This function sends the commands contained in the command_buffer to
the lua state from where it is executed. These could be commands passed
 via the command line or loaded elsewhere.*/
void
chi::Console::FlushConsole()
{
  try
  {
    for (auto& command : command_buffer_)
    {
      bool error = luaL_dostring(console_state_, command.c_str());
      if (error)
      {
        Chi::log.LogAll() << lua_tostring(console_state_, -1);
        lua_pop(console_state_, 1);
      }
    }
  }
  catch (const std::exception& e)
  {
    Chi::log.LogAllError() << e.what();
    Chi::Exit(EXIT_FAILURE);
  }
}
