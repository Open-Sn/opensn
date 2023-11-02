#include "lbts_transient_solver.h"

#include "chi_runtime.h"
#include "console/chi_console.h"

/**Provides a callback interface to lua after each timestep. Users,
 * can setup all sorts of controls like adaptive timestepping and
 * outputs.*/
void lbs::DiscOrdTransientSolver::PostStepCallBackFunction() const
{
  const std::string fname = "lbs::TransientSolver::PostStepCallBackFunction";

  if (transient_options_.console_call_back_function.empty()) return;

  auto& L = chi::console.GetConsoleState();
  const auto& lua_func_name = transient_options_.console_call_back_function;

  //============= Load lua function
  lua_getglobal(L, lua_func_name.c_str());

  //============= Error check lua function
  if (not lua_isfunction(L, -1))
    throw std::logic_error(fname + " attempted to access lua-function, " +
                           lua_func_name + ", but it seems the function"
                                           " could not be retrieved.");

  //============= Push arguments
  //There are no arguments

  //============= Call lua function
  //0 arguments, 0 result (double), 0=original error object
  double lua_return;
  if (lua_pcall(L,0,0,0) == 0)
  {
  }
  else
    throw std::logic_error(fname + " attempted to call lua-function, " +
                           lua_func_name + ", but the call failed.");

//  lua_pop(L,1); //pop the double, or error code
}