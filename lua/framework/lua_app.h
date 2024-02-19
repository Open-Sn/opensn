#pragma once

#include <string>
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;

namespace opensnlua
{

class LuaApp
{
public:
  LuaApp(const mpi::Communicator& comm);
  ~LuaApp();

  int Run(int argc, char** argv);

protected:
  /**Initializes PetSc for use by all entities.*/
  int InitPetSc(int argc, char** argv);

  /**
   * Parses input arguments.
   * \param argc int    Number of arguments supplied.
   * \param argv char** Array of strings representing each argument.
   */
  void ParseArguments(int argc, char** argv);

  int RunInteractive(int argc, char** argv);
  int RunBatch(int argc, char** argv);

private:
  // run_time quantities
  bool termination_posted_ = false;
  bool sim_option_interactive_ = true;
  bool allow_petsc_error_handler_ = false;
  bool supress_beg_end_timelog_ = false;
  bool dump_registry_ = false;
};

} // namespace opensnlua
