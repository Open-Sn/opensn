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

  ~LuaApp()
  {}

  int Run(int argc, char** argv);

protected:
  int InitPetSc(int argc, char** argv);

  int ParseArguments(int argc, char** argv);

  int RunInteractive(int argc, char** argv);

  int RunBatch(int argc, char** argv);

private:
  bool sim_option_interactive_;
  bool allow_petsc_error_handler_;
};

} // namespace opensnlua
