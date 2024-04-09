// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

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

  int Run(int argc, char** argv);

protected:
  int InitPetSc(int argc, char** argv);

  int ProcessArguments(int argc, char** argv);

  int RunInteractive(int argc, char** argv);

  int RunBatch(int argc, char** argv);

private:
  bool sim_option_interactive_;
  bool allow_petsc_error_handler_;
};

} // namespace opensnlua
