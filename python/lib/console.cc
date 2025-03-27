// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/console.h"
#include "framework/utils/utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iostream>
#include <functional>
#include <regex>

using namespace opensn;
namespace py = pybind11;

namespace opensnpy
{

Console& console = Console::GetInstance();

void
Console::BindModule(std::function<void(py::module&)> bind_function)
{
  py::module main = py::module::import("__main__");
  bind_function(main);
}

void
Console::BindBarrier(const mpi::Communicator& comm)
{
  // clang-format off
  py::module main = py::module::import("__main__");
  main.def(
    "MPIBarrier",
    [&comm](void)
    {
      comm.barrier();
    },
    "MPI barrier for console."
  );
  // clang-format on
}

void
Console::InitConsole()
{
  assert(Py_IsInitialized());

  for (const auto& command : command_buffer_)
  {
    try
    {
      py::exec(command);
      py::exec("import sys; sys.stdout.flush(); sys.stderr.flush();");
      std::cout << std::flush;
      std::cerr << std::flush;
    }
    catch (const py::error_already_set& e)
    {
      opensn::log.LogAllError() << e.what();
      opensn::mpi_comm.abort(EXIT_FAILURE);
    }
  }
}

void
Console::ExecuteFile(const std::string& input_filename) const
{
  try
  {
    py::eval_file(input_filename);
    py::exec("import sys; sys.stdout.flush(); sys.stderr.flush();");
    std::cout << std::flush;
    std::cerr << std::flush;
  }
  catch (const py::error_already_set& e)
  {
    opensn::log.LogAllError() << e.what();
    opensn::mpi_comm.abort(EXIT_FAILURE);
  }
}

} // namespace opensnpy
