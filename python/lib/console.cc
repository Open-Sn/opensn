// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/console.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include <iostream>
#include <functional>

using namespace opensn;
namespace py = pybind11;

namespace opensnpy
{

Console& console = Console::GetInstance();

Console&
Console::GetInstance() noexcept
{
  static Console singleton;
  return singleton;
}

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
Console::ExecutePythonCommand(const std::string& command) const
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
    if (opensn::mpi_comm.rank() == 0)
      std::cerr << e.what() << std::endl;
  }
}

void
Console::FlushConsole()
{
  assert(Py_IsInitialized());

  for (const auto& command : command_buffer_)
    ExecutePythonCommand(command);
}

void
Console::RunConsoleLoop() const
{
  assert(Py_IsInitialized());

  if (opensn::mpi_comm.rank() == 0)
    std::cout << "Console loop started. Type \"exit\" to quit.\n" << std::endl;

  while (true)
  {
    std::string command;

    opensn::mpi_comm.barrier();
    if (opensn::mpi_comm.rank() == 0)
      std::getline(std::cin, command);
    mpi_comm.broadcast(command, 0);
    if (command == "exit")
      break;

    ExecutePythonCommand(command);
  }

  if (opensn::mpi_comm.rank() == 0)
    std::cout << "Console loop stopped successfully." << std::endl;
  ;
}

int
Console::ExecuteFile(const std::string& fileName) const
{
  if (fileName.empty())
    return EXIT_FAILURE;

  try
  {
    py::eval_file(fileName);
    py::exec("import sys; sys.stdout.flush(); sys.stderr.flush();");
    std::cout << std::flush;
    std::cerr << std::flush;
  }
  catch (const py::error_already_set& e)
  {
    if (opensn::mpi_comm.rank() == 0)
      std::cerr << "Python Error while executing file: " << fileName << " - " << e.what()
                << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

} // namespace opensnpy
