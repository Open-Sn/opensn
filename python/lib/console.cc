// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/console.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include <iostream>

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
Console::BindModule(const std::string& module_name, void (*wrap_function)(py::module&))
{
  py::module main = py::module::import("__main__");
  py::module mod = main.def_submodule(module_name.c_str(), ("Submodule: " + module_name).c_str());
  wrap_function(mod);
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
  assert(PyIsInitialized());

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
