// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
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
    [comm]()
    {
      comm.barrier();
    },
    "MPI barrier for console."
  );
  // clang-format on
}

void
Console::BindAllReduce(const mpi::Communicator& comm)
{
  // clang-format off
  py::module main = py::module::import("__main__");
  main.def(
    "MPIAllReduce",
    [comm](double value)
    {
      double out = 0.0;
      comm.all_reduce(value, out, mpi::op::sum<double>());
      return out;
    },
    "MPI all-reduce sum for a scalar."
  );
  main.def(
    "MPIAllReduce",
    [comm](const std::vector<double>& values)
    {
      std::vector<double> out(values.size(), 0.0);
      if (not values.empty())
      {
        const auto count = static_cast<int>(values.size());
        comm.all_reduce(values.data(), count, out.data(), mpi::op::sum<double>());
      }
      return out;
    },
    "MPI all-reduce sum for a list of doubles."
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
