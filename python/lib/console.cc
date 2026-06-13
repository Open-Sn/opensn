// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/console.h"
#include "framework/utils/utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iostream>
#include <functional>
#include <regex>
#include <stdexcept>

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

  // Scalar double — op: "sum" (default), "max", "min"
  main.def(
    "MPIAllReduce",
    [comm](double value, const std::string& op) -> double
    {
      double out = 0.0;
      if (op == "sum")      comm.all_reduce(value, out, mpi::op::sum<double>());
      else if (op == "max") comm.all_reduce(value, out, mpi::op::max<double>());
      else if (op == "min") comm.all_reduce(value, out, mpi::op::min<double>());
      else throw std::invalid_argument("MPIAllReduce: unknown op '" + op + "'");
      return out;
    },
    py::arg("value"), py::arg("op") = "sum",
    "MPI all-reduce for a scalar double."
  );

  // Vector of doubles — op: "sum" (default), "max", "min"
  main.def(
    "MPIAllReduce",
    [comm](const std::vector<double>& values, const std::string& op) -> std::vector<double>
    {
      std::vector<double> out(values.size(), 0.0);
      if (not values.empty())
      {
        const auto n = static_cast<int>(values.size());
        if (op == "sum")      comm.all_reduce(values.data(), n, out.data(), mpi::op::sum<double>());
        else if (op == "max") comm.all_reduce(values.data(), n, out.data(), mpi::op::max<double>());
        else if (op == "min") comm.all_reduce(values.data(), n, out.data(), mpi::op::min<double>());
        else throw std::invalid_argument("MPIAllReduce: unknown op '" + op + "'");
      }
      return out;
    },
    py::arg("values"), py::arg("op") = "sum",
    "MPI all-reduce for a list of doubles."
  );

  // Scalar int — op: "sum" (default), "max", "min", "bor"
  main.def(
    "MPIAllReduce",
    [comm](int value, const std::string& op) -> int
    {
      int out = 0;
      if (op == "sum")       comm.all_reduce(value, out, mpi::op::sum<int>());
      else if (op == "max")  comm.all_reduce(value, out, mpi::op::max<int>());
      else if (op == "min")  comm.all_reduce(value, out, mpi::op::min<int>());
      else if (op == "bor")  MPI_Allreduce(&value, &out, 1, MPI_INT, MPI_BOR, comm);
      else throw std::invalid_argument("MPIAllReduce: unknown op '" + op + "'");
      return out;
    },
    py::arg("value"), py::arg("op") = "sum",
    "MPI all-reduce for a scalar integer."
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
