// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "mpicpp-lite/mpicpp-lite.h"
#include <vector>
#include <string>
#include <Python.h>
#include <pybind11/embed.h>

namespace mpi = mpicpp_lite;
namespace py = pybind11;

namespace opensnpy
{

class Console
{
public:
  static Console& GetInstance() noexcept;
  std::vector<std::string>& GetCommandBuffer() { return command_buffer_; }
  void FlushConsole();
  void RunConsoleLoop() const;
  int ExecuteFile(const std::string& fileName) const;
  void BindModule(std::function<void(py::module&)> bind_function);
  void BindBarrier(const mpi::Communicator& comm);

private:
  Console() noexcept = default;
  void ExecutePythonCommand(const std::string& command) const;

  std::vector<std::string> command_buffer_;
};

extern Console& console;

} // namespace opensnpy
