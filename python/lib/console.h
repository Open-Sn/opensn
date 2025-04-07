// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/utils/utils.h"
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
  static Console& GetInstance() noexcept
  {
    static Console singleton;
    return singleton;
  }

  std::vector<std::string>& GetCommandBuffer() { return command_buffer_; }
  void InitConsole();
  void ExecuteFile(const std::string& input_filename) const;
  void BindBarrier(const mpi::Communicator& comm);

  static void BindModule(std::function<void(py::module&)> bind_function);

private:
  Console() noexcept = default;
  Console(const Console&) = delete;
  Console& operator=(const Console&) = delete;

  std::vector<std::string> command_buffer_;
};

extern Console& console;

} // namespace opensnpy
