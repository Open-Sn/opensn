// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <memory>

namespace py = pybind11;

namespace opensn
{

/// Environment initializer and finalizer for OpenSn.
class PyEnv
{
public:
  PyEnv();

  ~PyEnv();

  /// Default environement.
  static std::unique_ptr<PyEnv> p_default_env;
};

} // namespace opensn
