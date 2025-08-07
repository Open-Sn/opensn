// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "python/lib/py_env.h"
#include "framework/runtime.h"
#include "config.h"
#include "petsc.h"

PYBIND11_MODULE(pyopensn, pyopensn)
{
  using namespace opensn;

  // Metadata
  pyopensn.doc() = "Python interface for OpenSn.";
  pyopensn.attr("__version__") = PROJECT_VERSION;
  // Environment
  if (PyEnv::p_default_env == nullptr)
  {
    PyEnv::p_default_env = new PyEnv();
    ::Py_AtExit(
      []() noexcept
      {
        delete PyEnv::p_default_env;
        PyEnv::p_default_env = nullptr;
      });
  }

  // Wrap libraries
  py_context(pyopensn);
  py_math(pyopensn);
  py_aquad(pyopensn);
  py_mesh(pyopensn);
  py_logvol(pyopensn);
  py_source(pyopensn);
  py_xs(pyopensn);
  py_ffunc(pyopensn);
  py_response(pyopensn);
  py_solver(pyopensn);
}
