// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include "petsc.h"
#include "framework/runtime.h"

#include "config.h"

PYBIND11_MODULE(pyopensn, pyopensn) {
    using namespace opensn;
    // metadata
    pyopensn.doc() = "Python interface for OpenSn.";
    pyopensn.attr("__version__") = PROJECT_VERSION;
    // environment
    if (PyEnv::p_default_env == nullptr) {
        PyEnv::p_default_env = new PyEnv();
    }
    ::Py_AtExit(
        [](void) noexcept {
            delete PyEnv::p_default_env;
            PyEnv::p_default_env = nullptr;
        }
    );
    // wrap libraries
    py_aquad(pyopensn);
    py_ffunc(pyopensn);
    py_logvol(pyopensn);
    py_mat(pyopensn);
    py_mesh(pyopensn);
    py_solver(pyopensn);
    py_source(pyopensn);
    py_xs(pyopensn);
}
