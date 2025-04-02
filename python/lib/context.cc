// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "python/lib/py_env.h"

namespace opensn
{

// Wrap Finalize
void
WrapFinalize(py::module& context)
{
  // clang-format off
  context.def(
    "Finalize",
    [](void)
    {
      delete PyEnv::p_default_env;
      PyEnv::p_default_env = nullptr;
    },
    R"(
    Finalize OpenSn context.

    This function can only be invoked once at the end of the runtime process. In standard Python
    module mode, the context is automatically finalized when the interpreter ends the lifecycle of
    variables. However, in environments like IPython within Jupyter, functions registered with
    `atexit` may not be executed. In such cases, users must explicitly call this finalize function
    followed by MPI Finalize to properly terminate the context; otherwise, an MPI error will occur.
    )"
  );
  // clang-format on
}

// Wrap the context components of OpenSn
void
py_context(py::module& pyopensn)
{
  py::module context = pyopensn.def_submodule("context", "Context manager module.");
  WrapFinalize(context);
}

} // namespace opensn
