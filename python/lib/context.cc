// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "python/lib/py_env.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include "petscsys.h"
#include <pybind11/eval.h>
#include <stdexcept>

namespace opensn
{

// Wrap Finalize
static void
WrapFinalize(py::module& context)
{
  // clang-format off
  // finalize
  context.def(
    "Finalize",
    [](void)
    {
      delete PyEnv::p_default_env;
      PyEnv::p_default_env = nullptr;
    },
    R"(
    Finalize OpenSn context.

    This function can only be invoked **once at the end of the runtime process**.

    In standard Python module mode, the context is automatically finalized when the interpreter ends
    the lifecycle of variables. However, in environments like IPython within Jupyter, functions
    registered with `atexit` may not be executed.

    In such cases, users must explicitly call this finalize function followed by MPI Finalize to
    properly terminate the context; otherwise, an MPI error will occur. This process can be
    registered to the ``post_execute`` event of IPython console as follows:

    .. code-block::

       from IPython import get_ipython

       def finalize_env():
           Finalize()
           MPI.Finalize()

       ipython_instance = get_ipython()
       if ipython_instance is not None:
           ipython_instance.events.register("post_execute", finalize_env)

    This approach allows the module to safely execute the Python script generated by ``nbconvert``.
    )"
  );
  // clang-format on
}

// Wrap settings
static void
WrapSettings(py::module& context)
{
  // clang-format off
  // log settings
  context.def(
    "SetVerbosityLevel",
    [](int level)
    {
      log.SetVerbosity(level);
    },
    "Set verbosity level (0 to 3). Default is 0.",
    py::arg("level")
  );
  context.def(
    "UseColor",
    [](bool cfg)
    {
      suppress_color = (!cfg);
    },
    "Enable/disable color output. Default is True.",
    py::arg("config") = true
  );
  // PETSc error handler
  context.def(
    "EnablePETScErrorHandler",
    []()
    {
      ::PetscOptionsSetValue(nullptr, "-no_signal_handler", "false");
    },
    "Allow PETSc error handler."
  );
  // Caliper reporting
  context.def(
    "SetCaliperConfig",
    [](const std::string& config)
    {
      if (use_caliper)
      {
        throw std::runtime_error("This function can only be called before enabling Cailper.");
      }
      cali_config = config;
    },
    R"(
    Set configuration to the Caliper manager.

    This function can only be called before using :py:func:`pyopensn.context.EnableCaliper`. Note
    that this function does not start the Caliper manager immediately.

    Parameters
    ----------
    config: str, default='runtime-report(calc.inclusive=true),max_column_width=80'
        Configuration.
    )",
    py::arg("config") = "runtime-report(calc.inclusive=true),max_column_width=80"
  );
  context.def(
    "EnableCaliper",
    []()
    {
      // check if caliper is already initialized
      if (use_caliper)
      {
        throw std::runtime_error("Caliper is already set.");
      }
      use_caliper = true;
      // initialize Caliper
      cali_mgr.add(cali_config.c_str());
      cali_set_global_string_byname("opensn.version", GetVersionStr().c_str());
      cali_set_global_string_byname("opensn.input", input_path.c_str());
      cali_mgr.start();
    },
    R"(
    Start the Caliper manager and mark the program begin.
    )"
  );
  // clang-format on
}

// Wrap sys.argv translator
static void
WrapSysArgv(py::module& context)
{
  // clang-format off
  context.def(
    "InitializeWithArgv",
    [](py::list sys_argv)
    {
      // looping over each argument of the list
      for (int i_arg = 0; i_arg < sys_argv.size(); ++i_arg)
      {
        // skip for the first argument (Python script name)
        if (i_arg == 0)
        {
          continue;
        }
        // cast argument to string
        std::string arg = sys_argv[i_arg].cast<std::string>();
        // color
        if (arg == "-c" || arg == "--suppress-color")
        {
          suppress_color = true;
          continue;
        }
        // verbosity
        if (arg == "-v" || arg == "--verbose")
        {
          std::string next_arg = sys_argv[++i_arg].cast<std::string>();
          log.SetVerbosity(std::atoi(next_arg.c_str()));
          continue;
        }
        // caliper
        if (arg == "--caliper")
        {
          std::string next_arg = sys_argv[++i_arg].cast<std::string>();
          use_caliper = true;
          cali_config = next_arg;
          cali_mgr.add(cali_config.c_str());
          cali_set_global_string_byname("opensn.version", GetVersionStr().c_str());
          cali_set_global_string_byname("opensn.input", input_path.c_str());
          cali_mgr.start();
          continue;
        }
        // PETSc handler
        if (arg == "--allow-petsc-error-handler")
        {
          ::PetscOptionsSetValue(nullptr, "-no_signal_handler", "false");
          continue;
        }
        // execute Python statement
        if (arg == "--py")
        {
          std::string next_arg = sys_argv[++i_arg].cast<std::string>();
          try
          {
            py::exec(next_arg);
          }
          catch (py::error_already_set& e)
          {
            py::object exc_type = e.type();
            std::string type_name = py::str(exc_type.attr("__name__"));
            log.LogAllError() << "Caught an exception of type: [" << type_name
              << "] when executing code of the sys.argv.\nException message: " << e.what() << "\n";
            e.restore();
          }
          continue;
        }
        // throw error for invalid argument
        throw std::runtime_error("Invalid argument: " + arg + "\n");
      }
    },
    R"(
    Overwrite OpenSn settings using ``sys.argv``.

    This is a temporary solution for dealing with command line argument mode of the console. Once
    the console is gone, this functionality will also be killed.

    Parameters
    ----------
    sys_argv: List[str], default=sys.argv
        Argument vector to be used. Default to ``sys.argv``.
    )",
    py::arg_v("sys_argv", py::module::import("sys").attr("argv").cast<py::list>(), "sys.argv")
  );
  // clang-format on
}

// Wrap the context components of OpenSn
void
py_context(py::module& pyopensn)
{
  py::module context = pyopensn.def_submodule("context", "Context manager module.");
  WrapFinalize(context);
  WrapSettings(context);
  WrapSysArgv(context);
}

} // namespace opensn
