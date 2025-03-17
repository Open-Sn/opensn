// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include "petscsys.h"
#include <memory>
#include <stdexcept>

namespace opensn
{

// Wrap settings
static void
WrapSettings(py::module& settings)
{
  // clang-format off
  // log settings
  settings.def(
    "SetVerbosityLevel",
    [](int level)
    {
      log.SetVerbosity(level);
    },
    "Set verbosity level (0 to 3). Default is 0.",
    py::arg("level")
  );
  settings.def(
    "UseColor",
    [](bool cfg)
    {
      suppress_color = (!cfg);
    },
    "Enable/disable color output. Default is True.",
    py::arg("config") = true
  );
  // PETSc error handler
  settings.def(
    "EnablePETScErrorHandler",
    []()
    {
      ::PetscOptionsInsertString(nullptr, "-no_signal_handler");
    },
    "Allow PETSc error handler."
  );
  // Caliper reporting
  settings.def(
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

    This function can only be called before using :py:func:`pyopensn.settings.EnableCaliper`. Note
    that this function does not start the Caliper manager immediately.

    Parameters
    ----------
    config: str, default='runtime-report(calc.inclusive=true),max_column_width=80'
        Configuration.
    )",
    py::arg("config") = "runtime-report(calc.inclusive=true),max_column_width=80"
  );
  settings.def(
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
      // CALI_MARK_BEGIN(opensn::program.c_str());
    },
    R"(
    Start the Caliper manager and mark the program begin.
    )"
  );
  // clang-format on
}

// Wrap the settings components of OpenSn
void
py_settings(py::module& pyopensn)
{
  py::module settings = pyopensn.def_submodule("settings", "Settings module.");
  WrapSettings(settings);
}

} // namespace opensn
