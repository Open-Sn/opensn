// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_app.h"
#include "python/lib/py_api.h"
#include "python/lib/console.h"
#include "python/lib/py_wrappers.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include "cxxopts/cxxopts.h"
#include "petsc.h"
#include <string>

using namespace opensn;
namespace py = pybind11;

namespace opensnpy
{

PyApp::PyApp(const mpi::Communicator& comm)
  : sim_option_interactive_(true), allow_petsc_error_handler_(false)
{
  opensn::mpi_comm = comm;

  py::module sys = py::module::import("sys");
  py::exec("rank = " + std::to_string(comm.rank()));
  py::exec("size = " + std::to_string(comm.size()));
  py::exec("opensn_console = True");

  console.BindBarrier(comm);

  console.BindModule(WrapYlm);
  console.BindModule(WrapVector3);
  console.BindModule(WrapFunctors);

  console.BindModule(WrapQuadrature);
  console.BindModule(WrapProductQuadrature);
  console.BindModule(WrapCurvilinearQuadrature);
  console.BindModule(WrapSLDFESQuadrature);

  console.BindModule(WrapMesh);
  console.BindModule(WrapMeshGenerator);
  console.BindModule(WrapGraphPartitioner);

  console.BindModule(WrapLogicalVolume);

  console.BindModule(WrapPointSource);
  console.BindModule(WrapVolumetricSource);

  console.BindModule(WrapMultiGroupXS);

  console.BindModule(WrapFieldFunction);
  console.BindModule(WrapFieldFunctionGridBased);
  console.BindModule(WrapFieldFunctionInterpolation);

  console.BindModule(WrapResEval);

  console.BindModule(WrapSolver);
  console.BindModule(WrapLBS);
  console.BindModule(WrapSteadyState);
  console.BindModule(WrapNLKEigen);
  console.BindModule(WrapPIteration);
  console.BindModule(WrapPRK);

  console.BindModule(WrapDiffusion);

  console.BindModule(WrapPostProcessor);
  console.BindModule(WrapPrinter);

}

int
PyApp::InitPetSc(int argc, char** argv)
{
  PetscOptionsInsertString(nullptr, "-error_output_stderr");
  if (!allow_petsc_error_handler_)
    PetscOptionsInsertString(nullptr, "-no_signal_handler");
  PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));
  return 0;
}

int
PyApp::Run(int argc, char** argv)
{
  if (opensn::mpi_comm.rank() == 0)
  {
    std::cout << opensn::program << " version " << GetVersionStr() << "\n"
              << Timer::GetLocalDateTimeString() << " Running " << opensn::program << " with "
              << opensn::mpi_comm.size() << " processes.\n"
              << opensn::program << " number of arguments supplied: " << argc - 1 << "\n"
              << std::endl;
  }

  int error_code = ProcessArguments(argc, argv);
  if (!error_code)
  {
    PetscOptionsSetValue(NULL, "-options_left", "0");
    InitPetSc(argc, argv);
    opensn::Initialize();
    console.FlushConsole();

    if (sim_option_interactive_)
      error_code = RunInteractive();
    else
      error_code = RunBatch();

    opensn::Finalize();
    PetscFinalize();

    if (opensn::mpi_comm.rank() == 0)
    {
      std::cout << "\nElapsed execution time: " << program_timer.GetTimeString() << "\n"
                << Timer::GetLocalDateTimeString() << " " << opensn::program
                << " finished execution." << std::endl;
    }
  }
  cali_mgr.flush();
  return error_code;
}

int
PyApp::ProcessArguments(int argc, char** argv)
{
  try
  {
    cxxopts::Options options(LowerCase(opensn::program), "");

    /* clang-format off */
    options.add_options("User")
    ("h,help",                      "Help message")
    ("c,suppress-color",            "Suppress color output")
    ("v,verbose",                   "Verbosity level (0 to 3). Default is 0.", cxxopts::value<int>())
    ("caliper",                     "Enable Caliper reporting",
      cxxopts::value<std::string>()->implicit_value("runtime-report(calc.inclusive=true),max_column_width=80"))
    ("i,input",                     "Input file", cxxopts::value<std::string>())
    ("p,py",                        "Python expression",
      cxxopts::value<std::vector<std::string>>());

    options.add_options("Dev")
      ("help-dev",                  "Developer options help")
      ("allow-petsc-error-handler", "Allow PETSc error handler");
    /* clang-format on */

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      if (opensn::mpi_comm.rank() == 0)
        std::cout << options.help({"User"}) << std::endl;
      return 1;
    }

    if (result.count("help-dev"))
    {
      if (opensn::mpi_comm.rank() == 0)
        std::cout << options.help({"Dev"}) << std::endl;
      return 1;
    }

    if (result.count("verbose"))
    {
      int verbosity = result["verbose"].as<int>();
      opensn::log.SetVerbosity(verbosity);
    }

    if (result.count("allow-petsc-error-handler"))
      allow_petsc_error_handler_ = true;

    if (result.count("suppress-color"))
      opensn::suppress_color = true;

    if (result.count("caliper"))
    {
      opensn::use_caliper = true;
      opensn::cali_config = result["caliper"].as<std::string>();
    }

    if (result.count("input"))
    {
      opensn::input_path = result["input"].as<std::string>();
      sim_option_interactive_ = false;
    }

    if (result.count("py"))
    {
      for (const auto& pyarg : result["py"].as<std::vector<std::string>>())
        console.GetCommandBuffer().push_back(pyarg);
    }
  }
  catch (const std::exception& e)
  {
    if (opensn::mpi_comm.rank() == 0)
      std::cerr << e.what() << std::endl;
    return 1;
  }
  return 0;
}

int
PyApp::RunInteractive()
{
  if (std::filesystem::exists(input_path))
    console.ExecuteFile(opensn::input_path.string());
  console.RunConsoleLoop();

  return 0;
}

int
PyApp::RunBatch()
{
  if (std::filesystem::exists(opensn::input_path))
    return console.ExecuteFile(opensn::input_path.string());

  if (opensn::mpi_comm.rank() == 0)
    std::cerr << "Could not open file " << opensn::input_path.string() << "." << std::endl;
  return EXIT_FAILURE;
}

} // namespace opensnpy
