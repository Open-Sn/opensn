// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_app.h"
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

PyApp::PyApp(const mpi::Communicator& comm) : allow_petsc_error_handler_(false)
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

  console.BindModule(WrapQuadraturePointPhiTheta);
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

  console.BindModule(WrapProblem);
  console.BindModule(WrapSolver);
  console.BindModule(WrapLBS);
  console.BindModule(WrapSteadyState);
  console.BindModule(WrapKEigen);
  console.BindModule(WrapNLKEigen);
  console.BindModule(WrapPIteration);
  console.BindModule(WrapDiscreteOrdinatesKEigenAcceleration);
}

int
PyApp::InitPETSc(int argc, char** argv)
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

  if (ProcessArguments(argc, argv))
  {
    PetscOptionsSetValue(NULL, "-options_left", "0");
    InitPETSc(argc, argv);

    opensn::Initialize();
    console.InitConsole();
    console.ExecuteFile(opensn::input_path.string());
    opensn::Finalize();

    PetscFinalize();

    if (opensn::mpi_comm.rank() == 0)
    {
      std::cout << "\nElapsed execution time: " << program_timer.GetTimeString() << "\n"
                << Timer::GetLocalDateTimeString() << " " << opensn::program
                << " finished execution." << std::endl;
    }
  }
  else
    return EXIT_FAILURE;

  cali_mgr.flush();
  return EXIT_SUCCESS;
}

bool
PyApp::ProcessArguments(int argc, char** argv)
{
  cxxopts::Options options(LowerCase(opensn::program), "");

  try
  {
    /* clang-format off */
    options.add_options("User")
    ("h,help",                      "Help message")
    ("c,suppress-color",            "Suppress color output")
    ("v,verbose",                   "Verbosity level (0 to 3). Default is 0.", cxxopts::value<int>())
    ("caliper",                     "Enable Caliper reporting",
      cxxopts::value<std::string>()->implicit_value("runtime-report(calc.inclusive=true),max_column_width=80"))
    ("i,input",                     "Input file", cxxopts::value<std::string>())
    ("allow-petsc-error-handler",   "Allow PETSc error handler")
    ("p,py",                        "Python expression", cxxopts::value<std::vector<std::string>>());
    /* clang-format on */

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      if (opensn::mpi_comm.rank() == 0)
        std::cout << options.help({"User"}) << std::endl;
      return false;
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

    if (result.count("py"))
    {
      for (const auto& pyarg : result["py"].as<std::vector<std::string>>())
        console.GetCommandBuffer().push_back(pyarg);
    }

    opensn::input_path = result["input"].as<std::string>();
    if (not std::filesystem::exists(input_path) or not std::filesystem::is_regular_file(input_path))
    {
      if (opensn::mpi_comm.rank() == 0)
        std::cerr << "Invalid input file: " << input_path.string() << "\n" << std::endl;
      return false;
    }
  }
  catch (const std::exception& e)
  {
    if (opensn::mpi_comm.rank() == 0)
      std::cerr << e.what() << "\n" << options.help({"User"}) << std::endl;
    return false;
  }

  return true;
}

} // namespace opensnpy
