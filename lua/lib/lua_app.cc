// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/lua_app.h"
#include "lua/lib/console.h"
#include "framework/event_system/event.h"
#include "framework/event_system/system_wide_event_publisher.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"
#include "caliper/cali.h"
#include "cxxopts/cxxopts.h"
#include "petsc.h"
#include <string>
#ifndef NDEBUG
#include <unistd.h>
#endif

using namespace opensn;

namespace opensnlua
{

LuaApp::LuaApp(const mpi::Communicator& comm)
  : sim_option_interactive_(true), allow_petsc_error_handler_(false)
{
  opensn::mpi_comm = comm;

  auto& L = Console::GetInstance().GetConsoleState();
  luabridge::getGlobalNamespace(L)
    .addVariable("location_id", opensn::mpi_comm.rank())
    .addVariable("number_of_processes", opensn::mpi_comm.size());
}

int
LuaApp::InitPetSc(int argc, char** argv)
{
  PetscOptionsInsertString(nullptr, "-error_output_stderr");

  if (not allow_petsc_error_handler_)
    PetscOptionsInsertString(nullptr, "-no_signal_handler");

  PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));

  return 0;
}

int
LuaApp::Run(int argc, char** argv)
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
      error_code = RunInteractive(argc, argv);
    else
      error_code = RunBatch(argc, argv);

    opensn::Finalize();
    PetscFinalize();

    if (opensn::mpi_comm.rank() == 0)
    {
      std::cout << "\n"
                << "Elapsed execution time: " << program_timer.GetTimeString() << "\n"
                << Timer::GetLocalDateTimeString() << " " << opensn::program
                << " finished execution." << std::endl;
    }
  }

  if (opensn::mpi_comm.rank() == 0)
    std::cout << std::endl;
  cali_mgr.flush();

  return error_code;
}

int
LuaApp::ProcessArguments(int argc, char** argv)
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
    ("l,lua",                       "Lua expression",
      cxxopts::value<std::vector<std::string>>());

    options.add_options("Dev")
      ("help-dev",                  "Developer options help")
      ("allow-petsc-error-handler", "Allow PETSc error handler")
      ("dump-object-registry",      "Dump object registry");
    /* clang-format on */

    auto result = options.parse(argc, argv);

    // Note that the order of evaluation of the command-line options is important!
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

    if (result.count("dump-object-registry"))
    {
      ObjectFactory::GetInstance().DumpRegister();
      console.DumpRegister();
      return 1;
    }

    if (result.count("caliper"))
    {
      opensn::use_caliper = true;
      auto config = result["caliper"].as<std::string>();
      if (!config.empty())
      {
        std::string error = cali_mgr.check(config.c_str());
        if (!error.empty())
          throw std::runtime_error("Invalid Caliper config: " + config);
        opensn::cali_config = config;
      }
    }

    if (result.count("input"))
    {
      opensn::input_path = result["input"].as<std::string>();
      sim_option_interactive_ = false;
    }

    if (result.count("lua"))
    {
      for (auto larg : result["lua"].as<std::vector<std::string>>())
        console.GetCommandBuffer().push_back(larg);
    }
  }
  catch (const cxxopts::exceptions::exception& e)
  {
    if (opensn::mpi_comm.rank() == 0)
      std::cerr << e.what() << std::endl;
    return 1;
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
LuaApp::RunInteractive(int argc, char** argv)
{
  if (std::filesystem::exists(input_path))
  {
    try
    {
      console.ExecuteFile(opensn::input_path.string(), argc, argv);
    }
    catch (const std::exception& excp)
    {
      opensn::log.LogAllError() << excp.what();
      // No quitting if file execution fails
    }
  }

  console.RunConsoleLoop();

  return 0;
}

int
LuaApp::RunBatch(int argc, char** argv)
{
  int error_code = 0;
  if (std::filesystem::exists(input_path))
  {
    try
    {
      error_code = console.ExecuteFile(opensn::input_path.string(), argc, argv);
    }
    catch (const std::exception& excp)
    {
      opensn::log.LogAllError() << excp.what();
      error_code = EXIT_FAILURE;
    }
  }
  else
  {
    opensn::log.Log0Error() << "Could not open file " << opensn::input_path.string() << ".";
    error_code = EXIT_FAILURE;
  }

  return error_code;
}

} // namespace opensnlua
