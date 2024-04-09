#include "lua/framework/lua_app.h"
#include "lua/framework/console/console.h"
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
    std::cout << opensn::name << " version " << GetVersionStr() << "\n"
              << Timer::GetLocalDateTimeString() << " Running " << opensn::name << " with "
              << opensn::mpi_comm.size() << " processes.\n"
              << opensn::name << " number of arguments supplied: " << argc - 1 << "\n"
              << std::endl;
  }

  int error_code = ProcessArguments(argc, argv);

  if (!error_code)
  {
    InitPetSc(argc, argv);
    opensn::Initialize();
    console.PostMPIInfo(opensn::mpi_comm.rank(), opensn::mpi_comm.size());
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
                << Timer::GetLocalDateTimeString() << " " << opensn::name << " finished execution."
                << std::endl;
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
    cxxopts::Options options(LowerCase(opensn::name), "");

    /* clang-format off */
    options.add_options("User")
    ("h,help",                      "Help message")
    ("c,suppress-color",            "Suppress color output")
    ("v,verbose",                   "Verbosity level (0 to 3). Default is 0.",
      cxxopts::value<int>()->implicit_value("0"))
    ("caliper",                     "Enable Caliper reporting",
      cxxopts::value<std::string>()->implicit_value("runtime-report(calc.inclusive=true),max_column_width=80"))
    ("positional",                  "Positional arugments", cxxopts::value<std::vector<std::string>>());

    options.add_options("Dev")
      ("help-dev",                  "Developer options help")
      ("allow-petsc-error-handler", "Allow PETSc error handler")
      ("dump-object-registry",      "Dump object registry");
    /* clang-format on */

    options.positional_help("[FILE] [LUA ASSIGNMENT]");
    options.parse_positional("positional");

    bool found_filename = false;
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

    if (result.count("verbose"))
    {
      int verbosity = result["verbose"].as<int>();
      opensn::log.SetVerbosity(verbosity);
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

    if (result.count("positional"))
    {
      auto args = result["positional"].as<std::vector<std::string>>();
      for (const auto& arg : args)
      {
        if (arg.find('=') != std::string::npos)
          console.GetCommandBuffer().push_back(arg);
        else
        {
          if (not found_filename)
          {
            opensn::input_path = arg;
            sim_option_interactive_ = false;
            found_filename = true;
          }
          else
            throw std::runtime_error("Invalid option " + arg);
        }
      }
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
#ifndef NDEBUG
  opensn::log.Log() << "Waiting...";
  if (opensn::mpi_comm.rank() == 0)
    for (int k = 0; k < 2; ++k)
    {
      usleep(1000000);
      opensn::log.Log() << k;
    }
  mpi_comm.barrier();
#endif

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
