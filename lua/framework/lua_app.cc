#include "lua/framework/lua_app.h"
#include "framework/logging/log.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"
#include "framework/utils/timer.h"
#include "framework/event_system/event.h"
#include "framework/event_system/system_wide_event_publisher.h"
#ifndef NDEBUG
#include <unistd.h>
#endif
#include "petsc.h"

using namespace opensn;

namespace opensnlua
{

const std::string command_line_help_string_ =
  "\nUsage: exe inputfile [options values]\n"
  "\n"
  "     -v                          Level of verbosity. Default 0.\n"
  "                                 Can be either 0, 1 or 2.\n"
  "     a=b                         Executes argument as a lua string. "
  "i.e. x=2 or y=[[\"string\"]]\n"
  "     --allow-petsc-error-handler Allows petsc error handler.\n"
  "     --suppress-color            Suppresses the printing of color.\n"
  "                                 useful for unit tests requiring a diff.\n"
  "     --dump-object-registry      Dumps the object registry.\n"
  "\n\n\n";

LuaApp::LuaApp(const mpi::Communicator& comm)
{
  opensn::mpi_comm = comm;
}

LuaApp::~LuaApp()
{
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
  InitPetSc(argc, argv);
  opensn::Initialize();
  console.PostMPIInfo(opensn::mpi_comm.rank(), opensn::mpi_comm.size());

  opensn::log.Log() << opensn::name << " version " << GetVersionStr();
  opensn::log.Log() << Timer::GetLocalDateTimeString() << " Running " << opensn::name
                    << " in interactive-mode with " << opensn::mpi_comm.size() << " processes.";
  opensn::log.Log() << opensn::name << " number of arguments supplied: " << argc - 1;
  opensn::log.LogAll();
  console.FlushConsole();

  int error_code = 0;
  ParseArguments(argc, argv);
  if (not termination_posted_)
  {
    if (sim_option_interactive_)
      error_code = RunInteractive(argc, argv);
    else
      error_code = RunBatch(argc, argv);
  }

  opensn::Finalize();
  PetscFinalize();

  opensn::log.Log() << "Elapsed execution time: " << program_timer.GetTimeString();
  opensn::log.Log() << Timer::GetLocalDateTimeString() << " " << opensn::name
                    << " finished execution.";

  return error_code;
}

void
LuaApp::ParseArguments(int argc, char** argv)
{
  bool input_file_found = false;
  for (int i = 1; i < argc; i++)
  {
    std::string argument(argv[i]);

    opensn::log.Log() << "Parsing argument " << i << " " << argument;

    if (argument.find("-h") != std::string::npos or argument.find("--help") != std::string::npos)
    {
      opensn::log.Log() << command_line_help_string_;
      termination_posted_ = true;
    }
    else if (argument.find("--allow-petsc-error-handler") != std::string::npos)
    {
      allow_petsc_error_handler_ = true;
    }
    else if (argument.find("--suppress-color") != std::string::npos)
    {
      opensn::suppress_color = true;
    }
    else if (argument.find("--dump-object-registry") != std::string::npos)
    {
      dump_registry_ = true;
      termination_posted_ = true;
    }
    // No-graphics option
    else if (argument.find("-b") != std::string::npos)
    {
      sim_option_interactive_ = false;
    } //-b
    // Verbosity
    else if (argument.find("-v") != std::string::npos)
    {
      if ((i + 1) >= argc)
      {
        std::cerr << "Invalid option used with command line argument "
                     "-v. Options are 0,1 or 2."
                  << std::endl;
        Exit(EXIT_FAILURE);
      }
      else
      {
        std::string v_option(argv[i + 1]);
        try
        {
          int level = std::stoi(v_option);
          opensn::log.SetVerbosity(level);
        }
        catch (const std::invalid_argument& e)
        {
          std::cerr << "Invalid option used with command line argument "
                       "-v. Options are 0,1 or 2."
                    << std::endl;
          Exit(EXIT_FAILURE);
        }
      }

    } //-v
    else if ((argument.find('=') == std::string::npos) and (not input_file_found))
    {
      opensn::input_path = argument;
      input_file_found = true;
      sim_option_interactive_ = false;
    } // no =
    else if (argument.find('=') != std::string::npos)
    {
      console.GetCommandBuffer().push_back(argument);
    }
  }

  if (dump_registry_)
  {
    ObjectFactory::GetInstance().DumpRegister();
    console.DumpRegister();
  }
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
  else
    opensn::log.Log0Error() << "Could not open file " << opensn::input_path.string() << ".";

  console.RunConsoleLoop();

  return 0;
}

int
LuaApp::RunBatch(int argc, char** argv)
{
  if (argc <= 1)
    opensn::log.Log() << command_line_help_string_;
  console.FlushConsole();

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
  if (std::filesystem::exists(input_path) and (not termination_posted_))
  {
    try
    {
      error_code = console.ExecuteFile(opensn::input_path.string(), argc, argv);
    }
    catch (const std::exception& excp)
    {
      opensn::log.LogAllError() << excp.what();
      Exit(EXIT_FAILURE);
    }
  }
  else
  {
    opensn::log.Log0Error() << "Could not open file " << opensn::input_path.string() << ".";
    Exit(EXIT_FAILURE);
  }

  return error_code;
}

} // namespace opensnlua
