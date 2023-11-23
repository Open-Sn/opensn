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
  "     --suppress-beg-end-timelog  Suppresses time logs at the \n"
  "                                 beginning and end of execution.\n"
  "     --suppress-color            Suppresses the printing of color.\n"
  "                                 useful for unit tests requiring a diff.\n"
  "     --dump-object-registry      Dumps the object registry.\n"
  "\n\n\n";

LuaApp::LuaApp(MPI_Comm comm)
{
  int location_id = 0, number_processes = 1;
  MPI_Comm_rank(comm, &location_id);
  MPI_Comm_size(comm, &number_processes);

  opensn::mpi.SetCommunicator(comm);
  opensn::mpi.SetLocationID(location_id);
  opensn::mpi.SetProcessCount(number_processes);
}

LuaApp::~LuaApp()
{
}

int
LuaApp::InitPetSc(int argc, char** argv)
{
  PETSC_COMM_WORLD = mpi.comm;
  PetscOptionsInsertString(nullptr, "-error_output_stderr");
  if (not allow_petsc_error_handler_) PetscOptionsInsertString(nullptr, "-no_signal_handler");

  PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));

  return 0;
}

int
LuaApp::Run(int argc, char** argv)
{
  InitPetSc(argc, argv);
  opensn::Initialize();
  console.PostMPIInfo(opensn::mpi.location_id, opensn::mpi.process_count);

  int error_code = 0;
  ParseArguments(argc, argv);
  if (not termination_posted_)
  {
    if (sim_option_interactive_) error_code = RunInteractive(argc, argv);
    else
      error_code = RunBatch(argc, argv);
  }

  opensn::Finalize();
  PetscFinalize();

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
    else if (argument.find("--suppress-beg-end-timelog") != std::string::npos)
    {
      supress_beg_end_timelog_ = true;
    }
    else if (argument.find("--allow-petsc-error-handler") != std::string::npos)
    {
      allow_petsc_error_handler_ = true;
    }
    else if (argument.find("--suppress-color") != std::string::npos)
    {
      opensn::Chi::suppress_color_ = true;
    }
    else if (argument.find("--dump-object-registry") != std::string::npos)
    {
      dump_registry_ = true;
      termination_posted_ = true;
    }
    // No-graphics option
    else if (argument.find("-b") != std::string::npos) { sim_option_interactive_ = false; } //-b
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
      input_file_name_ = argument;
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
  if (not supress_beg_end_timelog_)
  {
    opensn::log.Log() << Timer::GetLocalDateTimeString() << " Running " << opensn::name
                      << " in interactive-mode with " << mpi.process_count << " processes.";

    opensn::log.Log() << opensn::name << " version " << GetVersionStr();
  }

  opensn::log.Log() << opensn::name << " number of arguments supplied: " << argc - 1;

  opensn::log.LogAll();

  console.FlushConsole();

  const auto& input_fname = input_file_name_;

  if (not input_fname.empty())
  {
    try
    {
      console.ExecuteFile(input_fname, argc, argv);
    }
    catch (const std::exception& excp)
    {
      opensn::log.LogAllError() << excp.what();
      // No quitting if file execution fails
    }
  }

  console.RunConsoleLoop();

  if (not supress_beg_end_timelog_)
  {
    opensn::log.Log() << "Final program time " << program_timer.GetTimeString();
    opensn::log.Log() << Timer::GetLocalDateTimeString() << " " << opensn::name
                      << " finished execution.";
  }

  return 0;
}

int
LuaApp::RunBatch(int argc, char** argv)
{
  if (not supress_beg_end_timelog_)
  {
    opensn::log.Log() << Timer::GetLocalDateTimeString() << " Running " << opensn::name
                      << " in batch-mode with " << mpi.process_count << " processes.";

    opensn::log.Log() << opensn::name << " version " << GetVersionStr();
  }

  opensn::log.Log() << opensn::name << " number of arguments supplied: " << argc - 1;

  if (argc <= 1) opensn::log.Log() << command_line_help_string_;
  console.FlushConsole();

#ifndef NDEBUG
  opensn::log.Log() << "Waiting...";
  if (mpi.location_id == 0)
    for (int k = 0; k < 2; ++k)
    {
      usleep(1000000);
      opensn::log.Log() << k;
    }

  mpi.Barrier();
#endif

  const auto& input_fname = input_file_name_;
  int error_code = 0;

  if ((not input_fname.empty()) and (not termination_posted_))
  {
    try
    {
      error_code = console.ExecuteFile(input_fname, argc, argv);
    }
    catch (const std::exception& excp)
    {
      opensn::log.LogAllError() << excp.what();
      Exit(EXIT_FAILURE);
    }
  }

  if (not supress_beg_end_timelog_)
  {
    opensn::log.Log() << "\nFinal program time " << program_timer.GetTimeString();
    opensn::log.Log() << Timer::GetLocalDateTimeString() << " " << opensn::name
                      << " finished execution of " << input_file_name_;
  }

  return error_code;
}

} // namespace opensnlua
