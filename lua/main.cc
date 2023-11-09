#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/event_system/system_wide_event_publisher.h"
#include "framework/event_system/event.h"
#include "framework/event_system/event_codes.h"
#include "modules/modules_lua.h"

/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int
main(int argc, char** argv)
{
  int location_id = 0, number_processes = 1;
  MPI_Comm communicator = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(communicator, &location_id);
  MPI_Comm_size(communicator, &number_processes);

  Chi::mpi.SetCommunicator(communicator);
  Chi::mpi.SetLocationID(location_id);
  Chi::mpi.SetProcessCount(number_processes);

  chi_modules::lua_utils::LoadRegisteredLuaItems();
  Chi::console.PostMPIInfo(location_id, number_processes);

  Chi::run_time::ParseArguments(argc, argv);

  Chi::run_time::InitPetSc(argc, argv);

  auto& t_main = Chi::log.CreateTimingBlock("ChiTech");
  t_main.TimeSectionBegin();
  chi::SystemWideEventPublisher::GetInstance().PublishEvent(
    chi::Event("ProgramStart", chi::GetStandardEventCode("ProgramStart")));

  int error_code;
  if (Chi::run_time::sim_option_interactive_) error_code = Chi::RunInteractive(argc, argv);
  else
    error_code = Chi::RunBatch(argc, argv);

  Chi::Finalize();

  return error_code;
}
