#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/event_system/system_wide_event_publisher.h"
#include "framework/event_system/event.h"
#include "modules/modules_lua.h"

using namespace opensn;

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

  opensn::mpi.SetCommunicator(communicator);
  opensn::mpi.SetLocationID(location_id);
  opensn::mpi.SetProcessCount(number_processes);

  opensnlua::LoadRegisteredLuaItems();
  opensn::console.PostMPIInfo(location_id, number_processes);

  opensn::Chi::run_time::ParseArguments(argc, argv);

  opensn::Chi::run_time::InitPetSc(argc, argv);

  auto& t_main = opensn::log.CreateTimingBlock("ChiTech");
  t_main.TimeSectionBegin();
  opensn::SystemWideEventPublisher::GetInstance().PublishEvent(Event("ProgramStart"));

  int error_code;
  if (opensn::Chi::run_time::sim_option_interactive_)
    error_code = opensn::Chi::RunInteractive(argc, argv);
  else
    error_code = opensn::Chi::RunBatch(argc, argv);

  opensn::Finalize();

  return error_code;
}
