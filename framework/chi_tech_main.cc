#include "chi_runtime.h"

//######################################################### Program entry point
/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int
main(int argc, char** argv)
{
  Chi::Initialize(argc, argv, MPI_COMM_WORLD);

  int error_code;
  if (Chi::run_time::sim_option_interactive_) error_code = Chi::RunInteractive(argc, argv);
  else
    error_code = Chi::RunBatch(argc, argv);

  Chi::Finalize();

  return error_code;
}