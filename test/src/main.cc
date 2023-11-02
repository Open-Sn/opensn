#include "chi_runtime.h"

int main(int argc, char** argv)
{
  Chi::Initialize(argc, argv, Chi::mpi.comm);

  int error_code;
  if (Chi::run_time::sim_option_interactive_)
    error_code = Chi::RunInteractive(argc, argv);
  else
    error_code = Chi::RunBatch(argc, argv);

  Chi::Finalize();

  return error_code;
}