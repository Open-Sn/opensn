#include <gtest/gtest.h>
#include "framework/runtime.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include "petsc.h"

int
main(int argc, char** argv)
{
  mpicpp_lite::Environment mpi_env;
  opensn::mpi_comm = MPI_COMM_WORLD;

  PetscInitialize(nullptr, nullptr, nullptr, nullptr);

  testing::InitGoogleTest(&argc, argv);
  const auto exit_code = RUN_ALL_TESTS();

  PetscFinalize();

  return exit_code;
}
