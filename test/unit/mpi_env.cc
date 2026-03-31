#include "test/unit/mpi_env.h"
#include "framework/runtime.h"

void
MPIEnvironment::SetUp()
{
  opensn::Initialize();
}

void
MPIEnvironment::TearDown()
{
  opensn::Finalize();
}
