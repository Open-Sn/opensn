#include "test/unit/opensn_unit_test.h"
#include "framework/runtime.h"

void
OpenSnUnitTest::SetUp()
{
  opensn::Initialize();
}

void
OpenSnUnitTest::TearDown()
{
  opensn::Finalize();
}
