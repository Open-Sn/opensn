#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/lib/console.h"

using namespace opensn;

namespace unit_tests
{

void
TestCFunction()
{
  opensn::log.Log() << "Hello from a C function";
}

BIND_FUNCTION(unit_tests, TestCFunction);

} //  namespace unit_tests
