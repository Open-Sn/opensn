#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/parameters/input_parameters.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock TestCFunction(const InputParameters&);

RegisterWrapperFunctionInNamespace(unit_tests, TestCFunction, nullptr, TestCFunction);

ParameterBlock
TestCFunction(const InputParameters&)
{
  opensn::log.Log() << "Hello from a C function";
  return ParameterBlock();
}

} //  namespace unit_tests
