#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/parameters/input_parameters.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock ExampleTest(const InputParameters&);

RegisterWrapperFunction(chi_unit_tests, ExampleTest, nullptr, ExampleTest);

ParameterBlock
ExampleTest(const InputParameters&)
{
  opensn::log.Log() << "This is an example test";

  return ParameterBlock();
}

} //  namespace unit_tests
