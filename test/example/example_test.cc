#include "framework/console/console.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock ExampleTest(const InputParameters&);

RegisterWrapperFunction(chi_unit_tests, ExampleTest, nullptr, ExampleTest);

ParameterBlock
ExampleTest(const InputParameters&)
{
  opensn::Chi::log.Log() << "This is an example test";

  return ParameterBlock();
}

} //  namespace unit_tests
