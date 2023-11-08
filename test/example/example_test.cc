#include "framework/console/chi_console.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace chi_unit_tests
{

chi::ParameterBlock ExampleTest(const chi::InputParameters&);

RegisterWrapperFunction(chi_unit_tests, ExampleTest, nullptr, ExampleTest);

chi::ParameterBlock
ExampleTest(const chi::InputParameters&)
{
  Chi::log.Log() << "This is an example test";

  return chi::ParameterBlock();
}

} // namespace chi_unit_tests
