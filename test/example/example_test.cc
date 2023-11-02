#include "console/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_unit_tests
{

chi::ParameterBlock ExampleTest(const chi::InputParameters&);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/ExampleTest,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/ExampleTest);

chi::ParameterBlock ExampleTest(const chi::InputParameters&)
{
  Chi::log.Log() << "This is an example test";

  return chi::ParameterBlock();
}

} // namespace chi_unit_tests