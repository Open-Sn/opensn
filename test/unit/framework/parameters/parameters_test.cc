#include "test/unit/opensn_unit_test.h"
#include "framework/utils/utils.h"
#include "framework/parameters/input_parameters.h"
#include <gmock/gmock.h>

using namespace opensn;

namespace
{

InputParameters
GetInputParameters()
{
  InputParameters params;
  params.AddRequiredParameter<std::string>("name", "Name");
  return params;
}

} // namespace

class ParametersTest : public OpenSnUnitTest
{
};

TEST_F(ParametersTest, AssignAParameter)
{
  ParameterBlock name("name", "this_is_my_name");

  ParameterBlock input("params");
  input.AddParameter(name);

  auto params = GetInputParameters();
  params.AssignParameters(input);

  EXPECT_EQ(params.GetParamValue<std::string>("name"), "this_is_my_name");
}

TEST_F(ParametersTest, AssignNonexistentParameter)
{
  ParameterBlock non("non-existent", 42.);

  ParameterBlock input("params");
  input.AddParameter(non);

  auto params = GetInputParameters();
  EXPECT_THROW(params.AssignParameters(input), std::invalid_argument);
}

TEST_F(ParametersTest, AssignIncorrectlyTypedParameter)
{
  ParameterBlock name("name", 42.);

  ParameterBlock input("params");
  input.AddParameter(name);

  auto params = GetInputParameters();
  EXPECT_THROW(params.AssignParameters(input), std::invalid_argument);
}
