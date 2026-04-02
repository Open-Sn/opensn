#include "framework/utils/utils.h"
#include "framework/parameters/input_parameters.h"
#include "framework/data_types/allowable_range.h"
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

TEST(ParametersTest, AssignAParameter)
{
  ParameterBlock name("name", "this_is_my_name");

  ParameterBlock input("params");
  input.AddParameter(name);

  auto params = GetInputParameters();
  params.AssignParameters(input);

  EXPECT_EQ(params.GetParamValue<std::string>("name"), "this_is_my_name");
}

TEST(ParametersTest, AssignNonexistentParameter)
{
  ParameterBlock non("non-existent", 42.);

  ParameterBlock input("params");
  input.AddParameter(non);

  auto params = GetInputParameters();
  EXPECT_THROW(params.AssignParameters(input), std::invalid_argument);
}

TEST(ParametersTest, AssignIncorrectlyTypedParameter)
{
  ParameterBlock name("name", 42.);

  ParameterBlock input("params");
  input.AddParameter(name);

  auto params = GetInputParameters();
  EXPECT_THROW(params.AssignParameters(input), std::invalid_argument);
}

TEST(ParametersTest, RequiredMissingAndRenamed)
{
  InputParameters params;
  params.SetObjectType("UnitTest");
  params.AddRequiredParameter<int>("count", "Count");

  ParameterBlock input_missing("params");
  EXPECT_THROW(params.AssignParameters(input_missing), std::invalid_argument);

  params.MarkParameterRenamed("count", "Use new_count instead");
  ParameterBlock input_with_old("params");
  input_with_old.AddParameter(ParameterBlock("count", 3));
  EXPECT_THROW(params.AssignParameters(input_with_old), std::invalid_argument);
}

TEST(ParametersTest, DeprecationAndConstraints)
{
  InputParameters params;
  params.SetObjectType("UnitTest");
  params.AddOptionalParameter<int>("mode", 1, "Mode");
  params.AddOptionalParameter<double>("scale", 1.0, "Scale");

  params.MarkParameterDeprecatedError("mode", "Deprecated mode");
  params.ConstrainParameterRange("scale", AllowableRangeLowHighLimit::New(0.0, 2.0));

  ParameterBlock input("params");
  input.AddParameter(ParameterBlock("mode", 2));
  input.AddParameter(ParameterBlock("scale", 3.0));

  EXPECT_THROW(params.AssignParameters(input), std::runtime_error);

  ParameterBlock input_bad_scale("params");
  input_bad_scale.AddParameter(ParameterBlock("scale", 3.0));
  EXPECT_THROW(params.AssignParameters(input_bad_scale), std::invalid_argument);
}
