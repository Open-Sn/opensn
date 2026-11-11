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

TEST(ParametersTest, ArrayElementsKeepIndexOrder)
{
  // Appended in index order, past "9" where string and numeric order differ.
  ParameterBlock in_order("in_order");
  in_order.ChangeToArray();
  std::vector<int> expected;
  for (int k = 0; k < 25; ++k)
  {
    in_order.AddParameter(std::to_string(k), 100 + k);
    expected.push_back(100 + k);
  }
  EXPECT_EQ(in_order.GetVectorValue<int>(), expected);

  // From a std::vector.
  EXPECT_EQ(ParameterBlock("from_vector", expected).GetVectorValue<int>(), expected);

  // Added out of order: sorted by index.
  ParameterBlock out_of_order("out_of_order");
  out_of_order.ChangeToArray();
  out_of_order.AddParameter("1", 11);
  out_of_order.AddParameter("0", 10);
  out_of_order.AddParameter("2", 12);
  EXPECT_EQ(out_of_order.GetVectorValue<int>(), std::vector<int>({10, 11, 12}));

  // With a gap: the name equal to the size goes before the larger index.
  ParameterBlock with_gap("with_gap");
  with_gap.ChangeToArray();
  with_gap.AddParameter("0", 10);
  with_gap.AddParameter("5", 15);
  with_gap.AddParameter("2", 12);
  EXPECT_EQ(with_gap.GetVectorValue<int>(), std::vector<int>({10, 12, 15}));

  // A duplicate index is rejected.
  EXPECT_THROW(in_order.AddParameter("3", 0), std::invalid_argument);
  EXPECT_THROW(with_gap.AddParameter("5", 0), std::invalid_argument);

  // A block converted to an array is sorted by index, not by name. Here the names sorted as
  // strings would end in "9" (= size - 1), but sorted as indices they end in "10", so adding
  // "10" again takes the checked path and is rejected.
  ParameterBlock converted("converted");
  for (int k = 0; k <= 10; ++k)
    if (k != 7)
      converted.AddParameter(std::to_string(k), 100 + k);
  converted.ChangeToArray();
  EXPECT_EQ(converted.GetVectorValue<int>(),
            std::vector<int>({100, 101, 102, 103, 104, 105, 106, 108, 109, 110}));
  EXPECT_THROW(converted.AddParameter("10", 0), std::invalid_argument);
}
