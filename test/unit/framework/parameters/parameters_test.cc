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

  // A name that is not an index is rejected, and the block is left unchanged.
  ParameterBlock not_indexed("not_indexed");
  not_indexed.AddParameter("10", 110);
  not_indexed.AddParameter("9", 109);
  not_indexed.AddParameter("x", 0);
  EXPECT_THROW(not_indexed.ChangeToArray(), std::logic_error);
  EXPECT_EQ(not_indexed.GetType(), ParameterBlockType::BLOCK);
  std::vector<std::string> names;
  for (const auto& param : not_indexed.GetParameters())
    names.push_back(param.GetName());
  EXPECT_EQ(names, std::vector<std::string>({"10", "9", "x"}));
}

TEST(ParametersTest, IntegerForFloatParameter)
{
  InputParameters params;
  params.SetObjectType("UnitTest");
  params.AddOptionalParameter<double>("tolerance", 1.0e-6, "Tolerance");
  params.AddOptionalParameter<int>("count", 1, "Count");
  params.ConstrainParameterRange("tolerance", AllowableRangeLowHighLimit::New(0.0, 10.0));

  ParameterBlock input("params");
  input.AddParameter(ParameterBlock("tolerance", 2));
  params.AssignParameters(input);
  EXPECT_EQ(params.GetParamValue<double>("tolerance"), 2.0);

  // The range of the floating point parameter applies to the converted value.
  ParameterBlock out_of_range("params");
  out_of_range.AddParameter(ParameterBlock("tolerance", 20));
  EXPECT_THROW(params.AssignParameters(out_of_range), std::invalid_argument);

  // A floating point value is not converted to an integer.
  ParameterBlock float_for_integer("params");
  float_for_integer.AddParameter(ParameterBlock("count", 2.5));
  EXPECT_THROW(params.AssignParameters(float_for_integer), std::invalid_argument);
}

TEST(ParametersTest, ArrayOfIntegersAndFloats)
{
  ParameterBlock edges("edges");
  edges.ChangeToArray();
  edges.AddParameter("0", 0);
  edges.AddParameter("1", 1.5);
  edges.AddParameter("2", 3);
  EXPECT_EQ(edges.GetVectorValue<double>(), std::vector<double>({0.0, 1.5, 3.0}));
  // An integer array is not read as floating point values with a fractional entry.
  EXPECT_THROW(edges.GetVectorValue<int>(), std::logic_error);

  ParameterBlock with_string("with_string");
  with_string.ChangeToArray();
  with_string.AddParameter("0", 1.0);
  with_string.AddParameter("1", std::string("two"));
  EXPECT_THROW(with_string.GetVectorValue<double>(), std::logic_error);

  ParameterBlock negative("negative");
  negative.ChangeToArray();
  negative.AddParameter("0", 1);
  negative.AddParameter("1", -1);
  EXPECT_THROW(negative.GetVectorValue<unsigned int>(), std::logic_error);
}
