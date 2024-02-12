#include "params_test_01.h"

#include "framework/object_factory.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace unit_tests
{

OpenSnRegisterObjectInNamespace(unit_testsB, TestObject);

InputParameters
TestObject::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("\\defgroup unit_tests__TestObject unit_testsB.TestObject\n"
                               "\\ingroup DocUnitTests\n"
                               "General test object");

  params.AddOptionalParameter("solver_type", "A", "The solver type.");
  params.AddRequiredParameter<std::string>("coupled_field", "The text name of the coupled field.");
  params.AddRequiredParameterBlock("sub_obj1",
                                   "A block of parameters for unit_testsB::TestSubObject");

  ParameterBlock sub_obj2_param_block("sub_obj2");
  sub_obj2_param_block.AddParameter("num_groups", 99);
  params.AddOptionalParameterBlock(
    "sub_obj2", sub_obj2_param_block, "A block of parameters for unit_testsB::TestSubObject");

  params.AddOptionalParameter("limiter_type", 1, "Type of limiter to use in the solver");
  params.MarkParamaterDeprecatedWarning("limiter_type");

  params.AddOptionalParameter("scheme", "Zorba", "What scheme to use");
  params.MarkParamaterDeprecatedError("scheme");

  params.AddRequiredParameter<bool>("format", "What output format to use");
  params.MarkParamaterDeprecatedError("format");

  params.AddOptionalParameter("use_my_stuff", false, "Yeah please do");
  params.MarkParamaterRenamed("use_my_stuff", "Renamed to \"use_zaks_stuff\".");

  params.AddRequiredParameter<bool>("use_ragusas_stuff", "If you want");
  params.MarkParamaterRenamed("use_ragusas_stuff", "Renamed to \"use_complicated_stuff\".");

  params.AddOptionalParameter<int>(
    "groupset_num_subsets",
    1,
    "The number of subsets to apply to the set of groups in this set. This is "
    "useful for increasing pipeline size for parallel simulations");

  params.ConstrainParameterRange("groupset_num_subsets", AllowableRangeLowLimit::New<int>(1));
  return params;
}

TestObject::TestObject(const InputParameters& params)
  : solver_type_(params.GetParamValue<std::string>("solver_type")),
    sub_obj1_(InputParameters::MakeForObject<TestSubObject>(params.GetParam("sub_obj1"))),
    sub_obj2_(InputParameters::MakeForObject<TestSubObject>(params.GetParam("sub_obj2")))
{
  opensn::log.Log() << "TestObject created "
                    << "solver_type=" << solver_type_;
}

OpenSnRegisterObjectInNamespace(unit_testsB, TestSubObject);

InputParameters
TestSubObject::GetInputParameters()
{
  InputParameters params;

  params.SetGeneralDescription("\\defgroup unit_tests__TestSubObject unit_testsB.TestSubObject\n"
                               "\\ingroup DocUnitTests\n"
                               "General test sub-object");

  params.AddRequiredParameter<size_t>("num_groups", "Number of groups to use in the simulation");

  return params;
}

TestSubObject::TestSubObject(const InputParameters& params)
  : num_groups_(params.GetParamValue<size_t>("num_groups"))
{
  opensn::log.Log() << "TestSubObject created "
                    << "num_groups=" << num_groups_;
}

OpenSnRegisterObjectInNamespace(unit_testsB, ChildTestObject);

InputParameters
ChildTestObject::GetInputParameters()
{
  InputParameters params = TestObject::GetInputParameters();

  params.SetGeneralDescription(
    "\\defgroup unit_tests__ChildTestObject unit_testsB.ChildTestObject\n"
    "\\ingroup DocUnitTests\n"
    "General test child-object inheriting option from parent");

  params.ChangeExistingParamToOptional("coupled_field", "Q");
  params.ChangeExistingParamToRequired<std::string>("solver_type");

  params.AddOptionalParameter("num_sub_groups", 1, "Number of sub-groups to use in the simultion");

  return params;
}

ChildTestObject::ChildTestObject(const InputParameters& params)
  : TestObject(params), num_sub_groups_(params.GetParamValue<int>("num_sub_groups"))
{
  opensn::log.Log() << "ChildTestObject created "
                    << "num_sub_groups=" << num_sub_groups_;
}

} // namespace unit_tests
