#pragma once

#include "framework/object.h"
#include "framework/parameters/input_parameters.h"

using namespace opensn;

namespace unit_tests
{

class TestSubObject : public Object
{
private:
  const size_t num_groups_;

public:
  static InputParameters GetInputParameters();
  explicit TestSubObject(const InputParameters& params);
};

class TestObject : public Object
{
private:
  const std::string solver_type_;
  TestSubObject sub_obj1_;
  TestSubObject sub_obj2_;

public:
  static InputParameters GetInputParameters();
  explicit TestObject(const InputParameters& params);
};

class ChildTestObject : public TestObject
{
private:
  const int num_sub_groups_;

public:
  static InputParameters GetInputParameters();
  explicit ChildTestObject(const InputParameters& params);
};

} // namespace unit_tests
