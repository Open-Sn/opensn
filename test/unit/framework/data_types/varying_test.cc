#include "test/unit/opensn_unit_test.h"
#include "framework/data_types/varying.h"
#include <gmock/gmock.h>

using namespace opensn;

class VaryingTest : public OpenSnUnitTest
{
};

TEST_F(VaryingTest, Test1)
{
  Varying v(12);
  EXPECT_EQ(v.GetIntegerValue(), 12);

  v = true;
  EXPECT_TRUE(v.GetBoolValue());
  EXPECT_TRUE(v.GetValue<bool>());
}

TEST_F(VaryingTest, Test2)
{
  Varying v(12);
  EXPECT_EQ(v.GetIntegerValue(), 12);
  v = 12.0;
  EXPECT_DOUBLE_EQ(v.GetFloatValue(), 12.);
  EXPECT_DOUBLE_EQ(v.GetValue<double>(), 12.);
  EXPECT_DOUBLE_EQ(v.GetValue<float>(), 12.);
}

TEST_F(VaryingTest, Test3)
{
  Varying v(12.0);
  EXPECT_EQ(v.GetByteSize(), 8);
  EXPECT_DOUBLE_EQ(v.GetFloatValue(), 12.);
  v = 12;
  EXPECT_EQ(v.GetIntegerValue(), 12);
  EXPECT_EQ(v.GetValue<int>(), 12);
  EXPECT_EQ(v.GetValue<size_t>(), 12);
}

TEST_F(VaryingTest, Test4)
{
  Varying v(std::string("Hello"));
  EXPECT_EQ(v.GetStringValue(), "Hello");
  EXPECT_EQ(v.GetValue<std::string>(), "Hello");
}
