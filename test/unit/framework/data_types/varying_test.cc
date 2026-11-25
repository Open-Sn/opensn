#include "framework/data_types/varying.h"
#include <gmock/gmock.h>

using namespace opensn;

TEST(VaryingTest, Test1)
{
  Varying v(12);
  EXPECT_EQ(v.GetIntegerValue(), 12);

  v = true;
  EXPECT_TRUE(v.GetBoolValue());
  EXPECT_TRUE(v.GetValue<bool>());
}

TEST(VaryingTest, Test2)
{
  Varying v(12);
  EXPECT_EQ(v.GetIntegerValue(), 12);
  v = 12.0;
  EXPECT_DOUBLE_EQ(v.GetFloatValue(), 12.);
  EXPECT_DOUBLE_EQ(v.GetValue<double>(), 12.);
  EXPECT_DOUBLE_EQ(v.GetValue<float>(), 12.);
}

TEST(VaryingTest, Test3)
{
  Varying v(12.0);
  EXPECT_EQ(v.GetByteSize(), 8);
  EXPECT_DOUBLE_EQ(v.GetFloatValue(), 12.);
  v = 12;
  EXPECT_EQ(v.GetIntegerValue(), 12);
  EXPECT_EQ(v.GetValue<int>(), 12);
  EXPECT_EQ(v.GetValue<size_t>(), 12);
}

TEST(VaryingTest, Test4)
{
  Varying v(std::string("Hello"));
  EXPECT_EQ(v.GetStringValue(), "Hello");
  EXPECT_EQ(v.GetValue<std::string>(), "Hello");
}

TEST(VaryingTest, IntegerServesAsFloat)
{
  Varying v(-3);
  EXPECT_EQ(v.GetValue<double>(), -3.0);
  EXPECT_EQ(v.GetValue<float>(), -3.0F);
  // A floating point value does not serve as an integer, and a string does not serve as a float.
  EXPECT_THROW(Varying(2.5).GetValue<int>(), std::logic_error);
  EXPECT_THROW(Varying(std::string("1.0")).GetValue<double>(), std::logic_error);
}
