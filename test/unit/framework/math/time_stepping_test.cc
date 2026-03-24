#include "test/unit/opensn_unit_test.h"
#include "framework/math/math_time_stepping.h"
#include <gtest/gtest.h>

using namespace opensn;

class TimeSteppingTest : public OpenSnUnitTest
{
};

TEST_F(TimeSteppingTest, SteppingMethodStringName)
{
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::NONE), "none");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::EXPLICIT_EULER), "explicit_euler");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::IMPLICIT_EULER), "implicit_euler");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::CRANK_NICOLSON), "crank_nicholson");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::THETA_SCHEME), "theta_scheme");

  EXPECT_THROW(SteppingMethodStringName(static_cast<SteppingMethod>(999)), std::logic_error);
}
