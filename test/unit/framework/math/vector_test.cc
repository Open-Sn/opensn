#include "test/unit/opensn_unit_test.h"
#include "framework/math/vector.h"
#include <gtest/gtest.h>

using namespace opensn;

class VectorTest : public OpenSnUnitTest
{
};

TEST_F(VectorTest, Resize)
{
  Vector<double> a(3);
  EXPECT_EQ(a.Rows(), 3);
  a.Resize(4);
  EXPECT_EQ(a.Rows(), 4);
}

TEST_F(VectorTest, Set)
{
  Vector<double> a(3);
  a.Set(123.);
  EXPECT_DOUBLE_EQ(a(0), 123.);
  EXPECT_DOUBLE_EQ(a(1), 123.);
  EXPECT_DOUBLE_EQ(a(2), 123.);
}

TEST_F(VectorTest, Scale)
{
  Vector<double> a(3);
  a.Set(123.);
  a.Scale(0.1);
  EXPECT_DOUBLE_EQ(a(0), 12.3);
  EXPECT_DOUBLE_EQ(a(1), 12.3);
  EXPECT_DOUBLE_EQ(a(2), 12.3);
}

TEST_F(VectorTest, Normalized)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;

  auto b_normalized = b.Normalized();
  EXPECT_DOUBLE_EQ(b_normalized(0), 0.6);
  EXPECT_DOUBLE_EQ(b_normalized(1), 0.8);
}

TEST_F(VectorTest, Add2)
{
  Vector<double> b(2);
  b.Set(1.);

  Vector<double> added(2);
  added.Set(1.);
  added.Add(b);
  EXPECT_DOUBLE_EQ(added(0), 2.);
  EXPECT_DOUBLE_EQ(added(1), 2.);
}

TEST_F(VectorTest, Subtract2)
{
  Vector<double> b(2);
  b.Set(1.);

  Vector<double> subtracted(2);
  subtracted.Set(2.);
  subtracted.Subtract(b);
  EXPECT_DOUBLE_EQ(subtracted(0), 1.);
  EXPECT_DOUBLE_EQ(subtracted(1), 1.);
}

TEST_F(VectorTest, Dot)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;

  Vector<double> c(2);
  c(0) = 5;
  c(1) = 2;

  auto dp1 = b.Dot(c);
  EXPECT_DOUBLE_EQ(dp1, 23.);
}

// free functions

TEST_F(VectorTest, FVec2Norm)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;
  EXPECT_DOUBLE_EQ(Vec2Norm(b), 5.);
}

TEST_F(VectorTest, FAdd)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;

  Vector<double> c(2);
  c(0) = 5;
  c(1) = 2;

  auto d = Add(b, c);
  EXPECT_DOUBLE_EQ(d(0), 8.);
  EXPECT_DOUBLE_EQ(d(1), 6.);
}

TEST_F(VectorTest, FSubtract)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;

  Vector<double> c(2);
  c(0) = 5;
  c(1) = 2;

  auto e = Subtract(b, c);
  EXPECT_DOUBLE_EQ(e(0), -2.);
  EXPECT_DOUBLE_EQ(e(1), 2.);
}

TEST_F(VectorTest, FDot)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;

  Vector<double> c(2);
  c(0) = 5;
  c(1) = 2;

  auto dp0 = Dot(b, c);
  EXPECT_DOUBLE_EQ(dp0, 23.);
}

TEST_F(VectorTest, FScaled)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;

  auto b_scaled = Scaled(b, 4.);
  EXPECT_DOUBLE_EQ(b_scaled(0), 12.);
  EXPECT_DOUBLE_EQ(b_scaled(1), 16.);
}

TEST_F(VectorTest, FScale)
{
  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;

  Scale(b, 2.);
  EXPECT_DOUBLE_EQ(b(0), 6.);
  EXPECT_DOUBLE_EQ(b(1), 8.);
}
