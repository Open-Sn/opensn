#include "test/unit/opensn_unit_test.h"
#include "framework/mesh/mesh_vector.h"

using namespace opensn;

class Vector3Test : public OpenSnUnitTest
{
};

void
ExpectVectorEqual(const Vector3& v, const double x, const double y = 0.0, const double z = 0.0)
{
  EXPECT_DOUBLE_EQ(v.x, x);
  EXPECT_DOUBLE_EQ(v.y, y);
  EXPECT_DOUBLE_EQ(v.z, z);
}

void
ExpectVectorEqual(const Vector3& a, const Vector3& b)
{
  EXPECT_DOUBLE_EQ(a.x, b.x);
  EXPECT_DOUBLE_EQ(a.y, b.y);
  EXPECT_DOUBLE_EQ(a.z, b.z);
}

TEST_F(Vector3Test, Constructors)
{
  ExpectVectorEqual(Vector3(), 0);
  ExpectVectorEqual(Vector3(1), 1);
  ExpectVectorEqual(Vector3(1, 2), 1, 2);
  ExpectVectorEqual(Vector3(1, 2, 3), 1, 2, 3);
  ExpectVectorEqual(Vector3({}), 0);
  ExpectVectorEqual(Vector3({4}), 4);
  ExpectVectorEqual(Vector3({4, 5}), 4, 5);
  ExpectVectorEqual(Vector3({4, 5, 6}), 4, 5, 6);
  ExpectVectorEqual(Vector3(std::vector<double>{}), 0);
  ExpectVectorEqual(Vector3(std::vector<double>{7}), 7);
  ExpectVectorEqual(Vector3(std::vector<double>{7, 8}), 7, 8);
  ExpectVectorEqual(Vector3(std::vector<double>{7, 8, 9}), 7, 8, 9);
}

TEST_F(Vector3Test, Plus)
{
  const Vector3 a(1, 2, 3);
  const Vector3 b(4, 5, 6);
  ExpectVectorEqual(a + b, a.x + b.x, a.y + b.y, a.z + b.z);
}

TEST_F(Vector3Test, PlusEquals)
{
  Vector3 a(1, 2, 3);
  const Vector3 b(4, 5, 6);
  const Vector3 a_copy = a;
  a += b;
  ExpectVectorEqual(a, a_copy.x + b.x, a_copy.y + b.y, a_copy.z + b.z);
}

TEST_F(Vector3Test, Shifted)
{
  const Vector3 a(1, 2, 3);
  const double val = 2;
  const auto b = a.Shifted(val);
  ExpectVectorEqual(b, a.x + val, a.y + val, a.z + val);
}

TEST_F(Vector3Test, Shift)
{
  Vector3 a(1, 2, 3);
  const auto a_copy = a;
  const double val = 5;
  a.Shift(val);
  ExpectVectorEqual(a, a_copy.x + val, a_copy.y + val, a_copy.z + val);
}

TEST_F(Vector3Test, Minus)
{
  const Vector3 a(1, 2, 3);
  const Vector3 b(4, 5, 6);
  ExpectVectorEqual(a - b, a.x - b.x, a.y - b.y, a.z - b.z);
}

TEST_F(Vector3Test, MinsEquals)
{
  Vector3 a(1, 2, 3);
  const Vector3 b(4, 5, 6);
  const Vector3 a_copy = a;
  a -= b;
  ExpectVectorEqual(a, a_copy.x - b.x, a_copy.y - b.y, a_copy.z - b.z);
}

TEST_F(Vector3Test, Times)
{
  const Vector3 a(1, 2, 3);
  const double val = 2;
  ExpectVectorEqual(a * val, a.x * val, a.y * val, a.z * val);
}

TEST_F(Vector3Test, TimesEquals)
{
  Vector3 a(1, 2, 3);
  const Vector3 a_copy = a;
  const double val = 2;
  a *= val;
  ExpectVectorEqual(a, a_copy.x * val, a_copy.y * val, a_copy.z * val);
}

TEST_F(Vector3Test, Div)
{
  const Vector3 a(1, 2, 3);
  const double val = 2;
  ExpectVectorEqual(a / val, a.x / val, a.y / val, a.z / val);
}

TEST_F(Vector3Test, DivEquals)
{
  Vector3 a(1, 2, 3);
  const Vector3 a_copy = a;
  const double val = 2;
  a /= val;
  ExpectVectorEqual(a, a_copy.x / val, a_copy.y / val, a_copy.z / val);
}

TEST_F(Vector3Test, DivComponent)
{
  const Vector3 a(1, 2, 3);
  const Vector3 b(6, 5, 4);
  ExpectVectorEqual(a / b, a.x / b.x, a.y / b.y, a.z / b.z);
}

TEST_F(Vector3Test, DivComponentEquals)
{
  Vector3 a(1, 2, 3);
  const auto a_copy = a;
  const Vector3 b(6, 5, 4);
  a /= b;
  ExpectVectorEqual(a, a_copy.x / b.x, a_copy.y / b.y, a_copy.z / b.z);
}

TEST_F(Vector3Test, Cross)
{
  const Vector3 a(1, 2, 3);
  const Vector3 b(6, 5, 4);
  ExpectVectorEqual(
    a.Cross(b), a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

TEST_F(Vector3Test, Dot)
{
  const Vector3 a(1, 2, 3);
  const Vector3 b(6, 5, 4);
  EXPECT_DOUBLE_EQ(a.Dot(b), a.x * b.x + a.y * b.y + a.z * b.z);
}

TEST_F(Vector3Test, Norm)
{
  const Vector3 a(1, 2, 3);
  EXPECT_DOUBLE_EQ(a.Norm(), std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z));
}

TEST_F(Vector3Test, NormSquare)
{
  const Vector3 a(1, 2, 3);
  EXPECT_DOUBLE_EQ(a.NormSquare(), a.x * a.x + a.y * a.y + a.z * a.z);
}

TEST_F(Vector3Test, Normalize)
{
  const Vector3 a(1, 2, 3);
  const auto norm = a.Norm();
  Vector3 a_copy = a;
  a_copy.Normalize();
  ExpectVectorEqual(a / norm, a_copy);
}

TEST_F(Vector3Test, Normalized)
{
  const Vector3 a(1, 2, 3);
  const auto norm = a.Norm();
  ExpectVectorEqual(a.Normalized(), a / norm);
}

TEST_F(Vector3Test, InverseZeroIfSmaller)
{
  const double tol = 1;
  const double val = 2 * -tol;
  for (size_t i = 0; i < Vector3::Size(); ++i)
  {
    const Vector3 a(i == 0 ? val : 0, i == 1 ? val : 0, i == 2 ? val : 0);
    Vector3 b;
    b(i) = 1.0 / val;
    ExpectVectorEqual(a.InverseZeroIfSmaller(tol), b);
  }
}

TEST_F(Vector3Test, InverseOneIfSmaller)
{
  const double tol = 10;
  const double val = 2 * tol;
  for (size_t i = 0; i < Vector3::Size(); ++i)
  {
    const Vector3 a(i == 0 ? val : 0, i == 1 ? val : 0, i == 2 ? val : 0);
    Vector3 b(1, 1, 1);
    b(i) = 1.0 / val;
    ExpectVectorEqual(a.InverseOneIfSmaller(tol), b);
  }
}

TEST_F(Vector3Test, Inverse)
{
  const Vector3 a(1, 2, 3);
  ExpectVectorEqual(a.Inverse(), 1 / a.x, 1 / a.y, 1 / a.z);

  const Vector3 b(0, 1, 2);
  EXPECT_THROW(
    {
      try
      {
        b.Inverse();
      }
      catch (const std::runtime_error& e)
      {
        EXPECT_STREQ("Division by zero in Vector3::Inverse.", e.what());
        throw;
      }
    },
    std::runtime_error);
}

TEST_F(Vector3Test, Size)
{
  EXPECT_EQ(Vector3::Size(), std::size_t(3));
}
