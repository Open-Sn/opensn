#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/math/quadratures/spatial/line_quadrature.h"
#include "framework/math/quadratures/spatial/triangle_quadrature.h"
#include <gtest/gtest.h>

using namespace opensn;

TEST(QuadratureTest, GaussLegendre)
{
  GaussLegendreQuadrature quad(4);

  ASSERT_EQ(quad.qpoints.size(), 4);
  EXPECT_NEAR(quad.qpoints[0].x, -0.86113631159405, 1e-12);
  EXPECT_NEAR(quad.qpoints[0].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[0].z, 0., 1e-12);

  EXPECT_NEAR(quad.qpoints[1].x, -0.33998104358486, 1e-12);
  EXPECT_NEAR(quad.qpoints[1].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[1].z, 0., 1e-12);

  EXPECT_NEAR(quad.qpoints[2].x, 0.33998104358486, 1e-12);
  EXPECT_NEAR(quad.qpoints[2].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[2].z, 0., 1e-12);

  EXPECT_NEAR(quad.qpoints[3].x, 0.86113631159405, 1e-12);
  EXPECT_NEAR(quad.qpoints[3].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[3].z, 0., 1e-12);

  ASSERT_EQ(quad.weights.size(), 4);
  EXPECT_NEAR(quad.weights[0], 0.34785484513745, 1e-12);
  EXPECT_NEAR(quad.weights[1], 0.65214515486255, 1e-12);
  EXPECT_NEAR(quad.weights[2], 0.65214515486255, 1e-12);
  EXPECT_NEAR(quad.weights[3], 0.34785484513745, 1e-12);
}

TEST(QuadratureTest, GaussChebyshev)
{
  GaussChebyshevQuadrature quad(4);

  ASSERT_EQ(quad.qpoints.size(), 4);
  EXPECT_NEAR(quad.qpoints[0].x, -0.92387953251129, 1e-12);
  EXPECT_NEAR(quad.qpoints[0].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[0].z, 0., 1e-12);

  EXPECT_NEAR(quad.qpoints[1].x, -0.38268343236509, 1e-12);
  EXPECT_NEAR(quad.qpoints[1].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[1].z, 0., 1e-12);

  EXPECT_NEAR(quad.qpoints[2].x, 0.38268343236509, 1e-12);
  EXPECT_NEAR(quad.qpoints[2].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[2].z, 0., 1e-12);

  EXPECT_NEAR(quad.qpoints[3].x, 0.92387953251129, 1e-12);
  EXPECT_NEAR(quad.qpoints[3].y, 0., 1e-12);
  EXPECT_NEAR(quad.qpoints[3].z, 0., 1e-12);

  ASSERT_EQ(quad.weights.size(), 4);
  EXPECT_NEAR(quad.weights[0], 0.78539816339745, 1e-12);
  EXPECT_NEAR(quad.weights[1], 0.78539816339745, 1e-12);
  EXPECT_NEAR(quad.weights[2], 0.78539816339745, 1e-12);
  EXPECT_NEAR(quad.weights[3], 0.78539816339745, 1e-12);
}

TEST(QuadratureTest, Legendre)
{
  EXPECT_NEAR(Legendre(0, 0.25), 1., 1e-12);
  EXPECT_NEAR(Legendre(1, 0.25), 0.25, 1e-12);
}

TEST(QuadratureTest, dLegendredx)
{
  EXPECT_NEAR(dLegendredx(0, 0.25), 0., 1e-12);
  EXPECT_NEAR(dLegendredx(1, 0.25), 1., 1e-12);
}

TEST(QuadratureTest, Ylm)
{
  EXPECT_NEAR(Ylm(0, 0, 45 * M_PI / 180.0, 45 * M_PI / 180.0), 1.0, 1e-12);
  EXPECT_NEAR(Ylm(1, 0, 45 * M_PI / 180.0, 45 * M_PI / 180.0), 0.70710678118655, 1e-12);
}

TEST(QuadratureTest, GaussQuadratureSetRange)
{
  GaussLegendreQuadrature quad(QuadratureOrder::FIRST);
  quad.SetRange({2.0, 4.0});

  ASSERT_EQ(quad.qpoints.size(), 1);
  ASSERT_EQ(quad.weights.size(), 1);
  EXPECT_NEAR(quad.qpoints[0].x, 3.0, 1e-12);
  EXPECT_NEAR(quad.weights[0], 2.0, 1e-12);
}

TEST(QuadratureTest, LineQuadratureSetRange)
{
  LineQuadrature quad(QuadratureOrder::FIRST);

  ASSERT_EQ(quad.qpoints.size(), 1);
  ASSERT_EQ(quad.weights.size(), 1);
  EXPECT_NEAR(quad.qpoints[0].x, 0.5, 1e-12);
  EXPECT_NEAR(quad.weights[0], 1.0, 1e-12);

  quad.SetRange({-1.0, 1.0});
  EXPECT_NEAR(quad.qpoints[0].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.weights[0], 2.0, 1e-12);

  EXPECT_THROW(quad.SetRange({1.0, 1.0}), std::invalid_argument);
}

TEST(QuadratureTest, TriangleQuadrature)
{
  TriangleQuadrature q1(QuadratureOrder::FIRST);
  ASSERT_EQ(q1.qpoints.size(), 1);
  ASSERT_EQ(q1.weights.size(), 1);
  EXPECT_NEAR(q1.weights[0], 0.5, 1e-12);

  TriangleQuadrature q2(QuadratureOrder::SECOND);
  ASSERT_EQ(q2.qpoints.size(), 3);
  ASSERT_EQ(q2.weights.size(), 3);
  EXPECT_NEAR(q2.weights[0] + q2.weights[1] + q2.weights[2], 0.5, 1e-12);

  EXPECT_THROW(
    {
      TriangleQuadrature bad_quad(QuadratureOrder::FIFTH);
      (void)bad_quad;
    },
    std::invalid_argument);
}
