#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/math/quadratures/spatial/line_quadrature.h"
#include "framework/math/quadratures/spatial/triangle_quadrature.h"
#include "framework/math/quadratures/angular/lebedev_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/triangular_quadrature.h"
#include <gtest/gtest.h>
#include <numeric>

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

// Angular quadrature unit tests

TEST(QuadratureTest, LebedevQuadrature3DXYZ)
{
  // Order 3: 6-point set, scattering_order=1 -> 4 moments
  LebedevQuadrature3DXYZ quad(3, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 6);
  EXPECT_EQ(quad.GetDimension(), 3u);
  EXPECT_EQ(quad.GetNumMoments(), 4u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-10);

  // Check point-wise values
  const double w = 1.0 / 6.0;
  EXPECT_NEAR(quad.GetOmegas()[0].x, 1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].x, -1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].y, 1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].y, -1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[4].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[4].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[4].z, 1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[5].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[5].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[5].z, -1.0, 1e-12);

  // Check point-wise weights
  for (size_t i = 0; i < 6; ++i)
    EXPECT_NEAR(quad.GetWeights()[i], w, 1e-12);
}

TEST(QuadratureTest, LebedevQuadrature2DXY)
{
  // Order 3 upper hemisphere: 5 points, scattering_order=1 -> 3 moments
  LebedevQuadrature2DXY quad(3, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 5);
  EXPECT_EQ(quad.GetDimension(), 2u);
  EXPECT_EQ(quad.GetNumMoments(), 3u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-10);

  // Check point-wise values
  EXPECT_NEAR(quad.GetOmegas()[0].x, 1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].x, -1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].y, 1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].y, -1.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].z, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[4].x, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[4].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[4].z, 1.0, 1e-12);

  // Check point-wise weights
  const double w_eq = 1.0 / 6.0;
  const double w_pole = 1.0 / 3.0;
  EXPECT_NEAR(quad.GetWeights()[0], w_eq, 1e-12);
  EXPECT_NEAR(quad.GetWeights()[1], w_eq, 1e-12);
  EXPECT_NEAR(quad.GetWeights()[2], w_eq, 1e-12);
  EXPECT_NEAR(quad.GetWeights()[3], w_eq, 1e-12);
  EXPECT_NEAR(quad.GetWeights()[4], w_pole, 1e-12);
}

TEST(QuadratureTest, SLDFEsqQuadrature3DXYZ)
{
  // Level 0: 24 spherical quadrilaterals × 4 sub-points = 96 directions
  // scattering_order=1 -> 4 moments (3D)
  SLDFEsqQuadrature3DXYZ quad(0, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 96);
  EXPECT_EQ(quad.GetDimension(), 3u);
  EXPECT_EQ(quad.GetNumMoments(), 4u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-10);

  // Reference values for the first face (4 points)
  const double tol = 1e-10;
  const std::array<std::array<double, 4>, 4> ref = {{
    {0.0160210613009911, 0.9581267809141927, 0.2024760130176001, 0.2024760130176001},
    {0.0101028063271445, 0.7745966692417845, 0.6109051323703643, 0.1636915368706672},
    {0.0054399927113861, 0.6675548182224416, 0.5264838861104881, 0.5264838861104881},
    {0.0101028063271445, 0.7745966692417845, 0.1636915368706672, 0.6109051323703643},
  }};

  // Checks point-wise values and weights
  for (size_t i = 0; i < 4; ++i)
  {
    EXPECT_NEAR(quad.GetWeights()[i], ref[i][0], tol);
    EXPECT_NEAR(quad.GetOmegas()[i].x, ref[i][1], tol);
    EXPECT_NEAR(quad.GetOmegas()[i].y, ref[i][2], tol);
    EXPECT_NEAR(quad.GetOmegas()[i].z, ref[i][3], tol);
  }
}

TEST(QuadratureTest, GLProductQuadrature1DSlab)
{
  // n_polar=4: 4 directions along the slab axis (phi=0), scattering_order=1 -> 2 moments (1D)
  GLProductQuadrature1DSlab quad(4, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 4);
  EXPECT_EQ(quad.GetDimension(), 1u);
  EXPECT_EQ(quad.GetNumMoments(), 2u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-12);

  // Check point-wise values
  const double mu0 = 0.86113631159405;
  const double mu1 = 0.33998104358486;
  const double x0 = std::sqrt(1.0 - mu0 * mu0);
  const double x1 = std::sqrt(1.0 - mu1 * mu1);
  EXPECT_NEAR(quad.GetOmegas()[0].x, x0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].z, mu0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].x, x1, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].z, mu1, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].x, x1, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].z, -mu1, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].x, x0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].y, 0.0, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].z, -mu0, 1e-12);

  // Check point-wise weights
  const double w0 = 0.34785484513745 / 2.0;
  const double w1 = 0.65214515486255 / 2.0;
  EXPECT_NEAR(quad.GetWeights()[0], w0, 1e-12);
  EXPECT_NEAR(quad.GetWeights()[1], w1, 1e-12);
  EXPECT_NEAR(quad.GetWeights()[2], w1, 1e-12);
  EXPECT_NEAR(quad.GetWeights()[3], w0, 1e-12);
}

TEST(QuadratureTest, GLCProductQuadrature2DXY)
{
  // n_polar=2, n_azim=4: upper hemisphere only -> 1 polar × 4 azimuthal = 4 directions
  // scattering_order=1 -> 3 moments (2D)
  GLCProductQuadrature2DXY quad(2, 4, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 4);
  EXPECT_EQ(quad.GetDimension(), 2u);
  EXPECT_EQ(quad.GetNumMoments(), 3u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-12);

  // Check pointwise values
  const double v = 1.0 / std::sqrt(3.0);
  const double w = 0.25;
  EXPECT_NEAR(quad.GetOmegas()[0].x, v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].y, v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[0].z, v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].x, -v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].y, v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[1].z, v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].x, -v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].y, -v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[2].z, v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].x, v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].y, -v, 1e-12);
  EXPECT_NEAR(quad.GetOmegas()[3].z, v, 1e-12);

  // Check point-wise weights
  for (size_t i = 0; i < 4; ++i)
    EXPECT_NEAR(quad.GetWeights()[i], w, 1e-12);
}

TEST(QuadratureTest, GLCProductQuadrature3DXYZ)
{
  // n_polar=2, n_azim=4: full sphere -> 2 polar × 4 azimuthal = 8 directions
  // scattering_order=1 -> 4 moments (3D)
  GLCProductQuadrature3DXYZ quad(2, 4, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 8);
  EXPECT_EQ(quad.GetDimension(), 3u);
  EXPECT_EQ(quad.GetNumMoments(), 4u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-12);

  // Check point-wise values and weights
  const double v = 1.0 / std::sqrt(3.0);
  const double s = std::sqrt(2.0 / 3.0);
  const double w = 0.125;
  const std::array<double, 4> phi_vals = {
    M_PI / 4.0, 3.0 * M_PI / 4.0, 5.0 * M_PI / 4.0, 7.0 * M_PI / 4.0};
  for (size_t i = 0; i < 4; ++i)
  {
    const double cx = s * std::cos(phi_vals[i]);
    const double cy = s * std::sin(phi_vals[i]);
    EXPECT_NEAR(quad.GetOmegas()[2 * i + 0].x, cx, 1e-12);
    EXPECT_NEAR(quad.GetOmegas()[2 * i + 0].y, cy, 1e-12);
    EXPECT_NEAR(quad.GetOmegas()[2 * i + 0].z, v, 1e-12);
    EXPECT_NEAR(quad.GetOmegas()[2 * i + 1].x, cx, 1e-12);
    EXPECT_NEAR(quad.GetOmegas()[2 * i + 1].y, cy, 1e-12);
    EXPECT_NEAR(quad.GetOmegas()[2 * i + 1].z, -v, 1e-12);
    EXPECT_NEAR(quad.GetWeights()[2 * i + 0], w, 1e-12);
    EXPECT_NEAR(quad.GetWeights()[2 * i + 1], w, 1e-12);
  }
}

TEST(QuadratureTest, GLCTriangularQuadrature2DXY)
{
  // n_polar=4, upper hemisphere only: 4+8=12 directions
  // scattering_order=1 -> 3 moments (2D)
  GLCTriangularQuadrature2DXY quad(4, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 12);
  EXPECT_EQ(quad.GetDimension(), 2u);
  EXPECT_EQ(quad.GetNumMoments(), 3u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-12);

  // Reference values (W, X, Y, Z)
  const double tol = 5e-4;
  const std::array<std::array<double, 4>, 12> ref = {{
    {0.0870, 0.3595, 0.3595, 0.8611},
    {0.0870, -0.3595, 0.3595, 0.8611},
    {0.0870, -0.3595, -0.3595, 0.8611},
    {0.0870, 0.3595, -0.3595, 0.8611},
    {0.0815, 0.8688, 0.3599, 0.3400},
    {0.0815, 0.3599, 0.8688, 0.3400},
    {0.0815, -0.3599, 0.8688, 0.3400},
    {0.0815, -0.8688, 0.3599, 0.3400},
    {0.0815, -0.8688, -0.3599, 0.3400},
    {0.0815, -0.3599, -0.8688, 0.3400},
    {0.0815, 0.3599, -0.8688, 0.3400},
    {0.0815, 0.8688, -0.3599, 0.3400},
  }};

  // Check point-wise values and weights
  for (size_t i = 0; i < 12; ++i)
  {
    EXPECT_NEAR(quad.GetWeights()[i], ref[i][0], tol);
    EXPECT_NEAR(quad.GetOmegas()[i].x, ref[i][1], tol);
    EXPECT_NEAR(quad.GetOmegas()[i].y, ref[i][2], tol);
    EXPECT_NEAR(quad.GetOmegas()[i].z, ref[i][3], tol);
  }
}

TEST(QuadratureTest, GLCTriangularQuadrature3DXYZ)
{
  // n_polar=4, full sphere: 4+8+8+4=24 directions
  // scattering_order=1 -> 4 moments (3D)
  GLCTriangularQuadrature3DXYZ quad(4, 1);
  ASSERT_EQ(quad.GetOmegas().size(), 24);
  EXPECT_EQ(quad.GetDimension(), 3u);
  EXPECT_EQ(quad.GetNumMoments(), 4u);

  // Check weight sum = 1.0
  const double weight_sum =
    std::accumulate(quad.GetWeights().begin(), quad.GetWeights().end(), 0.0);
  EXPECT_NEAR(weight_sum, 1.0, 1e-12);

  // Reference values (W, X, Y, Z)
  const double tol = 5e-4;
  const std::array<std::array<double, 4>, 24> ref = {{
    {0.0435, 0.3595, 0.3595, 0.8611},    {0.0435, -0.3595, 0.3595, 0.8611},
    {0.0435, -0.3595, -0.3595, 0.8611},  {0.0435, 0.3595, -0.3595, 0.8611},
    {0.0408, 0.8688, 0.3599, 0.3400},    {0.0408, 0.3599, 0.8688, 0.3400},
    {0.0408, -0.3599, 0.8688, 0.3400},   {0.0408, -0.8688, 0.3599, 0.3400},
    {0.0408, -0.8688, -0.3599, 0.3400},  {0.0408, -0.3599, -0.8688, 0.3400},
    {0.0408, 0.3599, -0.8688, 0.3400},   {0.0408, 0.8688, -0.3599, 0.3400},
    {0.0408, 0.8688, 0.3599, -0.3400},   {0.0408, 0.3599, 0.8688, -0.3400},
    {0.0408, -0.3599, 0.8688, -0.3400},  {0.0408, -0.8688, 0.3599, -0.3400},
    {0.0408, -0.8688, -0.3599, -0.3400}, {0.0408, -0.3599, -0.8688, -0.3400},
    {0.0408, 0.3599, -0.8688, -0.3400},  {0.0408, 0.8688, -0.3599, -0.3400},
    {0.0435, 0.3595, 0.3595, -0.8611},   {0.0435, -0.3595, 0.3595, -0.8611},
    {0.0435, -0.3595, -0.3595, -0.8611}, {0.0435, 0.3595, -0.3595, -0.8611},
  }};

  // Check point-wise values and weights
  for (size_t i = 0; i < 24; ++i)
  {
    EXPECT_NEAR(quad.GetWeights()[i], ref[i][0], tol) << "weight mismatch at index " << i;
    EXPECT_NEAR(quad.GetOmegas()[i].x, ref[i][1], tol) << "omega.x mismatch at index " << i;
    EXPECT_NEAR(quad.GetOmegas()[i].y, ref[i][2], tol) << "omega.y mismatch at index " << i;
    EXPECT_NEAR(quad.GetOmegas()[i].z, ref[i][3], tol) << "omega.z mismatch at index " << i;
  }
}
