// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/reflecting_boundary.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace opensn;

namespace
{

// Two 2D directions that are mirror images across a constant-x plane, with given weights.
class MirrorPairQuadrature : public AngularQuadrature
{
public:
  MirrorPairQuadrature(double weight, double mirror_weight)
    : AngularQuadrature(AngularQuadratureType::PRODUCT_QUADRATURE, 2, 0)
  {
    omegas_ = {Vector3(0.6, 0.0, 0.8), Vector3(-0.6, 0.0, 0.8)};
    weights_ = {weight, mirror_weight};
  }

  std::string GetName() const override { return "MirrorPairQuadrature"; }
};

// The weight tolerance that MapReflectedDirections accepts.
double
WeightTolerance(double weight, double mirror_weight)
{
  return 64.0 * std::numeric_limits<double>::epsilon() *
         std::max({1.0, std::fabs(weight), std::fabs(mirror_weight)});
}

} // namespace

TEST(ReflectingBoundaryTest, MapsDirectionsToTheirMirrorImages)
{
  const GLCProductQuadrature3DXYZ quadrature(4, 8, 0);
  const std::vector<Vector3> normals = {Vector3(-1.0, 0.0, 0.0),
                                        Vector3(1.0, 0.0, 0.0),
                                        Vector3(0.0, -1.0, 0.0),
                                        Vector3(0.0, 1.0, 0.0),
                                        Vector3(0.0, 0.0, -1.0),
                                        Vector3(0.0, 0.0, 1.0)};
  for (const auto& normal : normals)
  {
    const auto reflected = ReflectingBoundary::MapReflectedDirections(
      quadrature, normal, CoordinateSystemType::CARTESIAN);
    ASSERT_EQ(reflected.size(), quadrature.GetNumAngles());
    for (size_t n = 0; n < reflected.size(); ++n)
    {
      const auto& omega = quadrature.GetOmega(n);
      const auto& mirror = quadrature.GetOmega(reflected[n]);
      const auto expected = omega - 2.0 * normal * omega.Dot(normal);
      EXPECT_NEAR((mirror - expected).Norm(), 0.0, 1.0e-12);
      EXPECT_EQ(reflected[reflected[n]], n);
      EXPECT_NEAR(quadrature.GetWeight(reflected[n]),
                  quadrature.GetWeight(n),
                  WeightTolerance(quadrature.GetWeight(n), quadrature.GetWeight(reflected[n])));
    }
  }
}

TEST(ReflectingBoundaryTest, MapsRZDirectionsOnTheAxisAndConstantZPlanes)
{
  const GLCProductQuadrature2DRZ quadrature(4, 8, 0);
  // The symmetry axis (rmin) and the constant-z planes (zmin, zmax).
  const std::vector<Vector3> normals = {
    Vector3(-1.0, 0.0, 0.0), Vector3(0.0, -1.0, 0.0), Vector3(0.0, 1.0, 0.0)};
  for (const auto& normal : normals)
  {
    const auto reflected = ReflectingBoundary::MapReflectedDirections(
      quadrature, normal, CoordinateSystemType::CYLINDRICAL);
    ASSERT_EQ(reflected.size(), quadrature.GetNumAngles());
    for (size_t n = 0; n < reflected.size(); ++n)
    {
      const auto& omega = quadrature.GetOmega(n);
      const auto expected = omega - 2.0 * normal * omega.Dot(normal);
      EXPECT_NEAR((quadrature.GetOmega(reflected[n]) - expected).Norm(), 0.0, 1.0e-12);
      EXPECT_EQ(reflected[reflected[n]], n);
    }
  }
}

TEST(ReflectingBoundaryTest, RejectsRZOuterRadialBoundary)
{
  const GLCProductQuadrature2DRZ quadrature(4, 8, 0);
  EXPECT_THROW(ReflectingBoundary::MapReflectedDirections(
                 quadrature, Vector3(1.0, 0.0, 0.0), CoordinateSystemType::CYLINDRICAL),
               std::invalid_argument);
}

TEST(ReflectingBoundaryTest, MapsSphericalMuToMinusMu)
{
  const GLProductQuadrature1DSpherical quadrature(8, 0);
  const auto reflected = ReflectingBoundary::MapReflectedDirections(
    quadrature, Vector3(0.0, 0.0, 1.0), CoordinateSystemType::SPHERICAL);
  ASSERT_EQ(reflected.size(), quadrature.GetNumAngles());
  for (size_t n = 0; n < reflected.size(); ++n)
  {
    const auto& omega = quadrature.GetOmega(n);
    const auto& mirror = quadrature.GetOmega(reflected[n]);
    EXPECT_NE(reflected[n], n);
    EXPECT_NEAR(mirror.z, -omega.z, 1.0e-12);
    EXPECT_NEAR(mirror.x, omega.x, 1.0e-12);
    EXPECT_EQ(reflected[reflected[n]], n);
  }
}

TEST(ReflectingBoundaryTest, RejectsMirrorDirectionsWithDifferentWeights)
{
  const Vector3 normal(-1.0, 0.0, 0.0);

  // Equal weights: the pair maps onto itself.
  const auto reflected = ReflectingBoundary::MapReflectedDirections(
    MirrorPairQuadrature(0.5, 0.5), normal, CoordinateSystemType::CARTESIAN);
  ASSERT_EQ(reflected.size(), 2U);
  EXPECT_EQ(reflected[0], 1U);
  EXPECT_EQ(reflected[1], 0U);

  // Different weights: the boundary would not conserve the reflected current.
  EXPECT_THROW(ReflectingBoundary::MapReflectedDirections(
                 MirrorPairQuadrature(0.5, 0.75), normal, CoordinateSystemType::CARTESIAN),
               std::logic_error);
}

TEST(ReflectingBoundaryTest, RejectsQuadratureWithoutMirrorDirection)
{
  // The mirror of (0.6, 0, 0.8) across a constant-z plane, (0.6, 0, -0.8), is not in the set.
  EXPECT_THROW(ReflectingBoundary::MapReflectedDirections(MirrorPairQuadrature(0.5, 0.5),
                                                          Vector3(0.0, 0.0, 1.0),
                                                          CoordinateSystemType::CARTESIAN),
               std::logic_error);
}
