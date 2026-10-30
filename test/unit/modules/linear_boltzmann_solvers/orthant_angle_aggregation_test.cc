// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_runtime_builder.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "gtest/gtest.h"
#include <array>
#include <set>

using namespace opensn;

namespace
{

void
ExpectOrthantGroups(const AngularQuadrature& quadrature,
                    const unsigned int dimension,
                    const size_t expected_groups)
{
  const auto groups =
    detail::GroupDirectionsByOrthant("OrthantAngleAggregationTest", quadrature, dimension);
  ASSERT_EQ(groups.size(), expected_groups);

  std::set<size_t> visited;
  for (const auto& group : groups)
  {
    ASSERT_FALSE(group.empty());
    const auto& reference = quadrature.GetOmega(group.front());
    const std::array<double, 3> reference_components = {reference.x, reference.y, reference.z};
    for (const size_t direction : group)
    {
      EXPECT_TRUE(visited.insert(direction).second);
      const auto& omega = quadrature.GetOmega(direction);
      const std::array<double, 3> components = {omega.x, omega.y, omega.z};
      for (unsigned int d = 0; d < dimension; ++d)
        EXPECT_EQ(components[d] > 0.0, reference_components[d] > 0.0);
    }
  }
  EXPECT_EQ(visited.size(), quadrature.GetNumAngles());
}

class AxisDirectionQuadrature : public AngularQuadrature
{
public:
  AxisDirectionQuadrature() : AngularQuadrature(AngularQuadratureType::PRODUCT_QUADRATURE, 2, 0)
  {
    omegas_ = {{1.0, 0.0, 0.0}};
    weights_ = {1.0};
  }

  std::string GetName() const override { return "AxisDirectionQuadrature"; }
};

} // namespace

TEST(OrthantAngleAggregationTest, Groups2DDirectionsIntoQuadrants)
{
  const GLCProductQuadrature2DXY quadrature(4, 8, 0);
  ExpectOrthantGroups(quadrature, 2, 4);
}

TEST(OrthantAngleAggregationTest, Groups3DDirectionsIntoOctants)
{
  const GLCProductQuadrature3DXYZ quadrature(4, 8, 0);
  ExpectOrthantGroups(quadrature, 3, 8);
}

TEST(OrthantAngleAggregationTest, RejectsAxisDirections)
{
  const AxisDirectionQuadrature quadrature;
  EXPECT_THROW(detail::GroupDirectionsByOrthant("OrthantAngleAggregationTest", quadrature, 2),
               std::logic_error);
}

TEST(OrthantAngleAggregationTest, OrthantGroupSPDSUsesMeanDirection)
{
  const GLCProductQuadrature3DXYZ quadrature(4, 8, 0);
  for (const auto& group : detail::GroupDirectionsByOrthant("Test", quadrature, 3))
  {
    Vector3 sum;
    for (const size_t direction : group)
      sum += quadrature.GetOmega(direction);
    const auto expected = sum.Normalized();
    const auto omega = detail::SweepOrderingOmega(quadrature, group, AngleAggregationType::OCTANT);
    EXPECT_NEAR(omega.x, expected.x, 1.0e-15);
    EXPECT_NEAR(omega.y, expected.y, 1.0e-15);
    EXPECT_NEAR(omega.z, expected.z, 1.0e-15);
    EXPECT_NEAR(omega.Norm(), 1.0, 1.0e-15);
    // Every member streams through each coordinate plane in the same direction as the mean.
    for (const size_t direction : group)
    {
      const auto& member = quadrature.GetOmega(direction);
      EXPECT_GT(member.x * omega.x, 0.0);
      EXPECT_GT(member.y * omega.y, 0.0);
      EXPECT_GT(member.z * omega.z, 0.0);
    }
  }
}

TEST(OrthantAngleAggregationTest, OtherGroupsUseFirstDirection)
{
  const GLCProductQuadrature3DXYZ quadrature(4, 8, 0);
  const DirIDs single = {3};
  const auto omega = detail::SweepOrderingOmega(quadrature, single, AngleAggregationType::POLAR);
  EXPECT_EQ(omega.x, quadrature.GetOmega(3).x);
  EXPECT_EQ(omega.y, quadrature.GetOmega(3).y);
  EXPECT_EQ(omega.z, quadrature.GetOmega(3).z);

  // Azimuthal (curvilinear) groups keep their first direction.
  const DirIDs pair = {0, 1};
  const auto azimuthal =
    detail::SweepOrderingOmega(quadrature, pair, AngleAggregationType::AZIMUTHAL);
  EXPECT_EQ(azimuthal.x, quadrature.GetOmega(0).x);
  EXPECT_EQ(azimuthal.y, quadrature.GetOmega(0).y);
  EXPECT_EQ(azimuthal.z, quadrature.GetOmega(0).z);
}
