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
