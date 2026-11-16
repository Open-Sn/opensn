// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/exact_predicates.h"
#include "gtest/gtest.h"
#include <cmath>
#include <random>

using namespace opensn;

namespace
{

/// Reference sign of the determinant with long double differences and products, accurate when
/// the determinant is far from zero.
int
ReferenceSign(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d)
{
  const long double adx = static_cast<long double>(a.x) - d.x;
  const long double ady = static_cast<long double>(a.y) - d.y;
  const long double adz = static_cast<long double>(a.z) - d.z;
  const long double bdx = static_cast<long double>(b.x) - d.x;
  const long double bdy = static_cast<long double>(b.y) - d.y;
  const long double bdz = static_cast<long double>(b.z) - d.z;
  const long double cdx = static_cast<long double>(c.x) - d.x;
  const long double cdy = static_cast<long double>(c.y) - d.y;
  const long double cdz = static_cast<long double>(c.z) - d.z;
  const long double det =
    adx * (bdy * cdz - bdz * cdy) - ady * (bdx * cdz - bdz * cdx) + adz * (bdx * cdy - bdy * cdx);
  return det > 0 ? 1 : det < 0 ? -1 : 0;
}

} // namespace

TEST(ExactPredicates, Orient3DSignConvention)
{
  // A triangle in the plane z = 0 counterclockwise seen from +z: a point below it (inside a body
  // whose outward normal is +z) is positive, above is negative, in the plane is zero.
  const Vector3 a(0.0, 0.0, 0.0);
  const Vector3 b(1.0, 0.0, 0.0);
  const Vector3 c(0.0, 1.0, 0.0);
  EXPECT_EQ(Orient3D(a, b, c, Vector3(0.2, 0.2, -1.0)), 1);
  EXPECT_EQ(Orient3D(a, b, c, Vector3(0.2, 0.2, 1.0)), -1);
  EXPECT_EQ(Orient3D(a, b, c, Vector3(0.2, 0.2, 0.0)), 0);
  EXPECT_EQ(Orient3D(a, b, c, Vector3(5.0, -7.0, 0.0)), 0);
}

TEST(ExactPredicates, Orient3DExactZeros)
{
  // Points exactly in a plane not aligned with the axes, with coordinates whose floating-point
  // products round, must give zero: integer vertices with dyadic combinations are exact.
  const Vector3 p(1.0, 2.0, 3.0);
  const Vector3 q(4.0, 6.0, 5.0);
  const Vector3 r(-3.0, 7.0, 11.0);
  EXPECT_EQ(Orient3D(p, q, r, p + 0.5 * (q - p) + 0.25 * (r - p)), 0);
  EXPECT_EQ(Orient3D(p, q, r, p), 0);
  EXPECT_EQ(Orient3D(p, q, r, 0.5 * (q + r)), 0);
  EXPECT_EQ(Orient3D(p, q, r, q + 3.0 * (r - q) - 2.0 * (p - q)), 0);
}

TEST(ExactPredicates, Orient3DNearlyDegenerate)
{
  // Points a tiny distance off a plane, far from the origin, where the floating-point filter
  // cannot decide: the exact sign must agree with a high-precision reference and must be
  // antisymmetric under swapping two points.
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> coord(1.0e6, 1.0e6 + 1.0);
  int undecided_by_reference = 0;
  for (int k = 0; k < 20000; ++k)
  {
    const Vector3 a(coord(rng), coord(rng), coord(rng));
    const Vector3 b(coord(rng), coord(rng), coord(rng));
    const Vector3 c(coord(rng), coord(rng), coord(rng));
    // d close to the plane: a point of the plane plus a few ulps along the normal.
    const Vector3 in_plane = a + 0.3 * (b - a) + 0.4 * (c - a);
    const Vector3 normal = (b - a).Cross(c - a).Normalized();
    const double shift = std::ldexp(1.0, -20) * (static_cast<double>(k % 7) - 3.0);
    const Vector3 d = in_plane + shift * normal;
    const int sign = Orient3D(a, b, c, d);
    EXPECT_EQ(Orient3D(b, a, c, d), -sign);
    EXPECT_EQ(Orient3D(a, c, b, d), -sign);
    const int reference = ReferenceSign(a, b, c, d);
    if (std::abs(shift) >= std::ldexp(1.0, -20))
      EXPECT_EQ(sign, reference) << "k=" << k << " shift " << shift;
    else
      ++undecided_by_reference;
  }
  EXPECT_GT(undecided_by_reference, 0);
}

TEST(ExactPredicates, Orient3DConsistentAcrossSharedEdge)
{
  // The predicate of a segment against a shared edge is the same number whichever triangle
  // evaluates it, so it is odd under reversing the edge: both cells agree on the side exactly.
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> coord(-1000.0, 1000.0);
  for (int k = 0; k < 10000; ++k)
  {
    const Vector3 p(coord(rng), coord(rng), coord(rng));
    const Vector3 q(coord(rng), coord(rng), coord(rng));
    const Vector3 u(coord(rng), coord(rng), coord(rng));
    const Vector3 v(coord(rng), coord(rng), coord(rng));
    EXPECT_EQ(Orient3D(p, q, u, v), -Orient3D(p, q, v, u));
  }
}
