// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"
#include "framework/mesh/logical_volume/rcc_logical_volume.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/mesh/logical_volume/sphere_logical_volume.h"
#include "framework/parameters/parameter_block.h"
#include <limits>
#include <map>
#include <string>

using namespace opensn;

namespace
{

template <typename T, typename V>
T
MakeVolume(const std::map<std::string, V>& values)
{
  ParameterBlock block("");
  for (const auto& [name, value] : values)
    block.AddParameter(ParameterBlock(name, value));
  auto params = T::GetInputParameters();
  params.AssignParameters(block);
  return T(params);
}

/// Checks that every point of a 61 x 61 x 61 lattice on [-4, 4]^3 that is inside the volume is
/// inside its box, and that some points are inside.
void
ExpectBoxContainsVolume(const LogicalVolume& volume)
{
  Vector3 lo, hi;
  ASSERT_TRUE(volume.GetBoundingBox(lo, hi));
  constexpr int num_points = 61;
  const auto coordinate = [](const int n) { return -4.0 + 8.0 * n / (num_points - 1); };
  int num_inside = 0;
  for (int i = 0; i < num_points; ++i)
    for (int j = 0; j < num_points; ++j)
      for (int k = 0; k < num_points; ++k)
      {
        const Vector3 p(coordinate(i), coordinate(j), coordinate(k));
        if (not volume.Inside(p))
          continue;
        ++num_inside;
        EXPECT_TRUE(p.x >= lo.x and p.x <= hi.x and p.y >= lo.y and p.y <= hi.y and p.z >= lo.z and
                    p.z <= hi.z)
          << "point (" << p.x << ", " << p.y << ", " << p.z
          << ") is inside the volume but not its box";
      }
  EXPECT_GT(num_inside, 100);
}

} // namespace

TEST(LogicalVolumeTest, RPPBoundingBox)
{
  const auto rpp = MakeVolume<RPPLogicalVolume, double>(
    {{"xmin", -1.0}, {"xmax", 2.0}, {"ymin", 0.5}, {"ymax", 0.75}, {"zmin", -3.0}, {"zmax", 1.0}});
  Vector3 lo, hi;
  ASSERT_TRUE(rpp.GetBoundingBox(lo, hi));
  EXPECT_EQ(lo.x, -1.0);
  EXPECT_EQ(hi.x, 2.0);
  EXPECT_EQ(lo.y, 0.5);
  EXPECT_EQ(hi.y, 0.75);
  EXPECT_EQ(lo.z, -3.0);
  EXPECT_EQ(hi.z, 1.0);
  ExpectBoxContainsVolume(rpp);

  const auto slab = MakeVolume<RPPLogicalVolume, bool>({{"infx", true}, {"infz", true}});
  ASSERT_TRUE(slab.GetBoundingBox(lo, hi));
  constexpr double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(lo.x, -inf);
  EXPECT_EQ(hi.x, inf);
  EXPECT_EQ(lo.y, 0.0);
  EXPECT_EQ(hi.y, 1.0);
  EXPECT_EQ(lo.z, -inf);
  EXPECT_EQ(hi.z, inf);
}

TEST(LogicalVolumeTest, SphereBoundingBox)
{
  const auto sphere =
    MakeVolume<SphereLogicalVolume, double>({{"r", 1.5}, {"x", 0.3}, {"y", -0.7}, {"z", 1.1}});
  ExpectBoxContainsVolume(sphere);

  // Points on the surface are inside, so they must be in the box.
  Vector3 lo, hi;
  ASSERT_TRUE(sphere.GetBoundingBox(lo, hi));
  ASSERT_TRUE(sphere.Inside(Vector3(1.8, -0.7, 1.1)));
  EXPECT_LE(hi.x, 1.8 + 1.0e-10);
  EXPECT_GE(hi.x, 1.8);
  ASSERT_TRUE(sphere.Inside(Vector3(0.3, -0.7, -0.4)));
  EXPECT_LE(lo.z, -0.4);
}

TEST(LogicalVolumeTest, RCCBoundingBox)
{
  // Along +z, along -z, nearly along z, and oblique.
  const std::vector<std::map<std::string, double>> cylinders = {
    {{"r", 0.8}, {"x0", 0.1}, {"y0", 0.2}, {"z0", -1.0}, {"vx", 0.0}, {"vy", 0.0}, {"vz", 2.5}},
    {{"r", 0.8}, {"x0", 0.1}, {"y0", 0.2}, {"z0", 1.5}, {"vx", 0.0}, {"vy", 0.0}, {"vz", -2.5}},
    {{"r", 0.6}, {"x0", -0.5}, {"y0", 0.4}, {"z0", -2.0}, {"vx", 1.0e-7}, {"vy", 0.0}, {"vz", 3.0}},
    {{"r", 0.7}, {"x0", -1.0}, {"y0", 1.0}, {"z0", -0.5}, {"vx", 2.0}, {"vy", -1.5}, {"vz", 1.0}},
    {{"r", 1.0}, {"x0", 1.0}, {"y0", -1.0}, {"z0", 0.0}, {"vx", -2.0}, {"vy", 1.0}, {"vz", 0.0}}};
  for (const auto& values : cylinders)
    ExpectBoxContainsVolume(MakeVolume<RCCLogicalVolume, double>(values));
}
