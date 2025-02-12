#include "test/unit/opensn_unit_test.h"
#include "framework/data_types/byte_array.h"
#include "framework/mesh/mesh_vector.h"
#include <gtest/gtest.h>

using namespace opensn;

class ByteArrayTest : public OpenSnUnitTest
{
};

TEST_F(ByteArrayTest, WriteRead)
{
  ByteArray barr;

  Vector3 v3a(3.0, 2.0, 1.0);
  barr.Write<double>(1.01234567890123456789);
  barr.Write<int>(-15600700);
  barr.Write<bool>(false);
  barr.Write<bool>(true);
  barr.Write<double>(v3a.x);
  barr.Write<double>(v3a.y);
  barr.Write<double>(v3a.z);

  EXPECT_DOUBLE_EQ(barr.Read<double>(), 1.01234567890123456789);
  EXPECT_EQ(barr.Read<int>(), -15600700);
  EXPECT_FALSE(barr.Read<bool>());
  EXPECT_TRUE(barr.Read<bool>());
  auto vec3_x = barr.Read<double>();
  auto vec3_y = barr.Read<double>();
  auto vec3_z = barr.Read<double>();
  Vector3 vec3b(vec3_x, vec3_y, vec3_z);
  EXPECT_DOUBLE_EQ(vec3_x, 3.0);
  EXPECT_DOUBLE_EQ(vec3_y, 2.0);
  EXPECT_DOUBLE_EQ(vec3_z, 1.0);
}

TEST_F(ByteArrayTest, Seek)
{
  ByteArray seeker;
  seeker.Write<bool>(false);
  seeker.Write<double>(1.01234567890123456789);

  EXPECT_FALSE(seeker.EndOfBuffer());
  EXPECT_EQ(seeker.Offset(), 0);
  seeker.Seek(seeker.Size() - sizeof(double));
  EXPECT_EQ(seeker.Offset(), sizeof(bool));
  EXPECT_DOUBLE_EQ(seeker.Read<double>(), 1.01234567890123456789);
}
