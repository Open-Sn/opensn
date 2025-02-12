#include "test/unit/opensn_unit_test.h"
#include "framework/data_types/ndarray.h"
#include <gtest/gtest.h>

using namespace opensn;

class NDArrayTest : public OpenSnUnitTest
{
};

TEST_F(NDArrayTest, CtorVector)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST_F(NDArrayTest, CtorArray)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST_F(NDArrayTest, CtorInitList)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST_F(NDArrayTest, CtorVectorWithInit)
{
  NDArray<double, 3> arr({2, 2, 2}, 0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST_F(NDArrayTest, CtorArrayWithInit)
{
  NDArray<double, 3> arr({2, 2, 2}, 0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST_F(NDArrayTest, CtorInitListWithInit)
{
  NDArray<double, 3> arr({2, 2, 2}, 0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST_F(NDArrayTest, CtorDefault)
{
  NDArray<double, 2> arr;
  EXPECT_EQ(arr.size(), 0);
}

TEST_F(NDArrayTest, Rank)
{
  NDArray<double, 2> arr2;
  EXPECT_EQ(arr2.rank(), 2);

  NDArray<double, 4> arr4;
  EXPECT_EQ(arr4.rank(), 4);
}

TEST_F(NDArrayTest, Iterators)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(1.0);

  const auto& const_arr = arr;
  for (auto i = const_arr.cbegin(); i != const_arr.cend(); ++i)
    EXPECT_DOUBLE_EQ(*i, 1.0);
}

TEST_F(NDArrayTest, Size)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(1.0);

  EXPECT_EQ(arr.size(), 8);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 1.0);
}

TEST_F(NDArrayTest, Empty)
{
  NDArray<double, 3> arr;
  EXPECT_TRUE(arr.empty());
}

TEST_F(NDArrayTest, NotEmpty)
{
  NDArray<double, 3> arr({2, 2, 2});
  EXPECT_FALSE(arr.empty());
}
