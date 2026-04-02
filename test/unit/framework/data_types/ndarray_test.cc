#include "framework/data_types/ndarray.h"
#include <gtest/gtest.h>

using namespace opensn;

TEST(NDArrayTest, CtorVector)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST(NDArrayTest, CtorArray)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST(NDArrayTest, CtorInitList)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST(NDArrayTest, CtorVectorWithInit)
{
  NDArray<double, 3> arr({2, 2, 2}, 0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST(NDArrayTest, CtorArrayWithInit)
{
  NDArray<double, 3> arr({2, 2, 2}, 0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST(NDArrayTest, CtorInitListWithInit)
{
  NDArray<double, 3> arr({2, 2, 2}, 0.0);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 0.);
}

TEST(NDArrayTest, CtorDefault)
{
  NDArray<double, 2> arr;
  EXPECT_EQ(arr.size(), 0);
}

TEST(NDArrayTest, Rank)
{
  NDArray<double, 2> arr2;
  EXPECT_EQ(arr2.rank(), 2);

  NDArray<double, 4> arr4;
  EXPECT_EQ(arr4.rank(), 4);
}

TEST(NDArrayTest, Iterators)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(1.0);

  const auto& const_arr = arr;
  for (auto i = const_arr.cbegin(); i != const_arr.cend(); ++i)
    EXPECT_DOUBLE_EQ(*i, 1.0);
}

TEST(NDArrayTest, Size)
{
  NDArray<double, 3> arr({2, 2, 2});
  arr.set(1.0);

  EXPECT_EQ(arr.size(), 8);
  for (auto val : arr)
    EXPECT_DOUBLE_EQ(val, 1.0);
}

TEST(NDArrayTest, Empty)
{
  NDArray<double, 3> arr;
  EXPECT_TRUE(arr.empty());
}

TEST(NDArrayTest, NotEmpty)
{
  NDArray<double, 3> arr({2, 2, 2});
  EXPECT_FALSE(arr.empty());
}
