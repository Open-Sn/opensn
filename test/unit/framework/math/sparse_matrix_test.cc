#include "test/unit/opensn_unit_test.h"
#include "framework/math/sparse_matrix/sparse_matrix.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace opensn;

class SparseMatrixTest : public OpenSnUnitTest
{
};

TEST_F(SparseMatrixTest, Insert)
{
  SparseMatrix matrix(4, 4);
  auto& mat = matrix;
  mat.Insert(0, 0, 1.0);
  mat.Insert(0, 1, 1.1);
  mat.Insert(0, 2, 1.2);
  mat.Insert(0, 3, 1.3);
  mat.Insert(1, 0, 1.9);
  mat.Insert(1, 1, 2.0);
  mat.Insert(1, 2, 2.1);
  mat.Insert(2, 1, 2.9);
  mat.Insert(2, 2, 3.0);
  mat.Insert(2, 3, 3.1);
  mat.Insert(3, 2, 3.9);
  mat.Insert(3, 3, 4.0);

  EXPECT_DOUBLE_EQ(mat.GetValueIJ(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(0, 1), 1.1);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(0, 2), 1.2);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(0, 3), 1.3);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(1, 0), 1.9);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(1, 2), 2.1);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(2, 1), 2.9);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(2, 2), 3.0);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(2, 3), 3.1);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(3, 2), 3.9);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(3, 3), 4.0);
}

TEST_F(SparseMatrixTest, Row)
{
  SparseMatrix m(3, 3);
  m.Insert(0, 0, 1.);
  m.Insert(0, 1, 2.);
  m.Insert(0, 2, 3.);
  m.Insert(1, 0, 4.);
  m.Insert(1, 1, 5.);
  m.Insert(1, 2, 6.);
  m.Insert(2, 0, 7.);
  m.Insert(2, 1, 8.);
  m.Insert(2, 2, 9.);

  std::vector<double> vals;
  for (const auto& entry : m.Row(2))
    vals.push_back(entry.value);
  EXPECT_THAT(vals, ::testing::ElementsAre(7., 8., 9.));

  for (const auto& [row_index, column_index, value] : m.Row(2))
    value *= 2;

  vals.clear();
  std::vector<size_t> rows;
  std::vector<size_t> cols;
  for (const auto& entry : m.Row(2))
  {
    rows.push_back(entry.row_index);
    cols.push_back(entry.column_index);
    vals.push_back(entry.value);
  }
  EXPECT_THAT(rows, ::testing::ElementsAre(2, 2, 2));
  EXPECT_THAT(cols, ::testing::ElementsAre(0, 1, 2));
  EXPECT_THAT(vals, ::testing::ElementsAre(14., 16., 18.));
}

TEST_F(SparseMatrixTest, Compress)
{
  SparseMatrix mat(4, 4);
  mat.Insert(0, 0, 1.0);
  mat.Insert(0, 1, 1.1);
  mat.Insert(0, 2, 1.2);
  mat.Insert(0, 3, 1.3);
  mat.Insert(1, 0, 1.9);
  mat.Insert(1, 1, 2.0);
  mat.Insert(1, 2, 2.1);
  mat.Insert(2, 1, 2.9);
  mat.Insert(2, 2, 3.0);
  mat.Insert(2, 3, 3.1);
  mat.Insert(3, 2, 3.9);
  mat.Insert(3, 3, 4.0);

  std::vector<size_t> rows;
  std::vector<size_t> cols;
  std::vector<double> vals;
  for (const auto& entry : mat)
  {
    rows.push_back(entry.row_index);
    cols.push_back(entry.column_index);
    vals.push_back(entry.value);
  }
  EXPECT_THAT(rows, ::testing::ElementsAre(0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3));
  EXPECT_THAT(cols, ::testing::ElementsAre(0, 1, 2, 3, 0, 1, 2, 1, 2, 3, 2, 3));
  EXPECT_THAT(vals,
              ::testing::ElementsAre(1.0, 1.1, 1.2, 1.3, 1.9, 2.0, 2.1, 2.9, 3.0, 3.1, 3.9, 4.0));
}

TEST_F(SparseMatrixTest, RowsCols)
{
  SparseMatrix mat(3, 4);
  EXPECT_EQ(mat.GetNumRows(), 3);
  EXPECT_EQ(mat.GetNumCols(), 4);
}

TEST_F(SparseMatrixTest, InsertAdd)
{
  SparseMatrix mat(3, 4);
  mat.Insert(0, 0, 2.);
  mat.InsertAdd(0, 0, 5.);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(0, 0), 7.);
}

TEST_F(SparseMatrixTest, SetDiagonal)
{
  SparseMatrix mat(3, 3);
  mat.SetDiagonal({1., 2., 3.});
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(0, 0), 1.);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(1, 1), 2.);
  EXPECT_DOUBLE_EQ(mat.GetValueIJ(2, 2), 3.);
}
