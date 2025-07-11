#include "test/unit/opensn_unit_test.h"
#include "framework/data_types/dense_matrix.h"
#include <gtest/gtest.h>

using namespace opensn;

class DenseMatrixTest : public OpenSnUnitTest
{
};

TEST_F(DenseMatrixTest, ColRows)
{
  DenseMatrix<double> a(2, 4);

  EXPECT_EQ(a.Rows(), 2);
  EXPECT_EQ(a.Columns(), 4);
}

TEST_F(DenseMatrixTest, OpCall)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;

  EXPECT_DOUBLE_EQ(a(0, 0), 1);
  EXPECT_DOUBLE_EQ(a(0, 1), 2);
  EXPECT_DOUBLE_EQ(a(0, 2), 3);
  EXPECT_DOUBLE_EQ(a(0, 3), 4);
  EXPECT_DOUBLE_EQ(a(1, 0), -1);
  EXPECT_DOUBLE_EQ(a(1, 1), -2);
  EXPECT_DOUBLE_EQ(a(1, 2), -3);
  EXPECT_DOUBLE_EQ(a(1, 3), -4);
}

TEST_F(DenseMatrixTest, SetRow)
{
  DenseMatrix<double> a(2, 4);
  Vector<double> r0(4);
  r0(0) = 10;
  r0(1) = 11;
  r0(2) = 12;
  r0(3) = 13;
  a.SetRow(0, r0);

  EXPECT_DOUBLE_EQ(a(0, 0), 10);
  EXPECT_DOUBLE_EQ(a(0, 1), 11);
  EXPECT_DOUBLE_EQ(a(0, 2), 12);
  EXPECT_DOUBLE_EQ(a(0, 3), 13);
}

TEST_F(DenseMatrixTest, FSwapRows)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  SwapRows(a, 0, 1);

  EXPECT_DOUBLE_EQ(a(0, 0), -1);
  EXPECT_DOUBLE_EQ(a(0, 1), -2);
  EXPECT_DOUBLE_EQ(a(0, 2), -3);
  EXPECT_DOUBLE_EQ(a(0, 3), -4);
  EXPECT_DOUBLE_EQ(a(1, 0), 1);
  EXPECT_DOUBLE_EQ(a(1, 1), 2);
  EXPECT_DOUBLE_EQ(a(1, 2), 3);
  EXPECT_DOUBLE_EQ(a(1, 3), 4);
}

TEST_F(DenseMatrixTest, SetDiagonal)
{
  DenseMatrix<double> diag(3, 3, 0.);
  diag.SetDiagonal(12.);

  EXPECT_DOUBLE_EQ(diag(0, 0), 12.);
  EXPECT_DOUBLE_EQ(diag(1, 1), 12.);
  EXPECT_DOUBLE_EQ(diag(2, 2), 12.);
}

TEST_F(DenseMatrixTest, FTranspose)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_trans = Transpose(a);

  EXPECT_DOUBLE_EQ(a_trans(0, 0), 1);
  EXPECT_DOUBLE_EQ(a_trans(1, 0), 2);
  EXPECT_DOUBLE_EQ(a_trans(2, 0), 3);
  EXPECT_DOUBLE_EQ(a_trans(3, 0), 4);
  EXPECT_DOUBLE_EQ(a_trans(0, 1), -1);
  EXPECT_DOUBLE_EQ(a_trans(1, 1), -2);
  EXPECT_DOUBLE_EQ(a_trans(2, 1), -3);
  EXPECT_DOUBLE_EQ(a_trans(3, 1), -4);
}

TEST_F(DenseMatrixTest, FMultScalar)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_mult = Mult(a, 2.);

  EXPECT_DOUBLE_EQ(a_mult(0, 0), 2);
  EXPECT_DOUBLE_EQ(a_mult(0, 1), 4);
  EXPECT_DOUBLE_EQ(a_mult(0, 2), 6);
  EXPECT_DOUBLE_EQ(a_mult(0, 3), 8);
  EXPECT_DOUBLE_EQ(a_mult(1, 0), -2);
  EXPECT_DOUBLE_EQ(a_mult(1, 1), -4);
  EXPECT_DOUBLE_EQ(a_mult(1, 2), -6);
  EXPECT_DOUBLE_EQ(a_mult(1, 3), -8);
}

TEST_F(DenseMatrixTest, FMultVec)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;

  Vector<double> v0(4);
  v0(0) = 1;
  v0(1) = -1;
  v0(2) = 3;
  v0(3) = -2;

  auto a_mult_vec = Mult(a, v0);

  EXPECT_DOUBLE_EQ(a_mult_vec(0), 0);
  EXPECT_DOUBLE_EQ(a_mult_vec(1), 0);
}

TEST_F(DenseMatrixTest, FMultMat)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;

  DenseMatrix<double> b(4, 3);
  b(0, 0) = -1;
  b(0, 1) = 0;
  b(0, 2) = 2;
  b(1, 0) = -4;
  b(1, 1) = 0;
  b(1, 2) = -2;
  b(2, 0) = 3;
  b(2, 1) = 5;
  b(2, 2) = 8;
  b(3, 0) = 1;
  b(3, 1) = -1;
  b(3, 2) = 0;

  auto ab = Mult(a, b);

  EXPECT_DOUBLE_EQ(ab(0, 0), 4);
  EXPECT_DOUBLE_EQ(ab(0, 1), 11);
  EXPECT_DOUBLE_EQ(ab(0, 2), 22);
  EXPECT_DOUBLE_EQ(ab(1, 0), -4);
  EXPECT_DOUBLE_EQ(ab(1, 1), -11);
  EXPECT_DOUBLE_EQ(ab(1, 2), -22);
}

TEST_F(DenseMatrixTest, FAdd)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_mult = Mult(a, 2.);
  auto apb = Add(a, a_mult);

  EXPECT_DOUBLE_EQ(apb(0, 0), 3);
  EXPECT_DOUBLE_EQ(apb(0, 1), 6);
  EXPECT_DOUBLE_EQ(apb(0, 2), 9);
  EXPECT_DOUBLE_EQ(apb(0, 3), 12);
  EXPECT_DOUBLE_EQ(apb(1, 0), -3);
  EXPECT_DOUBLE_EQ(apb(1, 1), -6);
  EXPECT_DOUBLE_EQ(apb(1, 2), -9);
  EXPECT_DOUBLE_EQ(apb(1, 3), -12);
}

TEST_F(DenseMatrixTest, FSubtract)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_mult = Mult(a, 2.);
  auto amb = Subtract(a, a_mult);

  EXPECT_DOUBLE_EQ(amb(0, 0), -1);
  EXPECT_DOUBLE_EQ(amb(0, 1), -2);
  EXPECT_DOUBLE_EQ(amb(0, 2), -3);
  EXPECT_DOUBLE_EQ(amb(0, 3), -4);
  EXPECT_DOUBLE_EQ(amb(1, 0), 1);
  EXPECT_DOUBLE_EQ(amb(1, 1), 2);
  EXPECT_DOUBLE_EQ(amb(1, 2), 3);
  EXPECT_DOUBLE_EQ(amb(1, 3), 4);
}

TEST_F(DenseMatrixTest, FScale)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_scaled = a;
  Scale(a_scaled, 0.5);

  EXPECT_DOUBLE_EQ(a_scaled(0, 0), 0.5);
  EXPECT_DOUBLE_EQ(a_scaled(0, 1), 1);
  EXPECT_DOUBLE_EQ(a_scaled(0, 2), 1.5);
  EXPECT_DOUBLE_EQ(a_scaled(0, 3), 2);
  EXPECT_DOUBLE_EQ(a_scaled(1, 0), -0.5);
  EXPECT_DOUBLE_EQ(a_scaled(1, 1), -1);
  EXPECT_DOUBLE_EQ(a_scaled(1, 2), -1.5);
  EXPECT_DOUBLE_EQ(a_scaled(1, 3), -2);
}

TEST_F(DenseMatrixTest, FScaled)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_scaled = Scaled(a, 0.5);

  EXPECT_DOUBLE_EQ(a_scaled(0, 0), 0.5);
  EXPECT_DOUBLE_EQ(a_scaled(0, 1), 1);
  EXPECT_DOUBLE_EQ(a_scaled(0, 2), 1.5);
  EXPECT_DOUBLE_EQ(a_scaled(0, 3), 2);
  EXPECT_DOUBLE_EQ(a_scaled(1, 0), -0.5);
  EXPECT_DOUBLE_EQ(a_scaled(1, 1), -1);
  EXPECT_DOUBLE_EQ(a_scaled(1, 2), -1.5);
  EXPECT_DOUBLE_EQ(a_scaled(1, 3), -2);
}

TEST_F(DenseMatrixTest, FDeterminant)
{
  DenseMatrix<double> diag(3, 3, 0.);
  diag.SetDiagonal(12.);
  auto d0 = Determinant(diag);

  EXPECT_DOUBLE_EQ(d0, 12. * 12. * 12.);
}

TEST_F(DenseMatrixTest, Add)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_mult = Mult(a, 2.);
  a.Add(a_mult);

  EXPECT_DOUBLE_EQ(a(0, 0), 3);
  EXPECT_DOUBLE_EQ(a(0, 1), 6);
  EXPECT_DOUBLE_EQ(a(0, 2), 9);
  EXPECT_DOUBLE_EQ(a(0, 3), 12);
  EXPECT_DOUBLE_EQ(a(1, 0), -3);
  EXPECT_DOUBLE_EQ(a(1, 1), -6);
  EXPECT_DOUBLE_EQ(a(1, 2), -9);
  EXPECT_DOUBLE_EQ(a(1, 3), -12);
}

TEST_F(DenseMatrixTest, Subtract)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_mult = Mult(a, 2.);
  a.Subtract(a_mult);

  EXPECT_DOUBLE_EQ(a(0, 0), -1);
  EXPECT_DOUBLE_EQ(a(0, 1), -2);
  EXPECT_DOUBLE_EQ(a(0, 2), -3);
  EXPECT_DOUBLE_EQ(a(0, 3), -4);
  EXPECT_DOUBLE_EQ(a(1, 0), 1);
  EXPECT_DOUBLE_EQ(a(1, 1), 2);
  EXPECT_DOUBLE_EQ(a(1, 2), 3);
  EXPECT_DOUBLE_EQ(a(1, 3), 4);
}

TEST_F(DenseMatrixTest, MultVec)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;

  Vector<double> v0(4);
  v0(0) = 1;
  v0(1) = -1;
  v0(2) = 3;
  v0(3) = -2;

  auto mv1 = a.Mult(v0);

  EXPECT_DOUBLE_EQ(mv1(0), 0);
  EXPECT_DOUBLE_EQ(mv1(1), 0);
}

TEST_F(DenseMatrixTest, MultMat)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;

  DenseMatrix<double> b(4, 3);
  b(0, 0) = -1;
  b(0, 1) = 0;
  b(0, 2) = 2;
  b(1, 0) = -4;
  b(1, 1) = 0;
  b(1, 2) = -2;
  b(2, 0) = 3;
  b(2, 1) = 5;
  b(2, 2) = 8;
  b(3, 0) = 1;
  b(3, 1) = -1;
  b(3, 2) = 0;
  auto mm1 = a.Mult(b);

  EXPECT_DOUBLE_EQ(mm1(0, 0), 4);
  EXPECT_DOUBLE_EQ(mm1(0, 1), 11);
  EXPECT_DOUBLE_EQ(mm1(0, 2), 22);
  EXPECT_DOUBLE_EQ(mm1(1, 0), -4);
  EXPECT_DOUBLE_EQ(mm1(1, 1), -11);
  EXPECT_DOUBLE_EQ(mm1(1, 2), -22);
}

TEST_F(DenseMatrixTest, Transposed)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  auto a_trans = a.Transposed();

  EXPECT_DOUBLE_EQ(a_trans(0, 0), 1);
  EXPECT_DOUBLE_EQ(a_trans(1, 0), 2);
  EXPECT_DOUBLE_EQ(a_trans(2, 0), 3);
  EXPECT_DOUBLE_EQ(a_trans(3, 0), 4);
  EXPECT_DOUBLE_EQ(a_trans(0, 1), -1);
  EXPECT_DOUBLE_EQ(a_trans(1, 1), -2);
  EXPECT_DOUBLE_EQ(a_trans(2, 1), -3);
  EXPECT_DOUBLE_EQ(a_trans(3, 1), -4);
}

TEST_F(DenseMatrixTest, Transpose)
{
  DenseMatrix<double> a(2, 4);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 3) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 3) = -4;
  a.Transpose();

  EXPECT_DOUBLE_EQ(a(0, 0), 1);
  EXPECT_DOUBLE_EQ(a(1, 0), 2);
  EXPECT_DOUBLE_EQ(a(2, 0), 3);
  EXPECT_DOUBLE_EQ(a(3, 0), 4);
  EXPECT_DOUBLE_EQ(a(0, 1), -1);
  EXPECT_DOUBLE_EQ(a(1, 1), -2);
  EXPECT_DOUBLE_EQ(a(2, 1), -3);
  EXPECT_DOUBLE_EQ(a(3, 1), -4);
}
