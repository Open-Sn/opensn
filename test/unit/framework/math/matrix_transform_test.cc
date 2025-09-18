// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/opensn_unit_test.h"
#include "framework/math/math.h"
#include <gtest/gtest.h>
#include <cmath>

using namespace opensn;

class MatrixTransformTest : public OpenSnUnitTest
{
};

// ---------------------------------------------------------------------------
// Transpose
// ---------------------------------------------------------------------------

TEST_F(MatrixTransformTest, TransposeNonSquare)
{
  // 2x3 input
  const std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  const auto AT = Transpose(A);

  ASSERT_EQ(AT.size(), 3u);
  for (const auto& row : AT)
    ASSERT_EQ(row.size(), 2u);

  EXPECT_DOUBLE_EQ(AT[0][0], 1.0);
  EXPECT_DOUBLE_EQ(AT[0][1], 4.0);
  EXPECT_DOUBLE_EQ(AT[1][0], 2.0);
  EXPECT_DOUBLE_EQ(AT[1][1], 5.0);
  EXPECT_DOUBLE_EQ(AT[2][0], 3.0);
  EXPECT_DOUBLE_EQ(AT[2][1], 6.0);
}

TEST_F(MatrixTransformTest, TransposeSquare)
{
  const std::vector<std::vector<double>> A = {{1.0, 2.0}, {3.0, 4.0}};
  const auto AT = Transpose(A);

  ASSERT_EQ(AT.size(), 2u);
  EXPECT_DOUBLE_EQ(AT[0][0], 1.0);
  EXPECT_DOUBLE_EQ(AT[0][1], 3.0);
  EXPECT_DOUBLE_EQ(AT[1][0], 2.0);
  EXPECT_DOUBLE_EQ(AT[1][1], 4.0);
}

TEST_F(MatrixTransformTest, TransposeOfTransposeIsOriginal)
{
  const std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  const auto ATA = Transpose(Transpose(A));

  ASSERT_EQ(ATA.size(), A.size());
  for (size_t i = 0; i < A.size(); ++i)
    for (size_t j = 0; j < A[i].size(); ++j)
      EXPECT_DOUBLE_EQ(ATA[i][j], A[i][j]);
}

TEST_F(MatrixTransformTest, TransposeEmptyThrows)
{
  EXPECT_THROW(Transpose({}), std::runtime_error);
}

// ---------------------------------------------------------------------------
// InvertMatrix
// ---------------------------------------------------------------------------

TEST_F(MatrixTransformTest, InvertIdentity2x2)
{
  const std::vector<std::vector<double>> I = {{1.0, 0.0}, {0.0, 1.0}};
  const auto Iinv = InvertMatrix(I);

  ASSERT_EQ(Iinv.size(), 2u);
  EXPECT_NEAR(Iinv[0][0], 1.0, 1e-12);
  EXPECT_NEAR(Iinv[0][1], 0.0, 1e-12);
  EXPECT_NEAR(Iinv[1][0], 0.0, 1e-12);
  EXPECT_NEAR(Iinv[1][1], 1.0, 1e-12);
}

TEST_F(MatrixTransformTest, InvertKnown2x2)
{
  // [[2, 1], [1, 1]]^{-1} = [[1, -1], [-1, 2]]
  const std::vector<std::vector<double>> A = {{2.0, 1.0}, {1.0, 1.0}};
  const auto Ainv = InvertMatrix(A);

  ASSERT_EQ(Ainv.size(), 2u);
  EXPECT_NEAR(Ainv[0][0], 1.0, 1e-12);
  EXPECT_NEAR(Ainv[0][1], -1.0, 1e-12);
  EXPECT_NEAR(Ainv[1][0], -1.0, 1e-12);
  EXPECT_NEAR(Ainv[1][1], 2.0, 1e-12);
}

TEST_F(MatrixTransformTest, InvertRoundTripIsIdentity)
{
  // A * A^{-1} should be the identity (within tolerance)
  const std::vector<std::vector<double>> A = {{4.0, 3.0}, {6.0, 5.0}};
  const auto Ainv = InvertMatrix(A);

  // Multiply A * Ainv row by row
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
    {
      double val = 0.0;
      for (size_t k = 0; k < 2; ++k)
        val += A[i][k] * Ainv[k][j];
      EXPECT_NEAR(val, (i == j) ? 1.0 : 0.0, 1e-10);
    }
}

TEST_F(MatrixTransformTest, InvertDiagonal3x3)
{
  const std::vector<std::vector<double>> D = {{2.0, 0.0, 0.0}, {0.0, 4.0, 0.0}, {0.0, 0.0, 0.5}};
  const auto Dinv = InvertMatrix(D);

  EXPECT_NEAR(Dinv[0][0], 0.5, 1e-12);
  EXPECT_NEAR(Dinv[1][1], 0.25, 1e-12);
  EXPECT_NEAR(Dinv[2][2], 2.0, 1e-12);
}

TEST_F(MatrixTransformTest, InvertEmptyThrows)
{
  EXPECT_THROW(InvertMatrix({}), std::runtime_error);
}

// ---------------------------------------------------------------------------
// OrthogonalizeMatrixSpan
// ---------------------------------------------------------------------------

TEST_F(MatrixTransformTest, OrthogonalizeProducesOrthonormalColumns)
{
  // Two non-orthogonal vectors with uniform weights
  const std::vector<std::vector<double>> A = {
    {1.0, 1.0},
    {0.0, 1.0},
    {0.0, 0.0},
  };
  const std::vector<double> w = {1.0, 1.0, 1.0};

  const auto Q = OrthogonalizeMatrixSpan(A, w);

  ASSERT_EQ(Q.size(), 3u);
  ASSERT_EQ(Q[0].size(), 2u);

  // Check norm = 1
  double n00 = 0.0;
  for (size_t r = 0; r < 3; ++r)
    n00 += w[r] * Q[r][0] * Q[r][0];
  EXPECT_NEAR(n00, 1.0, 1e-12);

  // Check norm = 1
  double n11 = 0.0;
  for (size_t r = 0; r < 3; ++r)
    n11 += w[r] * Q[r][1] * Q[r][1];
  EXPECT_NEAR(n11, 1.0, 1e-12);

  // Check inner-product = 1
  double n01 = 0.0;
  for (size_t r = 0; r < 3; ++r)
    n01 += w[r] * Q[r][0] * Q[r][1];
  EXPECT_NEAR(n01, 0.0, 1e-12);
}

TEST_F(MatrixTransformTest, OrthogonalizeSingleColumnNormalized)
{
  const std::vector<std::vector<double>> A = {{3.0}, {4.0}};
  const std::vector<double> w = {1.0, 1.0};

  const auto Q = OrthogonalizeMatrixSpan(A, w);

  // Weighted norm of the single column should be 1
  double norm_sq = w[0] * Q[0][0] * Q[0][0] + w[1] * Q[1][0] * Q[1][0];
  EXPECT_NEAR(norm_sq, 1.0, 1e-12);
}

TEST_F(MatrixTransformTest, OrthogonalizeNonUniformWeights)
{
  // Three-row two-column matrix with non-uniform weights
  const std::vector<std::vector<double>> A = {
    {1.0, 0.0},
    {1.0, 1.0},
    {1.0, 2.0},
  };
  const std::vector<double> w = {0.5, 1.0, 0.5};

  const auto Q = OrthogonalizeMatrixSpan(A, w);

  // Verify orthonormality under the given weights
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
    {
      double inner = 0.0;
      for (size_t r = 0; r < 3; ++r)
        inner += w[r] * Q[r][i] * Q[r][j];
      EXPECT_NEAR(inner, (i == j) ? 1.0 : 0.0, 1e-12);
    }
}

TEST_F(MatrixTransformTest, OrthogonalizeEmptyThrows)
{
  EXPECT_THROW(OrthogonalizeMatrixSpan({}, {}), std::runtime_error);
}
