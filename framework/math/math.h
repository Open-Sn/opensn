// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/data_types/vector.h"
#include "framework/data_types/dense_matrix.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <memory>

namespace opensn
{
class SparseMatrix;
class UnknownManager;
class SpatialDiscretization_FV;
class SpatialDiscretization_PWLD;
class SpatialDiscretization_PWLC;

using MatVec3 = std::vector<std::vector<Vector3>>;

/// Spatial discretization type.
enum class SpatialDiscretizationType
{
  UNDEFINED = 0,
  FINITE_VOLUME = 1,
  PIECEWISE_LINEAR_CONTINUOUS = 2,
  PIECEWISE_LINEAR_DISCONTINUOUS = 3,
  LAGRANGE_CONTINUOUS = 4,
  LAGRANGE_DISCONTINUOUS = 5
};

enum class NormType : int
{
  L1_NORM = 1,
  L2_NORM = 2,
  LINF_NORM = 3
};

/**Sample a Cumulative Distribution Function (CDF) given a probability.
 *
 * The supplied vector should contain the upper bin boundary for each bin and will return
 * the bin associated with the bin that brackets the supplied probability.
 */
int SampleCDF(double x, std::vector<double> cdf_bin);

/** Return the transpose of a matrix.
 *
 * For an m x n input A, returns an n x m matrix T where T[j][i] = A[i][j].
 * \throws std::runtime_error if the matrix is empty or has inconsistent row sizes.
 */
std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>>& matrix);

/** Invert a square matrix using SVD-based pseudo-inversion (LAPACK dgesvd).
 *
 * Singular values below machine-epsilon * n * sigma_max are treated as zero.
 * Emits a warning if the condition number exceeds 1e10.
 * \throws std::runtime_error if the matrix is empty, non-square, or SVD fails.
 */
std::vector<std::vector<double>> InvertMatrix(const std::vector<std::vector<double>>& matrix);

/** Orthogonalize the columns of a matrix under a weighted inner product.
 *
 * Applies double-pass modified Gram-Schmidt and normalizes each column so that
 * \f$\langle q_i, q_j \rangle_w = \delta_{ij}\f$.  Columns whose weighted norm
 * falls below 1e-14 are left as zero vectors (linearly dependent directions).
 * \throws std::runtime_error if the matrix is empty or has inconsistent row sizes.
 */
std::vector<std::vector<double>>
OrthogonalizeMatrixSpan(const std::vector<std::vector<double>>& matrix,
                        const std::vector<double>& weights);

/// Computes the factorial of an integer.
double Factorial(int x);

/** Return the azimuthal and polar angles for the given direction vector as a pair
 * [azimuthal-angle, polar-angle].
 */
std::pair<double, double> OmegaToPhiThetaSafe(const Vector3& omega);

/// Prints the Vector.
void PrintVector(const std::vector<double>& x);

/// Scales a vector in place by constant.
void Scale(std::vector<double>& x, const double& val);

/// Sets a constant value to a vector.
void Set(std::vector<double>& x, const double& val);

/// Multiplies the vector with a constant and returns result.
std::vector<double> Mult(const std::vector<double>& x, const double& val);

/** Return the 1-norm (Taxicab / Manhattan norm):
 * \f$\|\boldsymbol{x}\|_{1}=\sum_{i=1}^{n}\left|x_{i}\right|\f$
 */
double L1Norm(const std::vector<double>& x);

/** Return the 2-norm (Euclidean / Frobenius norm):
 * \f$\|\boldsymbol{x}\|_{2}=\sqrt{x_{1}^{2}+\cdots+x_{n}^{2}}\f$
 */
double L2Norm(const std::vector<double>& x);

/** Return the infinity-norm:
 * \f$\|\mathbf{x}\|_{\infty}=\max\left(\left|x_{1}\right|,\ldots,\left|x_{n}\right|\right)\f$
 */
double LInfNorm(const std::vector<double>& x);

/** Return the p-norm:
 * \f$\|\mathbf{x}\|_{p}=\left(\sum_{i=1}^{n}\left|x_{i}\right|^{p}\right)^{1/p}\f$
 */
double LpNorm(const std::vector<double>& x, const double& p);

/** Compute the dot product of two vectors:
 * \f$a \cdot b=\sum_{i=1}^{n} a_{i} b_{i}\f$
 */
double Dot(const std::vector<double>& x, const std::vector<double>& y);

/// Adds two vectors component-wise.
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

/// Subtracts two vectors component-wise.
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);

double ComputePointwiseChange(std::vector<double>& x, std::vector<double>& y);

double ComputeL2Change(std::vector<double>& x, std::vector<double>& y);

} // namespace opensn
