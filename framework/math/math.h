// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/math_incdef.h"
#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include <memory>

namespace opensn
{
class SparseMatrix;
class UnknownManager;
class CDFSampler;
class SpatialDiscretization;
class SpatialDiscretization_FV;
class SpatialDiscretization_PWLD;
class SpatialDiscretization_PWLC;

using MatDbl = std::vector<std::vector<double>>;
using MatVec3 = std::vector<std::vector<Vector3>>;

/// Coordinate system type.
enum class CoordinateSystemType
{
  UNDEFINED = 0,
  CARTESIAN = 1,
  CYLINDRICAL = 2,
  SPHERICAL = 3,
};

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

/**
 * Sample a Cumulative Distribution Function (CDF) given a probability.
 *
 * The supplied vector should contain the upper bin boundary for each
 * bin and will return the bin associated with the bin that brackets
 * the supplied probability.
 *
 * Example:
 * Suppose we sample bins 0-9. Suppose also that the probalities for each
 * bin is as follows:
 * - 0.1 bin 0
 * - 0.1 bin 1
 * - 0.5 bin 5
 * - 0.3 bin 8
 *
 * The CDF for this probability distribution will look like this
 * - bin 0 = 0.1
 * - bin 1 = 0.2
 * - bin 2 = 0.2
 * - bin 3 = 0.2
 * - bin 4 = 0.2
 * - bin 5 = 0.7
 * - bin 6 = 0.7
 * - bin 7 = 0.7
 * - bin 8 = 1.0
 * - bin 9 = 1.0
 *
 * Supplying a random number between 0 and 1 should indicate sampling one
 * of the bins 0,1,5 or 8. The most inefficient way to do this is to
 * linearly loop through the cdf and check \f$ cdf_{i-1} \ge \theta < cdf_i \f$.
 *  An optimized version of this sampling would be to perform a recursive
 *  block search which starts with a course view of the cdf and then gradually
 *  refines the view until the final linear search can be performed.*/
int SampleCDF(double x, std::vector<double> cdf_bin);

/// Computes the factorial of an integer.
double Factorial(int x);

/**
 * Determines the azimuthal- and polar-angle associated with the given direction vector.
 * Returns a pair = [azimuthal-angle,polar-angle].
 */
std::pair<double, double> OmegaToPhiThetaSafe(const Vector3& omega);

/// Prints the Vector.
void PrintVector(const std::vector<double>& x);

/// Scales a vector in place by constant.
void Scale(std::vector<double>& x, const double& val);

/// Sets a constant value to a vector.
void Set(std::vector<double>& x, const double& val);

/// Multiplies the vector with a constant and returns result.
std::vector<double> VecMul(const std::vector<double>& x, const double& val);

/**
 * Returns the 1-norm. Also known as the Taxicab or Manhattan norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{1}=\sum_{i=1}^{n}\left|x_{i}\right|
 * \f]
 */
double Vec1Norm(const std::vector<double>& x);

/**
 * Returns the 2-norm. Also known as the Euclidian or Frobenius norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{2}=\sqrt{x_{1}^{2}+\cdots+x_{n}^{2}}
 * \f]
 */
double Vec2Norm(const std::vector<double>& x);

/**
 * Returns the infinity-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{\infty}=\max \left(\left|x_{1}\right|,
 * \ldots,\left|x_{n}\right|\right) \f]
 */
double VecInfinityNorm(const std::vector<double>& x);

/**
 * Returns the p-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{p}=\left(\sum_{i=1}^{n}\left|x_{i}\right|^{p}\right)^{1 / p}
 * \f]
 */
double VecPNorm(const std::vector<double>& x, const double& p);

/**
 * Computes the dot product of two vectors.
 *
 * \f[
 * \mathrm{a} \cdot \mathrm{b}=\sum_{i=1}^{n} a_{i} b_{i}
 * \f]
 */
double Dot(const std::vector<double>& x, const std::vector<double>& y);

/// Adds two vectors component-wise.
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

/// Subtracts two vectors component-wise.
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);

/// Prints the contents of a matrix.
void PrintMatrix(const MatDbl& A);

/// Scales the matrix by a constant value.
void Scale(MatDbl& A, const double& val);

/// Sets all the entries of the matrix to a constant value.
void Set(MatDbl& A, const double& val);

/// Returns the transpose of a matrix.
MatDbl Transpose(const MatDbl& A);

/// Swaps two rows of a matrix.
void SwapRow(size_t r1, size_t r2, MatDbl& A);

/// Swaps two columns of a matrix.
void SwapColumn(size_t c1, size_t c2, MatDbl& A);

/// Multiply matrix with a constant and return result.
MatDbl MatMul(const MatDbl& A, const double c);

/// Multiply matrix with a vector and return resulting vector
std::vector<double> MatMul(const MatDbl& A, const std::vector<double>& x);

/// Mutliply two matrices and return result.
MatDbl MatMul(const MatDbl& A, const MatDbl& B);

/// Adds two matrices and returns the result.
MatDbl MatAdd(const MatDbl& A, const MatDbl& B);

/// Subtracts matrix A from B and returns the result.
MatDbl MatSubtract(const MatDbl& A, const MatDbl& B);

/// Computes the determinant of a matrix.
double Determinant(const MatDbl& A);

/// Returns a sub-matrix.
MatDbl SubMatrix(const size_t r, const size_t c, const MatDbl& A);

/// Gauss Elimination without pivoting.
void GaussElimination(MatDbl& A, std::vector<double>& b, int n);

/// Computes the inverse of a matrix using Gauss-Elimination with pivoting.
MatDbl InverseGEPivoting(const MatDbl& A);

/// Computes the inverse of a matrix.
MatDbl Inverse(const MatDbl& A);

/**
 * Performs power iteration to obtain the fundamental eigen mode. The eigen-value of the fundamental
 * mode is return whilst the eigen-vector is return via reference.
 */
double PowerIteration(const MatDbl& A,
                      std::vector<double>& e_vec,
                      int max_it = 2000,
                      double tol = 1.0e-13);

double ComputePointwiseChange(std::vector<double>& x, std::vector<double>& y);

double ComputeL2Change(std::vector<double>& x, std::vector<double>& y);

} // namespace opensn
