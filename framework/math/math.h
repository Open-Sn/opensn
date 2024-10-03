// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/math_incdef.h"
#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/vector.h"
#include "framework/math/dense_matrix.h"
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
std::vector<double> Mult(const std::vector<double>& x, const double& val);

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

double ComputePointwiseChange(std::vector<double>& x, std::vector<double>& y);

double ComputeL2Change(std::vector<double>& x, std::vector<double>& y);

} // namespace opensn
