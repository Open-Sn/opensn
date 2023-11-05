#pragma once

#include "opensn/framework/math/chi_math_incdef.h"

#include "opensn/framework/math/Quadratures/quadrature.h"
#include "opensn/framework/math/Quadratures/angular_quadrature_base.h"
#include "opensn/framework/math/UnknownManager/unknown_manager.h"

#include <memory>

typedef std::vector<double> VecDbl;
typedef std::vector<VecDbl> MatDbl;

namespace chi_math
{
class SparseMatrix;
class UnknownManager;
class CDFSampler;
class SpatialDiscretization;
class SpatialDiscretization_FV;
class SpatialDiscretization_PWLD;
class SpatialDiscretization_PWLC;

/**
 * Coordinate system type.
 */
enum class CoordinateSystemType
{
  UNDEFINED = 0,
  CARTESIAN = 1,
  CYLINDRICAL = 2,
  SPHERICAL = 3,
};

/**
 * Spatial discretization type.
 */
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

int SampleCDF(double x, std::vector<double> cdf_bin);

/**
 * Computes the factorial of an integer.
 */
double Factorial(int x);

/**
 * Determines the azimuthal- and polar-angle associated with the given direction vector.
 * Returns a pair = [azimuthal-angle,polar-angle].
 */
std::pair<double, double> OmegaToPhiThetaSafe(const chi_mesh::Vector3& omega);

/**
 * Prints the Vector.
 */
void PrintVector(const VecDbl& x);

/**
 * Scales a vector in place by constant.
 */
void Scale(VecDbl& x, const double& val);

/**
 * Sets a constant value to a vector.
 */
void Set(VecDbl& x, const double& val);

/**
 * Multiplies the vector with a constant and returns result.
 */
VecDbl VecMul(const VecDbl& x, const double& val);

/**
 * Returns the 1-norm. Also known as the Taxicab or Manhattan norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{1}=\sum_{i=1}^{n}\left|x_{i}\right|
 * \f]
 *
 */
double Vec1Norm(const VecDbl& x);

/**
 * Returns the 2-norm. Also known as the Euclidian or Frobenius norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{2}=\sqrt{x_{1}^{2}+\cdots+x_{n}^{2}}
 * \f]
 *
 */
double Vec2Norm(const VecDbl& x);

/** Returns the infinity-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{\infty}=\max \left(\left|x_{1}\right|,
 * \ldots,\left|x_{n}\right|\right) \f]
 *
 */
double VecInfinityNorm(const VecDbl& x);

/**
 * Returns the p-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{p}=\left(\sum_{i=1}^{n}\left|x_{i}\right|^{p}\right)^{1 / p}
 * \f]
 *
 */
double VecPNorm(const VecDbl& x, const double& p);

/** Computes the dot product of two vectors.
 *
 * \f[
 * \mathrm{a} \cdot \mathrm{b}=\sum_{i=1}^{n} a_{i} b_{i}
 * \f]
 */
double Dot(const VecDbl& x, const VecDbl& y);

/**
 * Adds two vectors component-wise.
 */
VecDbl operator+(const VecDbl& a, const VecDbl& b);

/**
 * Subtracts two vectors component-wise.
 */
VecDbl operator-(const VecDbl& a, const VecDbl& b);

/**
 * Prints the contents of a matrix.
 */
void PrintMatrix(const MatDbl& A);

/**
 * Scales the matrix by a constant value.
 */
void Scale(MatDbl& A, const double& val);

/**
 * Sets all the entries of the matrix to a constant value.
 */
void Set(MatDbl& A, const double& val);

/**
 * Returns the transpose of a matrix.
 */
MatDbl Transpose(const MatDbl& A);

/**
 * Swaps two rows of a matrix.
 */
void SwapRow(size_t r1, size_t r2, MatDbl& A);

/**
 * Swaps two columns of a matrix.
 */
void SwapColumn(size_t c1, size_t c2, MatDbl& A);

/**
 * Multiply matrix with a constant and return result.
 */
MatDbl MatMul(const MatDbl& A, const double c);

/**
 * Multiply matrix with a vector and return resulting vector
 */
VecDbl MatMul(const MatDbl& A, const VecDbl& x);

/**
 * Mutliply two matrices and return result.
 */
MatDbl MatMul(const MatDbl& A, const MatDbl& B);

/**
 * Adds two matrices and returns the result.
 */
MatDbl MatAdd(const MatDbl& A, const MatDbl& B);

/**
 * Subtracts matrix A from B and returns the result.
 */
MatDbl MatSubtract(const MatDbl& A, const MatDbl& B);

/**
 * Computes the determinant of a matrix.
 */
double Determinant(const MatDbl& A);

/**
 * Returns a sub-matrix.
 */
MatDbl SubMatrix(const size_t r, const size_t c, const MatDbl& A);

/**
 * Gauss Elimination without pivoting.
 */
void GaussElimination(MatDbl& A, VecDbl& b, int n);

/**
 * Computes the inverse of a matrix using Gauss-Elimination with pivoting.
 */
MatDbl InverseGEPivoting(const MatDbl& A);

/**
 * Computes the inverse of a matrix.
 */
MatDbl Inverse(const MatDbl& A);

/**
 * Performs power iteration to obtain the fundamental eigen mode. The eigen-value of the fundamental
 * mode is return whilst the eigen-vector is return via reference.
 */
double PowerIteration(const MatDbl& A, VecDbl& e_vec, int max_it = 2000, double tol = 1.0e-13);

} // namespace chi_math
