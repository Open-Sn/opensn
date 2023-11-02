#ifndef CHI_MATH_H
#define CHI_MATH_H

#include "chi_math_incdef.h"

#include "Quadratures/quadrature.h"
#include "Quadratures/angular_quadrature_base.h"
#include "UnknownManager/unknown_manager.h"

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

/**Coordinate system type.*/
enum class CoordinateSystemType
{
  UNDEFINED = 0,
  CARTESIAN = 1,
  CYLINDRICAL = 2,
  SPHERICAL = 3,
};
/**Spatial discretization type.*/
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

// 01 Utility
double Factorial(int x);

std::pair<double, double> OmegaToPhiThetaSafe(const chi_mesh::Vector3& omega);

// 02 Vector operations
void PrintVector(const VecDbl& x);
void Scale(VecDbl& x, const double& val);
void Set(VecDbl& x, const double& val);
VecDbl VecMul(const VecDbl& x, const double& val);
double Vec1Norm(const VecDbl& x);
double Vec2Norm(const VecDbl& x);
double VecInfinityNorm(const VecDbl& x);
double VecPNorm(const VecDbl& x, const double& p);
double Dot(const VecDbl& x, const VecDbl& y);
VecDbl operator+(const VecDbl& a, const VecDbl& b);
VecDbl operator-(const VecDbl& a, const VecDbl& b);

// 03 Matrix operations
void PrintMatrix(const MatDbl& A);
void Scale(MatDbl& A, const double& val);
void Set(MatDbl& A, const double& val);
MatDbl Transpose(const MatDbl& A);
void SwapRow(size_t r1, size_t r2, MatDbl& A);
void SwapColumn(size_t c1, size_t c2, MatDbl& A);
MatDbl MatMul(const MatDbl& A, const double c);
VecDbl MatMul(const MatDbl& A, const VecDbl& x);
MatDbl MatMul(const MatDbl& A, const MatDbl& B);
MatDbl MatAdd(const MatDbl& A, const MatDbl& B);
MatDbl MatSubtract(const MatDbl& A, const MatDbl& B);
double Determinant(const MatDbl& A);
MatDbl SubMatrix(const size_t r, const size_t c, const MatDbl& A);
void GaussElimination(MatDbl& A, VecDbl& b, int n);
MatDbl InverseGEPivoting(const MatDbl& A);
MatDbl Inverse(const MatDbl& A);

double PowerIteration(const MatDbl& A, VecDbl& e_vec, int max_it = 2000, double tol = 1.0e-13);

} // namespace chi_math

#endif
