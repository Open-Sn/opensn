// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/vector.h"
#include "framework/math/dense_matrix.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <vector>
#include <array>

namespace opensn
{
namespace SimplifiedLDFESQ
{
struct SphericalQuadrilateral;
struct FunctionWeightFromRho;
class Quadrature;

/**
 * Base Functor to inherit from to change the function to integrate in one of the integration
 * utilities.
 */
struct BaseFunctor
{
  virtual double operator()(double mu, double eta, double xi) { return 0.0; }
};
} // namespace SimplifiedLDFESQ

/// Serves as a general data structure for a spherical quadrilateral (SQ).
struct SimplifiedLDFESQ::SphericalQuadrilateral
{
  /// On square
  std::array<Vector3, 4> vertices_xy_tilde;
  /// On cube face
  std::array<Vector3, 4> vertices_xyz_prime;
  /// On unit sphere
  std::array<Vector3, 4> vertices_xyz;
  Vector3 centroid_xyz;

  Matrix3x3 rotation_matrix;
  Vector3 translation_vector;

  std::array<Vector3, 4> sub_sqr_points;
  std::array<double, 4> sub_sqr_weights;

  double area = 0.0;

  Vector3 octant_modifier;
};

/// Piecewise-linear Finite element quadrature using quadrilaterals.
class SimplifiedLDFESQ::Quadrature : public AngularQuadrature
{
public:
  enum class QuadraturePointOptimization
  {
    CENTROID,
    EMPIRICAL,
    ISOLATED,
    //    MULTI_VARIATE_SECANT
  };
  QuadraturePointOptimization qp_optimization_type = QuadraturePointOptimization::EMPIRICAL;
  std::string output_filename_prefix;

private:
  /// Inscribed cude side length
  static constexpr double a = 0.57735026919;
  int initial_level_ = 0;
  std::vector<Vector3> diagonal_vertices_;
  std::vector<SphericalQuadrilateral> initial_octant_SQs_;

public:
  std::vector<SphericalQuadrilateral> deployed_SQs;

private:
  std::vector<std::vector<SphericalQuadrilateral>> deployed_SQs_history_;

public:
  friend struct FunctionWeightFromRho;
  Quadrature() : AngularQuadrature(AngularQuadratureType::SLDFESQ) {}

  virtual ~Quadrature() {}

  /// Generates uniform spherical quadrilaterals from the subdivision of an inscribed cube.
  void GenerateInitialRefinement(int level);

private:
  /// Generates diagonal spacings.
  void GenerateDiagonalSpacings(int level);

  /// Generates the standard points on the reference face.
  void GenerateReferenceFaceVertices(const Matrix3x3& rotation_matrix,
                                     const Vector3& translation,
                                     int level);

  /// Develops LDFE quantities.
  void DevelopSQLDFEValues(SphericalQuadrilateral& sq, GaussLegendreQuadrature& legendre);

  /// Applies empirical quadrature point optimization.
  void EmpiricalQPOptimization(SphericalQuadrilateral& sq,
                               GaussLegendreQuadrature& legendre,
                               Vector3& sq_xy_tilde_centroid,
                               std::array<Vector3, 4>& radii_vectors_xy_tilde,
                               std::array<double, 4>& sub_sub_sqr_areas);

  void IsolatedQPOptimization(SphericalQuadrilateral& sq,
                              GaussLegendreQuadrature& legendre,
                              Vector3& sq_xy_tilde_centroid,
                              std::array<Vector3, 4>& radii_vectors_xy_tilde,
                              std::array<double, 4>& sub_sub_sqr_areas);

  /**
   * Computes the area of a cell. This routine uses Girard's theorem to get the area of a spherical
   * triangle using the spherical excess.
   */
  static double ComputeSphericalQuadrilateralArea(std::array<Vector3, 4>& vertices_xyz);

  /// Integrates shape functions to produce weights.
  static std::array<double, 4>
  IntegrateLDFEShapeFunctions(const SphericalQuadrilateral& sq,
                              std::array<Vector<double>, 4>& shape_coeffs,
                              const std::vector<Vector3>& legendre_qpoints,
                              const std::vector<double>& legendre_qweights);

  /// Deploys the current set of SQs to all octants.
  void CopyToAllOctants();

  /// Populates the quadrature abscissaes, weights and direction vectors.
  void PopulateQuadratureAbscissae();

private:
  /// Performs a simple Riemann integral of a base functor.
  double RiemannIntegral(BaseFunctor* F, int Ni = 20000);

  /// Performs a quadrature integral of a base functor using the supplied SQs.
  double QuadratureSSIntegral(BaseFunctor* F);

public:
  /// Performs a test integration of predefined cases.
  void TestIntegration(int test_case, double ref_solution, int RiemannN = 0);

public:
  /// Prints the quadrature to file.
  void PrintQuadratureToFile(const std::string& file_base);

public:
  /// Locally refines the cells.
  void LocallyRefine(const Vector3& ref_dir,
                     const double cone_size,
                     const bool dir_as_plane_normal = false);

private:
  /// Split a SQ.
  std::array<SphericalQuadrilateral, 4> SplitSQ(SphericalQuadrilateral& sq,
                                                GaussLegendreQuadrature& legendre);
};

/**
 * This is a utility function that encapsulates all the necessary functionality to determine shape
 * function coefficients and integrate accross a spherical quadrilateral.
 */
struct SimplifiedLDFESQ::FunctionWeightFromRho
{
  Quadrature& sldfesq;
  Vector3& centroid_xy_tilde;
  std::array<Vector3, 4>& radii_vectors_xy_tilde;
  SphericalQuadrilateral& sq;

  std::array<Vector<double>, 4> rhs;
  DenseMatrix<double> A;
  DenseMatrix<double> A_inv;
  std::array<Vector<double>, 4> c_coeffs;
  /// Legendre quadrature points
  std::vector<Vector3>& lqp;
  /// Legendre quadrature weights
  std::vector<double>& lqw;

  FunctionWeightFromRho(SimplifiedLDFESQ::Quadrature& sldfesq,
                        Vector3& centroid_xy_tilde,
                        std::array<Vector3, 4>& radii_vectors_xy_tilde,
                        SphericalQuadrilateral& sq,
                        GaussLegendreQuadrature& legendre_quadrature)
    : sldfesq(sldfesq),
      centroid_xy_tilde(centroid_xy_tilde),
      radii_vectors_xy_tilde(radii_vectors_xy_tilde),
      sq(sq),
      A(4, 4),
      A_inv(4, 4),
      lqp(legendre_quadrature.qpoints),
      lqw(legendre_quadrature.weights)
  {
    // Init RHS
    for (int i = 0; i < 4; ++i)
    {
      rhs[i] = Vector<double>(4);
      c_coeffs[i] = Vector<double>(4, 0.);
      for (int j = 0; j < 4; ++j)
        rhs[i](i) = 1.0;
    }
  }

  /**
   * Computes the quadrature point locations from rho, followed by the shape-function coefficients
   * and then the integral of the shape function to get the weights.
   */
  std::array<double, 4> operator()(const Vector<double>& rho)
  {
    // Determine qpoints from rho
    std::array<Vector3, 4> qpoints;
    for (int i = 0; i < 4; ++i)
    {
      auto xy_tilde = centroid_xy_tilde + rho(i) * radii_vectors_xy_tilde[i];
      auto xyz_prime = sq.rotation_matrix * xy_tilde + sq.translation_vector;
      qpoints[i] = xyz_prime.Normalized();
    }

    // Assemble A
    for (int i = 0; i < 4; ++i)
      A.SetRow(i, Vector<double>({1.0, qpoints[i][0], qpoints[i][1], qpoints[i][2]}));

    // Compute A-inverse
    A_inv = Inverse(A);

    // Compute coefficients
    for (int i = 0; i < 4; ++i)
      c_coeffs[i] = Mult(A_inv, rhs[i]);

    return sldfesq.IntegrateLDFEShapeFunctions(sq, c_coeffs, lqp, lqw);
  }
};

} // namespace opensn
