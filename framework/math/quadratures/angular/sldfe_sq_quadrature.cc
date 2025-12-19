// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/data_types/vector3.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace opensn
{

SLDFEsqQuadrature::SLDFEsqQuadrature(int level,
                                     unsigned int dimension,
                                     unsigned int scattering_order)
  : AngularQuadrature(AngularQuadratureType::SLDFEsq, dimension, scattering_order), level_(level)
{
  if (level < 0)
    throw std::invalid_argument("SLDFEsqQuadrature: level must be non-negative");
}

std::array<std::array<Vector3, 4>, 4>
SLDFEsqQuadrature::BuildSubSquares(const std::array<Vector3, 4>& vertices_xy_tilde,
                                   const Vector3& centroid)
{
  const auto mid01 = 0.5 * (vertices_xy_tilde[0] + vertices_xy_tilde[1]);
  const auto mid12 = 0.5 * (vertices_xy_tilde[1] + vertices_xy_tilde[2]);
  const auto mid23 = 0.5 * (vertices_xy_tilde[2] + vertices_xy_tilde[3]);
  const auto mid03 = 0.5 * (vertices_xy_tilde[0] + vertices_xy_tilde[3]);

  std::array<std::array<Vector3, 4>, 4> sub_squares{};
  sub_squares[0] = {vertices_xy_tilde[0], mid01, centroid, mid03};
  sub_squares[1] = {mid01, vertices_xy_tilde[1], mid12, centroid};
  sub_squares[2] = {centroid, mid12, vertices_xy_tilde[2], mid23};
  sub_squares[3] = {mid03, centroid, mid23, vertices_xy_tilde[3]};

  return sub_squares;
}

void
SLDFEsqQuadrature::GenerateRefinement()
{
  initial_octant_SQs_.clear();
  deployed_SQs.clear();

  // Define constants
  const Vector3 ihat = Vector3(1.0, 0.0, 0.0);
  const Vector3 jhat = Vector3(0.0, 1.0, 0.0);
  const Vector3 khat = Vector3(0.0, 0.0, 1.0);

  // Build rotation matrices for cube faces
  Matrix3x3 Rxface;
  Rxface.SetColJVec(0, jhat);
  Rxface.SetColJVec(1, khat);
  Rxface.SetColJVec(2, ihat);

  Matrix3x3 Ryface;
  Ryface.SetColJVec(0, ihat);
  Ryface.SetColJVec(1, khat);
  Ryface.SetColJVec(2, -1.0 * jhat);

  Matrix3x3 Rzface;
  Rzface.SetDiagonalVec(1.0, 1.0, 1.0);

  const std::array<std::pair<Matrix3x3, Vector3>, 3> face_transforms = {
    std::make_pair(Rxface, SLDFEsqQuadrature::a_ * ihat),
    std::make_pair(Ryface, SLDFEsqQuadrature::a_ * jhat),
    std::make_pair(Rzface, SLDFEsqQuadrature::a_ * khat)};

  // Generate general diagonal spacings in xy-tilde coordinates
  GenerateDiagonalSpacings();

  // Generate vertices for each face of inscribed cube
  for (const auto& [rotation, translation] : face_transforms)
    GenerateReferenceFaceVertices(rotation, translation);

  // Generate values for all octants
  CopyToAllOctants();

  // Populate quadriture points
  PopulateQuadratureAbscissae();

  // Build m2d and d2m operators
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

void
SLDFEsqQuadrature::GenerateDiagonalSpacings()
{
  // Define constants
  const int Ns = (level_ + 1); // Number of subdivisions
  const int Np = Ns + 1;       // Number of diagonal points

  const auto ihat = Vector3(1.0, 0.0, 0.0);

  // Define rotation matrix
  Matrix3x3 Rihat;
  auto n = Vector3(0.0, -1.0 / std::sqrt(2), 1.0 / std::sqrt(2));
  const auto& t = ihat;
  auto b = n.Cross(t).Normalized();

  Rihat.SetColJVec(0, t);
  Rihat.SetColJVec(1, b);
  Rihat.SetColJVec(2, n);

  // Generate p-points
  std::vector<Vector3> p_points(Np);
  double dphi = acos(a_) / Ns;
  double alpha = 0.10005;
  double beta = 1.0185;

  for (int i = 0; i < Np; ++i)
  {
    double phi = i * dphi * (1.0 + alpha * (cos(beta * M_PI_2 * i / Ns) - cos(beta * M_PI_2)));
    p_points[i] = Rihat * Vector3(cos(phi), sin(phi), 0.0);
  }

  // Compute tilde points
  std::vector<Vector3> tilde_points(Np);
  for (int i = 0; i < Np; ++i)
  {
    double R = SLDFEsqQuadrature::a_ / p_points[i].x;
    double x_tilde = p_points[i].y * R;
    double y_tilde = p_points[i].z * R;
    tilde_points[i] = Vector3(x_tilde, y_tilde, 0.0);
  }

  diagonal_vertices_ = tilde_points;
}

void
SLDFEsqQuadrature::GenerateReferenceFaceVertices(const Matrix3x3& rotation_matrix,
                                                 const Vector3& translation)
{
  int Ns = (level_ + 1); // Number of subdivisions
  int Np = Ns + 1;       // Number of diagonal points

  GaussLegendreQuadrature legendre(QuadratureOrder::THIRTYSECOND);

  // Generate xy_tilde values
  std::vector<std::vector<Vector3>> vertices_xy_tilde_ij;
  vertices_xy_tilde_ij.resize(Np, std::vector<Vector3>(Np));
  for (int i = 0; i < Np; ++i)
    for (int j = 0; j < Np; ++j)
      vertices_xy_tilde_ij[i][j] = Vector3(diagonal_vertices_[i].x, diagonal_vertices_[j].y, 0.0);

  // Generate SQs
  for (int i = 0; i < Ns; ++i)
  {
    for (int j = 0; j < Ns; ++j)
    {
      SphericalQuadrilateral sq;

      sq.rotation_matrix = rotation_matrix;
      sq.translation_vector = translation;

      // Set xy-tilde vertices
      sq.vertices_xy_tilde[0] = vertices_xy_tilde_ij[i][j];
      sq.vertices_xy_tilde[1] = vertices_xy_tilde_ij[i + 1][j];
      sq.vertices_xy_tilde[2] = vertices_xy_tilde_ij[i + 1][j + 1];
      sq.vertices_xy_tilde[3] = vertices_xy_tilde_ij[i][j + 1];
      auto& vxy = sq.vertices_xy_tilde;

      // Set xyz_prime vertices
      for (int v = 0; v < 4; ++v)
        sq.vertices_xyz_prime[v] = rotation_matrix * vxy[v] + translation;

      // Set xyz vertices
      for (int v = 0; v < 4; ++v)
        sq.vertices_xyz[v] = sq.vertices_xyz_prime[v].Normalized();

      // Compute SQ xyz-centroid
      for (auto& vertex : sq.vertices_xyz)
        sq.centroid_xyz += vertex;
      sq.centroid_xyz /= 4;
      sq.centroid_xyz.Normalize();

      auto v0 = sq.centroid_xyz.Normalized();
      auto v1 = sq.vertices_xyz[0];
      auto v2 = sq.vertices_xyz[1];

      // Correction orientation
      if ((v1 - v0).Cross(v2 - v0).Dot(v0) < 0.0)
      {
        std::reverse(sq.vertices_xy_tilde.begin(), sq.vertices_xy_tilde.end());
        std::reverse(sq.vertices_xyz_prime.begin(), sq.vertices_xyz_prime.end());
        std::reverse(sq.vertices_xyz.begin(), sq.vertices_xyz.end());
      }

      // Compute area
      sq.area = ComputeSphericalQuadrilateralArea(sq.vertices_xyz);

      // Set octant modifier
      sq.octant_modifier = Vector3(1.0, 1.0, 1.0);

      // Develop LDFE values
      DevelopSQLDFEValues(sq, legendre);

      initial_octant_SQs_.push_back(sq);
    }
  }
}

void
SLDFEsqQuadrature::EmpiricalQPOptimization(SphericalQuadrilateral& sq,
                                           GaussLegendreQuadrature& legendre,
                                           Vector3& sq_xy_tilde_centroid,
                                           std::array<Vector3, 4>& radii_vectors_xy_tilde,
                                           std::array<double, 4>& sub_sub_sqr_areas)
{
  auto ComputeWeights =
    FunctionWeightFromRho(*this, sq_xy_tilde_centroid, radii_vectors_xy_tilde, sq, legendre);
  double d = 1.0 / std::sqrt(3.0);
  Vector<double> rho({d, d, d, d});

  auto weights = ComputeWeights(rho);

  for (int i = 0; i < 4; ++i)
  {
    auto xy_tilde = sq_xy_tilde_centroid + rho(i) * radii_vectors_xy_tilde[i];
    auto xyz_prime = sq.rotation_matrix * xy_tilde + sq.translation_vector;
    sq.sub_sqr_points[i] = xyz_prime.Normalized();
    sq.sub_sqr_weights[i] = weights[i];
  }
}

void
SLDFEsqQuadrature::IsolatedQPOptimization(SphericalQuadrilateral& sq,
                                          GaussLegendreQuadrature& legendre,
                                          Vector3& sq_xy_tilde_centroid,
                                          std::array<Vector3, 4>& radii_vectors_xy_tilde,
                                          std::array<double, 4>& sub_sub_sqr_areas)
{
  auto& SA_i = sub_sub_sqr_areas;
  auto ComputeWeights =
    FunctionWeightFromRho(*this, sq_xy_tilde_centroid, radii_vectors_xy_tilde, sq, legendre);
  double d = 1.0 / std::sqrt(3.0);
  Vector<double> rho({d, d, d, d});
  double epsilon = 1.0e-1;
  Vector<double> delta({epsilon, epsilon, epsilon, epsilon});
  Vector<double> drho_df({0.0, 0.0, 0.0, 0.0});

  // Compute initial weights
  auto weights = ComputeWeights(rho);

  for (int k = 0; k < 150; ++k) // iteration
  {
    auto weights_offset_pos = ComputeWeights(Add(rho, delta));
    auto weights_offset_neg = ComputeWeights(Subtract(rho, delta));

    double rho_change_total = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      double slope = 0.0;
      slope += 0.5 * (weights_offset_pos[i] - weights[i]);
      slope -= 0.5 * (weights_offset_neg[i] - weights[i]);
      drho_df(i) = delta(i) / slope;

      double delta_rho = 1.0 * drho_df(i) * (SA_i[i] - weights[i]);
      rho(i) += delta_rho;
      rho(i) = std::fmax(0.0, rho(i));
      rho(i) = std::fmin(1.0, rho(i));
      rho_change_total -= 1.0 * drho_df(i) * (weights[i] - SA_i[i]);
    }

    // Update weights
    weights = ComputeWeights(rho);

    if (rho_change_total < 1.0e-2)
      break;
  }
  weights = ComputeWeights(rho);

  for (int i = 0; i < 4; ++i)
  {
    auto xy_tilde = sq_xy_tilde_centroid + rho(i) * radii_vectors_xy_tilde[i];
    auto xyz_prime = sq.rotation_matrix * xy_tilde + sq.translation_vector;
    sq.sub_sqr_points[i] = xyz_prime.Normalized();
    sq.sub_sqr_weights[i] = weights[i];
  }
}

void
SLDFEsqQuadrature::DevelopSQLDFEValues(SphericalQuadrilateral& sq,
                                       GaussLegendreQuadrature& legendre)
{
  // Determine sq tilde center
  Vector3 sq_tilde_center;
  for (const auto& v : sq.vertices_xy_tilde)
    sq_tilde_center += v;
  sq_tilde_center /= 4;

  // Determine off-set vectors
  std::array<Vector3, 4> vctoi;
  for (int v = 0; v < 4; ++v)
    vctoi[v] = sq.vertices_xy_tilde[v] - sq_tilde_center;

  const auto sub_sub_square_xy_tilde = BuildSubSquares(sq.vertices_xy_tilde, sq_tilde_center);

  // Determine sub-sub-square xyz
  std::array<std::array<Vector3, 4>, 4> sub_sub_square_xyz;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      sub_sub_square_xyz[i][j] =
        (sq.rotation_matrix * sub_sub_square_xy_tilde[i][j] + sq.translation_vector).Normalized();

  // Compute sub-sub-square area
  std::array<double, 4> SA_i = {0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < 4; ++i)
    SA_i[i] = ComputeSphericalQuadrilateralArea(sub_sub_square_xyz[i]);

  // Apply optimization
  if (qp_optimization_type == QuadraturePointOptimization::CENTROID)
  {
    for (int i = 0; i < 4; ++i)
    {
      for (int j = 0; j < 4; ++j)
        sq.sub_sqr_points[i] += sub_sub_square_xyz[i][j];
      sq.sub_sqr_points[i] /= 4.0;
      sq.sub_sqr_points[i].Normalize();
      sq.sub_sqr_weights[i] = SA_i[i];
    }
  }
  else if (qp_optimization_type == QuadraturePointOptimization::EMPIRICAL)
    EmpiricalQPOptimization(sq, legendre, sq_tilde_center, vctoi, SA_i);
  else if (qp_optimization_type == QuadraturePointOptimization::ISOLATED)
    IsolatedQPOptimization(sq, legendre, sq_tilde_center, vctoi, SA_i);
}

double
SLDFEsqQuadrature::ComputeSphericalQuadrilateralArea(std::array<Vector3, 4>& vertices_xyz)
{
  const auto num_verts = 4;

  // Compute centroid
  Vector3 centroid_xyz;
  for (auto& v : vertices_xyz)
    centroid_xyz += v;
  centroid_xyz /= num_verts;

  // Compute area via spherical excess
  double area = 0.0;
  auto v0 = centroid_xyz.Normalized();
  for (int v = 0; v < num_verts; ++v)
  {
    auto v1 = vertices_xyz[v];
    auto v2 = vertices_xyz[(v < (num_verts - 1)) ? v + 1 : 0];

    if ((v1 - v0).Cross(v2 - v0).Dot(v0) < 0.0)
      std::swap(v1, v2);

    // Lambda for spherical excess
    auto GetSphericalExcess = [](const Vector3& vA, const Vector3& vB, const Vector3& vC)
    {
      const auto& n = vA;
      auto vAB = vB - vA;
      auto vAC = vC - vA;
      auto tAB = vAB.Cross(n).Normalized();
      auto tAC = vAC.Cross(n).Normalized();
      auto bAB = n.Cross(tAB).Normalized();
      auto bAC = n.Cross(tAC).Normalized();
      double mu = std::max(-1.0, std::fmin(1.0, bAB.Dot(bAC)));

      return std::fabs(acos(mu));
    };

    double excess = GetSphericalExcess(v0, v1, v2) + GetSphericalExcess(v1, v2, v0) +
                    GetSphericalExcess(v2, v0, v1);

    area += excess - M_PI;
  }

  return area;
}

std::array<double, 4>
SLDFEsqQuadrature::IntegrateLDFEShapeFunctions(const SphericalQuadrilateral& sq,
                                               std::array<Vector<double>, 4>& shape_coeffs,
                                               const std::vector<Vector3>& legendre_qpoints,
                                               const std::vector<double>& legendre_qweights)
{
  // Lambda to evaluate LDFE shape func
  auto EvaluateShapeFunction = [](Vector<double>& shape_coeffs, Vector3& mu_eta_xi)
  {
    return shape_coeffs(0) + shape_coeffs(1) * mu_eta_xi(0) + shape_coeffs(2) * mu_eta_xi(1) +
           shape_coeffs(3) * mu_eta_xi(2);
  };

  // Determine integration bounds
  double x_tilde_max = 0.0;
  double x_tilde_min = 1.0;
  double y_tilde_max = 0.0;
  double y_tilde_min = 1.0;

  for (const auto& v : sq.vertices_xy_tilde)
  {
    x_tilde_max = std::fmax(x_tilde_max, v.x);
    x_tilde_min = std::fmin(x_tilde_min, v.x);
    y_tilde_max = std::fmax(y_tilde_max, v.y);
    y_tilde_min = std::fmin(y_tilde_min, v.y);
  }

  // Integrate Legendre Quadrature
  std::array<double, 4> integral = {0.0, 0.0, 0.0, 0.0};
  auto Nq = legendre_qpoints.size();
  double dx_tilde = (x_tilde_max - x_tilde_min);
  double dy_tilde = (y_tilde_max - y_tilde_min);

  for (std::size_t i = 0; i < Nq; ++i)
  {
    for (std::size_t j = 0; j < Nq; ++j)
    {
      // Determine xy_tilde
      double x_tilde = x_tilde_min + (1.0 + legendre_qpoints[j][0]) * dx_tilde / 2.0;
      double y_tilde = y_tilde_min + (1.0 + legendre_qpoints[i][0]) * dy_tilde / 2.0;
      Vector3 xy_tilde(x_tilde, y_tilde, 0.0);

      // Map to xyz
      auto xyz = (sq.rotation_matrix * xy_tilde + sq.translation_vector).Normalized();

      // Determine Jacobian
      double r = std::sqrt(x_tilde * x_tilde + y_tilde * y_tilde +
                           SLDFEsqQuadrature::a_ * SLDFEsqQuadrature::a_);
      double detJ = (SLDFEsqQuadrature::a_ / (r * r * r)) * dx_tilde * dy_tilde / 4.0;

      // Evaluate shape funcs and add to integral
      for (int k = 0; k < 4; ++k)
        integral[k] += EvaluateShapeFunction(shape_coeffs[k], xyz) * detJ * legendre_qweights[i] *
                       legendre_qweights[j];
    }
  }

  return integral;
}

void
SLDFEsqQuadrature::CopyToAllOctants()
{
  deployed_SQs.clear();
  deployed_SQs.reserve(initial_octant_SQs_.size() * 8);

  // clang-format off
  const std::array<Vector3, 8> octant_modifiers = {Vector3( 1.0,  1.0,  1.0),
                                                   Vector3(-1.0,  1.0,  1.0),
                                                   Vector3(-1.0, -1.0,  1.0),
                                                   Vector3( 1.0, -1.0,  1.0),
                                                   Vector3( 1.0,  1.0, -1.0),
                                                   Vector3(-1.0,  1.0, -1.0),
                                                   Vector3(-1.0, -1.0, -1.0),
                                                   Vector3( 1.0, -1.0, -1.0)};
  // clang-format on

  auto ApplyModifier = [](SphericalQuadrilateral sq, const Vector3& modifier)
  {
    for (auto& xyz : sq.vertices_xyz)
      xyz = xyz * modifier;
    sq.centroid_xyz = sq.centroid_xyz * modifier;
    for (auto& xyz : sq.sub_sqr_points)
      xyz = xyz * modifier;
    sq.octant_modifier = modifier;
    return sq;
  };

  for (const auto& modifier : octant_modifiers)
    for (const auto& sq : initial_octant_SQs_)
      deployed_SQs.push_back(ApplyModifier(sq, modifier));
}

void
SLDFEsqQuadrature::PopulateQuadratureAbscissae()
{
  abscissae.clear();
  weights.clear();
  omegas.clear();
  raw_weight_sum_ = 0.0;

  for (const auto& sq : deployed_SQs)
  {
    for (int i = 0; i < 4; ++i)
    {
      const auto& omega = sq.sub_sqr_points[i];
      const double weight = sq.sub_sqr_weights[i];
      const double theta = std::acos(omega.z);
      double phi = std::acos(omega.x / std::sin(theta));

      if (omega.y / std::sin(theta) < 0.0)
        phi = 2.0 * M_PI - phi;

      if (not FilterQuadraturePoint(omega))
        continue;

      raw_weight_sum_ += weight;
      abscissae.emplace_back(phi, theta);
      weights.push_back(weight);
      omegas.push_back(omega);
    }
  }

  // Normalize weights to 1.0
  for (auto& w : weights)
    w /= raw_weight_sum_;
}

bool
SLDFEsqQuadrature::FilterQuadraturePoint(const Vector3& omega) const
{
  return true;
}

SLDFEsqQuadrature3DXYZ::SLDFEsqQuadrature3DXYZ(int level, unsigned int scattering_order)
  : SLDFEsqQuadrature(level, 3, scattering_order)
{
  GenerateRefinement();
}

SLDFEsqQuadrature2DXY::SLDFEsqQuadrature2DXY(int level, unsigned int scattering_order)
  : SLDFEsqQuadrature(level, 2, scattering_order)
{
  GenerateRefinement();
}

bool
SLDFEsqQuadrature2DXY::FilterQuadraturePoint(const Vector3& omega) const
{
  return omega.z > 0.0;
}

double
SLDFEsqQuadrature::RiemannIntegral(BaseFunctor* F, int Ni)
{
  double dangle = M_PI_2 / Ni;
  double dtheta = dangle;
  double dphi = dangle;

  double I_riemann = 0.0;
  for (int i = 0; i < Ni; ++i)
  {
    double theta = (0.5 + i) * dtheta;
    for (int j = 0; j < Ni; ++j)
    {
      double phi = (0.5 + j) * dphi;
      double mu_r = cos(phi) * sin(theta);
      double eta_r = sin(phi) * sin(theta);
      double xi_r = cos(theta);
      double fval = (*F)(mu_r, eta_r, xi_r);
      I_riemann += fval * sin(theta) * dtheta * dphi;
    }
  }

  return I_riemann;
}

double
SLDFEsqQuadrature::QuadratureSSIntegral(BaseFunctor* F)
{
  double I_quadrature = 0.0;
  for (const auto& sq : initial_octant_SQs_)
  {
    for (int i = 0; i < 4; ++i)
    {
      double mu = sq.sub_sqr_points[i][2];
      double theta = acos(mu);
      double phi = acos(sq.sub_sqr_points[i][0] / sin(theta));
      double mu_r = cos(phi) * sin(theta);
      double eta_r = sin(phi) * sin(theta);
      double xi_r = cos(theta);
      double fval = (*F)(mu_r, eta_r, xi_r);
      I_quadrature += sq.sub_sqr_weights[i] * fval;
    }
  }

  return I_quadrature;
}

void
SLDFEsqQuadrature::TestIntegration(int test_case, double ref_solution, int RiemannN)
{
  struct Case1 : public BaseFunctor
  {
    double operator()(double mu, double eta, double xi) override
    {
      return pow(mu, 1.0) * pow(eta, 1.0) * pow(xi, 0.0);
    }
  };

  struct Case2 : public BaseFunctor
  {
    double operator()(double mu, double eta, double xi) override
    {
      return pow(mu, 3.0) * pow(eta, 1.0) * pow(xi, 1.0);
    }
  };

  struct Case3 : public BaseFunctor
  {
    double operator()(double mu, double eta, double xi) override
    {
      return pow(mu, 3.0) * pow(eta, 6.0) * pow(xi, 15.0);
    }
  };

  struct SphericalHarmonicF : public BaseFunctor
  {
    double operator()(double mu, double eta, double xi) override
    {
      double theta = acos(xi);
      double phi = acos(mu / sin(theta));
      return Ylm(15, 3, phi, theta);
    }
  };

  Case1 case1;
  Case2 case2;
  Case3 case3;
  SphericalHarmonicF SphF;
  BaseFunctor* F = nullptr;
  switch (test_case)
  {
    case 1:
      F = &case1;
      break;
    case 2:
      F = &case2;
      break;
    case 3:
      F = &case3;
      break;
    case 4:
      F = &SphF;
      break;
    default:
      F = &case1;
  }

  const auto Nd = initial_octant_SQs_.size() * 4;
  const int NR = RiemannN;
  double h = 1.0 / std::sqrt(8.0 * static_cast<double>(Nd));
  double I_riemann = ref_solution;
  if (NR > 0)
    I_riemann = std::fabs(RiemannIntegral(F, NR));
  double I_quadrature = std::fabs(QuadratureSSIntegral(F));

  log.Log() << "Riemann integral: " << std::setprecision(20) << std::scientific << I_riemann;
  log.Log() << "Quadrature integral: " << std::setprecision(10) << std::scientific << I_quadrature;
  log.Log() << "Error_RQ" << std::setw(5) << std::setfill('0') << Nd << '_' << std::setw(6)
            << std::setfill('0') << Nd * 8 << ": " << std::setw(2) << std::setfill(' ') << level_
            << ' ' << h << ' ' << std::scientific
            << std::fabs((I_riemann - I_quadrature) / ref_solution);
}

void
SLDFEsqQuadrature::PrintQuadratureToFile(const std::string& file_base)
{
  log.Log() << "Writing SLDFEsq quadrature to file";

  std::ofstream vert_file, cell_file, points_file;

  // Vertex file generation for plotting - is written to the path above
  vert_file.open(file_base + "_verts.csv");
  {
    vert_file << "x,y,z\n";
    for (const auto& sq : deployed_SQs)
    {
      for (int v = 0; v < 4; ++v)
      {
        const auto& v0 = sq.vertices_xyz_prime[v];
        const auto& v1 = sq.vertices_xyz_prime[(v < 3) ? v + 1 : 0];

        for (int d = 0; d <= 10; ++d)
        {
          auto vert = (1.0 - d / 10.0) * v0 + (d / 10.0) * v1;
          vert = vert * sq.octant_modifier;
          vert.Normalize();
          vert_file << vert.x << "," << vert.y << "," << vert.z << "\n";
        }
      }
    }
  }
  vert_file.close();

  // Indexing file for polygons for plotting
  cell_file.open(file_base + "_cells.csv");
  {
    cell_file << "Cell Index\n";

    int vi = 0;
    for (const auto& sq : deployed_SQs)
    {
      for (const auto& vert : sq.vertices_xyz)
        for (int d = 0; d <= 10; ++d)
          cell_file << vi++ << ",";
      cell_file << "\n";
    }
  }
  cell_file.close();

  // Formatted cell index file for each polygon
  points_file.open(file_base + "_points.csv");
  {
    points_file << "x,y,z,weights\n";

    const bool has_weights = weights.size() == omegas.size();
    for (size_t i = 0; i < omegas.size(); ++i)
    {
      const auto& point = omegas[i];
      const double normalized_weight = has_weights ? weights[i] : 0.0;
      const double raw_weight = normalized_weight * raw_weight_sum_;
      points_file << point[0] << "," << point[1] << "," << point[2] << "," << raw_weight << "\n";
    }
  }
  points_file.close();
}

std::array<SLDFEsqQuadrature::SphericalQuadrilateral, 4>
SLDFEsqQuadrature::SplitSQ(SphericalQuadrilateral& sq, GaussLegendreQuadrature& legendre)
{
  std::array<SphericalQuadrilateral, 4> new_sqs;

  // Determine sq tilde center
  Vector3 sq_tilde_center;
  for (const auto& v : sq.vertices_xy_tilde)
    sq_tilde_center += v;
  sq_tilde_center /= 4;

  const auto sub_sub_square_xy_tilde = BuildSubSquares(sq.vertices_xy_tilde, sq_tilde_center);

  // Compute SQ xyz-centroid, R,T,area, ldfe
  for (int i = 0; i < 4; ++i)
  {
    auto& child = new_sqs[i];
    child.vertices_xy_tilde = sub_sub_square_xy_tilde[i];
    child.rotation_matrix = sq.rotation_matrix;
    child.translation_vector = sq.translation_vector;
    child.octant_modifier = sq.octant_modifier;

    child.centroid_xyz = Vector3();
    for (int v = 0; v < 4; ++v)
    {
      child.vertices_xyz_prime[v] =
        child.rotation_matrix * child.vertices_xy_tilde[v] + child.translation_vector;
      child.vertices_xyz[v] = child.vertices_xyz_prime[v].Normalized();
      child.centroid_xyz += child.vertices_xyz[v];
    }

    child.centroid_xyz /= 4;
    child.centroid_xyz = child.centroid_xyz.Normalized() * sq.octant_modifier;
    child.area = ComputeSphericalQuadrilateralArea(child.vertices_xyz);
    DevelopSQLDFEValues(child, legendre);

    for (int v = 0; v < 4; ++v)
    {
      child.vertices_xyz[v] = child.vertices_xyz[v] * sq.octant_modifier;
      child.sub_sqr_points[v] = child.sub_sqr_points[v] * sq.octant_modifier;
    }
  }

  return new_sqs;
}

void
SLDFEsqQuadrature::LocallyRefine(const Vector3& ref_dir,
                                 const double cone_size,
                                 const bool dir_as_plane_normal)
{
  auto ref_dir_n = ref_dir.Normalized();
  double mu_cone = cos(cone_size);
  std::vector<SphericalQuadrilateral> new_deployment;
  new_deployment.reserve(deployed_SQs.size());

  GaussLegendreQuadrature legendre(QuadratureOrder::THIRTYSECOND);

  int num_refined = 0;
  for (auto& sq : deployed_SQs)
  {
    bool sq_to_be_split = false;

    if (not dir_as_plane_normal)
      sq_to_be_split = sq.centroid_xyz.Dot(ref_dir_n) > mu_cone;
    else
      sq_to_be_split = std::fabs(sq.centroid_xyz.Dot(ref_dir_n)) < (sin(cone_size));

    if (not sq_to_be_split)
      new_deployment.push_back(sq);
    else
    {
      auto new_sqs = SplitSQ(sq, legendre);
      for (auto& nsq : new_sqs)
        new_deployment.push_back(nsq);
      ++num_refined;
    }
  }

  deployed_SQs.clear();
  deployed_SQs = new_deployment;

  PopulateQuadratureAbscissae();
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  log.Log() << "SLDFEsq refined " << num_refined << " SQs.";
}

} // namespace opensn
