// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/mesh/mesh_vector.h"
#include "framework/runtime.h"
#include <algorithm>
#include <map>

namespace opensn
{

void
SimplifiedLDFESQ::Quadrature::GenerateInitialRefinement(int level)
{
  Timer timer;
  timer.Reset();
  initial_level_ = level;

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

  // Set translation vectors for cube faces
  auto txface = a * ihat;
  auto tyface = a * jhat;
  auto tzface = a * khat;

  // Generate general diagonal spacings in xy-tilde coordinates
  GenerateDiagonalSpacings(level);

  // Generate vertices for each face of inscribed cube
  GenerateReferenceFaceVertices(Rxface, txface, level);
  GenerateReferenceFaceVertices(Ryface, tyface, level);
  GenerateReferenceFaceVertices(Rzface, tzface, level);

  // Compute areas
  double total_area = 0.0;
  double area_max = -100.0;
  double area_min = 100.0;
  bool negative_weights_found = false;
  for (auto& sq : initial_octant_SQs_)
  {
    double area = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      area += sq.sub_sqr_weights[i];
      if (area < 0.0)
        negative_weights_found = true;
    }
    total_area += area;
    area_max = std::fmax(area_max, area);
    area_min = std::fmin(area_min, area);
  }
  double area_avg = total_area / initial_octant_SQs_.size();

  if (negative_weights_found)
    log.Log0Warning() << "SLDFESQ Quadrature detected negative weights.";

  // Print Statistics
  double time = timer.GetTime() / 1000.0;
  log.Log0Verbose1() << "Number of dirs/octant: " << initial_octant_SQs_.size();
  log.Log0Verbose1() << "Total weight         : " << total_area;
  log.Log0Verbose1() << "Total weight/(pi/2)  : " << total_area / M_PI_2;
  log.Log0Verbose1() << "Area Max/Min         : " << area_max / area_min;
  log.Log0Verbose1() << "Area Max/Avg         : " << area_max / area_avg;

  CopyToAllOctants();

  // Populate quadriture points
  PopulateQuadratureAbscissae();

  log.Log0Verbose1() << "Time taken           : " << time;
}

void
SimplifiedLDFESQ::Quadrature::GenerateDiagonalSpacings(int level)
{
  // Define constants
  const int Ns = (level + 1); // Number of subdivisions
  const int Np = Ns + 1;      // Number of diagonal points

  const auto ihat = Vector3(1.0, 0.0, 0.0);

  // Define rotation matrix
  Matrix3x3 Rihat;
  auto n = Vector3(0.0, -1.0 / sqrt(2), 1.0 / sqrt(2));
  auto& t = ihat;
  auto b = n.Cross(t).Normalized();

  Rihat.SetColJVec(0, t);
  Rihat.SetColJVec(1, b);
  Rihat.SetColJVec(2, n);

  // Generate p-points
  std::vector<Vector3> p_points(Np);
  double dphi = acos(a) / Ns;
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
    double R = a / p_points[i].x;
    double x_tilde = p_points[i].y * R;
    double y_tilde = p_points[i].z * R;

    tilde_points[i] = Vector3(x_tilde, y_tilde, 0.0);
  }

  diagonal_vertices_ = tilde_points;
}

void
SimplifiedLDFESQ::Quadrature::GenerateReferenceFaceVertices(const Matrix3x3& rotation_matrix,
                                                            const Vector3& translation,
                                                            int level)
{
  int Ns = (level + 1); // Number of subdivisions
  int Np = Ns + 1;      // Number of diagonal points

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
    } // for j
  }   // for i
}

void
SimplifiedLDFESQ::Quadrature::EmpiricalQPOptimization(
  SphericalQuadrilateral& sq,
  GaussLegendreQuadrature& legendre,
  Vector3& sq_xy_tilde_centroid,
  std::array<Vector3, 4>& radii_vectors_xy_tilde,
  std::array<double, 4>& sub_sub_sqr_areas)
{
  FunctionWeightFromRho ComputeWeights(
    *this, sq_xy_tilde_centroid, radii_vectors_xy_tilde, sq, legendre);
  double d = 1.0 / sqrt(3.0);
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
SimplifiedLDFESQ::Quadrature::IsolatedQPOptimization(SphericalQuadrilateral& sq,
                                                     GaussLegendreQuadrature& legendre,
                                                     Vector3& sq_xy_tilde_centroid,
                                                     std::array<Vector3, 4>& radii_vectors_xy_tilde,
                                                     std::array<double, 4>& sub_sub_sqr_areas)
{
  auto& SA_i = sub_sub_sqr_areas;

  // Declare algorithm utilities
  FunctionWeightFromRho ComputeWeights(
    *this, sq_xy_tilde_centroid, radii_vectors_xy_tilde, sq, legendre);
  double d = 1.0 / sqrt(3.0);
  Vector<double> rho({d, d, d, d});
  double epsilon = 1.0e-1;
  Vector<double> delta({epsilon, epsilon, epsilon, epsilon});
  Vector<double> drho_df({0.0, 0.0, 0.0, 0.0});

  // Compute initial weights
  auto weights = ComputeWeights(rho);

  // Apply algorithm
  log.Log() << "=================================================== ";
  for (int k = 0; k < 150; ++k) // iteration
  {
    //    constexpr int N = 4;
    //    double fac = 1.0/N;
    //    std::array<std::array<double,4>,N> weights_offset;

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

      //      delta = {0.0,0.0,0.0,0.0}; delta[i] = epsilon;
      //      auto weights_offset_pos = ComputeWeights(rho + delta);
      //      double slope = (weights_offset_pos[i]-weights[i])/epsilon;
      //
      //      double delta_rho = 10.0*slope*(SA_i[i]-weights[i]);

      rho(i) += delta_rho;
      rho(i) = std::fmax(0.0, rho(i));
      rho(i) = std::fmin(1.0, rho(i));
      rho_change_total -= 1.0 * drho_df(i) * (weights[i] - SA_i[i]);
    }

    // Update weights
    weights = ComputeWeights(rho);
    double change = 0.0;
    for (int i = 0; i < 4; ++i)
      change = std::fabs((weights[i] - SA_i[i]) / weights[i]);

    log.Log() << "Weights: " << weights[0] << " " << weights[1] << " " << weights[2] << " "
              << weights[3] << " ";
    log.Log() << "Areas: " << SA_i[0] << " " << SA_i[1] << " " << SA_i[2] << " " << SA_i[3] << "\n";
    log.Log() << "rhos: " << rho(0) << " " << rho(1) << " " << rho(2) << " " << rho(3) << "\n";
    log.Log() << k << " " << std::fabs(change);
    log.Log() << "  ";

    if (rho_change_total < 1.0e-2)
      break;
    //    if (std::fabs(change) < 1.0e-2) break;
  }
  //  log.Log() << "rhos: "
  //                     << rho[0]/(1.0/sqrt(3.0)) << " "
  //                     << rho[1]/(1.0/sqrt(3.0)) << " "
  //                     << rho[2]/(1.0/sqrt(3.0)) << " "
  //                     << rho[3]/(1.0/sqrt(3.0)) << "\n";

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
SimplifiedLDFESQ::Quadrature::DevelopSQLDFEValues(SphericalQuadrilateral& sq,
                                                  GaussLegendreQuadrature& legendre)
{
  // Determine sq tilde center
  Vector3 sq_tilde_center;
  for (const auto& v : sq.vertices_xy_tilde)
    sq_tilde_center += v;
  sq_tilde_center /= 4;

  // Determine off-set vectors
  auto& vc = sq_tilde_center;
  std::array<Vector3, 4> vctoi;
  for (int v = 0; v < 4; ++v)
    vctoi[v] = sq.vertices_xy_tilde[v] - vc;

  // Determine sub-sub-squares
  std::array<std::array<Vector3, 4>, 4> sub_sub_square_xy_tilde;
  std::map<std::string, Vector3> vm;

  for (int v = 0; v < 4; ++v)
    vm[std::to_string(v)] = sq.vertices_xy_tilde[v];

  vm["01"] = 0.5 * (sq.vertices_xy_tilde[0] + sq.vertices_xy_tilde[1]);
  vm["12"] = 0.5 * (sq.vertices_xy_tilde[1] + sq.vertices_xy_tilde[2]);
  vm["23"] = 0.5 * (sq.vertices_xy_tilde[2] + sq.vertices_xy_tilde[3]);
  vm["03"] = 0.5 * (sq.vertices_xy_tilde[0] + sq.vertices_xy_tilde[3]);
  vm["c"] = sq_tilde_center;

  auto& sst = sub_sub_square_xy_tilde;
  sst[0] = {vm["0"], vm["01"], vm["c"], vm["03"]};
  sst[1] = {vm["01"], vm["1"], vm["12"], vm["c"]};
  sst[2] = {vm["c"], vm["12"], vm["2"], vm["23"]};
  sst[3] = {vm["03"], vm["c"], vm["23"], vm["3"]};

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
    for (int i = 0; i < 4; ++i)
    {
      for (int j = 0; j < 4; ++j)
        sq.sub_sqr_points[i] += sub_sub_square_xyz[i][j];
      sq.sub_sqr_points[i] /= 4.0;
      sq.sub_sqr_points[i].Normalize();

      sq.sub_sqr_weights[i] = SA_i[i];
    } // for i
  else if (qp_optimization_type == QuadraturePointOptimization::EMPIRICAL)
    EmpiricalQPOptimization(sq, legendre, vc, vctoi, SA_i);
  else if (qp_optimization_type == QuadraturePointOptimization::ISOLATED)
    IsolatedQPOptimization(sq, legendre, vc, vctoi, SA_i);
}

double
SimplifiedLDFESQ::Quadrature::ComputeSphericalQuadrilateralArea(
  std::array<Vector3, 4>& vertices_xyz)
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
SimplifiedLDFESQ::Quadrature::IntegrateLDFEShapeFunctions(
  const SphericalQuadrilateral& sq,
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

  for (auto& v : sq.vertices_xy_tilde)
  {
    x_tilde_max = std::fmax(x_tilde_max, v.x);
    x_tilde_min = std::fmin(x_tilde_min, v.x);
    y_tilde_max = std::fmax(y_tilde_max, v.y);
    y_tilde_min = std::fmin(y_tilde_min, v.y);
  }

  // Integrate Legendre Quadrature
  std::array<double, 4> integral = {0.0, 0.0, 0.0, 0.0};
  int Nq = legendre_qpoints.size();
  double dx_tilde = (x_tilde_max - x_tilde_min);
  double dy_tilde = (y_tilde_max - y_tilde_min);

  for (int i = 0; i < Nq; ++i)
  {
    for (int j = 0; j < Nq; ++j)
    {
      // Determine xy_tilde
      double x_tilde = x_tilde_min + (1.0 + legendre_qpoints[j][0]) * dx_tilde / 2.0;
      double y_tilde = y_tilde_min + (1.0 + legendre_qpoints[i][0]) * dy_tilde / 2.0;
      Vector3 xy_tilde(x_tilde, y_tilde, 0.0);

      // Map to xyz
      auto xyz = (sq.rotation_matrix * xy_tilde + sq.translation_vector).Normalized();

      // Determine Jacobian
      double r = sqrt(x_tilde * x_tilde + y_tilde * y_tilde + a * a);
      double detJ = (a / (r * r * r)) * dx_tilde * dy_tilde / 4.0;

      // Evaluate shape funcs and add to integral
      for (int k = 0; k < 4; ++k)
        integral[k] += EvaluateShapeFunction(shape_coeffs[k], xyz) * detJ * legendre_qweights[i] *
                       legendre_qweights[j];
    } // for j
  }   // for i

  return integral;
}

void
SimplifiedLDFESQ::Quadrature::CopyToAllOctants()
{
  deployed_SQs.clear(); // just to be sure
  deployed_SQs.reserve(initial_octant_SQs_.size() * 8);

  // Define modifying variables
  Vector3 octant_mod(1.0, 1.0, 1.0);

  // Top NE octant, no change
  for (auto& sq : initial_octant_SQs_)
    deployed_SQs.push_back(sq);

  // Top NW octant
  octant_mod = Vector3(-1.0, 1.0, 1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)
      xyz = xyz * octant_mod;
    auto& vcc = new_sq.centroid_xyz;
    vcc = vcc * octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points)
      xyz = xyz * octant_mod;
    new_sq.octant_modifier = octant_mod;

    deployed_SQs.push_back(new_sq);
  }

  // Top SW octant
  octant_mod = Vector3(-1.0, -1.0, 1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)
      xyz = xyz * octant_mod;
    auto& vcc = new_sq.centroid_xyz;
    vcc = vcc * octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points)
      xyz = xyz * octant_mod;
    new_sq.octant_modifier = octant_mod;

    deployed_SQs.push_back(new_sq);
  }

  // Top SE octant
  octant_mod = Vector3(1.0, -1.0, 1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)
      xyz = xyz * octant_mod;
    auto& vcc = new_sq.centroid_xyz;
    vcc = vcc * octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points)
      xyz = xyz * octant_mod;
    new_sq.octant_modifier = octant_mod;

    deployed_SQs.push_back(new_sq);
  }

  // Bot NE octant
  octant_mod = Vector3(1.0, 1.0, -1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)
      xyz = xyz * octant_mod;
    auto& vcc = new_sq.centroid_xyz;
    vcc = vcc * octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points)
      xyz = xyz * octant_mod;
    new_sq.octant_modifier = octant_mod;

    deployed_SQs.push_back(new_sq);
  }

  // Bot NW octant
  octant_mod = Vector3(-1.0, 1.0, -1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)
      xyz = xyz * octant_mod;
    auto& vcc = new_sq.centroid_xyz;
    vcc = vcc * octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points)
      xyz = xyz * octant_mod;
    new_sq.octant_modifier = octant_mod;

    deployed_SQs.push_back(new_sq);
  }

  // Bot SW octant
  octant_mod = Vector3(-1.0, -1.0, -1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)
      xyz = xyz * octant_mod;
    auto& vcc = new_sq.centroid_xyz;
    vcc = vcc * octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points)
      xyz = xyz * octant_mod;
    new_sq.octant_modifier = octant_mod;

    deployed_SQs.push_back(new_sq);
  }

  // Bot SE octant
  octant_mod = Vector3(1.0, -1.0, -1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)
      xyz = xyz * octant_mod;
    auto& vcc = new_sq.centroid_xyz;
    vcc = vcc * octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points)
      xyz = xyz * octant_mod;
    new_sq.octant_modifier = octant_mod;

    deployed_SQs.push_back(new_sq);
  }

  // Make history entry
  deployed_SQs_history_.push_back(deployed_SQs);
}

void
SimplifiedLDFESQ::Quadrature::PopulateQuadratureAbscissae()
{
  abscissae.clear();
  weights.clear();
  omegas.clear();

  for (const auto& sq : deployed_SQs)
  {
    for (int i = 0; i < 4; ++i)
    {
      const auto& omega = sq.sub_sqr_points[i];
      const double weight = sq.sub_sqr_weights[i];

      double theta = acos(omega.z);
      double phi = acos(omega.x / sin(theta));

      if (omega.y / sin(theta) < 0.0)
        phi = 2.0 * M_PI - phi;

      const auto abscissa = QuadraturePointPhiTheta(phi, theta);

      abscissae.push_back(abscissa);
      weights.push_back(weight);
      omegas.push_back(omega);
    }
  }
}

double
SimplifiedLDFESQ::Quadrature::RiemannIntegral(BaseFunctor* F, int Ni)
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
    } // for j
  }   // for i

  return I_riemann;
}

double
SimplifiedLDFESQ::Quadrature::QuadratureSSIntegral(BaseFunctor* F)
{
  double I_quadrature = 0.0;
  for (const auto& sq : initial_octant_SQs_)
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

  return I_quadrature;
}

void
SimplifiedLDFESQ::Quadrature::TestIntegration(int test_case, double ref_solution, int RiemannN)
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
  BaseFunctor* F = &case1;
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

  const int Nd = initial_octant_SQs_.size() * 4;
  const int NR = RiemannN;

  double h = 1.0 / sqrt(8.0 * Nd);
  double I_riemann = ref_solution;
  if (NR > 0)
    I_riemann = std::fabs(RiemannIntegral(F, NR));

  double I_quadrature = std::fabs(QuadratureSSIntegral(F));

  char buff0[200], buff1[200], buff2[200];
  snprintf(buff0, 200, "Riemann integral: %.20e\n", I_riemann);
  snprintf(buff1, 200, "Quadrature integral: %.10e\n", I_quadrature);
  snprintf(buff2,
           200,
           "Error_RQ%05d_%06d: %2d %f %e\n",
           Nd,
           Nd * 8,
           initial_level_,
           h,
           std::fabs((I_riemann - I_quadrature) / ref_solution));

  log.Log() << buff0;
  log.Log() << buff1;
  log.Log() << buff2;
}

void
SimplifiedLDFESQ::Quadrature::PrintQuadratureToFile(const std::string& file_base)
{
  log.Log() << "Printing SLDFE-Quadrature to file.";

  std::ofstream vert_file, cell_file, points_file, python_file;
  vert_file.open(file_base + "verts.txt");
  {
    for (const auto& sq : deployed_SQs)
      for (int v = 0; v < 4; ++v)
      {
        auto& v0 = sq.vertices_xyz_prime[v];
        auto& v1 = sq.vertices_xyz_prime[(v < 3) ? v + 1 : 0];

        for (int d = 0; d <= 10; ++d)
        {
          auto vert = (1.0 - d / 10.0) * v0 + (d / 10.0) * v1;
          vert = vert * sq.octant_modifier;
          vert.Normalize();
          // std::cout << "Verts " << vert.x << " " << vert.y << " " << vert.z << std::endl;
          vert_file << vert.x << " " << vert.y << " " << vert.z << "\n";
        }
      }
  }
  vert_file.close();

  cell_file.open(file_base + "cells.txt");
  {
    int vi = 0;
    for (const auto& sq : deployed_SQs)
    {
      for (const auto& vert : sq.vertices_xyz)
      {
        for (int d = 0; d <= 10; ++d)
          cell_file << vi++ << " ";
      }
      cell_file << "\n";
    }
  }
  cell_file.close();

  points_file.open(file_base + "points.txt");
  {
    for (auto& sq : deployed_SQs)
    {
      int ss = -1;
      for (const auto& point : sq.sub_sqr_points)
      {
        ++ss;
        for (int i = 0; i < 3; ++i)
          points_file << point[i] << " ";
        points_file << sq.sub_sqr_weights[ss];
        points_file << "\n";
      }
    }
  }
  points_file.close();

  python_file.open(file_base + "python.py");
  python_file << "import matplotlib.pyplot as plt\n"
                 "from mpl_toolkits import mplot3d\n"
                 "import mpl_toolkits.mplot3d.art3d as art3d\n"
                 "import mpl_toolkits.mplot3d as ax3\n"
                 "import matplotlib.transforms as mpltransform\n"
                 "\n"
                 "import numpy as np\n"
                 "import math\n"
                 "\n"
                 "#====================================== Read vertices\n"
                 "verts = []\n"
                 "verts_file = open(\""
              << file_base
              << "verts.txt\")\n"
                 "for line in verts_file:\n"
                 "    words = line.split()\n"
                 "    "
                 "verts.append(np.array([float(words[0]),float(words[1]),float("
                 "words[2])]))\n"
                 "verts_file.close()\n"
                 "\n"
                 "#====================================== Read cells\n"
                 "cells = []\n"
                 "cells_file = open(\""
              << file_base
              << "cells.txt\")\n"
                 "for line in cells_file:\n"
                 "    words = line.split()\n"
                 "    cell = []\n"
                 "    for word in words:\n"
                 "        cell.append(int(word))\n"
                 "    cells.append(cell)\n"
                 "cells_file.close()\n"
                 "\n"
                 "#====================================== Read points\n"
                 "points = []\n"
                 "weightsum=0.0\n"
                 "points_file = open(\""
              << file_base
              << "points.txt\")\n"
                 "for line in points_file:\n"
                 "    words = line.split()\n"
                 "    point = []\n"
                 "    for word in words:\n"
                 "        point.append(float(word))\n"
                 "    points.append(point)\n"
                 "    weightsum += point[3]\n"
                 "points_file.close()\n"
                 "\n"
                 "print(\"Weightsum check: \",weightsum,weightsum/4/math.pi)\n"
                 "\n"
                 "points_array = np.array(points)\n"
                 "\n"
                 "#====================================== Generate polygons\n"
                 "patches = []\n"
                 "for cell in cells:\n"
                 "\n"
                 "    vertex_list = []\n"
                 "    for index in cell:\n"
                 "        vertex_list.append(verts[index])\n"
                 "\n"
                 "    polygon = art3d.Poly3DCollection([vertex_list])\n"
                 "    polygon.set_color([1.0,1.0,1.0,1.0])\n"
                 "    polygon.set_edgecolor([0.0,0.0,0.0,1.0])\n"
                 "    patches.append(polygon)\n"
                 "\n"
                 "#====================================== Plot polygons\n"
                 "fig = plt.figure(figsize=(10,8.5))\n"
                 "ax = fig.add_subplot(111, projection='3d')\n"
                 "\n"
                 "ax.view_init(20,45)\n"
                 "limit = 1\n"
                 "\n"
                 "for poly in patches:\n"
                 "    ax.add_collection3d(poly)\n"
                 "\n"
                 "if limit==8:\n"
                 "    ax.set_xlim([0.0,1.0])\n"
                 "    ax.set_ylim([0.0,1.0])\n"
                 "    ax.set_zlim([0.0,1.0])\n"
                 "else:\n"
                 "    ax.set_xlim([-1.0,1.0])\n"
                 "    ax.set_ylim([-1.0,1.0])\n"
                 "    ax.set_zlim([-1.0,1.0])\n"
                 "\n"
                 "ax.margins(0.5)\n"
                 "# ax.set_xlabel(r\"$\\mu$\")\n"
                 "# ax.set_ylabel(r\"$\\eta$\")\n"
                 "# ax.set_zlabel(r\"$\\xi$\")\n"
                 "ax.set_xlabel(r\"$x$\")\n"
                 "ax.set_ylabel(r\"$y$\")\n"
                 "ax.set_zlabel(r\"$z$\")\n"
                 "plt.savefig(\"CQuad.png\")\n"
                 "plt.show()\n";
  python_file.close();

  log.Log() << "Done printing SLDFE-Quadrature to file.";
}

std::array<SimplifiedLDFESQ::SphericalQuadrilateral, 4>
SimplifiedLDFESQ::Quadrature::SplitSQ(SphericalQuadrilateral& sq, GaussLegendreQuadrature& legendre)
{
  std::array<SphericalQuadrilateral, 4> new_sqs;

  // Determine sq tilde center
  Vector3 sq_tilde_center;
  for (const auto& v : sq.vertices_xy_tilde)
    sq_tilde_center += v;
  sq_tilde_center /= 4;

  // Determine sub-sub-squares Tilde coordinates
  std::array<std::array<Vector3, 4>, 4> sub_sub_square_xy_tilde;
  std::map<std::string, Vector3> vm;

  for (int v = 0; v < 4; ++v)
    vm[std::to_string(v)] = sq.vertices_xy_tilde[v];

  vm["01"] = 0.5 * (sq.vertices_xy_tilde[0] + sq.vertices_xy_tilde[1]);
  vm["12"] = 0.5 * (sq.vertices_xy_tilde[1] + sq.vertices_xy_tilde[2]);
  vm["23"] = 0.5 * (sq.vertices_xy_tilde[2] + sq.vertices_xy_tilde[3]);
  vm["03"] = 0.5 * (sq.vertices_xy_tilde[0] + sq.vertices_xy_tilde[3]);
  vm["c"] = sq_tilde_center;

  auto& sst = sub_sub_square_xy_tilde;
  sst[0] = {vm["0"], vm["01"], vm["c"], vm["03"]};
  sst[1] = {vm["01"], vm["1"], vm["12"], vm["c"]};
  sst[2] = {vm["c"], vm["12"], vm["2"], vm["23"]};
  sst[3] = {vm["03"], vm["c"], vm["23"], vm["3"]};

  for (int i = 0; i < 4; ++i)
    new_sqs[i].vertices_xy_tilde = sst[i];

  // Determine xyz-prime
  for (int i = 0; i < 4; ++i)
    for (int v = 0; v < 4; ++v)
      new_sqs[i].vertices_xyz_prime[v] = sq.rotation_matrix * sst[i][v] + sq.translation_vector;

  // Compute xyz
  for (int i = 0; i < 4; ++i)
    for (int v = 0; v < 4; ++v)
      new_sqs[i].vertices_xyz[v] = new_sqs[i].vertices_xyz_prime[v].Normalized();

  // Compute SQ xyz-centroid, R,T,area, ldfe
  for (int i = 0; i < 4; ++i)
  {
    for (int v = 0; v < 4; ++v)
      new_sqs[i].centroid_xyz += new_sqs[i].vertices_xyz[v];
    new_sqs[i].centroid_xyz /= 4;
    new_sqs[i].centroid_xyz = new_sqs[i].centroid_xyz.Normalized() * sq.octant_modifier;

    new_sqs[i].rotation_matrix = sq.rotation_matrix;
    new_sqs[i].translation_vector = sq.translation_vector;

    new_sqs[i].area = ComputeSphericalQuadrilateralArea(new_sqs[i].vertices_xyz);
    DevelopSQLDFEValues(new_sqs[i], legendre);
    new_sqs[i].octant_modifier = sq.octant_modifier;

    for (int v = 0; v < 4; ++v)
    {
      new_sqs[i].vertices_xyz[v] = new_sqs[i].vertices_xyz[v] * sq.octant_modifier;
      new_sqs[i].sub_sqr_points[v] = new_sqs[i].sub_sqr_points[v] * sq.octant_modifier;
    }
  }

  return new_sqs;
}

void
SimplifiedLDFESQ::Quadrature::LocallyRefine(const Vector3& ref_dir,
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
  deployed_SQs_history_.push_back(new_deployment);

  PopulateQuadratureAbscissae();

  log.Log() << "SLDFESQ refined " << num_refined << " SQs.";
}

} // namespace opensn
