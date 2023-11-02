#include "sldfe_sq.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/***/
void
chi_math::SimplifiedLDFESQ::Quadrature::IsolatedQPOptimization(
  SphericalQuadrilateral& sq,
  chi_math::QuadratureGaussLegendre& legendre,
  chi_mesh::Vertex& sq_xy_tilde_centroid,
  std::array<chi_mesh::Vector3, 4>& radii_vectors_xy_tilde,
  std::array<double, 4>& sub_sub_sqr_areas)
{
  auto& SA_i = sub_sub_sqr_areas;

  //============================================= Declare algorithm utilities
  FUNCTION_WEIGHT_FROM_RHO ComputeWeights(
    *this, sq_xy_tilde_centroid, radii_vectors_xy_tilde, sq, legendre);
  double d = 1.0 / sqrt(3.0);
  chi_math::DynamicVector<double> rho = {d, d, d, d};
  double epsilon = 1.0e-1;
  chi_math::DynamicVector<double> delta = {epsilon, epsilon, epsilon, epsilon};
  chi_math::DynamicVector<double> drho_df = {0.0, 0.0, 0.0, 0.0};

  //============================================= Compute initial weights
  auto weights = ComputeWeights(rho);

  //============================================= Apply algorithm
  Chi::log.Log() << "=================================================== ";
  for (int k = 0; k < 150; ++k) // iteration
  {
    //    constexpr int N = 4;
    //    double fac = 1.0/N;
    //    std::array<std::array<double,4>,N> weights_offset;

    auto weights_offset_pos = ComputeWeights(rho + delta);
    auto weights_offset_neg = ComputeWeights(rho - delta);

    double rho_change_total = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      double slope = 0.0;
      slope += 0.5 * (weights_offset_pos[i] - weights[i]);
      slope -= 0.5 * (weights_offset_neg[i] - weights[i]);
      drho_df[i] = delta[i] / slope;

      double delta_rho = 1.0 * drho_df[i] * (SA_i[i] - weights[i]);

      //      delta = {0.0,0.0,0.0,0.0}; delta[i] = epsilon;
      //      auto weights_offset_pos = ComputeWeights(rho + delta);
      //      double slope = (weights_offset_pos[i]-weights[i])/epsilon;
      //
      //      double delta_rho = 10.0*slope*(SA_i[i]-weights[i]);

      rho[i] += delta_rho;
      rho[i] = std::fmax(0.0, rho[i]);
      rho[i] = std::fmin(1.0, rho[i]);
      rho_change_total -= 1.0 * drho_df[i] * (weights[i] - SA_i[i]);
    }

    //================================= Update weights
    weights = ComputeWeights(rho);
    double change = 0.0;
    for (int i = 0; i < 4; ++i)
      change = std::fabs((weights[i] - SA_i[i]) / weights[i]);

    Chi::log.Log() << "Weights: " << weights[0] << " " << weights[1] << " " << weights[2] << " "
                   << weights[3] << " ";
    Chi::log.Log() << "Areas: " << SA_i[0] << " " << SA_i[1] << " " << SA_i[2] << " " << SA_i[3]
                   << "\n";
    Chi::log.Log() << "rhos: " << rho[0] << " " << rho[1] << " " << rho[2] << " " << rho[3] << "\n";
    Chi::log.Log() << k << " " << std::fabs(change);
    Chi::log.Log() << "  ";

    if (rho_change_total < 1.0e-2) break;
    //    if (std::fabs(change) < 1.0e-2) break;
  }
  //  chi::log.Log() << "rhos: "
  //                     << rho[0]/(1.0/sqrt(3.0)) << " "
  //                     << rho[1]/(1.0/sqrt(3.0)) << " "
  //                     << rho[2]/(1.0/sqrt(3.0)) << " "
  //                     << rho[3]/(1.0/sqrt(3.0)) << "\n";

  weights = ComputeWeights(rho);

  for (int i = 0; i < 4; ++i)
  {
    auto xy_tilde = sq_xy_tilde_centroid + rho[i] * radii_vectors_xy_tilde[i];
    auto xyz_prime = sq.rotation_matrix * xy_tilde + sq.translation_vector;
    sq.sub_sqr_points[i] = xyz_prime.Normalized();
    sq.sub_sqr_weights[i] = weights[i];
  }
}