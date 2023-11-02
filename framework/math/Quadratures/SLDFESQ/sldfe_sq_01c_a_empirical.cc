#include "sldfe_sq.h"

//###################################################################
/**Applies empirical quadrature point optimization.*/
void chi_math::SimplifiedLDFESQ::Quadrature::
  EmpiricalQPOptimization(
    SphericalQuadrilateral &sq,
    chi_math::QuadratureGaussLegendre &legendre,
    chi_mesh::Vertex& sq_xy_tilde_centroid,
    std::array<chi_mesh::Vector3, 4>& radii_vectors_xy_tilde,
    std::array<double,4>& sub_sub_sqr_areas)
{
  FUNCTION_WEIGHT_FROM_RHO ComputeWeights(*this,
                                          sq_xy_tilde_centroid,
                                          radii_vectors_xy_tilde, sq, legendre);
  double d = 1.0/sqrt(3.0);
  chi_math::DynamicVector<double> rho = {d, d, d, d};

  auto weights = ComputeWeights(rho);

  for (int i=0; i<4; ++i)
  {
    auto xy_tilde  = sq_xy_tilde_centroid +
                     rho[i]*radii_vectors_xy_tilde[i];
    auto xyz_prime = sq.rotation_matrix*xy_tilde + sq.translation_vector;

    sq.sub_sqr_points[i] = xyz_prime.Normalized();
    sq.sub_sqr_weights[i] = weights[i];
  }
}