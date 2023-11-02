#include "chi_math.h"

//###################################################################
/**Computes the factorial of an integer.*/
double
chi_math::Factorial(const int x)
{
  double factorial_value = 1.0;
  for (int i = 2; i <= x; ++i)
    factorial_value *= i;

  return factorial_value;
}

//###################################################################
/**Determines the azimuthal- and polar-angle associated with
 * the given direction vector.
 * Returns a pair = [azimuthal-angle,polar-angle].*/
std::pair<double, double>
chi_math::OmegaToPhiThetaSafe(const chi_mesh::Vector3& omega)
{
  // Notes: asin maps [-1,+1] to [-pi/2,+pi/2]
  //        acos maps [-1,+1] to [0,pi]
  // This mapping requires some logic for determining the azimuthal angle.
  //
  const auto omega_hat = omega.Normalized();

  double mu = omega_hat.z;
  mu = std::min(mu, 1.0);
  mu = std::max(mu, -1.0);

  double theta = acos(mu);

  //===== Handling omega aligned to k_hat
  if (std::fabs(omega_hat.z) < 1.0e-16) return {0.0, theta};

  double cos_phi = omega_hat.x / sin(theta);
  cos_phi = std::min(cos_phi, 1.0);
  cos_phi = std::max(cos_phi, -1.0);

  //===== Computing varphi for NE and NW quadrant
  if (omega_hat.y >= 0.0) return {acos(cos_phi), theta};
  //===== Computing varphi for SE and SW quadrant
  else
    return {2.0 * M_PI - acos(cos_phi), theta};
}
