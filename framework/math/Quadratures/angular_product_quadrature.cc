#include "opensn/framework/math/Quadratures/angular_product_quadrature.h"
#include "opensn/framework/math/Quadratures/quadrature_gausslegendre.h"
#include "opensn/framework/math/Quadratures/quadrature_gausschebyshev.h"

#include <cmath>
#include <sstream>

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"

//###################################################################
/**Initializes the quadrature with custom angles and weights.*/
void
chi_math::ProductQuadrature::AssembleCosines(const std::vector<double>& azimuthal,
                                             const std::vector<double>& polar,
                                             const std::vector<double>& in_weights,
                                             bool verbose)
{
  size_t Na = azimuthal.size();
  size_t Np = polar.size();
  size_t Nw = in_weights.size();

  if (Nw != Na * Np)
  {
    Chi::log.LogAllError() << "Product Quadrature, InitializeWithCustom: mismatch in the amount "
                              "angles and weights. Number of azimuthal angles times number "
                              "polar angles must equal the amount of weights.";
    Chi::Exit(EXIT_FAILURE);
  }

  azimu_ang_ = azimuthal;
  polar_ang_ = polar;

  if (verbose)
  {
    Chi::log.Log() << "Azimuthal angles:";
    for (const auto& ang : azimu_ang_)
      Chi::log.Log() << ang;

    Chi::log.Log() << "Polar angles:";
    for (const auto& ang : polar_ang_)
      Chi::log.Log() << ang;
  }

  //================================================== Create angle pairs
  map_directions_.clear();
  for (unsigned int j = 0; j < Np; ++j)
    map_directions_.emplace(j, std::vector<unsigned int>());

  abscissae_.clear();
  weights_.clear();
  std::stringstream ostr;
  double weight_sum = 0.0;
  for (unsigned int i = 0; i < Na; ++i)
  {
    for (unsigned int j = 0; j < Np; ++j)
    {
      map_directions_[j].emplace_back(i * Np + j);

      const auto abscissa = chi_math::QuadraturePointPhiTheta(azimu_ang_[i], polar_ang_[j]);

      abscissae_.emplace_back(abscissa);

      const double weight = in_weights[i * Np + j];
      weights_.emplace_back(weight);
      weight_sum += weight;

      if (verbose)
      {
        char buf[200];
        snprintf(buf,
                 200,
                 "Varphi=%.2f Theta=%.2f Weight=%.3e\n",
                 abscissa.phi * 180.0 / M_PI,
                 abscissa.theta * 180.0 / M_PI,
                 weight);
        ostr << buf;
      }
    }
  }

  //================================================== Create omega list
  omegas_.clear();
  for (const auto& qpoint : abscissae_)
  {
    chi_mesh::Vector3 new_omega;
    new_omega.x = sin(qpoint.theta) * cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta) * sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

    omegas_.emplace_back(new_omega);

    if (verbose) Chi::log.Log() << "Quadrature angle=" << new_omega.PrintS();
  }

  if (verbose)
  {
    Chi::log.Log() << ostr.str() << "\n"
                   << "Weight sum=" << weight_sum;
  }
}

//###################################################################
/**Optimizes the angular quadrature for polar symmetry by removing
 * all the direction with downward pointing polar angles.
 *
 * \param normalization float. (Optional) The default is a negative number
 *                             which does not apply any normalization. If a
 *                             positive number is provided, the weights will be
 *                             normalized to sum to this number.*/
void
chi_math::ProductQuadrature::OptimizeForPolarSymmetry(const double normalization)
{
  std::vector<chi_math::QuadraturePointPhiTheta> new_abscissae;
  std::vector<double> new_weights;
  std::vector<chi_mesh::Vector3> new_omegas;
  std::vector<double> new_polar_ang;
  std::vector<double> new_azimu_ang;

  const size_t num_pol = polar_ang_.size();
  const size_t num_azi = azimu_ang_.size();

  std::vector<unsigned int> new_polar_map;
  for (size_t p = 0; p < num_pol; ++p)
    if (polar_ang_[p] < M_PI_2)
    {
      new_polar_ang.push_back(polar_ang_[p]);
      new_polar_map.push_back(p);
    }
  new_azimu_ang = azimu_ang_;

  const size_t new_num_pol = new_polar_ang.size();
  double weight_sum = 0.0;
  for (size_t a = 0; a < num_azi; ++a)
    for (size_t p = 0; p < new_num_pol; ++p)
    {
      const auto pmap = new_polar_map[p];
      const auto dmap = GetAngleNum(pmap, a);
      new_weights.push_back(weights_[dmap]);
      weight_sum += weights_[dmap];
    }

  if (normalization > 0.0)
    for (double& w : new_weights)
      w *= normalization / weight_sum;

  AssembleCosines(new_azimu_ang, new_polar_ang, new_weights, false);
  polar_ang_ = new_polar_ang;
  azimu_ang_ = new_azimu_ang;
}

//###################################################################
/**Constructor for Angular Gauss-Legendre.*/
chi_math::AngularQuadratureProdGL::AngularQuadratureProdGL(int Nphemi, bool verbose)
  : chi_math::ProductQuadrature()
{
  chi_math::QuadratureGaussLegendre gl_polar(Nphemi * 2);

  //================================================= Create azimuthal angles
  azimu_ang_.clear();
  azimu_ang_.emplace_back(0.0);

  //================================================== Create polar angles
  polar_ang_.clear();
  for (unsigned int j = 0; j < (Nphemi * 2); ++j)
    polar_ang_.emplace_back(M_PI - acos(gl_polar.qpoints_[j][0]));

  //================================================== Create combined weights
  auto& weights = gl_polar.weights_;

  //================================================== Initialize
  AssembleCosines(azimu_ang_, polar_ang_, weights, verbose);
}

//###################################################################
/**Constructor for Angular Gauss-Legendre-Legendre.*/
chi_math::AngularQuadratureProdGLL::AngularQuadratureProdGLL(int Na, int Np, bool verbose)
{
  chi_math::QuadratureGaussLegendre gl_polar(Np * 2);
  chi_math::QuadratureGaussLegendre gl_azimu(Na * 4);

  //================================================= Create azimuthal angles
  azimu_ang_.clear();
  for (unsigned int i = 0; i < (Na * 4); ++i)
    azimu_ang_.emplace_back(M_PI * gl_azimu.qpoints_[i][0] + M_PI);

  //================================================== Create polar angles
  polar_ang_.clear();
  for (unsigned int j = 0; j < (Np * 2); ++j)
    polar_ang_.emplace_back(M_PI - acos(gl_polar.qpoints_[j][0]));

  //================================================== Create combined weights
  std::vector<double> weights;
  for (unsigned int i = 0; i < azimu_ang_.size(); ++i)
    for (unsigned int j = 0; j < polar_ang_.size(); ++j)
      weights.emplace_back(M_PI * gl_azimu.weights_[i] * gl_polar.weights_[j]);

  //================================================== Initialize
  AssembleCosines(azimu_ang_, polar_ang_, weights, verbose);
}

//###################################################################
/**Constructor for Angular Gauss-Legendre-Chebyshev.*/
chi_math::AngularQuadratureProdGLC::AngularQuadratureProdGLC(int Na, int Np, bool verbose)
{
  chi_math::QuadratureGaussLegendre gl_polar(Np * 2);
  chi_math::QuadratureGaussChebyshev gc_azimu(Na * 4);

  //================================================= Create azimuthal angles
  azimu_ang_.clear();
  for (unsigned int i = 0; i < (Na * 4); ++i)
    azimu_ang_.emplace_back(M_PI * (2 * (i + 1) - 1) / (Na * 4));

  //================================================== Create polar angles
  polar_ang_.clear();
  for (unsigned int j = 0; j < (Np * 2); ++j)
    polar_ang_.emplace_back(M_PI - acos(gl_polar.qpoints_[j][0]));

  //================================================== Create combined weights
  std::vector<double> weights;
  for (unsigned int i = 0; i < azimu_ang_.size(); ++i)
    for (unsigned int j = 0; j < polar_ang_.size(); ++j)
      weights.emplace_back(2 * gc_azimu.weights_[i] * gl_polar.weights_[j]);

  //================================================== Initialize
  AssembleCosines(azimu_ang_, polar_ang_, weights, verbose);
}

//###################################################################
/**Constructor for Custom Angular Product Quadrature.*/
chi_math::AngularQuadratureProdCustom::AngularQuadratureProdCustom(
  const std::vector<double>& azimuthal,
  const std::vector<double>& polar,
  const std::vector<double>& in_weights,
  bool verbose)
{
  size_t Na = azimuthal.size();
  size_t Np = polar.size();
  size_t Nw = in_weights.size();

  if (Nw != Na * Np)
  {
    Chi::log.LogAllError() << "Product Quadrature, InitializeWithCustom: mismatch in the amount "
                              "angles and weights. Number of azimuthal angles times number "
                              "polar angles must equal the amount of weights.";
    Chi::Exit(EXIT_FAILURE);
  }

  AssembleCosines(azimuthal, polar, weights_, verbose);
}
