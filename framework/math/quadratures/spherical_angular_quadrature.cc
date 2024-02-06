#include "framework/math/quadratures/spherical_angular_quadrature.h"

#include <algorithm>
#include <limits>
#include <numeric>

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

SphericalAngularQuadrature::SphericalAngularQuadrature(const Quadrature& quad_polar,
                                                       const bool verbose)
  : CurvilinearAngularQuadrature()
{
  Initialize(quad_polar, verbose);
}

void
SphericalAngularQuadrature::Initialize(const Quadrature& quad_polar, const bool verbose)
{
  //  copies of input quadratures
  auto polar_quad(quad_polar);

  //  --------------------------------------------------------------------------
  //  verifications and corrections (if possible)
  //  --------------------------------------------------------------------------
  const auto eps = std::numeric_limits<double>::epsilon();

  if (polar_quad.weights_.size() == 0)
    throw std::invalid_argument("SphericalAngularQuadrature::Initialize : "
                                "invalid polar quadrature size = " +
                                std::to_string(polar_quad.weights_.size()));

  //  --------------------------------------------------------------------------
  //  verifications on polar quadrature
  //  --------------------------------------------------------------------------
  const double polar_quad_sum_weights = 2;
  const auto polar_quad_span = std::pair<double, double>(-1, +1);

  //  weights sum to 2
  const auto integral_weights =
    std::accumulate(polar_quad.weights_.begin(), polar_quad.weights_.end(), 0.0);
  if (std::abs(integral_weights) > 0)
  {
    const auto fac = polar_quad_sum_weights / integral_weights;
    if (std::abs(fac - 1) > eps)
      for (auto& w : polar_quad.weights_)
        w *= fac;
  }
  else
    throw std::invalid_argument("SphericalAngularQuadrature::Initialize : "
                                "polar quadrature weights sum to zero.");

  //  defined on range [-1;+1]
  if (std::abs(polar_quad.GetRange().first - polar_quad_span.first) > eps or
      std::abs(polar_quad.GetRange().second - polar_quad_span.second) > eps)
    polar_quad.SetRange(polar_quad_span);

  //  abscissae sorted in ascending order
  auto lt_qp = [](const QuadraturePointXYZ& qp0, const QuadraturePointXYZ& qp1)
  { return qp0[0] < qp1[0]; };
  if (!std::is_sorted(polar_quad.qpoints_.begin(), polar_quad.qpoints_.end(), lt_qp))
    throw std::invalid_argument("SphericalAngularQuadrature::Initialize : "
                                "polar quadrature abscissae not in ascending order.");

  //  existence of zero-weight abscissae at the start and at the end of the
  //  interval
  if (std::abs(polar_quad.weights_.front()) > eps and
      std::abs(polar_quad.qpoints_.front()[0] - polar_quad_span.first) > eps)
  {
    polar_quad.weights_.emplace(polar_quad.weights_.begin(), 0);
    polar_quad.qpoints_.emplace(polar_quad.qpoints_.begin(), polar_quad_span.first);
  }
  if (std::abs(polar_quad.weights_.back()) > eps and
      std::abs(polar_quad.qpoints_.back()[0] - polar_quad_span.second) > eps)
  {
    polar_quad.weights_.emplace(polar_quad.weights_.end(), 0);
    polar_quad.qpoints_.emplace(polar_quad.qpoints_.end(), polar_quad_span.second);
  }

  //  --------------------------------------------------------------------------
  //  product quadrature : initialisation
  //  --------------------------------------------------------------------------

  //  compute weights, abscissae $(0, \vartheta_{p})$ and direction vectors
  //  $\omega_{p} := ((1-\mu_{p}^{2})^{1/2}, 0, \mu_{p})$
  weights_.clear();
  abscissae_.clear();
  omegas_.clear();
  for (size_t p = 0; p < polar_quad.weights_.size(); ++p)
  {
    const auto pol_wei = polar_quad.weights_[p];
    const auto pol_abs = polar_quad.qpoints_[p][0];
    const auto pol_com = std::sqrt(1 - pol_abs * pol_abs);

    const auto weight = pol_wei;
    const auto abscissa = QuadraturePointPhiTheta(0, std::acos(pol_abs));
    const auto omega = Vector3(pol_com, 0, pol_abs);

    weights_.emplace_back(weight);
    abscissae_.emplace_back(abscissa);
    omegas_.emplace_back(omega);
  }
  weights_.shrink_to_fit();
  abscissae_.shrink_to_fit();
  omegas_.shrink_to_fit();

  //  map of direction indices
  map_directions_.clear();
  for (size_t p = 0; p < polar_quad.weights_.size(); ++p)
  {
    std::vector<unsigned int> vec_directions_p;
    vec_directions_p.emplace_back(p);
    map_directions_.emplace(p, vec_directions_p);
  }

  //  --------------------------------------------------------------------------
  //  curvilinear product quadrature : compute additional parametrising factors
  //  --------------------------------------------------------------------------
  InitializeParameters();

  //  --------------------------------------------------------------------------
  //  print
  //  --------------------------------------------------------------------------
  if (verbose)
  {
    log.Log() << "map_directions" << std::endl;
    for (const auto& dir : map_directions_)
    {
      log.Log() << "polar level " << dir.first << " : ";
      for (const auto& q : dir.second)
        log.Log() << q << ", ";
      log.Log() << std::endl;
    }
    log.Log() << "curvilinear product quadrature : spherical" << std::endl;
    for (size_t k = 0; k < weights_.size(); ++k)
      log.Log() << "angle index " << k << ": weight = " << weights_[k] << ", (phi, theta) = ("
                << abscissae_[k].phi << ", " << abscissae_[k].theta << ")"
                << ", omega = " << omegas_[k].PrintStr()
                << ", fac_diamond_difference = " << fac_diamond_difference_[k]
                << ", fac_streaming_operator = " << fac_streaming_operator_[k] << std::endl;
    const auto sum_weights = std::accumulate(weights_.begin(), weights_.end(), 0.0);
    log.Log() << "sum(weights) = " << sum_weights << std::endl;
  }
}

void
SphericalAngularQuadrature::InitializeParameters()
{
  fac_diamond_difference_.resize(weights_.size(), 1);
  fac_streaming_operator_.resize(weights_.size(), 0);

  //  interface quantities initialised to starting direction values
  double alpha_interface = 0;
  std::vector<double> mu_interface(2, omegas_[map_directions_[0].front()].z);

  //  initialisation permits to forego start direction and final direction
  for (size_t p = 1; p < map_directions_.size() - 1; ++p)
  {
    const auto k = map_directions_[p][0];
    const auto w_p = weights_[k];
    const auto mu_p = omegas_[k].z;

    alpha_interface -= w_p * mu_p;

    mu_interface[0] = mu_interface[1];
    mu_interface[1] += w_p;

    const auto tau = (mu_p - mu_interface[0]) / (mu_interface[1] - mu_interface[0]);

    fac_diamond_difference_[k] = tau;
    fac_streaming_operator_[k] = alpha_interface / (w_p * tau) + mu_p;
    fac_streaming_operator_[k] *= 2;
  }
}

void
SphericalAngularQuadrature::MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  if (m_to_ell_em_map_.empty())
  {
    if (dimension == 1)
      for (unsigned int l = 0; l <= scattering_order; ++l)
        m_to_ell_em_map_.emplace_back(l, 0);
    else
      throw std::invalid_argument("SphericalAngularQuadrature::MakeHarmonicIndices : "
                                  "invalid dimension.");
  }
}

} // namespace opensn
