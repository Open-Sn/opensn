#include "framework/math/quadratures/angular/cylindrical_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <limits>
#include <numeric>

namespace opensn
{

CylindricalQuadrature::CylindricalQuadrature(const SpatialQuadrature& quad_polar,
                                             const SpatialQuadrature& quad_azimu,
                                             const bool verbose)
  : CurvilinearQuadrature()
{
  const auto np = quad_polar.weights_.size();
  std::vector<SpatialQuadrature> quad_azimu_vec(np, quad_azimu);
  Initialize(quad_polar, quad_azimu_vec, verbose);
}

CylindricalQuadrature::CylindricalQuadrature(const SpatialQuadrature& quad_polar,
                                             const std::vector<SpatialQuadrature>& quad_azimu_vec,
                                             const bool verbose)
  : CurvilinearQuadrature()
{
  Initialize(quad_polar, quad_azimu_vec, verbose);
}

void
CylindricalQuadrature::Initialize(const SpatialQuadrature& quad_polar,
                                  const std::vector<SpatialQuadrature>& quad_azimu_vec,
                                  const bool verbose)
{
  auto polar_quad(quad_polar);
  auto azimu_quad_vec(quad_azimu_vec);

  // Verifications and corrections (if possible)
  const auto eps = std::numeric_limits<double>::epsilon();

  //  consistency among polar quadrature and azimuthal quadratures
  if (polar_quad.weights_.size() != azimu_quad_vec.size())
    throw std::invalid_argument("CylindricalQuadrature::Initialize : "
                                "number of azimuthal quadratures does not correspond "
                                "to number of polar points of the polar quadrature.");

  // At present, this class does not handle correctly reduced geometries
  if (polar_quad.weights_.size() == 0)
    throw std::invalid_argument("CylindricalQuadrature::Initialize : "
                                "invalid polar quadrature size = " +
                                std::to_string(polar_quad.weights_.size()));

  for (const auto& azimu_quad : azimu_quad_vec)
    if (azimu_quad.weights_.size() == 0)
      throw std::invalid_argument("CylindricalQuadrature::Initialize : "
                                  "invalid azimuthal quadrature size = " +
                                  std::to_string(azimu_quad.weights_.size()));

  // Verifications on polar quadrature
  const double polar_quad_sum_weights = 2;
  const auto polar_quad_span = std::pair<double, double>(-1, +1);

  // Weights sum to 2
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
    throw std::invalid_argument("CylindricalQuadrature::Initialize : "
                                "polar quadrature weights sum to zero.");

  // Defined on range [-1;+1]
  if (std::abs(polar_quad.GetRange().first - polar_quad_span.first) > eps or
      std::abs(polar_quad.GetRange().second - polar_quad_span.second) > eps)
    polar_quad.SetRange(polar_quad_span);

  // Verifications on azimuthal quadrature
  const double azimu_quad_sum_weights = M_PI;
  const auto azimu_quad_span = std::pair<double, double>(-1, +1);

  for (auto& azimu_quad : azimu_quad_vec)
  {
    // Weights sum to $\pi$
    const auto integral_weights =
      std::accumulate(azimu_quad.weights_.begin(), azimu_quad.weights_.end(), 0.0);
    if (std::abs(integral_weights) > 0)
    {
      const auto fac = azimu_quad_sum_weights / integral_weights;
      if (std::abs(fac - 1) > eps)
        for (auto& w : azimu_quad.weights_)
          w *= fac;
    }
    else
      throw std::invalid_argument("CylindricalQuadrature::Initialize : "
                                  "azimuthal quadrature weights sum to zero.");

    // Defined on range [-1;+1]
    if (std::abs(azimu_quad.GetRange().first - azimu_quad_span.first) > eps or
        std::abs(azimu_quad.GetRange().second - azimu_quad_span.second) > eps)
      azimu_quad.SetRange(azimu_quad_span);

    // Abscissae sorted in ascending order
    auto lt_qp = [](const Vector3& qp0, const Vector3& qp1) { return qp0[0] < qp1[0]; };
    if (!std::is_sorted(azimu_quad.qpoints_.begin(), azimu_quad.qpoints_.end(), lt_qp))
      throw std::invalid_argument("CylindricalQuadrature::Initialize : "
                                  "azimuthal quadrature abscissae not in ascending order.");

    // Existence of zero-weight abscissae at the start and at the end of the interval
    if (std::abs(azimu_quad.weights_.front()) > eps and
        std::abs(azimu_quad.qpoints_.front()[0] - azimu_quad_span.first) > eps)
    {
      azimu_quad.weights_.emplace(azimu_quad.weights_.begin(), 0);
      azimu_quad.qpoints_.emplace(azimu_quad.qpoints_.begin(), azimu_quad_span.first);
    }
    if (std::abs(azimu_quad.weights_.back()) > eps and
        std::abs(azimu_quad.qpoints_.back()[0] - azimu_quad_span.second) > eps)
    {
      azimu_quad.weights_.emplace(azimu_quad.weights_.end(), 0);
      azimu_quad.qpoints_.emplace(azimu_quad.qpoints_.end(), azimu_quad_span.second);
    }
  }

  // Product quadrature initialization
  // Compute weights, abscissae $(\varphi, \vartheta)$ and direction vectors
  // $\omega_{pq} := (\mu_{pq}, \xi_{p}, \eta_{pq})$
  weights_.clear();
  abscissae_.clear();
  omegas_.clear();
  for (size_t p = 0; p < azimu_quad_vec.size(); ++p)
  {
    const auto pol_wei = polar_quad.weights_[p];
    const auto pol_abs = polar_quad.qpoints_[p][0];
    const auto pol_com = std::sqrt(1 - pol_abs * pol_abs);

    for (size_t q = 0; q < azimu_quad_vec[p].weights_.size(); ++q)
    {
      const auto& azimu_quad = azimu_quad_vec[p];

      const auto azi_wei = azimu_quad.weights_[q];
      const auto azi_abs = azimu_quad.qpoints_[q][0];
      const auto azi_com = std::sqrt(1 - azi_abs * azi_abs);

      const auto weight = pol_wei * azi_wei;
      const auto abscissa = QuadraturePointPhiTheta(std::acos(azi_abs), std::acos(pol_abs));
      const auto omega = Vector3(pol_com * azi_abs, pol_abs, pol_com * azi_com);

      weights_.emplace_back(weight);
      abscissae_.emplace_back(abscissa);
      omegas_.emplace_back(omega);
    }
  }
  weights_.shrink_to_fit();
  abscissae_.shrink_to_fit();
  omegas_.shrink_to_fit();

  // Map of direction indices
  unsigned int ind0 = 0;
  map_directions_.clear();
  for (size_t p = 0; p < azimu_quad_vec.size(); ++p)
  {
    std::vector<unsigned int> vec_directions_p;
    for (size_t q = 0; q < azimu_quad_vec[p].weights_.size(); ++q)
      vec_directions_p.emplace_back(ind0 + q);
    map_directions_.emplace(p, vec_directions_p);
    ind0 += azimu_quad_vec[p].weights_.size();
  }

  // Curvilinear product quadrature
  // Compute additional parametrising factors
  InitializeParameters();

  // Print
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
    log.Log() << "curvilinear product quadrature : cylindrical" << std::endl;
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
CylindricalQuadrature::InitializeParameters()
{
  fac_diamond_difference_.resize(weights_.size(), 1);
  fac_streaming_operator_.resize(weights_.size(), 0);
  for (size_t p = 0; p < map_directions_.size(); ++p)
  {
    double sum_q_weights = 0;
    for (size_t q = 0; q < map_directions_[p].size(); ++q)
      sum_q_weights += weights_[map_directions_[p][q]];
    const auto pi_sum_q_weights = M_PI / sum_q_weights;

    // Interface quantities initialised to starting direction values
    double alpha_interface = 0;
    double phi_interface = abscissae_[map_directions_[p].front()].phi;
    std::vector<double> mu_interface(2, std::cos(phi_interface));

    // Initialization permits to forego start direction and final direction
    for (size_t q = 1; q < map_directions_[p].size() - 1; ++q)
    {
      const auto k = map_directions_[p][q];
      const auto w_pq = weights_[k];
      const auto mu_pq = omegas_[k].x;
      const auto phi_pq = abscissae_[k].phi;

      alpha_interface -= w_pq * mu_pq;

      phi_interface -= w_pq * pi_sum_q_weights;
      mu_interface[0] = mu_interface[1];
      mu_interface[1] = std::cos(phi_interface);

      const auto mu = std::cos(phi_pq);
      const auto tau = (mu - mu_interface[0]) / (mu_interface[1] - mu_interface[0]);

      fac_diamond_difference_[k] = tau;
      fac_streaming_operator_[k] = alpha_interface / (w_pq * tau) + mu_pq;
    }
  }
}

void
CylindricalQuadrature::MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  if (m_to_ell_em_map_.empty())
  {
    if (dimension == 1)
    {
      for (unsigned int l = 0; l <= scattering_order; ++l)
        for (int m = 0; m <= l; ++m)
          if ((l + m) % 2 == 0)
            m_to_ell_em_map_.emplace_back(l, m);
    }
    else if (dimension == 2)
      for (unsigned int l = 0; l <= scattering_order; ++l)
        for (int m = 0; m <= l; ++m)
          m_to_ell_em_map_.emplace_back(l, m);
    else
      throw std::invalid_argument("CylindricalQuadrature::MakeHarmonicIndices : "
                                  "invalid dimension.");
  }
}

} // namespace opensn
