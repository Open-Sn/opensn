// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/spherical_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <limits>
#include <numeric>

namespace opensn
{

SphericalQuadrature::SphericalQuadrature(const GaussQuadrature& quad_polar, const bool verbose)
  : CurvilinearQuadrature()
{
  Initialize(quad_polar, verbose);
}

void
SphericalQuadrature::Initialize(const GaussQuadrature& quad_polar, const bool verbose)
{
  auto polar_quad(quad_polar);

  // Verifications and corrections (if possible)
  const auto eps = std::numeric_limits<double>::epsilon();

  if (polar_quad.weights.size() == 0)
    throw std::invalid_argument("SphericalQuadrature::Initialize : "
                                "invalid polar quadrature size = " +
                                std::to_string(polar_quad.weights.size()));

  // Verifications on polar quadrature
  const double polar_quad_sum_weights = 2;
  const auto polar_quad_span = std::pair<double, double>(-1, +1);

  //  weights sum to 2
  const auto integral_weights =
    std::accumulate(polar_quad.weights.begin(), polar_quad.weights.end(), 0.0);
  if (std::abs(integral_weights) > 0)
  {
    const auto fac = polar_quad_sum_weights / integral_weights;
    if (std::abs(fac - 1) > eps)
      for (auto& w : polar_quad.weights)
        w *= fac;
  }
  else
    throw std::invalid_argument("Sphericaluadrature::Initialize : "
                                "polar quadrature weights sum to zero.");

  // Defined on range [-1;+1]
  if (std::abs(polar_quad.GetRange().first - polar_quad_span.first) > eps or
      std::abs(polar_quad.GetRange().second - polar_quad_span.second) > eps)
    polar_quad.SetRange(polar_quad_span);

  // Abscissae sorted in ascending order
  auto lt_qp = [](const Vector3& qp0, const Vector3& qp1) { return qp0[0] < qp1[0]; };
  if (not std::is_sorted(polar_quad.qpoints.begin(), polar_quad.qpoints.end(), lt_qp))
    throw std::invalid_argument("SphericalQuadrature::Initialize : "
                                "polar quadrature abscissae not in ascending order.");

  // Existence of zero-weight abscissae at the start and at the end of the interval
  if (std::abs(polar_quad.weights.front()) > eps and
      std::abs(polar_quad.qpoints.front()[0] - polar_quad_span.first) > eps)
  {
    polar_quad.weights.emplace(polar_quad.weights.begin(), 0);
    polar_quad.qpoints.emplace(polar_quad.qpoints.begin(), polar_quad_span.first);
  }
  if (std::abs(polar_quad.weights.back()) > eps and
      std::abs(polar_quad.qpoints.back()[0] - polar_quad_span.second) > eps)
  {
    polar_quad.weights.emplace(polar_quad.weights.end(), 0);
    polar_quad.qpoints.emplace(polar_quad.qpoints.end(), polar_quad_span.second);
  }

  // Product quadrature initialization
  // Compute weights, abscissae $(0, \vartheta_{p})$ and direction vectors
  // $\omega_{p} := ((1-\mu_{p}^{2})^{1/2}, 0, \mu_{p})$
  weights.clear();
  abscissae.clear();
  omegas.clear();
  for (size_t p = 0; p < polar_quad.weights.size(); ++p)
  {
    const auto pol_wei = polar_quad.weights[p];
    const auto pol_abs = polar_quad.qpoints[p][0];
    const auto pol_com = std::sqrt(1 - pol_abs * pol_abs);

    const auto weight = pol_wei;
    const auto abscissa = QuadraturePointPhiTheta(0, std::acos(pol_abs));
    const auto omega = Vector3(pol_com, 0, pol_abs);

    weights.emplace_back(weight);
    abscissae.emplace_back(abscissa);
    omegas.emplace_back(omega);
  }
  weights.shrink_to_fit();
  abscissae.shrink_to_fit();
  omegas.shrink_to_fit();

  // Map of direction indices
  map_directions_.clear();
  for (size_t p = 0; p < polar_quad.weights.size(); ++p)
  {
    std::vector<unsigned int> vec_directions_p;
    vec_directions_p.emplace_back(p);
    map_directions_.emplace(p, vec_directions_p);
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
    log.Log() << "curvilinear product quadrature : spherical" << std::endl;
    for (size_t k = 0; k < weights.size(); ++k)
      log.Log() << "angle index " << k << ": weight = " << weights[k] << ", (phi, theta) = ("
                << abscissae[k].phi << ", " << abscissae[k].theta << ")"
                << ", omega = " << omegas[k].PrintStr()
                << ", fac_diamond_difference = " << fac_diamond_difference_[k]
                << ", fac_streaming_operator = " << fac_streaming_operator_[k] << std::endl;
    const auto sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    log.Log() << "sum(weights) = " << sum_weights << std::endl;
  }
}

void
SphericalQuadrature::InitializeParameters()
{
  fac_diamond_difference_.resize(weights.size(), 1);
  fac_streaming_operator_.resize(weights.size(), 0);

  // Interface quantities initialised to starting direction values
  double alpha_interface = 0;
  std::vector<double> mu_interface(2, omegas[map_directions_[0].front()].z);

  // Initialization permits to forego start direction and final direction
  for (size_t p = 1; p < map_directions_.size() - 1; ++p)
  {
    const auto k = map_directions_[p][0];
    const auto w_p = weights[k];
    const auto mu_p = omegas[k].z;

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
SphericalQuadrature::MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  if (m_to_ell_em_map_.empty())
  {
    if (dimension == 1)
      for (unsigned int l = 0; l <= scattering_order; ++l)
        m_to_ell_em_map_.emplace_back(l, 0);
    else
      throw std::invalid_argument("SphericalQuadrature::MakeHarmonicIndices : "
                                  "invalid dimension.");
  }
}

} // namespace opensn
