// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <algorithm>
#include <limits>
#include <numeric>

namespace opensn
{

GLProductQuadrature1DSpherical::GLProductQuadrature1DSpherical(int Npolar,
                                                               int scattering_order,
                                                               bool verbose)
  : CurvilinearProductQuadrature(1, scattering_order)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLProductQuadrature1DSpherical: Npolar must be even.");

  Initialize(Npolar, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

void
GLProductQuadrature1DSpherical::Initialize(int Npolar, const bool verbose)
{
  const auto quad_polar = GaussLegendreQuadrature(Npolar, verbose);
  auto polar_quad(quad_polar);

  // Verifications and corrections (if possible)
  const auto eps = std::numeric_limits<double>::epsilon();

  if (polar_quad.weights.empty())
    throw std::invalid_argument("GLProductQuadrature1DSpherical: "
                                "Invalid polar quadrature size = " +
                                std::to_string(polar_quad.weights.size()));

  // Verifications on polar quadrature
  const double polar_quad_sum_weights = 1.0;
  const auto polar_quad_span = std::pair<double, double>(-1, +1);

  //  weights sum to 1.0
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
    throw std::invalid_argument("GLProductQuadrature1DSpherical: "
                                "Polar quadrature weights sum to zero.");

  // Defined on range [-1;+1]
  if (std::abs(polar_quad.GetRange().first - polar_quad_span.first) > eps or
      std::abs(polar_quad.GetRange().second - polar_quad_span.second) > eps)
    polar_quad.SetRange(polar_quad_span);

  // Abscissae sorted in ascending order
  auto lt_qp = [](const Vector3& qp0, const Vector3& qp1) { return qp0[0] < qp1[0]; };
  if (not std::is_sorted(polar_quad.qpoints.begin(), polar_quad.qpoints.end(), lt_qp))
    throw std::invalid_argument("GLProductQuadrature1DSpherical: "
                                "Polar quadrature abscissae not in ascending order.");

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
GLProductQuadrature1DSpherical::InitializeParameters()
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
GLProductQuadrature1DSpherical::MakeHarmonicIndices()
{
  if (m_to_ell_em_map_.empty())
  {
    for (unsigned int l = 0; l <= scattering_order_; ++l)
      m_to_ell_em_map_.emplace_back(l, 0);
  }
}

GLCProductQuadrature2DRZ::GLCProductQuadrature2DRZ(int Npolar,
                                                   int Nazimuthal,
                                                   int scattering_order,
                                                   bool verbose)
  : CurvilinearProductQuadrature(2, scattering_order)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCProductQuadraturee2DRZ: Npolar must be even.");

  if (Nazimuthal % 4 != 0)
    throw std::invalid_argument("GLCProductQuadraturee2DRZ: Nazimuthal must be a multiple of 4.");

  const auto quad_polar = GaussLegendreQuadrature(Npolar, verbose);
  std::vector<GaussQuadrature> quad_azimuthal;
  for (auto n = 0; n < Npolar; ++n)
    quad_azimuthal.emplace_back(GaussChebyshevQuadrature(Nazimuthal, verbose));
  Initialize(quad_polar, quad_azimuthal, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

void
GLCProductQuadrature2DRZ::Initialize(const GaussQuadrature& quad_polar,
                                     const std::vector<GaussQuadrature>& quad_azimu_vec,
                                     const bool verbose)
{
  auto polar_quad(quad_polar);
  auto azimu_quad_vec(quad_azimu_vec);

  // Verifications and corrections (if possible)
  const auto eps = std::numeric_limits<double>::epsilon();

  //  consistency among polar quadrature and azimuthal quadratures
  if (polar_quad.weights.size() != azimu_quad_vec.size())
    throw std::invalid_argument("GLCProductQuadrature2DRZ: "
                                "Number of azimuthal quadratures does not correspond "
                                "to number of polar points of the polar quadrature.");

  // At present, this class does not handle correctly reduced geometries
  if (polar_quad.weights.empty())
    throw std::invalid_argument("GLCProductQuadrature2DRZ: "
                                "Invalid polar quadrature size = " +
                                std::to_string(polar_quad.weights.size()));

  for (const auto& azimu_quad : azimu_quad_vec)
    if (azimu_quad.weights.empty())
      throw std::invalid_argument("GLCProductQuadrature2DRZ: "
                                  "Invalid azimuthal quadrature size = " +
                                  std::to_string(azimu_quad.weights.size()));

  // Verifications on polar quadrature
  const double polar_quad_sum_weights = 2;
  const auto polar_quad_span = std::pair<double, double>(-1, +1);

  // Weights sum to 2
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
    throw std::invalid_argument("GLCProductQuadrature2DRZ: "
                                "Polar quadrature weights sum to zero.");

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
      std::accumulate(azimu_quad.weights.begin(), azimu_quad.weights.end(), 0.0);
    if (std::abs(integral_weights) > 0)
    {
      const auto fac = azimu_quad_sum_weights / integral_weights;
      if (std::abs(fac - 1) > eps)
        for (auto& w : azimu_quad.weights)
          w *= fac;
    }
    else
      throw std::invalid_argument("GLCProductQuadrature2DRZ: "
                                  "Azimuthal quadrature weights sum to zero.");

    // Defined on range [-1;+1]
    if (std::abs(azimu_quad.GetRange().first - azimu_quad_span.first) > eps or
        std::abs(azimu_quad.GetRange().second - azimu_quad_span.second) > eps)
      azimu_quad.SetRange(azimu_quad_span);

    // Abscissae sorted in ascending order
    auto lt_qp = [](const Vector3& qp0, const Vector3& qp1) { return qp0[0] < qp1[0]; };
    if (!std::is_sorted(azimu_quad.qpoints.begin(), azimu_quad.qpoints.end(), lt_qp))
      throw std::invalid_argument("GLCProductQuadrature2DRZ: "
                                  "Azimuthal quadrature abscissae not in ascending order.");

    // Existence of zero-weight abscissae at the start and at the end of the interval
    if (std::abs(azimu_quad.weights.front()) > eps and
        std::abs(azimu_quad.qpoints.front()[0] - azimu_quad_span.first) > eps)
    {
      azimu_quad.weights.emplace(azimu_quad.weights.begin(), 0);
      azimu_quad.qpoints.emplace(azimu_quad.qpoints.begin(), azimu_quad_span.first);
    }
    if (std::abs(azimu_quad.weights.back()) > eps and
        std::abs(azimu_quad.qpoints.back()[0] - azimu_quad_span.second) > eps)
    {
      azimu_quad.weights.emplace(azimu_quad.weights.end(), 0);
      azimu_quad.qpoints.emplace(azimu_quad.qpoints.end(), azimu_quad_span.second);
    }
  }

  // Product quadrature initialization
  // Compute weights, abscissae $(\varphi, \vartheta)$ and direction vectors
  // $\omega_{pq} := (\mu_{pq}, \xi_{p}, \eta_{pq})$
  weights.clear();
  abscissae.clear();
  omegas.clear();
  for (size_t p = 0; p < azimu_quad_vec.size(); ++p)
  {
    const auto pol_wei = polar_quad.weights[p];
    const auto pol_abs = polar_quad.qpoints[p][0];
    const auto pol_com = std::sqrt(1 - pol_abs * pol_abs);

    for (size_t q = 0; q < azimu_quad_vec[p].weights.size(); ++q)
    {
      const auto& azimu_quad = azimu_quad_vec[p];

      const auto azi_wei = azimu_quad.weights[q];
      const auto azi_abs = azimu_quad.qpoints[q][0];
      const auto azi_com = std::sqrt(1 - azi_abs * azi_abs);

      const auto weight = pol_wei * azi_wei;
      const auto abscissa = QuadraturePointPhiTheta(std::acos(azi_abs), std::acos(pol_abs));
      const auto omega = Vector3(pol_com * azi_abs, pol_abs, pol_com * azi_com);

      weights.emplace_back(weight);
      abscissae.emplace_back(abscissa);
      omegas.emplace_back(omega);
    }
  }
  weights.shrink_to_fit();
  abscissae.shrink_to_fit();
  omegas.shrink_to_fit();

  // Map of direction indices
  unsigned int ind0 = 0;
  map_directions_.clear();
  for (size_t p = 0; p < azimu_quad_vec.size(); ++p)
  {
    std::vector<unsigned int> vec_directions_p;
    for (size_t q = 0; q < azimu_quad_vec[p].weights.size(); ++q)
      vec_directions_p.emplace_back(ind0 + q);
    map_directions_.emplace(p, vec_directions_p);
    ind0 += azimu_quad_vec[p].weights.size();
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
GLCProductQuadrature2DRZ::InitializeParameters()
{
  fac_diamond_difference_.resize(weights.size(), 1);
  fac_streaming_operator_.resize(weights.size(), 0);
  for (size_t p = 0; p < map_directions_.size(); ++p)
  {
    double sum_q_weights = 0;
    for (size_t q = 0; q < map_directions_[p].size(); ++q)
      sum_q_weights += weights[map_directions_[p][q]];
    const auto pi_sum_q_weights = M_PI / sum_q_weights;

    // Interface quantities initialised to starting direction values
    double alpha_interface = 0;
    double phi_interface = abscissae[map_directions_[p].front()].phi;
    std::vector<double> mu_interface(2, std::cos(phi_interface));

    // Initialization permits to forego start direction and final direction
    for (size_t q = 1; q < map_directions_[p].size() - 1; ++q)
    {
      const auto k = map_directions_[p][q];
      const auto w_pq = weights[k];
      const auto mu_pq = omegas[k].x;
      const auto phi_pq = abscissae[k].phi;

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
GLCProductQuadrature2DRZ::MakeHarmonicIndices()
{
  m_to_ell_em_map_.clear();

  for (auto l = 0; l <= scattering_order_; ++l)
    for (auto m = 0; m <= l; ++m)
      m_to_ell_em_map_.emplace_back(l, m);
}

} // namespace opensn
