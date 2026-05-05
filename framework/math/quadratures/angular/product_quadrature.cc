// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <cmath>
#include <sstream>
#include <cassert>

namespace opensn
{

void
ProductQuadrature::AssembleCosines(const std::vector<double>& azimuthal,
                                   const std::vector<double>& polar,
                                   const std::vector<double>& wts,
                                   bool verbose)
{
  size_t Na = azimuthal.size();
  size_t Np = polar.size();

  azimu_ang_ = azimuthal;
  polar_ang_ = polar;

  if (verbose)
  {
    log.Log() << "Azimuthal angles:";
    for (const auto& ang : azimu_ang_)
      log.Log() << ang;

    log.Log() << "Polar angles:";
    for (const auto& ang : polar_ang_)
      log.Log() << ang;
  }

  // Create angle pairs
  map_directions_.clear();
  for (unsigned int j = 0; j < Np; ++j)
    map_directions_.emplace(j, std::vector<unsigned int>());

  abscissae_.clear();
  weights_.clear();
  std::stringstream ostr;
  weight_sum_ = 0.0;
  for (unsigned int i = 0; i < Na; ++i)
  {
    for (unsigned int j = 0; j < Np; ++j)
    {
      map_directions_[j].emplace_back(i * Np + j);

      const auto abscissa = QuadraturePointPhiTheta(azimu_ang_[i], polar_ang_[j]);

      abscissae_.emplace_back(abscissa);

      const double weight = wts[i * Np + j];
      weights_.emplace_back(weight);
      weight_sum_ += weight;

      if (verbose)
      {
        ostr << "Varphi=" << std::fixed << std::setprecision(2) << abscissa.phi * 180.0 / M_PI
             << " Theta=" << std::fixed << std::setprecision(2) << abscissa.theta * 180.0 / M_PI
             << " Weight=" << std::scientific << std::setprecision(3) << weight << '\n';
      }
    }
  }

  // Create omega list
  omegas_.clear();
  for (const auto& qpoint : abscissae_)
  {
    Vector3 new_omega;
    new_omega.x = sin(qpoint.theta) * cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta) * sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

    omegas_.emplace_back(new_omega);

    if (verbose)
      log.Log() << "Quadrature angle=" << new_omega.PrintStr();
  }

  // Normalize weights to 1.0
  for (auto& w : weights_)
    w /= weight_sum_;

  // Compute and store sum of weights after normalization
  weight_sum_ = 0.0;
  for (auto& w : weights_)
    weight_sum_ += w;
}

GLProductQuadrature1DSlab::GLProductQuadrature1DSlab(unsigned int Npolar,
                                                     unsigned int scattering_order,
                                                     bool verbose,
                                                     OperatorConstructionMethod method)
  : ProductQuadrature(1, scattering_order, method)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLProductQuadrature1DSlab: Npolar must be even.");

  n_polar_ = Npolar;
  n_azimuthal_ = 1;

  GaussLegendreQuadrature gl_polar(Npolar);

  // Create azimuthal angles
  azimu_ang_.clear();
  azimu_ang_.emplace_back(0.0);

  // Create polar angles
  polar_ang_.clear();
  for (unsigned int j = 0; j < Npolar; ++j)
    polar_ang_.emplace_back(M_PI - acos(gl_polar.qpoints[j][0]));

  // Create combined weights
  auto& weights = gl_polar.weights;

  // Initialize
  AssembleCosines(azimu_ang_, polar_ang_, weights, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

GLCProductQuadrature2DXY::GLCProductQuadrature2DXY(unsigned int Npolar,
                                                   unsigned int Nazimuthal,
                                                   unsigned int scattering_order,
                                                   bool verbose,
                                                   OperatorConstructionMethod method)
  : ProductQuadrature(2, scattering_order, method)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCProductQuadraturee2DXY: Npolar must be even.");

  if (Nazimuthal % 4 != 0)
    throw std::invalid_argument("GLCProductQuadraturee2DXY: Nazimuthal must be a multiple of 4.");

  n_polar_ = Npolar;
  n_azimuthal_ = Nazimuthal;

  GaussLegendreQuadrature gl_polar(Npolar);
  GaussChebyshevQuadrature gc_azimu(Nazimuthal);

  // Create azimuthal angles
  azimu_ang_.clear();
  for (unsigned int i = 0; i < Nazimuthal; ++i)
    azimu_ang_.emplace_back(M_PI * (2 * (i + 1) - 1) / Nazimuthal);

  // Create polar angles (keep the positive polar cosines)
  const unsigned int half = Npolar / 2;
  polar_ang_.resize(half);
  for (unsigned int j = 0; j < half; ++j)
    polar_ang_[j] = M_PI - std::acos(gl_polar.qpoints[j][0]);

  // Create combined weights
  std::vector<double> weights;
  for (auto i = 0; i < azimu_ang_.size(); ++i)
    for (auto j = 0; j < polar_ang_.size(); ++j)
      weights.emplace_back(2.0 * gc_azimu.weights[i] * gl_polar.weights[j]);

  // Initialize
  AssembleCosines(azimu_ang_, polar_ang_, weights, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

GLCProductQuadrature3DXYZ::GLCProductQuadrature3DXYZ(unsigned int Npolar,
                                                     unsigned int Nazimuthal,
                                                     unsigned int scattering_order,
                                                     bool verbose,
                                                     OperatorConstructionMethod method)
  : ProductQuadrature(3, scattering_order, method)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCProductQuadraturee3DXYZ: Npolar must be even.");

  if (Nazimuthal % 4 != 0)
    throw std::invalid_argument("GLCProductQuadraturee3DXYZ: Nazimuthal must be a multiple of 4.");

  n_polar_ = Npolar;
  n_azimuthal_ = Nazimuthal;

  GaussLegendreQuadrature gl_polar(Npolar);
  GaussChebyshevQuadrature gc_azimu(Nazimuthal);

  // Create azimuthal angles
  azimu_ang_.clear();
  for (unsigned int i = 0; i < Nazimuthal; ++i)
    azimu_ang_.emplace_back(M_PI * (2 * (i + 1) - 1) / Nazimuthal);

  // Create polar angles
  polar_ang_.clear();
  for (unsigned int j = 0; j < Npolar; ++j)
    polar_ang_.emplace_back(M_PI - acos(gl_polar.qpoints[j][0]));

  // Create combined weights
  std::vector<double> weights;
  for (auto i = 0; i < azimu_ang_.size(); ++i)
    for (auto j = 0; j < polar_ang_.size(); ++j)
      weights.emplace_back(2 * gc_azimu.weights[i] * gl_polar.weights[j]);

  // Initialize
  AssembleCosines(azimu_ang_, polar_ang_, weights, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();
}

} // namespace opensn
