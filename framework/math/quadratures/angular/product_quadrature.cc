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
  size_t Nw = wts.size();

  azimu_ang = azimuthal;
  polar_ang = polar;

  if (verbose)
  {
    log.Log() << "Azimuthal angles:";
    for (const auto& ang : azimu_ang)
      log.Log() << ang;

    log.Log() << "Polar angles:";
    for (const auto& ang : polar_ang)
      log.Log() << ang;
  }

  // Create angle pairs
  map_directions_.clear();
  for (unsigned int j = 0; j < Np; ++j)
    map_directions_.emplace(j, std::vector<unsigned int>());

  abscissae.clear();
  weights.clear();
  std::stringstream ostr;
  weight_sum_ = 0.0;
  for (unsigned int i = 0; i < Na; ++i)
  {
    for (unsigned int j = 0; j < Np; ++j)
    {
      map_directions_[j].emplace_back(i * Np + j);

      const auto abscissa = QuadraturePointPhiTheta(azimu_ang[i], polar_ang[j]);

      abscissae.emplace_back(abscissa);

      const double weight = wts[i * Np + j];
      weights.emplace_back(weight);
      weight_sum_ += weight;

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

  // Create omega list
  omegas.clear();
  for (const auto& qpoint : abscissae)
  {
    Vector3 new_omega;
    new_omega.x = sin(qpoint.theta) * cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta) * sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

    omegas.emplace_back(new_omega);

    if (verbose)
      log.Log() << "Quadrature angle=" << new_omega.PrintStr();
  }

  // Normalize weights to 1.0
  for (auto& w : weights)
    w /= weight_sum_;

  // Compute and store sum of weights after normalization
  weight_sum_ = 0.0;
  for (auto& w : weights)
    weight_sum_ += w;
}

GLProductQuadrature1DSlab::GLProductQuadrature1DSlab(int Npolar, int scattering_order, bool verbose)
  : ProductQuadrature(1, scattering_order)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLProductQuadrature1DSlab: Npolar must be even.");

  GaussLegendreQuadrature gl_polar(Npolar);

  // Create azimuthal angles
  azimu_ang.clear();
  azimu_ang.emplace_back(0.0);

  // Create polar angles
  polar_ang.clear();
  for (auto j = 0; j < Npolar; ++j)
    polar_ang.emplace_back(M_PI - acos(gl_polar.qpoints[j][0]));

  // Create combined weights
  auto& weights = gl_polar.weights;

  // Initialize
  AssembleCosines(azimu_ang, polar_ang, weights, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  log.Log() << "Using 1D Slab Gauss–Legendre product quadrature with " << omegas.size()
            << " angles and weight sum of " << std::fixed << std::setprecision(2) << weight_sum_
            << std::endl;
}

GLCProductQuadrature2DXY::GLCProductQuadrature2DXY(int Npolar,
                                                   int Nazimuthal,
                                                   int scattering_order,
                                                   bool verbose)
  : ProductQuadrature(2, scattering_order)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCProductQuadraturee2DXY: Npolar must be even.");

  if (Nazimuthal % 4 != 0)
    throw std::invalid_argument("GLCProductQuadraturee2DXY: Nazimuthal must be a multiple of 4.");

  GaussLegendreQuadrature gl_polar(Npolar);
  GaussChebyshevQuadrature gc_azimu(Nazimuthal);

  // Create azimuthal angles
  azimu_ang.clear();
  for (auto i = 0; i < Nazimuthal; ++i)
    azimu_ang.emplace_back(M_PI * (2 * (i + 1) - 1) / Nazimuthal);

  // Create polar angles (only take the half of the GL nodes < M_PI/2)
  const int half = Npolar / 2;
  polar_ang.resize(half);
  for (int j = 0; j < half; ++j)
    polar_ang[j] = M_PI - std::acos(gl_polar.qpoints[j][0]);

  // Create combined weights
  std::vector<double> weights;
  for (auto i = 0; i < azimu_ang.size(); ++i)
    for (auto j = 0; j < polar_ang.size(); ++j)
      weights.emplace_back(2.0 * gc_azimu.weights[i] * gl_polar.weights[j]);

  // Initialize
  AssembleCosines(azimu_ang, polar_ang, weights, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  log.Log() << "Using 2D XY Gauss–Legendre/Chebyshev product quadrature with " << omegas.size()
            << " angles and weight sum of " << std::fixed << std::setprecision(2) << weight_sum_
            << std::endl;
}

GLCProductQuadrature3DXYZ::GLCProductQuadrature3DXYZ(int Npolar,
                                                     int Nazimuthal,
                                                     int scattering_order,
                                                     bool verbose)
  : ProductQuadrature(3, scattering_order)
{
  if (Npolar % 2 != 0)
    throw std::invalid_argument("GLCProductQuadraturee3DXYZ: Npolar must be even.");

  if (Nazimuthal % 4 != 0)
    throw std::invalid_argument("GLCProductQuadraturee3DXYZ: Nazimuthal must be a multiple of 4.");

  GaussLegendreQuadrature gl_polar(Npolar);
  GaussChebyshevQuadrature gc_azimu(Nazimuthal);

  // Create azimuthal angles
  azimu_ang.clear();
  for (auto i = 0; i < Nazimuthal; ++i)
    azimu_ang.emplace_back(M_PI * (2 * (i + 1) - 1) / Nazimuthal);

  // Create polar angles
  polar_ang.clear();
  for (auto j = 0; j < Npolar; ++j)
    polar_ang.emplace_back(M_PI - acos(gl_polar.qpoints[j][0]));

  // Create combined weights
  std::vector<double> weights;
  for (auto i = 0; i < azimu_ang.size(); ++i)
    for (auto j = 0; j < polar_ang.size(); ++j)
      weights.emplace_back(2 * gc_azimu.weights[i] * gl_polar.weights[j]);

  // Initialize
  AssembleCosines(azimu_ang, polar_ang, weights, verbose);
  MakeHarmonicIndices();
  BuildDiscreteToMomentOperator();
  BuildMomentToDiscreteOperator();

  log.Log() << "Using 3D XYZ Gauss–Legendre/Chebyshev product quadrature with " << omegas.size()
            << " angles and weight sum of " << std::fixed << std::setprecision(2) << weight_sum_
            << std::endl;
}

} // namespace opensn
