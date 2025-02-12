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
  double weight_sum = 0.0;
  for (unsigned int i = 0; i < Na; ++i)
  {
    for (unsigned int j = 0; j < Np; ++j)
    {
      map_directions_[j].emplace_back(i * Np + j);

      const auto abscissa = QuadraturePointPhiTheta(azimu_ang[i], polar_ang[j]);

      abscissae.emplace_back(abscissa);

      const double weight = wts[i * Np + j];
      weights.emplace_back(weight);
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

  if (verbose)
  {
    log.Log() << ostr.str() << "\n"
              << "Weight sum=" << weight_sum;
  }
}

void
ProductQuadrature::OptimizeForPolarSymmetry(const double normalization)
{
  std::vector<QuadraturePointPhiTheta> new_abscissae;
  std::vector<double> new_weights;
  std::vector<Vector3> new_omegas;
  std::vector<double> new_polar_ang;
  std::vector<double> new_azimu_ang;

  const size_t num_pol = polar_ang.size();
  const size_t num_azi = azimu_ang.size();

  std::vector<unsigned int> new_polar_map;
  for (size_t p = 0; p < num_pol; ++p)
    if (polar_ang[p] < M_PI_2)
    {
      new_polar_ang.push_back(polar_ang[p]);
      new_polar_map.push_back(p);
    }
  new_azimu_ang = azimu_ang;

  const size_t new_num_pol = new_polar_ang.size();
  double weight_sum = 0.0;
  for (size_t a = 0; a < num_azi; ++a)
    for (size_t p = 0; p < new_num_pol; ++p)
    {
      const auto pmap = new_polar_map[p];
      const auto dmap = GetAngleNum(pmap, a);
      new_weights.push_back(weights[dmap]);
      weight_sum += weights[dmap];
    }

  if (normalization > 0.0)
    for (double& w : new_weights)
      w *= normalization / weight_sum;

  AssembleCosines(new_azimu_ang, new_polar_ang, new_weights, false);
  polar_ang = new_polar_ang;
  azimu_ang = new_azimu_ang;
}

GLProductQuadrature1DSlab::GLProductQuadrature1DSlab(int Npolar, bool verbose)
  : ProductQuadrature(1)
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
}

GLCProductQuadrature2DXY::GLCProductQuadrature2DXY(int Npolar, int Nazimuthal, bool verbose)
  : ProductQuadrature(2)
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

  OptimizeForPolarSymmetry(4.0 * M_PI);
}

GLCProductQuadrature3DXYZ::GLCProductQuadrature3DXYZ(int Npolar, int Nazimuthal, bool verbose)
  : ProductQuadrature(3)
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
}

} // namespace opensn
