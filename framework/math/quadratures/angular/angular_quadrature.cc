// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>
#include <numeric>

namespace opensn
{

void
AngularQuadrature::OptimizeForPolarSymmetry(const double normalization)
{
  std::vector<QuadraturePointPhiTheta> new_abscissae;
  std::vector<double> new_weights;
  std::vector<Vector3> new_omegas;

  const size_t num_dirs = omegas.size();
  double weight_sum = 0.0;
  for (size_t d = 0; d < num_dirs; ++d)
    if (omegas[d].z > 0.0)
    {
      new_abscissae.emplace_back(abscissae[d]);
      new_weights.emplace_back(weights[d]);
      new_omegas.emplace_back(omegas[d]);
      weight_sum += weights[d];
    }

  if (normalization > 0.0)
    for (double& w : new_weights)
      w *= normalization / weight_sum;

  abscissae = std::move(new_abscissae);
  weights = std::move(new_weights);
  omegas = std::move(new_omegas);
}

void
AngularQuadrature::MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  m_to_ell_em_map_.clear();

  if (dimension == 1)
    for (int ell = 0; ell <= scattering_order; ++ell)
      m_to_ell_em_map_.emplace_back(ell, 0);
  else if (dimension == 2)
    for (int ell = 0; ell <= scattering_order; ++ell)
      for (int m = -ell; m <= ell; m += 2)
        m_to_ell_em_map_.emplace_back(ell, m);
  else if (dimension == 3)
    for (int ell = 0; ell <= scattering_order; ++ell)
      for (int m = -ell; m <= ell; ++m)
        m_to_ell_em_map_.emplace_back(ell, m);
}

void
AngularQuadrature::BuildDiscreteToMomentOperator(unsigned int scattering_order, int dimension)
{
  if (d2m_op_built_)
    return;

  d2m_op_.clear();
  MakeHarmonicIndices(scattering_order, dimension);

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  for (const auto& ell_em : m_to_ell_em_map_)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n = 0; n < num_angles; ++n)
    {
      const auto& cur_angle = abscissae[n];
      double value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
      double w = weights[n];
      cur_mom.push_back(value * w);
    }

    d2m_op_.push_back(cur_mom);
  }
  d2m_op_built_ = true;

  // Verbose printout
  std::stringstream outs;
  outs << "\nQuadrature d2m operator:\n";
  for (int n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (int m = 0; m < num_moms; ++m)
    {
      outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << d2m_op_[m][n]
           << " ";
    }
    outs << "\n";
  }

  log.Log0Verbose1() << outs.str();
}

void
AngularQuadrature::BuildMomentToDiscreteOperator(unsigned int scattering_order, int dimension)
{
  if (m2d_op_built_)
    return;

  m2d_op_.clear();
  MakeHarmonicIndices(scattering_order, dimension);

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  const auto normalization = std::accumulate(weights.begin(), weights.end(), 0.0);

  for (const auto& ell_em : m_to_ell_em_map_)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n = 0; n < num_angles; ++n)
    {
      const auto& cur_angle = abscissae[n];
      double value = ((2.0 * ell_em.ell + 1.0) / normalization) *
                     Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
      cur_mom.push_back(value);
    }

    m2d_op_.push_back(cur_mom);
  } // for m
  m2d_op_built_ = true;

  // Verbose printout
  std::stringstream outs;

  outs << "\nQuadrature m2d operator:\n";
  for (int n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (int m = 0; m < num_moms; ++m)
    {
      outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << m2d_op_[m][n]
           << " ";
    }
    outs << "\n";
  }

  log.Log0Verbose1() << outs.str();
}

std::vector<std::vector<double>> const&
AngularQuadrature::GetDiscreteToMomentOperator() const
{
  const std::string fname = __FUNCTION__;
  if (not d2m_op_built_)
    throw std::logic_error(fname +
                           ": Called but D2M operator not yet built. "
                           "Make a call to BuildDiscreteToMomentOperator before using this.");
  return d2m_op_;
}

std::vector<std::vector<double>> const&
AngularQuadrature::GetMomentToDiscreteOperator() const
{
  const std::string fname = __FUNCTION__;
  if (not m2d_op_built_)
    throw std::logic_error(fname +
                           ": Called but M2D operator not yet built. "
                           "Make a call to BuildMomentToDiscreteOperator before using this.");
  return m2d_op_;
}

const std::vector<AngularQuadrature::HarmonicIndices>&
AngularQuadrature::GetMomentToHarmonicsIndexMap() const
{
  const std::string fname = __FUNCTION__;
  if (not(d2m_op_built_ or m2d_op_built_))
    throw std::logic_error(fname + ": Called but map not yet built. "
                                   "Make a call to either BuildDiscreteToMomentOperator or"
                                   "BuildMomentToDiscreteOperator before using this.");
  return m_to_ell_em_map_;
}

AngularQuadratureCustom::AngularQuadratureCustom(std::vector<double>& azimuthal,
                                                 std::vector<double>& polar,
                                                 std::vector<double>& weights,
                                                 bool verbose)
{
  size_t Na = azimuthal.size();
  size_t Np = polar.size();
  size_t Nw = weights.size();

  if ((Na - Np != 0) or (Na - Nw != 0))
  {
    log.LogAllError() << "AngularQuadrature::InitializeWithCustom: supplied"
                         " vectors need to be of equal length.";
    Exit(EXIT_FAILURE);
  }

  // Create angle pairs
  std::stringstream ostr;
  double weight_sum = 0.0;

  for (unsigned int i = 0; i < Na; ++i)
  {
    const auto abscissa = QuadraturePointPhiTheta(azimuthal[i], polar[i]);

    abscissae.push_back(abscissa);

    const double weight = weights[i];
    weights.push_back(weight);
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

  // Create omega list
  for (const auto& qpoint : abscissae)
  {
    Vector3 new_omega;
    new_omega.x = sin(qpoint.theta) * cos(qpoint.phi);
    new_omega.y = sin(qpoint.theta) * sin(qpoint.phi);
    new_omega.z = cos(qpoint.theta);

    omegas.push_back(new_omega);
  }

  if (verbose)
  {
    log.Log() << ostr.str() << "\n"
              << "Weight sum=" << weight_sum;
  }
}

} // namespace opensn
