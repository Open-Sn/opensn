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
AngularQuadrature::MakeHarmonicIndices()
{
  m_to_ell_em_map_.clear();

  if (dimension_ == 1)
    for (auto ell = 0; ell <= scattering_order_; ++ell)
      m_to_ell_em_map_.emplace_back(ell, 0);
  else if (dimension_ == 2)
    for (auto ell = 0; ell <= scattering_order_; ++ell)
      for (auto m = -ell; m <= ell; m += 2)
        m_to_ell_em_map_.emplace_back(ell, m);
  else if (dimension_ == 3)
    for (auto ell = 0; ell <= scattering_order_; ++ell)
      for (auto m = -ell; m <= ell; ++m)
        m_to_ell_em_map_.emplace_back(ell, m);
}

void
AngularQuadrature::BuildDiscreteToMomentOperator()
{
  d2m_op_.clear();

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  for (const auto& ell_em : m_to_ell_em_map_)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (auto n = 0; n < num_angles; ++n)
    {
      const auto& cur_angle = abscissae[n];
      double value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
      double w = weights[n];
      cur_mom.push_back(value * w);
    }

    d2m_op_.push_back(cur_mom);
  }

  // Verbose printout
  std::stringstream outs;
  outs << "\nQuadrature d2m operator:\n";
  for (auto n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (auto m = 0; m < num_moms; ++m)
    {
      outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << d2m_op_[m][n]
           << " ";
    }
    outs << "\n";
  }

  log.Log0Verbose1() << outs.str();
}

void
AngularQuadrature::BuildMomentToDiscreteOperator()
{
  m2d_op_.clear();

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  const auto normalization = std::accumulate(weights.begin(), weights.end(), 0.0);

  for (const auto& ell_em : m_to_ell_em_map_)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (auto n = 0; n < num_angles; ++n)
    {
      const auto& cur_angle = abscissae[n];
      double value = ((2.0 * ell_em.ell + 1.0) / normalization) *
                     Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
      cur_mom.push_back(value);
    }

    m2d_op_.push_back(cur_mom);
  }

  // Verbose printout
  std::stringstream outs;

  outs << "\nQuadrature m2d operator:\n";
  for (auto n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (auto m = 0; m < num_moms; ++m)
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
  return d2m_op_;
}

std::vector<std::vector<double>> const&
AngularQuadrature::GetMomentToDiscreteOperator() const
{
  return m2d_op_;
}

const std::vector<AngularQuadrature::HarmonicIndices>&
AngularQuadrature::GetMomentToHarmonicsIndexMap() const
{
  return m_to_ell_em_map_;
}

} // namespace opensn
