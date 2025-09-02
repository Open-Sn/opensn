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
    for (auto ell = 0; ell <= scattering_order_; ++ell) // NOLINT
      m_to_ell_em_map_.emplace_back(ell, 0);
  else if (dimension_ == 2)
    for (auto ell = 0; ell <= scattering_order_; ++ell) // NOLINT
      for (auto m = -ell; m <= ell; m += 2)
        m_to_ell_em_map_.emplace_back(ell, m);
  else if (dimension_ == 3)
    for (auto ell = 0; ell <= scattering_order_; ++ell) // NOLINT
      for (auto m = -ell; m <= ell; ++m)
        m_to_ell_em_map_.emplace_back(ell, m);
}

void
AngularQuadrature::BuildDiscreteToMomentOperator()
{
  d2m_op_.clear();

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  d2m_op_.assign(num_angles, std::vector<double>(num_moms, 0.0));

  for (size_t n = 0; n < num_angles; ++n)
  {
    const auto& ang = abscissae[n];
    const double w = weights[n];
    auto& row = d2m_op_[n];

    for (size_t m = 0; m < num_moms; ++m)
    {
      const auto& ell_em = m_to_ell_em_map_[m];
      const double val = Ylm(ell_em.ell, ell_em.m, ang.phi, ang.theta);
      row[m] = val * w;
    }
  }

  std::stringstream outs;
  outs << "\nQuadrature d2m operator (angle rows, moment columns):\n";
  for (size_t n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (size_t m = 0; m < num_moms; ++m)
    {
      outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << d2m_op_[n][m]
           << ' ';
    }
  }
  log.Log0Verbose1() << outs.str();
}

void
AngularQuadrature::BuildMomentToDiscreteOperator()
{
  m2d_op_.clear();

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  const double normalization = std::accumulate(weights.begin(), weights.end(), 0.0);

  m2d_op_.assign(num_angles, std::vector<double>(num_moms, 0.0));

  for (size_t n = 0; n < num_angles; ++n)
  {
    const auto& ang = abscissae[n];
    auto& row = m2d_op_[n];

    for (size_t m = 0; m < num_moms; ++m)
    {
      const auto& ell_em = m_to_ell_em_map_[m];
      const double val =
        ((2.0 * ell_em.ell + 1.0) / normalization) * Ylm(ell_em.ell, ell_em.m, ang.phi, ang.theta);
      row[m] = val;
    }
  }

  std::stringstream outs;
  outs << "\nQuadrature m2d operator (angle rows, moment columns):\n";
  for (size_t n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (size_t m = 0; m < num_moms; ++m)
    {
      outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << m2d_op_[n][m]
           << ' ';
    }
    outs << '\n';
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
