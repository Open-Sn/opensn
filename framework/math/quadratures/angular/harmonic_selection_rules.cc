// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/harmonic_selection_rules.h"
#include "framework/logging/log.h"
#include <algorithm>
#include <sstream>

namespace opensn
{

std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::SelectHarmonics(const SelectionParameters& params)
{
  // For Standard method use all harmonics
  if (params.construction_method == OperatorConstructionMethod::STANDARD)
  {
    return SelectStandard(params);
  }

  // Handle Cartesian Product quadratures for Galerkin methods
  if (params.quadrature_type == AngularQuadratureType::ProductQuadrature)
  {
    if (params.dimension == 2)
      return Select2DCartesianProduct(params);
    else if (params.dimension == 3)
      return Select3DCartesianProduct(params);
  }

  // If no rules match, throw a runtime error
  throw std::runtime_error(
    "No specific Galerkin rules implemented for this quadrature type and dimension. "
    "Harmonic selection failed.");
}

// 2D Product Quadrature Rules
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::Select2DCartesianProduct(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;

  const unsigned int N = params.N;
  const unsigned int max_ell = 2 * (N - 1); // Rule 1: ℓ = 0, ..., 2(N-1)

  opensn::log.Log0Verbose1() << "2D Cartesian Product Galerkin rules:";
  opensn::log.Log0Verbose1() << "  N = " << N;
  opensn::log.Log0Verbose1() << "  max_ell = " << max_ell;

  for (unsigned int ell = 0; ell <= max_ell; ++ell)
  {
    // Rule 2: Loop over m = -ℓ, ..., ℓ
    for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); ++m)
    {
      if (CheckCartesian2DRules(ell, m, N))
      {
        harmonics.emplace_back(ell, m);
        opensn::log.Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
      }
      else
      {
        opensn::log.Log0Verbose2() << "    Rejected: ℓ=" << ell << ", m=" << m;
      }
    }
  }

  opensn::log.Log0Verbose1() << "2D Cartesian Product: Selected " << harmonics.size()
                             << " harmonics out of " << ((max_ell + 1) * (max_ell + 1))
                             << " possible";

  return harmonics;
}

bool
HarmonicSelectionRules::CheckCartesian2DRules(unsigned int ell, int m, unsigned int N)
{
  const unsigned int abs_m = std::abs(m);

  // Rule 3: If ℓ + |m| odd OR (ℓ ≥ N and m > 2N - ℓ - 2), reject
  if ((ell + abs_m) % 2 == 1)
  {
    return false; // ℓ + |m| is odd
  }

  if (ell >= N && m > static_cast<int>(2 * N - ell - 2))
  {
    return false; // ℓ ≥ N and m > 2N - ℓ - 2
  }

  // Rule 4: If ℓ ≥ N and ℓ even and (m = 0 or m < -2N + ℓ), reject
  if (ell >= N && ell % 2 == 0)
  {
    if (m == 0 || m < static_cast<int>(-2 * N + ell))
    {
      return false;
    }
  }

  // Rule 5: If ℓ ≥ N and ℓ odd and m < -2N + ℓ + 2, reject
  if (ell >= N && ell % 2 == 1)
  {
    if (m < static_cast<int>(-2 * N + ell + 2))
    {
      return false;
    }
  }

  // Rule 6: Otherwise, accept
  return true;
}

// 3D Product Quadrature rules
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::Select3DCartesianProduct(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;

  const unsigned int N = params.N;
  const unsigned int max_ell = 2 * (N - 1); // Rule 1: ℓ = 0, ..., 2(N-1)

  opensn::log.Log0Verbose1() << "3D Cartesian Product Galerkin rules:";
  opensn::log.Log0Verbose1() << "  N = " << N;
  opensn::log.Log0Verbose1() << "  max_ell = " << max_ell;

  for (unsigned int ell = 0; ell <= max_ell; ++ell)
  {
    // Rule 2: Loop over m = -ℓ, ..., ℓ
    for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); ++m)
    {
      if (CheckCartesian3DRules(ell, m, N))
      {
        harmonics.emplace_back(ell, m);
        opensn::log.Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
      }
      else
      {
        opensn::log.Log0Verbose2() << "    Rejected: ℓ=" << ell << ", m=" << m;
      }
    }
  }

  opensn::log.Log0Verbose1() << "3D Cartesian Product: Selected " << harmonics.size()
                             << " harmonics out of " << ((max_ell + 1) * (max_ell + 1))
                             << " possible";

  return harmonics;
}

bool
HarmonicSelectionRules::CheckCartesian3DRules(unsigned int ell, int m, unsigned int N)
{
  // Rule 3: If ℓ ≥ N and m ≤ 0 and (m < -N or m ≥ N - ℓ), reject
  if (ell >= N && m <= 0)
  {
    if (m < -static_cast<int>(N) || m >= static_cast<int>(N - ell))
    {
      return false;
    }
  }

  // Rule 4: If ℓ ≥ N and m > 0 and (m ≥ N or m ≤ ℓ - N), reject
  if (ell >= N && m > 0)
  {
    if (m >= static_cast<int>(N) || m <= static_cast<int>(ell - N))
    {
      return false;
    }
  }

  // Rule 5: Otherwise, accept
  return true;
}

// SLDFESQ + Lebedev Implementation Go Here

// Standard method (original behavior)
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::SelectStandard(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;

  if (params.dimension == 1)
  {
    for (unsigned int ell = 0; ell <= params.scattering_order; ++ell)
      harmonics.emplace_back(ell, 0);
  }
  else if (params.dimension == 2)
  {
    for (unsigned int ell = 0; ell <= params.scattering_order; ++ell)
      for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); m += 2)
        harmonics.emplace_back(ell, m);
  }
  else if (params.dimension == 3)
  {
    for (unsigned int ell = 0; ell <= params.scattering_order; ++ell)
      for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); ++m)
        harmonics.emplace_back(ell, m);
  }

  return harmonics;
}

} // namespace opensn
