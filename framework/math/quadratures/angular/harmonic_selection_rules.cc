// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/harmonic_selection_rules.h"
#include "framework/math/quadratures/angular/harmonic_selection_rules_emperical.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
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

  // Handle different quadrature types for Galerkin methods
  switch (params.quadrature_type)
  {
    case AngularQuadratureType::ProductQuadrature:
      if (params.dimension == 2)
        return Select2DCartesianProduct(params);
      else if (params.dimension == 3)
        return Select3DCartesianProduct(params);
      break;

    case AngularQuadratureType::LebedevQuadrature:
      return SelectLebedev(params);

    case AngularQuadratureType::SLDFESQ:
      return SelectSLDFESQ(params);

    default:
      break;
  }

  // If no rules match, fall back to Standard Method with warning
  Logger::GetInstance().Log0Warning()
    << "No specific Galerkin rules implemented for this quadrature type and dimension. "
    << "Falling back to standard harmonic selection.";
  return SelectStandard(params);
}

// 2D Product Quadrature Rules
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::Select2DCartesianProduct(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;

  const unsigned int N = params.quadrature_order > 0 ? params.quadrature_order : params.N;
  const unsigned int max_ell = 2 * (N - 1); // Rule 1: ℓ = 0, ..., 2(N-1)

  Logger::GetInstance().Log0Verbose1() << "2D Cartesian Product Galerkin rules:";
  Logger::GetInstance().Log0Verbose1() << "  N = " << N;
  Logger::GetInstance().Log0Verbose1() << "  max_ell = " << max_ell;

  for (unsigned int ell = 0; ell <= max_ell; ++ell)
  {
    // Rule 2: Loop over m = -ℓ, ..., ℓ
    for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); ++m)
    {
      if (CheckCartesian2DRules(ell, m, N))
      {
        harmonics.emplace_back(ell, m);
        Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
      }
      else
      {
        Logger::GetInstance().Log0Verbose2() << "    Rejected: ℓ=" << ell << ", m=" << m;
      }
    }
  }

  Logger::GetInstance().Log0Verbose1()
    << "2D Cartesian Product: Selected " << harmonics.size() << " harmonics out of "
    << ((max_ell + 1) * (max_ell + 1)) << " possible";

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

  const unsigned int N = params.quadrature_order > 0 ? params.quadrature_order : params.N;
  const unsigned int max_ell = 2 * (N - 1); // Rule 1: ℓ = 0, ..., 2(N-1)

  Logger::GetInstance().Log0Verbose1() << "3D Cartesian Product Galerkin rules:";
  Logger::GetInstance().Log0Verbose1() << "  N = " << N;
  Logger::GetInstance().Log0Verbose1() << "  max_ell = " << max_ell;

  for (unsigned int ell = 0; ell <= max_ell; ++ell)
  {
    // Rule 2: Loop over m = -ℓ, ..., ℓ
    for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); ++m)
    {
      if (CheckCartesian3DRules(ell, m, N))
      {
        harmonics.emplace_back(ell, m);
        Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
      }
      else
      {
        Logger::GetInstance().Log0Verbose2() << "    Rejected: ℓ=" << ell << ", m=" << m;
      }
    }
  }

  Logger::GetInstance().Log0Verbose1()
    << "3D Cartesian Product: Selected " << harmonics.size() << " harmonics out of "
    << ((max_ell + 1) * (max_ell + 1)) << " possible";

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

// Lebedev Quadrature Rules
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::SelectLebedev(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;

  // Determine the Lebedev order
  int lebedev_order = params.quadrature_order;
  if (lebedev_order == 0)
  {
    // Obtain from Lookup Table
    try
    {
      lebedev_order = DetermineLebedevOrder(params.num_angles);
    }
    // If not in Lookup Table - Error
    catch (const std::runtime_error& e)
    {
      Logger::GetInstance().Log0Warning()
        << "Could not determine Lebedev order from number of angles (" << params.num_angles << "). "
        << e.what() << " Falling back to standard harmonic selection.";
      return SelectStandard(params);
    }
  }

  // Look up empirical harmonics
  auto order_num = EmpiricalHarmonicRules::LebedevHarmonics.find(lebedev_order);
  if (order_num != EmpiricalHarmonicRules::LebedevHarmonics.end())
  {
    Logger::GetInstance().Log0Verbose1()
      << "Lebedev Galerkin rules (order " << lebedev_order << "):";

    for (const auto& [ell, m] : order_num->second.harmonics)
    {
      harmonics.emplace_back(ell, m);
      Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
    }

    Logger::GetInstance().Log0Verbose1() << "Lebedev Galerkin: Selected " << harmonics.size()
                                         << " harmonics for order " << lebedev_order;
  }
  else
  {
    Logger::GetInstance().Log0Warning()
      << "No empirical harmonic rules found for Lebedev order " << lebedev_order
      << ". Falling back to standard harmonic selection.";
    return SelectStandard(params);
  }

  return harmonics;
}

// Helper function to map number of angles to Lebedev order
int
HarmonicSelectionRules::DetermineLebedevOrder(size_t num_angles)
{
  // Lebedev sets have specific numbers of points
  // Map num_angles to the corresponding order
  static const std::map<size_t, int> angle_to_order =
    {
      {6, 3},   // Order 3
      {14, 5},  // Order 5
      {26, 7},  // Order 7
      {38, 9},  // Order 9
      {50, 11}, // Order etc...
      {74, 13}, // FIX ME ---- ONLY 62 / 74 ORTHOGONAL VECTORS FOUND
      {86, 15},   {110, 17},  {146, 19},  {170, 21},  {194, 23},
      {230, 25}, // FIX ME ---- ONLY 224 / 230 ORTHOGONAL VECTORS FOUND
      {266, 27}, // FIX ME ---- ONLY 246 / 266 ORTHOGONAL VECTORS FOUND
      {302, 29},  {350, 31},  {434, 35},  {590, 41},  {770, 47},
      {974, 53},  {1202, 59}, {1454, 65}, {1730, 71}, {2030, 77},
      {2354, 83}, {2702, 89}, {3074, 95}, {3470, 101} // End of calculated Lebedev Quadrature
    };

  auto order_num = angle_to_order.find(num_angles);
  if (order_num != angle_to_order.end())
    return order_num->second;
  else
    throw std::runtime_error("Unknown Lebedev quadrature with " + std::to_string(num_angles) +
                             " angles");
}

// SLDFESQ Quadrature Rules
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::SelectSLDFESQ(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;

  // Determine the refinement level
  int refinement_level = params.quadrature_order;
  if (refinement_level == 0)
  {
    // Try to infer from number of angles
    try
    {
      refinement_level = DetermineSLDFELevel(params.num_angles);
    }
    catch (const std::runtime_error& e)
    {
      Logger::GetInstance().Log0Warning()
        << "Could not determine SLDFESQ refinement level from number of angles ("
        << params.num_angles << "). " << e.what()
        << " Falling back to standard harmonic selection.";
      return SelectStandard(params);
    }
  }

  // Look up empirical harmonics
  auto sldfe_harmonics = EmpiricalHarmonicRules::SLDFEHarmonics.find(refinement_level);
  if (sldfe_harmonics != EmpiricalHarmonicRules::SLDFEHarmonics.end())
  {
    Logger::GetInstance().Log0Verbose1()
      << "SLDFESQ Galerkin rules (level " << refinement_level << "):";

    for (const auto& [ell, m] : sldfe_harmonics->second.harmonics)
    {
      harmonics.emplace_back(ell, m);
      Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
    }

    Logger::GetInstance().Log0Verbose1() << "SLDFESQ Galerkin: Selected " << harmonics.size()
                                         << " harmonics for refinement level " << refinement_level;
  }
  else
  {
    Logger::GetInstance().Log0Warning()
      << "No empirical harmonic rules found for SLDFESQ refinement level " << refinement_level
      << ". Falling back to standard harmonic selection.";
    return SelectStandard(params);
  }

  return harmonics;
}

// Helper function to map number of angles to SLDFESQ refinement level
int
HarmonicSelectionRules::DetermineSLDFELevel(size_t num_angles)
{
  // SLDFESQ refinement levels have specific numbers of points
  // These are approximations - adjust based on actual SLDFESQ implementation
  static const std::map<size_t, int> angle_to_level = {
    {96, 0},   // Level 0
    {384, 1},  // Level 1
    {1536, 2}, // Level 2
    {6144, 3}, // Level 3
  };

  auto refinement_level = angle_to_level.find(num_angles);
  if (refinement_level != angle_to_level.end())
    return refinement_level->second;
  else
    throw std::runtime_error("Unknown SLDFESQ quadrature with " + std::to_string(num_angles) +
                             " angles");
}

// Standard method
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