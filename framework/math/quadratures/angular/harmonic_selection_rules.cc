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

    case AngularQuadratureType::SLDFEsq:
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

  const unsigned int n_polar = params.n_polar;
  const unsigned int n_azimuthal = params.n_azimuthal;
  const size_t num_dir = params.num_angles;

  const unsigned int num_azimu_90 = n_azimuthal / 4;
  const unsigned int L_crit = 2 * n_polar - 1;
  const unsigned int M_crit = 2 * num_azimu_90;

  Logger::GetInstance().Log0Verbose1() << "  npolar = " << n_polar;
  Logger::GetInstance().Log0Verbose1() << "  nb_phi_90 = " << num_azimu_90;
  Logger::GetInstance().Log0Verbose1() << "  L_crit = " << L_crit;
  Logger::GetInstance().Log0Verbose1() << "  M_crit = " << M_crit;

  unsigned int ind_L = 2; // parity toggle for L
  const unsigned int max_L = std::min(2 * (n_polar + num_azimu_90 - 1), params.scattering_order);

  for (unsigned int L = 0; L <= max_L; ++L)
  {
    if (harmonics.size() >= num_dir)
      break;

    // M=0 cosine only when L <= L_crit and even L
    if (L <= L_crit && ind_L == 2)
    {
      harmonics.emplace_back(L, 0);
      Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << L << ", m=" << 0;
      if (harmonics.size() == num_dir)
        break;
    }

    // Determine min_M
    unsigned int min_M;
    if (L <= L_crit)
      min_M = ind_L;
    else
      min_M = L + 1 - L_crit;

    // Determine max_M
    const unsigned int max_M = std::min(L, M_crit);

    // Loop over M with step of 2
    for (unsigned int M = min_M; M <= max_M; M += 2)
    {
      if (harmonics.size() >= num_dir)
        break;

      // Cosine only if M < M_crit
      if (M < M_crit)
      {
        harmonics.emplace_back(L, static_cast<int>(M));
        Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << L << ", m=" << M;
        if (harmonics.size() == num_dir)
          break;
      }

      // Sine always
      harmonics.emplace_back(L, -static_cast<int>(M));
      Logger::GetInstance().Log0Verbose2()
        << "    Accepted: ℓ=" << L << ", m=" << -static_cast<int>(M);
    }

    // Toggle L parity
    ind_L = 3 - ind_L;
  }

  Logger::GetInstance().Log0Verbose1()
    << "2D Product (General): Selected " << harmonics.size() << " harmonics";

  return harmonics;
}

// 3D Product Quadrature rules
// Note that (K, L) transcribes to (L, M) in this section
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::Select3DCartesianProduct(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;
  harmonics.reserve(params.num_angles);

  const unsigned int n_polar = params.n_polar;
  const unsigned int n_azimuthal = params.n_azimuthal;
  const size_t nb_dir = params.num_angles;

  const unsigned int nb_phi_90 = n_azimuthal / 4;
  const unsigned int npolar = n_polar / 2;

  const unsigned int K_crit = 2 * npolar - 1;
  const unsigned int L_crit = 2 * nb_phi_90;
  const unsigned int min_L_default = 1;

  Logger::GetInstance().Log0Verbose1() << "3D Product Galerkin (General) rules:";
  Logger::GetInstance().Log0Verbose1() << "  npolar       = " << npolar;
  Logger::GetInstance().Log0Verbose1() << "  nb_phi_90    = " << nb_phi_90;
  Logger::GetInstance().Log0Verbose1() << "  K_crit       = " << K_crit;
  Logger::GetInstance().Log0Verbose1() << "  L_crit       = " << L_crit;

  const unsigned int max_K = std::min(2 * (npolar + nb_phi_90) - 1, params.scattering_order);

  for (unsigned int K = 0; K <= max_K; ++K)
  {
    if (harmonics.size() >= nb_dir)
      break;

    unsigned int min_L;
    if (K <= K_crit)
    {
      // L = 0 cosine (append (K,0) once for K <= K_crit)
      harmonics.emplace_back(K, 0);
      Logger::GetInstance().Log0Verbose2() << "    Accepted: K=" << K << ", L=" << 0;
      if (harmonics.size() == nb_dir)
        break;
      min_L = min_L_default; // start inner L loop at 1
    }
    else
    {
      // For K > K_crit the min L is K - K_crit
      min_L = K - K_crit;
    }

    // max_L = min(K, L_crit)
    const unsigned int max_L = std::min(K, L_crit);

    // inner loop: L = min_L .. max_L (inclusive)
    for (unsigned int L = min_L; L <= max_L; ++L)
    {
      if (harmonics.size() >= nb_dir)
        break;

      // cosine only if L < L_crit
      if (L < L_crit)
      {
        harmonics.emplace_back(K, static_cast<int>(L));
        Logger::GetInstance().Log0Verbose2() << "    Accepted: K=" << K << ", L=" << L;
        if (harmonics.size() == nb_dir)
          break;
      }

      // sine always: append (K, -L)
      harmonics.emplace_back(K, -static_cast<int>(L));
      Logger::GetInstance().Log0Verbose2()
        << "    Accepted: K=" << K << ", L=" << -static_cast<int>(L);
    }
  }

  Logger::GetInstance().Log0Verbose1()
    << "3D Product (General): Selected " << harmonics.size() << " harmonics";

  return harmonics;
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
      // Only include harmonics with ℓ <= scattering_order
      if (ell <= params.scattering_order)
      {
        harmonics.emplace_back(ell, m);
        Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
      }
      else
      {
        Logger::GetInstance().Log0Verbose2()
          << "    Rejected (ℓ > scattering_order): ℓ=" << ell << ", m=" << m;
      }
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
      // Only include harmonics with ℓ <= scattering_order
      if (ell <= params.scattering_order)
      {
        harmonics.emplace_back(ell, m);
        Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << ell << ", m=" << m;
      }
      else
      {
        Logger::GetInstance().Log0Verbose2()
          << "    Rejected (ℓ > scattering_order): ℓ=" << ell << ", m=" << m;
      }
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