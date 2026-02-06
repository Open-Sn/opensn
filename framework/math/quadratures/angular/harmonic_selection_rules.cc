// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/harmonic_selection_rules.h"
#include "framework/math/quadratures/angular/harmonic_selection_rules_emperical.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/logging/log.h"
#include <algorithm>
#include <map>
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

    case AngularQuadratureType::TriangularQuadrature:
      if (params.dimension == 2)
        return Select2DTriangular(params);
      else if (params.dimension == 3)
        return Select3DTriangular(params);
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

  const unsigned int n_polar = params.n_polar / 2;
  const unsigned int n_azimuthal = params.n_azimuthal;
  const size_t num_dir = params.num_angles;

  const unsigned int num_azimu_90 = n_azimuthal / 4;
  const unsigned int L_crit = 2 * n_polar - 1;
  const unsigned int M_crit = 2 * num_azimu_90;

  Logger::GetInstance().Log0()
    << "HarmonicSelectionRules::Select2DCartesianProduct - GALERKIN method";
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

    // Determine min_M
    unsigned int min_M;
    if (L <= L_crit)
      min_M = ind_L;
    else
      min_M = L + 1 - L_crit;

    // Determine max_M
    const unsigned int max_M = std::min(L, M_crit);

    // Loop over M with step of 2 - add sine first, then cosine
    for (unsigned int M = min_M; M <= max_M; M += 2)
    {
      if (harmonics.size() >= num_dir)
        break;

      // Sine always (add negative m first)
      harmonics.emplace_back(L, -static_cast<int>(M));
      Logger::GetInstance().Log0() << "    Accepted: ℓ=" << L << ", m=-" << M;
      Logger::GetInstance().Log0Verbose2()
        << "    Accepted: ℓ=" << L << ", m=" << -static_cast<int>(M);
      if (harmonics.size() == num_dir)
        break;

      // Cosine only if M < M_crit (add positive m second)
      if (M < M_crit)
      {
        harmonics.emplace_back(L, static_cast<int>(M));
        Logger::GetInstance().Log0() << "    Accepted: ℓ=" << L << ", m=+" << M;
        Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << L << ", m=" << M;
        if (harmonics.size() == num_dir)
          break;
      }
    }

    if (harmonics.size() >= num_dir)
      break;

    // M=0 cosine only when L <= L_crit and even L (add last)
    if (L <= L_crit && ind_L == 2)
    {
      harmonics.emplace_back(L, 0);
      Logger::GetInstance().Log0() << "    Accepted: ℓ=" << L << ", m=" << 0;
      Logger::GetInstance().Log0Verbose2() << "    Accepted: ℓ=" << L << ", m=" << 0;
      if (harmonics.size() == num_dir)
        break;
    }

    // Toggle L parity
    ind_L = 3 - ind_L;
  }

  Logger::GetInstance().Log0() << "2D Product (General): Selected " << harmonics.size()
                               << " harmonics (needed " << num_dir << ")";
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

// 2D Triangular Quadrature Rules
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::Select2DTriangular(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;
  const size_t nb_dir = params.num_angles;

  // Compute directions per octant
  // For 2D: nb_dir = ndir_oct * 4 (4 quadrants)
  const size_t ndir_oct = nb_dir / 4;

  Logger::GetInstance().Log0Verbose1() << "2D Triangular Galerkin rules:";
  Logger::GetInstance().Log0Verbose1() << "  nb_dir     = " << nb_dir;
  Logger::GetInstance().Log0Verbose1() << "  ndir_oct   = " << ndir_oct;

  // Parity mode indices: ee=1, eo=2, oe=3, oo=4
  enum ParityMode
  {
    ee = 1,
    eo = 2,
    oe = 3,
    oo = 4
  };

  std::map<int, size_t> nb_C = {{ee, 0}, {eo, 0}, {oe, 0}, {oo, 0}}; // cosine counts
  std::map<int, size_t> nb_S = {{ee, 0}, {eo, 0}, {oe, 0}, {oo, 0}}; // sine counts

  int ind_K = 2; // 2 for even K, 1 for odd K (toggled each K)
  int K = -1;

  while (harmonics.size() < nb_dir)
  {
    K++;
    // mode = 7 - 3*ind_K -> 1 (ee) when K even, 4 (oo) when K odd
    int mode = 7 - 3 * ind_K;

    // cosine (K,0) only when mode==ee and we still have capacity
    if (mode == ee && nb_C[ee] < ndir_oct)
    {
      harmonics.emplace_back(K, 0);
      nb_C[ee]++;
      Logger::GetInstance().Log0Verbose2() << "    Accepted: K=" << K << ", L=0 (cosine)";
      if (harmonics.size() == nb_dir)
        break;
    }

    // Then L with the same parity as K: L = ind_K, ind_K+2, ..., K
    int Lcur = ind_K;
    while (Lcur <= K && harmonics.size() < nb_dir)
    {
      // cosine (K, +L)
      if (nb_C[mode] < ndir_oct)
      {
        harmonics.emplace_back(K, Lcur);
        nb_C[mode]++;
        Logger::GetInstance().Log0Verbose2()
          << "    Accepted: K=" << K << ", L=" << Lcur << " (cosine)";
        if (harmonics.size() == nb_dir)
          break;
      }
      // sine (K, -L)
      if (nb_S[mode] < ndir_oct && harmonics.size() < nb_dir)
      {
        harmonics.emplace_back(K, -Lcur);
        nb_S[mode]++;
        Logger::GetInstance().Log0Verbose2()
          << "    Accepted: K=" << K << ", L=" << -Lcur << " (sine)";
      }
      Lcur += 2;
    }

    // Toggle K parity
    ind_K = 3 - ind_K;
  }

  Logger::GetInstance().Log0Verbose1()
    << "2D Triangular: Selected " << harmonics.size() << " harmonics";

  return harmonics;
}

// 3D Triangular Quadrature Rules
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::Select3DTriangular(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;
  const size_t nb_dir = params.num_angles;

  // Compute directions per octant
  // For 3D: nb_dir = ndir_oct * 8 (8 octants)
  const size_t ndir_oct = nb_dir / 8;

  Logger::GetInstance().Log0Verbose1() << "3D Triangular Galerkin rules:";
  Logger::GetInstance().Log0Verbose1() << "  nb_dir     = " << nb_dir;
  Logger::GetInstance().Log0Verbose1() << "  ndir_oct   = " << ndir_oct;

  // Parity mode indices: ee=1, eo=2, oe=3, oo=4
  enum ParityMode
  {
    ee = 1,
    eo = 2,
    oe = 3,
    oo = 4
  };

  std::map<int, size_t> nb_C = {{ee, 0}, {eo, 0}, {oe, 0}, {oo, 0}}; // cosine counts
  std::map<int, size_t> nb_S = {{ee, 0}, {eo, 0}, {oe, 0}, {oo, 0}}; // sine counts

  int ind_K = 2; // 2 for even K, 1 for odd K (toggled each K)
  int K = -1;

  while (harmonics.size() < nb_dir)
  {
    K++;
    int ind_L = 1;             // 1=even L, 2=odd L (toggled as L increases)
    int id_ = 2 * (2 - ind_K); // id=0 for even K, 2 for odd K

    // L=0 cosine first
    int mode = ind_L + id_; // ee if even K, oe if odd K
    if (nb_C[mode] < ndir_oct)
    {
      harmonics.emplace_back(K, 0);
      nb_C[mode]++;
      Logger::GetInstance().Log0Verbose2() << "    Accepted: K=" << K << ", L=0 (cosine)";
      if (harmonics.size() == nb_dir)
        break;
    }

    // L = 1..K with toggling L parity
    for (int L = 1; L <= K; ++L)
    {
      ind_L = 3 - ind_L; // toggle L parity
      mode = ind_L + id_;

      // cosine (K, +L)
      if (nb_C[mode] < ndir_oct)
      {
        harmonics.emplace_back(K, L);
        nb_C[mode]++;
        Logger::GetInstance().Log0Verbose2()
          << "    Accepted: K=" << K << ", L=" << L << " (cosine)";
        if (harmonics.size() == nb_dir)
          break;
      }

      // sine (K, -L)
      if (nb_S[mode] < ndir_oct && harmonics.size() < nb_dir)
      {
        harmonics.emplace_back(K, -L);
        nb_S[mode]++;
        Logger::GetInstance().Log0Verbose2()
          << "    Accepted: K=" << K << ", L=" << -L << " (sine)";
      }
    }

    // Toggle K parity
    ind_K = 3 - ind_K;
  }

  Logger::GetInstance().Log0Verbose1()
    << "3D Triangular: Selected " << harmonics.size() << " harmonics";

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
      lebedev_order = DetermineLebedevOrder(params.num_angles, params.dimension);
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

  // Select the appropriate empirical harmonics map based on dimension
  const auto& lebedev_harmonics = (params.dimension == 2)
                                    ? EmpiricalHarmonicRules::LebedevHarmonics2D
                                    : EmpiricalHarmonicRules::LebedevHarmonics3D;

  // Look up empirical harmonics
  auto order_num = lebedev_harmonics.find(lebedev_order);
  if (order_num != lebedev_harmonics.end())
  {
    Logger::GetInstance().Log0Verbose1()
      << "Lebedev " << params.dimension << "D Galerkin rules (order " << lebedev_order << "):";

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

    Logger::GetInstance().Log0Verbose1()
      << "Lebedev " << params.dimension << "D Galerkin: Selected " << harmonics.size()
      << " harmonics for order " << lebedev_order;
  }
  else
  {
    Logger::GetInstance().Log0Warning()
      << "No empirical harmonic rules found for " << params.dimension << "D Lebedev order "
      << lebedev_order << ". Falling back to standard harmonic selection.";
    return SelectStandard(params);
  }

  return harmonics;
}

// Helper function to map number of angles to Lebedev order
int
HarmonicSelectionRules::DetermineLebedevOrder(size_t num_angles, unsigned int dimension)
{
  // 3D Lebedev sets have specific numbers of points (full sphere)
  static const std::map<size_t, int> angle_to_order_3d =
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

  // 2D Lebedev sets have fewer points (upper hemisphere only)
  // These values depend on how many points have z >= 0 for each order
  static const std::map<size_t, int> angle_to_order_2d = {
    {5, 3},   // Order 3
    {9, 5},   // Order 5
    {17, 7},  // Order 7
    {25, 9},  // Order 9
    {29, 11}, // Order etc...
    {45, 13},  {49, 15},   {61, 17},   {77, 19},   {93, 21},   {105, 23},  {125, 25}, {141, 27},
    {161, 29}, {185, 31},  {229, 35},  {309, 41},  {401, 47},  {505, 53},  {621, 59}, {749, 65},
    {889, 71}, {1041, 77}, {1205, 83}, {1381, 89}, {1569, 95}, {1769, 101}};

  const auto& angle_to_order = (dimension == 2) ? angle_to_order_2d : angle_to_order_3d;

  auto order_num = angle_to_order.find(num_angles);
  if (order_num != angle_to_order.end())
    return order_num->second;
  else
    throw std::runtime_error("Unknown " + std::to_string(dimension) + "D Lebedev quadrature with " +
                             std::to_string(num_angles) + " angles");
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
      refinement_level = DetermineSLDFELevel(params.num_angles, params.dimension);
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

  // Select the appropriate empirical harmonics map based on dimension
  const auto& sldfe_harmonics_map = (params.dimension == 2)
                                      ? EmpiricalHarmonicRules::SLDFEHarmonics2D
                                      : EmpiricalHarmonicRules::SLDFEHarmonics;

  // Look up empirical harmonics
  auto sldfe_harmonics = sldfe_harmonics_map.find(refinement_level);
  if (sldfe_harmonics != sldfe_harmonics_map.end())
  {
    Logger::GetInstance().Log0Verbose1()
      << "SLDFESQ " << params.dimension << "D Galerkin rules (level " << refinement_level << "):";

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

    Logger::GetInstance().Log0Verbose1()
      << "SLDFESQ " << params.dimension << "D Galerkin: Selected " << harmonics.size()
      << " harmonics for refinement level " << refinement_level;
  }
  else
  {
    Logger::GetInstance().Log0Warning()
      << "No empirical harmonic rules found for " << params.dimension << "D SLDFESQ refinement level "
      << refinement_level << ". Falling back to standard harmonic selection.";
    return SelectStandard(params);
  }

  return harmonics;
}

// Helper function to map number of angles to SLDFESQ refinement level
int
HarmonicSelectionRules::DetermineSLDFELevel(size_t num_angles, unsigned int dimension)
{
  // 3D SLDFESQ refinement levels have specific numbers of points (full sphere)
  // Formula: 96 × (level+1)² directions
  static const std::map<size_t, int> angle_to_level_3d = {
    {96, 0},   // Level 0: 96 × 1² = 96
    {384, 1},  // Level 1: 96 × 2² = 384
    {864, 2},  // Level 2: 96 × 3² = 864
    {1536, 3}, // Level 3: 96 × 4² = 1536
  };

  // 2D SLDFESQ sets have fewer points (upper hemisphere only, z > 0)
  // Formula: 48 × (level+1)² directions
  static const std::map<size_t, int> angle_to_level_2d = {
    {48, 0},   // Level 0: 48 × 1² = 48
    {192, 1},  // Level 1: 48 × 2² = 192
    {432, 2},  // Level 2: 48 × 3² = 432
    {768, 3},  // Level 3: 48 × 4² = 768
  };

  const auto& angle_to_level = (dimension == 2) ? angle_to_level_2d : angle_to_level_3d;

  auto refinement_level = angle_to_level.find(num_angles);
  if (refinement_level != angle_to_level.end())
    return refinement_level->second;
  else
    throw std::runtime_error("Unknown " + std::to_string(dimension) + "D SLDFESQ quadrature with " +
                             std::to_string(num_angles) + " angles");
}

// Standard method
std::vector<AngularQuadrature::HarmonicIndices>
HarmonicSelectionRules::SelectStandard(const SelectionParameters& params)
{
  std::vector<AngularQuadrature::HarmonicIndices> harmonics;

  Logger::GetInstance().Log0() << "HarmonicSelectionRules::SelectStandard - STANDARD method";
  Logger::GetInstance().Log0() << "  Dimension: " << params.dimension
                               << ", Scattering order: " << params.scattering_order;

  if (params.dimension == 1)
  {
    for (unsigned int ell = 0; ell <= params.scattering_order; ++ell)
      harmonics.emplace_back(ell, 0);
  }
  else if (params.dimension == 2)
  {
    Logger::GetInstance().Log0() << "  2D Harmonic selection:";
    for (unsigned int ell = 0; ell <= params.scattering_order; ++ell)
    {
      std::stringstream ss;
      ss << "    ℓ=" << ell << ": m = [";
      bool first = true;
      for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); m += 2)
      {
        harmonics.emplace_back(ell, m);
        if (!first)
          ss << ", ";
        ss << (m >= 0 ? "+" : "") << m;
        first = false;
      }
      ss << "]";
      Logger::GetInstance().Log0() << ss.str();
    }
  }
  else if (params.dimension == 3)
  {
    for (unsigned int ell = 0; ell <= params.scattering_order; ++ell)
      for (int m = -static_cast<int>(ell); m <= static_cast<int>(ell); ++m)
        harmonics.emplace_back(ell, m);
  }

  Logger::GetInstance().Log0() << "  Total harmonics selected: " << harmonics.size();

  return harmonics;
}

} // namespace opensn