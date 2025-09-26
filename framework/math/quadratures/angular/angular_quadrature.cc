// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/math.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <cmath>

namespace opensn
{

void
AngularQuadrature::MakeHarmonicIndices()
{
  m_to_ell_em_map_.clear();

  // For Galerkin methods with Product Quadrature, use rule-based selection
  if ((construction_method_ == OperatorConstructionMethod::GALERKIN_ONE ||
       construction_method_ == OperatorConstructionMethod::GALERKIN_THREE) &&
      type_ == AngularQuadratureType::ProductQuadrature)
  {
    // Set up parameters for rule-based selection
    HarmonicSelectionRules::SelectionParameters params;
    params.quadrature_type = type_;
    params.construction_method = construction_method_;
    params.dimension = dimension_;
    params.scattering_order = scattering_order_;
    params.num_angles = abscissae.size();
    params.N = galerkin_N_;

    // Get rule-selected harmonics
    m_to_ell_em_map_ = HarmonicSelectionRules::SelectHarmonics(params);
  }
  else
  {
    // Standard method
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
}

void
AngularQuadrature::BuildDiscreteToMomentOperator()
{
  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  switch (construction_method_)
  {
    case OperatorConstructionMethod::STANDARD:
    {
      d2m_op_.clear();

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
      outs << "\nQuadrature d2m operator (Standard Method):\n";
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
      break;
    }

    case OperatorConstructionMethod::GALERKIN_ONE:
    {
      log.Log0Verbose1() << "Building D2M operator by inverting M2D operator using PETSc";

      // Build M2D first if not already built
      if (m2d_op_.empty())
      {
        BuildMomentToDiscreteOperator();
      }

      // Check dimensions
      if (num_angles != num_moms)
      {
        throw std::runtime_error("Cannot invert M2D operator: number of directions (" +
                                 std::to_string(num_angles) + ") != number of moments (" +
                                 std::to_string(num_moms) + ")");
      }

      try
      {
        d2m_op_ = InvertMatrix(m2d_op_);
        log.Log0Verbose1() << "D2M operator successfully computed as inverse of M2D using PETSc";
      }
      catch (const std::exception& e)
      {
        log.LogAllError() << "Failed to invert M2D operator with PETSc: " << e.what();
        throw std::runtime_error("Galerkin Method One failed: unable to invert M2D operator");
      }
      break;
    }

    case OperatorConstructionMethod::GALERKIN_THREE:
    {
      d2m_op_.clear();

      // Build matrix of exact spherical harmonic values
      std::vector<std::vector<double>> exact_harmonics(num_angles, std::vector<double>(num_moms));

      for (size_t m = 0; m < num_moms; ++m)
      {
        const auto& ell_em = m_to_ell_em_map_[m];
        for (size_t n = 0; n < num_angles; ++n)
        {
          const auto& cur_angle = abscissae[n];
          exact_harmonics[n][m] = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        }
      }

      // Orthogonalize using Householder method
      std::vector<std::vector<double>> approx_harmonics = OrthogonalizeHouseholder(exact_harmonics);

      // Renormalize to match standard spherical harmonic normalization
      for (size_t m = 0; m < num_moms; ++m)
      {
        // Compute the norm under quadrature integration
        double norm_squared = 0.0;
        for (size_t n = 0; n < num_angles; ++n)
        {
          norm_squared += approx_harmonics[n][m] * approx_harmonics[n][m] * weights[n];
        }

        // The desired normalization for spherical harmonics is 4π/(2ℓ+1)
        const auto& ell_em = m_to_ell_em_map_[m];
        double desired_norm_squared = 4.0 * M_PI / (2.0 * ell_em.ell + 1.0);

        // Scale factor to achieve desired normalization
        double scale = std::sqrt(desired_norm_squared / norm_squared);

        // Apply scaling
        for (size_t n = 0; n < num_angles; ++n)
        {
          approx_harmonics[n][m] *= scale;
        }
      }

      // Build D2M operator using approximate harmonics
      for (size_t m = 0; m < num_moms; ++m)
      {
        std::vector<double> cur_mom(num_angles);
        for (size_t n = 0; n < num_angles; ++n)
        {
          cur_mom[n] = approx_harmonics[n][m] * weights[n];
        }
        d2m_op_.push_back(cur_mom);
      }

      // Verbose printout
      std::stringstream outs;
      outs << "\nQuadrature d2m operator (Galerkin Method 3):\n";
      for (size_t n = 0; n < num_angles; ++n)
      {
        outs << std::setw(5) << n;
        for (size_t m = 0; m < num_moms; ++m)
        {
          outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << d2m_op_[m][n]
               << " ";
        }
        outs << "\n";
      }

      log.Log0Verbose1() << outs.str();
      break;
    }
  }
}

void
AngularQuadrature::BuildMomentToDiscreteOperator()
{
  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  switch (construction_method_)
  {
    case OperatorConstructionMethod::STANDARD:
    case OperatorConstructionMethod::GALERKIN_ONE:
    {
      // Both STANDARD and GALERKIN_ONE use the same M2D construction
      m2d_op_.clear();

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
      outs << "\nQuadrature m2d operator (Standard Method):\n";
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
      break;
    }

    case OperatorConstructionMethod::GALERKIN_THREE:
    {
      m2d_op_.clear();

      // Build matrix of exact spherical harmonic values
      std::vector<std::vector<double>> exact_harmonics(num_angles, std::vector<double>(num_moms));

      for (size_t m = 0; m < num_moms; ++m)
      {
        const auto& ell_em = m_to_ell_em_map_[m];
        for (size_t n = 0; n < num_angles; ++n)
        {
          const auto& cur_angle = abscissae[n];
          exact_harmonics[n][m] = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        }
      }

      // Orthogonalize using Householder method
      std::vector<std::vector<double>> approx_harmonics = OrthogonalizeHouseholder(exact_harmonics);

      // Renormalize to match standard spherical harmonic normalization
      for (size_t m = 0; m < num_moms; ++m)
      {
        // Compute the norm under quadrature integration
        double norm_squared = 0.0;
        for (size_t n = 0; n < num_angles; ++n)
        {
          norm_squared += approx_harmonics[n][m] * approx_harmonics[n][m] * weights[n];
        }

        const auto& ell_em = m_to_ell_em_map_[m];
        double desired_norm_squared = 4.0 * M_PI / (2.0 * ell_em.ell + 1.0);
        double scale = std::sqrt(desired_norm_squared / norm_squared);

        for (size_t n = 0; n < num_angles; ++n)
        {
          approx_harmonics[n][m] *= scale;
        }
      }

      // Build M2D operator using approximate harmonics
      const auto normalization = std::accumulate(weights.begin(), weights.end(), 0.0);

      for (size_t m = 0; m < num_moms; ++m)
      {
        std::vector<double> cur_mom(num_angles);
        const auto& ell_em = m_to_ell_em_map_[m];

        for (size_t n = 0; n < num_angles; ++n)
        {
          cur_mom[n] = ((2.0 * ell_em.ell + 1.0) / normalization) * approx_harmonics[n][m];
        }
        m2d_op_.push_back(cur_mom);
      }

      // Verbose printout
      std::stringstream outs;
      outs << "\nQuadrature m2d operator (Galerkin Method 3):\n";
      for (size_t n = 0; n < num_angles; ++n)
      {
        outs << std::setw(5) << n;
        for (size_t m = 0; m < num_moms; ++m)
        {
          outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << m2d_op_[m][n]
               << " ";
        }
        outs << "\n";
      }

      log.Log0Verbose1() << outs.str();
      break;
    }
  }
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
