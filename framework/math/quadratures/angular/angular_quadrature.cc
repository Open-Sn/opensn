// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/quadratures/angular/harmonic_selection_rules.h"
#include "framework/math/math.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <limits>

namespace opensn
{

double
AngularQuadrature::GetWeightSum() const
{
  return std::accumulate(weights_.begin(), weights_.end(), 0.0);
}

void
AngularQuadrature::MakeHarmonicIndices()
{
  m_to_ell_em_map_.clear();

  // For Galerkin methods with Product Quadrature, use rule-based selection
  if ((construction_method_ == OperatorConstructionMethod::GALERKIN_ONE or
       construction_method_ == OperatorConstructionMethod::GALERKIN_THREE))
  {
    // Set up parameters for rule-based selection
    HarmonicSelectionRules::SelectionParameters params;
    params.quadrature_type = type_;
    params.construction_method = construction_method_;
    params.dimension = dimension_;
    // For GALERKIN_ONE, the scattering order is determined by the selection rules
    // so that the resulting operator is a square matrix (num_angles == num_moms).
    // The user-provided scattering_order is ignored.
    params.scattering_order = (construction_method_ == OperatorConstructionMethod::GALERKIN_ONE)
                                ? std::numeric_limits<unsigned int>::max()
                                : scattering_order_;
    params.num_angles = GetNumAngles();
    params.quadrature_order = GetQuadratureOrder();
    params.n_polar = GetNumPolarAngles();
    params.n_azimuthal = GetNumAzimuthalAngles();

    // Get rule-selected harmonics
    m_to_ell_em_map_ = HarmonicSelectionRules::SelectHarmonics(params);

    // For GALERKIN_ONE, update scattering_order_ to the maximum ell in the selected harmonics
    if (construction_method_ == OperatorConstructionMethod::GALERKIN_ONE and
        not m_to_ell_em_map_.empty())
    {
      unsigned int max_ell = 0;
      for (const auto& hi : m_to_ell_em_map_)
        max_ell = std::max(max_ell, hi.ell);
      scattering_order_ = max_ell;
    }
  }
  else
  {
    // Standard method
    log.Log0Verbose1() << "AngularQuadrature::MakeHarmonicIndices - STANDARD method";
    log.Log0Verbose1() << "  Dimension: " << dimension_
                       << ", Scattering order: " << scattering_order_;

    if (dimension_ == 1)
    {
      for (auto ell = 0U; ell <= scattering_order_; ++ell)
        m_to_ell_em_map_.emplace_back(ell, 0);
    }
    else if (dimension_ == 2)
    {
      log.Log0Verbose1() << "  2D Harmonic selection (STANDARD):";
      const auto max_ell = static_cast<int>(scattering_order_);
      for (int ell = 0; ell <= max_ell; ++ell)
      {
        std::stringstream ss;
        ss << "    ell=" << ell << ": m = [";
        bool first = true;
        for (int m = -ell; m <= ell; m += 2)
        {
          m_to_ell_em_map_.emplace_back(static_cast<unsigned int>(ell), m);
          if (!first)
            ss << ", ";
          ss << (m >= 0 ? "+" : "") << m;
          first = false;
        }
        ss << "]";
        log.Log0Verbose1() << ss.str();
      }
    }
    else if (dimension_ == 3)
    {
      const auto max_ell = static_cast<int>(scattering_order_);
      for (int ell = 0; ell <= max_ell; ++ell)
        for (int m = -ell; m <= ell; ++m)
          m_to_ell_em_map_.emplace_back(static_cast<unsigned int>(ell), m);
    }

    log.Log0Verbose1() << "  Total harmonics selected: " << m_to_ell_em_map_.size();
  }
}

void
AngularQuadrature::BuildDiscreteToMomentOperator()
{
  const size_t num_angles = GetNumAngles();
  const auto num_moms = m_to_ell_em_map_.size();

  switch (construction_method_)
  {
    case OperatorConstructionMethod::STANDARD:
    {
      d2m_op_.resize({num_angles, num_moms});
      std::fill(d2m_op_.begin(), d2m_op_.end(), 0.0);

      for (size_t n = 0; n < num_angles; ++n)
      {
        const auto& cur_angle = abscissae_[n];
        auto w = weights_[n];

        for (size_t m = 0; m < num_moms; ++m)
        {
          const auto& ell_em = m_to_ell_em_map_[m];
          auto value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
          d2m_op_(n, m) = value * w;
        }
      }

      // Verbose printout
      std::stringstream outs;
      outs << "\nQuadrature d2m operator (Standard Method):\n";
      for (size_t n = 0; n < num_angles; ++n)
      {
        outs << std::setw(5) << n;
        for (size_t m = 0; m < num_moms; ++m)
        {
          outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << d2m_op_(n, m)
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

        // OpenSn stores the TRANSPOSE of the D2M Matrix, rather than it directly
        d2m_op_ = Transpose(d2m_op_);
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
      // Build matrix of exact spherical harmonic values
      NDArray<double, 2> exact_harmonics({num_angles, num_moms}, 0.0);

      for (size_t m = 0; m < num_moms; ++m)
      {
        const auto& ell_em = m_to_ell_em_map_[m];
        for (size_t n = 0; n < num_angles; ++n)
        {
          const auto& cur_angle = abscissae_[n];
          exact_harmonics(n, m) = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        }
      }

      // Orthogonalize Matrix Span
      NDArray<double, 2> approx_harmonics = OrthogonalizeMatrixSpan(exact_harmonics, weights_);

      // Remove zero vectors (harmonics that are linearly dependent on the quadrature set).
      // This happens when no Galerkin rule exists and standard moments are used — some
      // harmonics cannot be resolved by the quadrature and collapse to zero after
      // orthogonalization.
      {
        const double zero_tol = 1e-14;
        std::vector<size_t> valid_cols;
        for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
        {
          double norm_sq = 0.0;
          for (size_t n = 0; n < num_angles; ++n)
            norm_sq += weights_[n] * approx_harmonics(n, m) * approx_harmonics(n, m);
          if (std::sqrt(norm_sq) > zero_tol)
            valid_cols.push_back(m);
        }

        const size_t num_pruned = m_to_ell_em_map_.size() - valid_cols.size();
        if (num_pruned > 0)
        {
          log.Log() << "GALERKIN_THREE D2M: Removing " << num_pruned
                    << " zero vector(s) from orthogonalized harmonic set.";

          std::vector<HarmonicIndices> new_map;
          NDArray<double, 2> new_harmonics({num_angles, valid_cols.size()}, 0.0);
          for (size_t i = 0; i < valid_cols.size(); ++i)
          {
            new_map.push_back(m_to_ell_em_map_[valid_cols[i]]);
            for (size_t n = 0; n < num_angles; ++n)
              new_harmonics(n, i) = approx_harmonics(n, valid_cols[i]);
          }
          m_to_ell_em_map_ = std::move(new_map);
          approx_harmonics = std::move(new_harmonics);
        }
      }

      // Renormalize columns to match standard spherical harmonic normalization
      for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
      {
        // Compute the norm under quadrature integration
        double norm_squared = 0.0;
        for (size_t n = 0; n < num_angles; ++n)
        {
          norm_squared += approx_harmonics(n, m) * approx_harmonics(n, m) * weights_[n];
        }

        // The desired normalization for spherical harmonics is 1/(2*ell+1)
        const auto& ell_em = m_to_ell_em_map_[m];
        double desired_norm_squared = 1 / (2.0 * ell_em.ell + 1.0);

        // Scale factor to achieve desired normalization
        auto scale = std::sqrt(desired_norm_squared / norm_squared);

        // Apply scaling
        for (size_t n = 0; n < num_angles; ++n)
        {
          approx_harmonics(n, m) *= scale;
        }
      }

      // Build D2M operator using approximate harmonics.
      // OpenSn storage is indexed [angle][moment].
      d2m_op_.resize({num_angles, m_to_ell_em_map_.size()});
      std::fill(d2m_op_.begin(), d2m_op_.end(), 0.0);
      for (size_t n = 0; n < num_angles; ++n)
        for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
          d2m_op_(n, m) = approx_harmonics(n, m) * weights_[n];

      // Verbose printout
      std::stringstream outs;
      outs << "\nQuadrature d2m operator (Galerkin Method 3):\n";
      for (size_t n = 0; n < num_angles; ++n)
      {
        outs << std::setw(5) << n;
        for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
        {
          outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << d2m_op_(n, m)
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
  const size_t num_angles = GetNumAngles();
  const auto num_moms = m_to_ell_em_map_.size();

  switch (construction_method_)
  {
    case OperatorConstructionMethod::STANDARD:
    case OperatorConstructionMethod::GALERKIN_ONE:
    {
      // Both STANDARD and GALERKIN_ONE use the same M2D construction
      m2d_op_.resize({num_angles, num_moms});
      std::fill(m2d_op_.begin(), m2d_op_.end(), 0.0);

      const auto normalization = std::accumulate(weights_.begin(), weights_.end(), 0.0);

      for (size_t n = 0; n < num_angles; ++n)
      {
        for (size_t m = 0; m < num_moms; ++m)
        {
          const auto& ell_em = m_to_ell_em_map_[m];
          const auto& cur_angle = abscissae_[n];
          double value = ((2.0 * ell_em.ell + 1.0) / normalization) *
                         Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
          m2d_op_(n, m) = value;
        }
      }

      // Verbose printout
      std::stringstream outs;
      outs << "\nQuadrature m2d operator (Standard/Galerkin-Quadrature 1 Method):\n";
      for (size_t n = 0; n < num_angles; ++n)
      {
        outs << std::setw(5) << n;
        for (size_t m = 0; m < num_moms; ++m)
        {
          outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << m2d_op_(n, m)
               << " ";
        }
        outs << "\n";
      }

      log.Log0Verbose1() << outs.str();
      break;
    }

    case OperatorConstructionMethod::GALERKIN_THREE:
    {
      // Build matrix of exact spherical harmonic values
      NDArray<double, 2> exact_harmonics({num_angles, num_moms}, 0.0);

      for (size_t m = 0; m < num_moms; ++m)
      {
        const auto& ell_em = m_to_ell_em_map_[m];
        for (size_t n = 0; n < num_angles; ++n)
        {
          const auto& cur_angle = abscissae_[n];
          exact_harmonics(n, m) = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        }
      }

      // Orthogonalize
      NDArray<double, 2> approx_harmonics = OrthogonalizeMatrixSpan(exact_harmonics, weights_);

      // Remove zero vectors (harmonics that are linearly dependent on the quadrature set).
      // This happens when no Galerkin rule exists and standard moments are used — some
      // harmonics cannot be resolved by the quadrature and collapse to zero after
      // orthogonalization.
      {
        const double zero_tol = 1e-14;
        std::vector<size_t> valid_cols;
        for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
        {
          double norm_sq = 0.0;
          for (size_t n = 0; n < num_angles; ++n)
            norm_sq += weights_[n] * approx_harmonics(n, m) * approx_harmonics(n, m);
          if (std::sqrt(norm_sq) > zero_tol)
            valid_cols.push_back(m);
        }

        const size_t num_pruned = m_to_ell_em_map_.size() - valid_cols.size();
        if (num_pruned > 0)
        {
          log.Log() << "GALERKIN_THREE M2D: Removing " << num_pruned
                    << " zero vector(s) from orthogonalized harmonic set.";

          std::vector<HarmonicIndices> new_map;
          NDArray<double, 2> new_harmonics({num_angles, valid_cols.size()}, 0.0);
          for (size_t i = 0; i < valid_cols.size(); ++i)
          {
            new_map.push_back(m_to_ell_em_map_[valid_cols[i]]);
            for (size_t n = 0; n < num_angles; ++n)
              new_harmonics(n, i) = approx_harmonics(n, valid_cols[i]);
          }
          m_to_ell_em_map_ = std::move(new_map);
          approx_harmonics = std::move(new_harmonics);
        }
      }

      // Renormalize columns to match standard spherical harmonic normalization
      for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
      {
        // Compute the norm under quadrature integration
        double norm_squared = 0.0;
        for (size_t n = 0; n < num_angles; ++n)
        {
          norm_squared += approx_harmonics(n, m) * approx_harmonics(n, m) * weights_[n];
        }

        // The desired normalization for spherical harmonics is 1/(2*ell+1)
        const auto& ell_em = m_to_ell_em_map_[m];
        double desired_norm_squared = 1.0 / (2.0 * ell_em.ell + 1.0);
        auto scale = std::sqrt(desired_norm_squared / norm_squared);

        for (size_t n = 0; n < num_angles; ++n)
        {
          approx_harmonics(n, m) *= scale;
        }
      }

      // Build M2D operator using approximate harmonics
      const auto normalization = std::accumulate(weights_.begin(), weights_.end(), 0.0);
      m2d_op_.resize({num_angles, m_to_ell_em_map_.size()});
      std::fill(m2d_op_.begin(), m2d_op_.end(), 0.0);
      for (size_t n = 0; n < num_angles; ++n)
      {
        for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
        {
          const auto& ell_em = m_to_ell_em_map_[m];
          m2d_op_(n, m) = ((2.0 * ell_em.ell + 1.0) / normalization) * approx_harmonics(n, m);
        }
      }

      // Verbose printout
      std::stringstream outs;
      outs << "\nQuadrature m2d operator (Galerkin Method 3):\n";
      for (size_t n = 0; n < num_angles; ++n)
      {
        outs << std::setw(5) << n;
        for (size_t m = 0; m < m_to_ell_em_map_.size(); ++m)
        {
          outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << m2d_op_(n, m)
               << " ";
        }
        outs << "\n";
      }

      log.Log0Verbose1() << outs.str();
      break;
    }
  }
}

NDArray<double, 2> const&
AngularQuadrature::GetDiscreteToMomentOperator() const
{
  return d2m_op_;
}

NDArray<double, 2> const&
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
