// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <vector>
#include <limits>
#include <cmath>
#include <iostream>

namespace opensn
{

/// Quarter-range isotropic incident fluxes on a boundary.
/// For x- boundary: constant flux for directions with μ>0 AND η>0, zero otherwise
/// For y- boundary: constant flux for directions with η>0 AND μ>0, zero otherwise
class QuarterRangeIsotropicBoundary : public SweepBoundary
{
private:
  std::vector<std::vector<double>> boundary_flux_; // [group][angle]
  Vector3 boundary_normal_;

public:
  QuarterRangeIsotropicBoundary(size_t num_groups,
                                std::vector<double> group_strength,
                                const Vector3& boundary_normal,
                                const std::vector<Vector3>& omega_directions,
                                CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(LBSBoundaryType::QUARTER_RANGE_ISOTROPIC, num_groups, coord_type),
      boundary_normal_(boundary_normal)
  {
    const size_t num_angles = omega_directions.size();
    boundary_flux_.resize(num_groups, std::vector<double>(num_angles, 0.0));

    // For each group, compute quarter-range isotropic flux
    for (size_t g = 0; g < num_groups; ++g)
    {
      const double strength = group_strength[g];

      // Quarter-range isotropic flux
      // For quarter-range (μ>0, η>0), the solid angle is π steradians
      const double quarter_flux = strength / M_PI;

      size_t num_in_quarter = 0; // Diagnostic counter
      for (size_t n = 0; n < num_angles; ++n)
      {
        const auto& omega = omega_directions[n];
        bool in_quarter_range = false;

        // Determine if direction is in quarter-range based on boundary normal
        // The quarter-range is defined as μ ∈ [0,1], η ∈ [0,1] (first quadrant in μ-η plane)

        const double tol = 1.0e-8;

        // x- face (normal ≈ [-1, 0, 0]): incoming particles have μ>0, quarter-range: μ>0 AND η>0
        if (std::abs(boundary_normal_.x + 1.0) < tol && std::abs(boundary_normal_.y) < tol &&
            std::abs(boundary_normal_.z) < tol)
        {
          in_quarter_range = (omega.x > 0.0 && omega.y > 0.0);
        }
        // x+ face (normal ≈ [1, 0, 0]): incoming particles have μ<0, quarter-range: μ<0 AND η>0
        else if (std::abs(boundary_normal_.x - 1.0) < tol && std::abs(boundary_normal_.y) < tol &&
                 std::abs(boundary_normal_.z) < tol)
        {
          in_quarter_range = (omega.x < 0.0 && omega.y > 0.0);
        }
        // y- face (normal ≈ [0, -1, 0]): incoming particles have η>0, quarter-range: μ>0 AND η>0
        else if (std::abs(boundary_normal_.y + 1.0) < tol && std::abs(boundary_normal_.x) < tol &&
                 std::abs(boundary_normal_.z) < tol)
        {
          in_quarter_range = (omega.x > 0.0 && omega.y > 0.0);
        }
        // y+ face (normal ≈ [0, 1, 0]): incoming particles have η<0, quarter-range: μ>0 AND η<0
        else if (std::abs(boundary_normal_.y - 1.0) < tol && std::abs(boundary_normal_.x) < tol &&
                 std::abs(boundary_normal_.z) < tol)
        {
          in_quarter_range = (omega.x > 0.0 && omega.y < 0.0);
        }
        // z- face (normal ≈ [0, 0, -1]): incoming particles have ξ>0, quarter-range: μ>0 AND η>0
        else if (std::abs(boundary_normal_.z + 1.0) < tol && std::abs(boundary_normal_.x) < tol &&
                 std::abs(boundary_normal_.y) < tol)
        {
          in_quarter_range = (omega.x > 0.0 && omega.y > 0.0 && omega.z > 0.0);
        }
        // z+ face (normal ≈ [0, 0, 1]): incoming particles have ξ<0, quarter-range: μ>0 AND η>0
        else if (std::abs(boundary_normal_.z - 1.0) < tol && std::abs(boundary_normal_.x) < tol &&
                 std::abs(boundary_normal_.y) < tol)
        {
          in_quarter_range = (omega.x > 0.0 && omega.y > 0.0 && omega.z < 0.0);
        }

        // Set boundary flux value
        boundary_flux_[g][n] = in_quarter_range ? quarter_flux : 0.0;
        if (in_quarter_range)
          ++num_in_quarter;
      }

      // Diagnostic output for first group only
      if (g == 0)
      {
        std::cout << "QuarterRangeIsotropic: " << num_in_quarter << " out of " << num_angles
                  << " directions in quarter-range" << std::endl;
        std::cout << "  Boundary normal: [" << boundary_normal_.x << ", " << boundary_normal_.y
                  << ", " << boundary_normal_.z << "]" << std::endl;
        std::cout << "  Angular flux value: " << quarter_flux << std::endl;
      }
    }
  }

  double* PsiIncoming(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int group_num) override
  {
    return &boundary_flux_[group_num][angle_num];
  }
};

} // namespace opensn