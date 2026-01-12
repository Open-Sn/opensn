// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/math/functions/function.h"
#include "framework/data_types/ndarray.h"
#include <memory>

namespace opensn
{

/// Arbitrary incident angular flux specified by a user-supplied function.
///
/// The AngularFluxFunction is expected to return the incoming angular flux for a given
/// energy group index and quadrature direction index. The returned values are uniform for
/// all faces with this boundary id.
class ArbitraryBoundary : public SweepBoundary
{
public:
  ArbitraryBoundary(size_t num_groups,
                    std::shared_ptr<AngularFluxFunction> angular_flux_function,
                    CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(LBSBoundaryType::ARBITRARY, num_groups, coord_type),
      angular_flux_function_(std::move(angular_flux_function)),
      boundary_flux_()
  {
  }

  double* PsiIncoming(std::uint32_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int group_num) override
  {
    AllocateSpace(angle_num);
    return &boundary_flux_(angle_num, group_num);
  }

private:
  std::shared_ptr<AngularFluxFunction> angular_flux_function_;
  NDArray<double, 2> boundary_flux_;
  std::size_t num_angles_ = 0;

  void AllocateSpace(unsigned int angle_num)
  {
    const auto required_angles = angle_num + 1;
    if (required_angles <= num_angles_ or not angular_flux_function_)
      return;

    num_angles_ = required_angles;
    boundary_flux_.resize(std::array<size_t, 2>{num_angles_, num_groups_});

    for (std::size_t angle = 0; angle < num_angles_; ++angle)
    {
      for (std::size_t group = 0; group < num_groups_; ++group)
      {
        const double value =
          (*angular_flux_function_)(static_cast<int>(group), static_cast<int>(angle));
        boundary_flux_(angle, group) = value;
      }
    }
  }
};

} // namespace opensn
