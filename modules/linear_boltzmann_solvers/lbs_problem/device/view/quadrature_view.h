// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/inline_macro.h"
#include <array>
#include <cstdint>
#include <cstddef>

namespace opensn
{

// NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)

/// Direction view.
struct DirectionView
{
  OPENSN_INLINE_HOST_DEV DirectionView() : omega{}, weight(0.0), m2d(nullptr), d2m(nullptr) {}
  OPENSN_INLINE_HOST_DEV void Update(const double* direction_data, const std::uint32_t& num_moments)
  {
    omega[0] = direction_data[0];
    omega[1] = direction_data[1];
    omega[2] = direction_data[2];
    weight = direction_data[3];
    m2d = direction_data + 4;
    d2m = m2d + num_moments;
  }

  std::array<double, 3> omega;
  double weight;
  const double* m2d;
  const double* d2m;
};

/// Quadrature view.
struct QuadratureView
{
  OPENSN_INLINE_HOST_DEV QuadratureView(const char* quad_data)
    : num_angles(0),
      num_moments(0),
      direction_data(
        reinterpret_cast<const double*>(reinterpret_cast<const std::uint32_t*>(quad_data) + 2))
  {
    // number of angles and number of moments
    const auto* num_angles_and_moments_data = reinterpret_cast<const std::uint32_t*>(quad_data);
    num_angles = num_angles_and_moments_data[0];
    num_moments = num_angles_and_moments_data[1];
  }

  OPENSN_INLINE_HOST_DEV void GetDirectionView(DirectionView& direction,
                                               const std::uint32_t& index) const
  {
    direction.Update(direction_data + static_cast<std::size_t>(index) * (4 + 2 * num_moments),
                     num_moments);
  }

  std::uint32_t num_angles;
  std::uint32_t num_moments;
  const double* direction_data;
};

// NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

} // namespace opensn
