// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/interpolator/linear_interpolator.h"
#include "framework/simd/simd.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

namespace opensn
{

namespace
{

void
vectorized_lerp(
  const double* lhs, const double* rhs, const double rhs_factor, double* out, const std::size_t n)
{
  const double lhs_factor = 1.0 - rhs_factor;
  const Simd lhs_factor_vector(1.0 - rhs_factor);
  const Simd rhs_factor_vector(rhs_factor);

  std::size_t i = 0;
  const std::size_t vectorized_size = (n / Simd::size) * Simd::size;
  for (; i < vectorized_size; i += Simd::size)
  {
    const Simd lhs_values(lhs + i);
    const Simd rhs_values(rhs + i);
    (lhs_factor_vector * lhs_values + rhs_factor_vector * rhs_values).StoreUnaligned(out + i);
  }

  for (; i < n; ++i)
    out[i] = lhs_factor * lhs[i] + rhs_factor * rhs[i];
}

} // namespace

LinearInterpolator::LinearInterpolator(const CartesianGrid& grid,
                                       const std::vector<std::shared_ptr<MultiGroupXS>>& xs,
                                       std::uint64_t flag)
  : Interpolator(grid, xs, flag)
{
}

std::vector<double>
LinearInterpolator::EvaluateContiguous(std::span<double> state_point)
{
  const auto shape = grid_.GetShape();
  const auto ndim = shape.size();

  if (ndim == 0)
  {
    const auto values = xs_data_.Get(std::span<const std::uint32_t>());
    return std::vector<double>{values.begin(), values.end()};
  }

  std::array<std::uint32_t, MAX_DIM> lower_index{};
  std::array<double, MAX_DIM> interpolation_factors{};
  for (std::size_t d = 0; d < ndim; ++d)
  {
    const auto grid_d = grid_.GetGrid(d);

    const double state_value = state_point[d];
    if (state_value == grid_d.back())
    {
      lower_index[d] = shape[d] - 2;
      interpolation_factors[d] = 1.0;
      continue;
    }

    const auto upper_it = std::upper_bound(grid_d.begin(), grid_d.end(), state_value);
    const auto upper_idx = static_cast<std::uint32_t>(upper_it - grid_d.begin());
    lower_index[d] = upper_idx - 1;

    const double lower_value = grid_d[lower_index[d]];
    const double upper_value = grid_d[upper_idx];
    interpolation_factors[d] = (state_value - lower_value) / (upper_value - lower_value);
  }

  const auto lower_index_span = std::span<const std::uint32_t>(lower_index.data(), ndim);
  const auto itemsize = xs_data_.Get(lower_index_span).size();
  const auto num_corners = std::size_t{1} << ndim;

  std::vector<double> hypercube_values(num_corners * itemsize);
  std::array<std::uint32_t, MAX_DIM> corner_index{};
  for (std::size_t corner = 0; corner < num_corners; ++corner)
  {
    std::copy_n(lower_index.begin(), ndim, corner_index.begin());
    for (std::size_t d = 0; d < ndim; ++d)
      corner_index[d] += (corner >> (ndim - 1 - d)) & 1U;

    const auto corner_values =
      xs_data_.Get(std::span<const std::uint32_t>(corner_index.data(), ndim));
    std::copy(corner_values.begin(),
              corner_values.end(),
              hypercube_values.begin() + static_cast<std::ptrdiff_t>(corner * itemsize));
  }

  std::size_t active_corner_count = num_corners;
  for (std::size_t d = 0; d < ndim; ++d)
  {
    const double factor = interpolation_factors[d];
    const auto pair_count = active_corner_count / 2;
    const auto total_values = pair_count * itemsize;
    double* out = hypercube_values.data();
    const double* lhs = hypercube_values.data();
    const double* rhs = lhs + total_values;
    vectorized_lerp(lhs, rhs, factor, out, total_values);

    active_corner_count = pair_count;
  }

  hypercube_values.resize(itemsize);
  return hypercube_values;
}

} // namespace opensn
