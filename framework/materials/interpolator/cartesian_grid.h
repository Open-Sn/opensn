// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/interpolator/max_dim.h"
#include <array>
#include <vector>
#include <span>

namespace opensn
{

/**
 * Rectilinear grid used as the domain for cross-section interpolation.
 * Each dimension stores a stricly ascending sequence of coordinates. The grid owns copies of
 * those coordinates and supports up to \ref MAX_DIM dimensions.
 */
class CartesianGrid
{
public:
  /// Construct an empty grid.
  CartesianGrid() = default;

  /**
   * Construct a grid from per-dimension coordinate arrays.
   * \param ndim Number of grid dimensions.
   * \param shape Number of coordinate points in each dimension.
   * \param points Coordinate arrays, one for each dimension. The coordinates have to be unique and
   * sorted in strictly increasing order.
   * \throws std::runtime_error If the dimensionality or coordinate data are invalid.
   */
  CartesianGrid(std::uint32_t ndim, const std::uint32_t* shape, const double* points[])
  {
    Construct(ndim, shape, points);
  }

  /**
   * Construct a grid from per-dimension coordinate vectors.
   * \param grid_data Coordinate points for each dimension.
   */
  CartesianGrid(const std::vector<std::vector<double>>& grid_data);

  /// Get the coordinates for a grid dimension.
  constexpr std::span<double> GetGrid(std::uint32_t d) noexcept
  {
    return std::span<double>{points_.data() + offset_[d], shape_[d]};
  }

  /// Get the coordinates for a grid dimension.
  constexpr std::span<const double> GetGrid(std::uint32_t d) const noexcept
  {
    return std::span<const double>{points_.data() + offset_[d], shape_[d]};
  }

  /// Get the total number of grid points across all dimensions.
  constexpr std::size_t GetSize() const noexcept { return size_; }

  /// Get the number of coordinate points in each grid dimension.
  constexpr std::span<const std::uint32_t> GetShape() const noexcept
  {
    return std::span<const std::uint32_t>{shape_.data(), ndim_};
  }

protected:
  std::uint32_t ndim_ = 0;
  std::array<std::uint32_t, MAX_DIM> shape_{};
  std::array<std::uint32_t, MAX_DIM> offset_{};
  std::vector<double> points_;
  std::size_t size_ = 0;

private:
  void Construct(std::uint32_t ndim, const std::uint32_t* shape, const double* const points[]);
};

} // namespace opensn
