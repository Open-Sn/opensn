// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/interpolator/cartesian_grid.h"
#include "framework/materials/interpolator/max_dim.h"

#include <algorithm>
#include <array>
#include <numeric>
#include <stdexcept>

namespace opensn
{

CartesianGrid::CartesianGrid(const std::vector<std::vector<double>>& grid_data)
{
  if (grid_data.size() > MAX_DIM)
    throw std::runtime_error("ndim > MAX_DIM.");

  std::array<std::uint32_t, MAX_DIM> shape{};
  std::array<const double*, MAX_DIM> points{};
  for (std::size_t d = 0; d < grid_data.size(); ++d)
  {
    shape[d] = static_cast<std::uint32_t>(grid_data[d].size());
    points[d] = grid_data[d].data();
  }

  Construct(grid_data.size(), shape.data(), points.data());
}

void
CartesianGrid::Construct(std::uint32_t ndim,
                         const std::uint32_t* shape,
                         const double* const points[])
{
  if (ndim == 0)
    return;

  if (ndim > MAX_DIM)
    throw std::runtime_error("ndim > MAX_DIM.");
  if (shape == nullptr)
    throw std::runtime_error("ndim > 0, but shape is a nulptr.");
  if (std::any_of(shape, shape + ndim, [](const std::uint32_t size) { return size == 0; }))
    throw std::runtime_error("shape cannot be zero.");
  if (points == nullptr)
    throw std::runtime_error("ndim > 0, but points is a nullptr.");

  if (std::any_of(shape, shape + ndim, [](std::uint32_t s) { return s < 2; }))
    throw std::runtime_error("Interpolation requires at least 2 points per dimension.");

  ndim_ = ndim;
  shape_.fill(0);
  offset_.fill(0);
  std::copy_n(shape, ndim, shape_.begin());
  std::exclusive_scan(shape_.begin(), shape_.begin() + ndim, offset_.begin(), std::uint32_t{0});

  size_ = std::accumulate(shape, shape + ndim, std::size_t{1}, std::multiplies<>{});
  points_.clear();
  points_.reserve(std::accumulate(shape, shape + ndim, std::uint32_t{0}));

  for (std::uint32_t d = 0; d < ndim; ++d)
  {
    const double* point_data = points[d];
    if (point_data == nullptr)
      throw std::runtime_error("Got nullptr.");

    for (std::uint32_t i = 1; i < shape_[d]; ++i)
    {
      if (point_data[i - 1] == point_data[i])
        throw std::runtime_error("Grid points must be unique in each dimension.");
      if (point_data[i - 1] > point_data[i])
        throw std::runtime_error("Grid points per-dimension must be sorted in ascending order.");
    }

    points_.insert(points_.end(), point_data, point_data + shape_[d]);
  }
}

} // namespace opensn
