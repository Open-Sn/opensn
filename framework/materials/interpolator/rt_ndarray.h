// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/interpolator/max_dim.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace opensn
{

/// Multi-dimensional array of which number of dimensions can be determined at run-time.
class RT_NDArray
{
public:
  RT_NDArray() = default;
  RT_NDArray(std::span<const std::uint32_t> shape, std::size_t itemsize);

  RT_NDArray(const std::vector<std::uint32_t>& shape, const std::size_t itemsize)
    : RT_NDArray(std::span<const std::uint32_t>(shape.data(), shape.size()), itemsize)
  {
  }

  RT_NDArray(const std::initializer_list<std::uint32_t>& shape, const std::size_t itemsize)
    : RT_NDArray(std::span<const std::uint32_t>(shape.begin(), shape.size()), itemsize)
  {
  }

  template <size_t N>
  RT_NDArray(const std::array<std::uint32_t, N>& shape, const std::size_t itemsize)
    : RT_NDArray(std::span<const std::uint32_t>(shape.data(), shape.size()), itemsize)
  {
  }

  /// @name Get (get without check)
  /// @{
  constexpr std::span<double> Get(std::span<const std::uint32_t> index) noexcept
  {
    const auto item_offset = ComputeItemOffsetUnchecked(index);
    return std::span<double>{storage_.data() + item_offset * itemsize_, itemsize_};
  }
  constexpr std::span<const double> Get(std::span<const std::uint32_t> index) const noexcept
  {
    const auto item_offset = ComputeItemOffsetUnchecked(index);
    return std::span<const double>{storage_.data() + item_offset * itemsize_, itemsize_};
  }
  template <size_t N>
  std::span<double> Get(const std::array<std::uint32_t, N>& index)
  {
    return Get(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  template <size_t N>
  std::span<const double> Get(const std::array<std::uint32_t, N>& index) const
  {
    return Get(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  std::span<double> Get(const std::vector<std::uint32_t>& index)
  {
    return Get(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  std::span<const double> Get(const std::vector<std::uint32_t>& index) const
  {
    return Get(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  std::span<double> Get(const std::initializer_list<std::uint32_t>& index)
  {
    return Get(std::span<const std::uint32_t>(index.begin(), index.size()));
  }
  std::span<const double> Get(const std::initializer_list<std::uint32_t>& index) const
  {
    return Get(std::span<const std::uint32_t>(index.begin(), index.size()));
  }
  /// @}

  /// @name At (get with check)
  /// @{
  std::span<double> At(std::span<const std::uint32_t> index);
  std::span<const double> At(std::span<const std::uint32_t> index) const;
  std::span<double> At(const std::vector<std::uint32_t>& index)
  {
    return At(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  std::span<const double> At(const std::vector<std::uint32_t>& index) const
  {
    return At(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  std::span<double> At(const std::initializer_list<std::uint32_t>& index)
  {
    return At(std::span<const std::uint32_t>(index.begin(), index.size()));
  }
  std::span<const double> At(const std::initializer_list<std::uint32_t>& index) const
  {
    return At(std::span<const std::uint32_t>(index.begin(), index.size()));
  }
  template <size_t N>
  std::span<double> At(const std::array<std::uint32_t, N>& index)
  {
    return At(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  template <size_t N>
  std::span<const double> At(const std::array<std::uint32_t, N>& index) const
  {
    return At(std::span<const std::uint32_t>(index.data(), index.size()));
  }
  /// @}

private:
  void Construct(std::span<const std::uint32_t> shape, std::size_t itemsize);
  constexpr std::size_t
  ComputeItemOffsetUnchecked(std::span<const std::uint32_t> index) const noexcept
  {
    std::size_t item_offset = 0;
    for (std::uint32_t d = 0; d < ndim_; ++d)
      item_offset += index[d] * strides_[d];

    return item_offset;
  }
  std::size_t ComputeItemOffsetChecked(std::span<const std::uint32_t> index) const;

private:
  std::uint32_t ndim_ = 0;
  std::array<std::uint32_t, MAX_DIM> shape_{};
  std::array<std::size_t, MAX_DIM> strides_{};
  std::size_t itemsize_ = 0;
  std::vector<double> storage_{};
};

} // namespace opensn
